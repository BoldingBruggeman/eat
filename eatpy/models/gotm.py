#!/usr/bin/env python

import datetime
import re
import collections
import logging
import copy
import os
import shutil
from typing import Optional, List, Tuple, Mapping, Any, Iterable, Iterator

import numpy as np
from mpi4py import MPI
import yaml
import netCDF4

from .. import shared


# Hack into yaml parser to preserve order of yaml nodes,
# represent NULL by empty string, skip interpretation of on/off as Boolean
def dict_representer(dumper, data):
    return dumper.represent_mapping(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items()
    )


def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))


def none_representer(self, _):
    return self.represent_scalar("tag:yaml.org,2002:null", "")


yaml_loader = yaml.SafeLoader
yaml_dumper = yaml.SafeDumper

# Do not convert on/off to bool
# [done by pyyaml according to YAML 1.1, dropped from YAML 1.2]
del yaml_loader.yaml_implicit_resolvers["o"]
del yaml_loader.yaml_implicit_resolvers["O"]

yaml.add_representer(collections.OrderedDict, dict_representer, Dumper=yaml_dumper)
yaml.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor, Loader=yaml_loader
)
yaml.add_representer(type(None), none_representer, Dumper=yaml_dumper)

# Regular expression for ISO 8601 datetimes
datetimere = re.compile(r"(\d\d\d\d).(\d\d).(\d\d) (\d\d).(\d\d).(\d\d)\s*")


def parse_memory_map(path: str):
    """Parse file describing layout of state in memory. This file is written by GOTM."""
    with open(path, "r") as f:
        for line in f:
            name, long_name, units, category, dims, start, length = line.split("\t")
            dim2length = collections.OrderedDict()
            for dim in dims.split(","):
                dimname, dimlength = dim.split("=")
                dimlength = int(dimlength)
                dim2length[dimname] = None if dimlength == -1 else dimlength
            yield name, {
                "long_name": long_name,
                "units": units,
                "dimensions": dim2length,
                "start": int(start) - 1,  # 1-based in Fortran, convert to 0-based
                "length": int(length),
            }


class File:
    def __init__(
        self, path: str, is_1d: bool = False, sd: Optional[float] = None,
    ):
        self.path = path
        self.f = open(path, "r")
        self.is_1d = is_1d
        self.iline = 0
        self.sd = sd
        self.read(first=True)

    def read(self, first: bool = False):
        if not first and self.next is None:
            # end-of-file reached previously
            return

        line = None
        while line is None:
            line = self.f.readline()
            if line == "":
                # EOF
                self.next = None
                return
            self.iline += 1
            line = line.split("#", 1)[0].strip()
            if line == "":
                # skip commented lines
                line = None

        datematch = datetimere.match(line)
        if datematch is None:
            raise Exception(
                "%s: line %i does not start with time (yyyy-mm-dd hh:mm:ss)."
                " Line contents: %s" % (self.path, self.iline, line)
            )
        refvals = map(
            int, datematch.group(1, 2, 3, 4, 5, 6)
        )  # Convert matched strings into integers
        try:
            curtime = datetime.datetime(*refvals)
        except ValueError:
            raise Exception(
                "%s line %i: %s is not a valid time"
                % (self.path, self.iline, line[:19])
            )
        data = line[datematch.end() :].rstrip("\n").split()
        if self.is_1d:
            # Depth-explicit variable (on each line: time, depth, value)
            if len(data) not in (2, 3):
                raise Exception(
                    "%s: line %i must contain two values (depth, observation) or three"
                    " values (depth, observation, sd) after the date + time, but it"
                    " contains %i values." % (self.path, self.iline, len(data))
                )
            if len(data) == 2 and self.sd is None:
                raise Exception(
                    "%s: line %i must contain three values (depth, observation, sd)"
                    " after the date + time, since the standard deviation of"
                    " observations has not been prescribed separately. However, the"
                    " line contains %i values." % (self.path, self.iline, len(data))
                )
            z = float(data[0])
            if not np.isfinite(z):
                raise Exception(
                    "%s: depth on line %i is not a valid number: %s."
                    % (self.path, self.iline, data[0])
                )
            value = float(data[1])
            sd = self.sd if self.sd is not None else float(data[2])
            self.next = curtime, z, value, sd
        else:
            # Depth-independent variable (on each line: time, value)
            if len(data) not in (1, 2):
                raise Exception(
                    "%s: line %i must contain one value (observation) or two values"
                    " (observation, sd) after the date + time, but it contains %i"
                    " values." % (self.path, self.iline, len(data))
                )
            if len(data) == 1 and self.sd is None:
                raise Exception(
                    "%s: line %i must contain two values (observation, sd)"
                    " after the date + time, since the standard deviation of"
                    " observations has not been prescribed separately. However, the"
                    " line contains %i values." % (self.path, self.iline, len(data))
                )
            value = float(data[0])
            sd = self.sd if self.sd is not None else float(data[1])
            self.next = curtime, None, value, sd


class GOTM(shared.Experiment):
    def __init__(
        self,
        start: Optional[datetime.datetime] = None,
        stop: Optional[datetime.datetime] = None,
        diagnostics_in_state: Iterable[str] = (),
        fabm_parameters_in_state: Iterable[str] = (),
        log_level: int = logging.INFO,
    ):
        self.diagnostics_in_state = diagnostics_in_state
        self.fabm_parameters_in_state = fabm_parameters_in_state

        # Note: the __init__ of the superclass will call configure_model,
        # so any variables needed there (diagnostics_in_state, fabm_parameters_in_state)
        # need to have been set before.
        super().__init__(log_level=log_level)

        if start is not None:
            self.start_time = max(self.start_time, start)
        if stop is not None:
            self.stop_time = min(self.stop_time, stop)
        self.time = None

        self.datasets: List[Tuple[str, File, int]] = []

    def configure_model(self):
        # Send names of diagnostics and FABM parameters that need to be added
        # to the model state
        parbuffers = [
            n.ljust(64).encode("ascii") for n in self.fabm_parameters_in_state
        ]
        diagbuffers = [n.ljust(64).encode("ascii") for n in self.diagnostics_in_state]
        ndiag = np.array(len(diagbuffers), dtype=np.intc)
        npar = np.array(len(parbuffers), dtype=np.intc)
        reqs = []
        for imodel in range(self.nmodel):
            reqs.append(self.comm_model.Isend(ndiag, dest=imodel + 1))
            for buf in diagbuffers:
                reqs.append(self.comm_model.Isend(buf, dest=imodel + 1))
            reqs.append(self.comm_model.Isend(npar, dest=imodel + 1))
            for buf in parbuffers:
                reqs.append(self.comm_model.Isend(buf, dest=imodel + 1))
        MPI.Request.Waitall(reqs)

        # Receive simulation start and stop
        timebuf = bytearray(19)
        self.comm_model.Recv(timebuf, source=1)
        self.start_time = datetime.datetime.strptime(
            timebuf.decode("ascii"), "%Y-%m-%d %H:%M:%S"
        )
        self.comm_model.Recv(timebuf, source=1)
        self.stop_time = datetime.datetime.strptime(
            timebuf.decode("ascii"), "%Y-%m-%d %H:%M:%S"
        )
        self.logger.info(
            "Model simulated period: %s - %s"
            % (self.start_time.isoformat(" "), self.stop_time.isoformat(" "))
        )

    def get_model_variables(self) -> Mapping[str, Any]:
        return collections.OrderedDict(parse_memory_map("da_variables.dat"))

    def add_observations(self, variable_name: str, path: str):
        """
        Args:
            variable_name: name of the model variable that is observed
                If this is a depth-explicit variable, the variable can optionally
                be followed by a depth index between square brackets
                (e.g., 0 for the bottom layer, -1 for the surface layer)
            path: path of the file with observations
                This must be a text file with 3 columns (datetime, value, sd) for
                depth-INdependent variables, or 4 columns (datetime, depth, value, sd)
                for depth-dependent variables. Values on a single line must be
                whitespace-separated. The datetime must be in ISO 8601 format without
                timezone designation, e.g., 1979-05-08 00:00:00.
        """
        # Ensure plugins are initialized, so that the final set of variables
        # that can be observed is known
        self.initialize_plugins()

        # Extract depth index from variable name, if any
        depth_index = None
        full_variable_name = variable_name
        if "[" in variable_name and variable_name[-1] == "]":
            variable_name, depth_index = variable_name[:-1].split("[")

        if variable_name not in self.variables:
            raise Exception(
                "Observed variable %s is not present in model state (after"
                " processing by plugins, if any). Available: %s"
                % (variable_name, ", ".join(sorted(self.variables)))
            )

        # Determine offset into state array
        is_1d = self.variables[variable_name]["length"] > 1
        offset = self.variables[variable_name]["start"]
        if depth_index is not None:
            offset += int(depth_index) % self.variables[variable_name]["length"]
            is_1d = False

        self.datasets.append((variable_name, File(path, is_1d=is_1d, sd=None), offset))

        self.logger.info(
            "Observations in %s mapped to model variable %s (%iD)."
            " Offset in state array: %i"
            % (path, full_variable_name, 1 if is_1d else 0, offset)
        )

    def observations(self) -> Iterator[datetime.datetime]:
        """Iterate over all observation times"""
        while True:
            # Find time of next observation
            self.time = None
            for _, obsfile, _ in self.datasets:
                if obsfile.next is not None:
                    current_time = obsfile.next[0]
                    if self.time is None or current_time < self.time:
                        self.time = current_time

            if self.time is None:
                # No more observations
                break

            if (self.start_time is None or self.time >= self.start_time) and (
                self.stop_time is None or self.time <= self.stop_time
            ):
                self.logger.debug(
                    "next observation time: %s" % self.time.isoformat(" ")
                )
                yield self.time
            else:
                self.logger.debug(
                    "skipping observations at time %s" % self.time.isoformat(" ")
                )
                for _, obsfile, _ in self.datasets:
                    while obsfile.next is not None and obsfile.next[0] == self.time:
                        obsfile.read()

    def collect_observations(self):
        """Collect observations applicable to the current time
        (the last time yielded by :meth:`observations`)."""
        self.depth_map = []
        values = []
        sds = []
        offsets = []
        for variable, obsfile, offset in self.datasets:
            zs = []
            while obsfile.next is not None and obsfile.next[0] == self.time:
                _, z, value, sd = obsfile.next
                self.logger.debug(
                    "- %s%s = %s (sd = %s)"
                    % (variable, "" if z is None else (" @ %.2f m" % z), value, sd)
                )
                values.append(value)
                sds.append(sd)
                offsets.append(offset)
                if obsfile.is_1d:
                    zs.append(z)
                obsfile.read()
            if zs:
                self.depth_map.append(
                    (variable, slice(len(values) - len(zs), len(values)), np.array(zs),)
                )
        self.values = np.array(values, dtype=float)
        self.sds = np.array(sds, dtype=float)
        self.offsets = np.array(offsets, dtype=np.intc)

    def get_observations(
        self, variables: Mapping[str, Any]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Map observations collected by :meth:`collect_observations` to the model
        state space and return them."""
        z_model = variables["z"]["model_data"][0, :]
        for model_variable, obs_slice, zs in self.depth_map:
            # Determine the index of layer that is nearest to the observation depth
            # (depth above mean sea level, typically negative)
            izs = np.abs(z_model[np.newaxis, :] - zs[:, np.newaxis]).argmin(axis=1)
            self.offsets[obs_slice] += izs
            self.logger.debug(
                "observed %s: %s (sd %s) @ %s"
                % (
                    model_variable,
                    self.values[obs_slice],
                    self.sds[obs_slice],
                    self.offsets[obs_slice],
                )
            )
        return self.offsets, self.values, self.sds


class YAMLFile:
    def __init__(self, path: str):
        with open(path) as f:
            self.root = yaml.load(f, yaml_loader)
        self.path = path

    def _find(self, variable_path: str, optional: bool = False) -> Tuple[dict, str]:
        root = self.root
        comps = variable_path.split("/")
        for comp in comps[:-1]:
            if comp not in root:
                raise KeyError("%s not found in %s" % (variable_path, self.path))
            root = root[comp]
        if comps[-1] not in root and not optional:
            raise KeyError("%s not found in %s" % (variable_path, self.path))
        return root, comps[-1]

    def __getitem__(self, key: str):
        parent, name = self._find(key)
        return parent[name]

    def __setitem__(self, key: str, value):
        parent, name = self._find(key, optional=True)
        parent[name] = value

    def __contains__(self, key: str):
        try:
            self[key]
        except KeyError:
            return False
        return True

    def get(self, key: str, default: Any=None):
        if key in self:
            return self[key]
        return default

class Ensemble:
    def __init__(self, n: int):
        self.n = n
        self.variable2values = collections.OrderedDict()
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger()

    def __getitem__(self, key: str):
        if key in self.variable2values:
            return self.variable2values[key]
        return self.template[key]


class YAMLEnsemble(Ensemble):
    def __init__(self, template_path: str, n: int, postfix: str = "_%04i"):
        super().__init__(n)
        self.template = YAMLFile(template_path)
        self.postfix = postfix

    def __setitem__(self, name: str, values):
        assert len(values) == self.n
        if name not in self.template:
            self.logger.warning(f"{name} not present in template {self.template.path}")
        if isinstance(values, np.ndarray):
            values = values.tolist()
        self.variable2values[name] = values

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        root_backup = copy.deepcopy(self.template.root)
        name, ext = os.path.splitext(self.template.path)
        self.logger.info(f"Using template {self.template.path}...")
        for variable in self.variable2values:
            value = self.template.get(variable, "NOT PRESENT")
            self.logger.info(f"  {variable}: {value}")
        for i in range(self.n):
            outpath = name + self.postfix % (i + 1) + ext
            self.logger.info(f"Writing {outpath}...")
            for variable, values in self.variable2values.items():
                self.logger.info(f"  {variable}: {values[i]}")
                self.template[variable] = values[i]
            with open(outpath, "w") as f:
                yaml.dump(self.template.root, f, yaml_dumper)
        self.template.root = root_backup


class RestartEnsemble(Ensemble):
    def __init__(self, template_path: str, n: int, postfix: str = "_%04i"):
        super().__init__(n)
        self.template_nc = netCDF4.Dataset(template_path)
        self.template = self.template_nc.variables
        self.template_path = template_path
        self.postfix = postfix

    def __setitem__(self, name: str, values):
        values = np.array(values)
        assert values.shape[0] == self.n
        if name not in self.template:
            raise Exception(f"{name} not present in template {self.template_path}")
        self.variable2values[name] = values

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.logger.info(f"Using template {self.template_path}...")
        for variable in self.variable2values:
            original = self.template[variable]
            self.logger.info(f"  {variable}: {np.min(original)} - {np.max(original)}")
        self.template_nc.close()
        name, ext = os.path.splitext(self.template_path)
        for i in range(self.n):
            outpath = name + self.postfix % (i + 1) + ext
            self.logger.info(f"Writing {outpath}...")
            shutil.copyfile(self.template_path, outpath)
            with netCDF4.Dataset(outpath, "r+") as nc:
                for variable, values in self.variable2values.items():
                    nc.variables[variable] = values[i, ...]
                    self.logger.info(
                        f"  {variable}: {values[i, ...].min()} - {values[i, ...].max()}"
                    )

