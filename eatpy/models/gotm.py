#!/usr/bin/env python

import datetime
import re
import collections
import logging
from typing import Optional, List, Tuple, Mapping, Any, Iterable, Iterator

import numpy as np
from mpi4py import MPI

from .. import shared

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
        curtime = datetime.datetime(*refvals)
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
