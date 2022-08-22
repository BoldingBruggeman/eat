import collections
from typing import Iterable, Optional
import argparse
import sys
import importlib
import logging

import numpy as np
from mpi4py import MPI

from . import shared
from . import output
from .gotm import obs
from . import _eat_filter_pdaf


class CvtHandler(shared.Plugin):
    def cvt(self, iter: int, state: np.ndarray, v_p: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt called but not implemented; state shape = %s, v_p shape = %s"
            % (state.shape, v_p.shape,)
        )

    def cvt_adj(self, iter: int, state: np.ndarray, Vv_p: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt_adj called but not implemented; state shape = %s, Vv_p shape = %s"
            % (state.shape, Vv_p.shape,)
        )

    def cvt_ens(self, iter: int, state: np.ndarray, v_p: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt_ens called but not implemented; state shape = %s, v_p shape = %s"
            % (state.shape, v_p.shape)
        )

    def cvt_adj_ens(self, iter: int, state: np.ndarray, Vv_p: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt_adj_ens called but not implemented; state shape = %s, Vv_p shape = %s"
            % (state.shape, Vv_p.shape)
        )


class PDAF(shared.Filter):
    """Filter class that wraps PDAF."""

    def __init__(
        self,
        comm: MPI.Comm,
        state_size: int,
        ensemble_size: int,
        cvt_handler: Optional[CvtHandler] = None,
    ):
        model_states = _eat_filter_pdaf.initialize(
            comm, state_size, ensemble_size, cvt_handler or CvtHandler()
        )
        super().__init__(model_states)

    def assimilate(self, iobs: np.ndarray, obs: np.ndarray, sds: np.ndarray):
        _eat_filter_pdaf.assimilate(iobs, obs, sds)

    def finalize(self):
        _eat_filter_pdaf.finalize()


def main(parse_args: bool = True, plugins: Iterable[shared.Plugin] = ()):
    # Enable appending to plugins
    plugins = list(plugins)

    # Ensure logging goes to console
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("filter")

    obs_handler = obs.ObservationHandler(verbose=True)

    # Parse command line arguments if requested.
    # This may append plugins, e.g., for output.
    if parse_args:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--debug", action="store_true", help="Show debug log messages"
        )
        parser.add_argument(
            "--output", help="NetCDF file to write forecast/analysis state to"
        )
        parser.add_argument(
            "-p",
            "--plugin",
            help=(
                "Name of Python plug-in class."
                " This must inherit from eatpy.shared.Plugin."
            ),
            action="append",
            default=[],
        )
        obs_handler.add_arguments(parser)
        args = parser.parse_args()
        if args.debug:
            logger.setLevel(logging.DEBUG)

        # Load user-specified plugins
        sys.path.insert(0, ".")
        for plugin in args.plugin:
            if "." not in plugin:
                raise Exception(
                    "Plugins must be be prefixed with a module name, followed by a dot,"
                    " for instance, eatpy.output.NetCDF."
                )
            smod = plugin.split("(", 1)[0].rsplit(".", 1)[0]
            mod = importlib.import_module(smod)
            p = eval(plugin[len(smod) + 1 :], mod.__dict__)
            if not isinstance(p, shared.Plugin) and not issubclass(p, shared.Plugin):
                raise Exception(
                    "%s is not an instance or subclass of shared.Plugin" % (plugin,)
                )
            if not isinstance(p, shared.Plugin):
                p = p()
            plugins.append(p)
        sys.path.pop(0)

        # Add a plugin for NetCDF output if --output was specified
        if args.output:
            plugins.insert(0, output.NetCDF(args.output))

    # Set up parallel communication
    comm, comm_model = shared.setup_mpi(shared.COLOR_FILTER)
    nmodel = comm_model.size - comm.size
    if not nmodel:
        logger.info("Running without models")

    # Receive size of state vector from lowest-ranking model
    state_size = np.array(0, dtype="i4")
    if nmodel:
        comm_model.Recv(state_size, source=1, tag=MPI.ANY_TAG)
    model_states = np.empty((nmodel, state_size))

    # Get state layout from the model
    variables = collections.OrderedDict(shared.parse_memory_map("da_variables.dat"))

    # Provide all variables with pointers to the model data
    # and rename "start" key to "model_start"
    for info in variables.values():
        start = info["start"]
        stop = start + info["length"]
        info["data"] = info["model_data"] = model_states[:, start:stop]
        info["model_start"] = start
        del info["start"]

    # Initialize plugins.
    # These may manipulate the state layout (add or remove variables);
    # keep a copy of all possible variables, needed by the observation handler
    all_variables = variables.copy()
    for plugin in plugins:
        plugin.initialize(variables, nmodel)
    all_variables.update(variables)  # add any new variables created by plugins

    # Determine final state size (plugins may have added or removed variables)
    if not variables:
        raise Exception(
            "No variables left for the DA filter to act upon"
            " (possibly due to plugins filtering out all)"
        )
    logger.info(
        "The filter will affect the following %i state variables:" % len(variables)
    )
    state_size = 0
    for name, info in variables.items():
        model_loc = "not in model"
        if "model_start" in info:
            start = info["model_start"]
            stop = start + info["length"]
            model_loc = "%i:%i in model" % (start, stop,)
        info["start"] = state_size
        state_size += info["length"]
        info["stop"] = state_size
        logger.info(
            "- %s = %s (%s) [%i values, %s]"
            % (name, info["long_name"], info["units"], info["length"], model_loc)
        )

    # Let the observation handler verify whether the model state and the set of
    # observed variables are compatible
    obs_handler.initialize(args, variables)

    # Create filter.
    # This will allocate an array to hold the state of all members: f.model_states
    cvt_handler = None
    for plugin in plugins:
        if isinstance(plugin, CvtHandler):
            cvt_handler = plugin
    f = PDAF(comm, state_size, nmodel, cvt_handler)

    # Create a mapping between state of the model and state as seen by the filter
    model2filter_state_map = []
    for info in variables.values():
        info["filter_data"] = f.model_states[:, info["start"]:info["stop"]]
        info["data"] = info["filter_data"]
        if "model_data" in info:
            model2filter_state_map.append((info["model_data"], info["filter_data"]))

    for time in obs_handler.observations():
        strtime = time.strftime("%Y-%m-%d %H:%M:%S").encode("ascii")
        logger.info("(-> model) {}".format(strtime))

        # Send time of next observation to all ensemble members,
        # so they can integrate to that point
        MPI.Request.Waitall(
            [
                comm_model.Issend(
                    [strtime, MPI.CHARACTER], dest=imodel + 1, tag=shared.TAG_TIMESTR
                )
                for imodel in range(nmodel)
            ]
        )

        # Start receiving state from models
        reqs = [
            comm_model.Irecv(
                model_states[imodel, :], source=imodel + 1, tag=shared.TAG_FORECAST
            )
            for imodel in range(nmodel)
        ]

        # While we are waiting for the model, collect all observations
        # for the current time
        obs_handler.collect_observations()

        # Ensure the model state has been received
        MPI.Request.Waitall(reqs)

        logger.info('Ensemble spread (root-mean-square differences):')
        for name, info in all_variables.items():
            d = info.get('model_data')
            if d is not None:
                rms = np.sqrt(np.mean(np.var(d, axis=0)))
                logger.info('  %s: %.3g %s' % (name, rms, info['units']))

        # Select the parts of the model state that the filter will act upon
        for model_state, filter_state in model2filter_state_map:
            filter_state[...] = model_state

        # Get observations mapped to to state indices
        obs_indices, obs_values, obs_sds = obs_handler.get_observations(all_variables)
        logger.debug(
            "Observations: %s (sd %s) @ %s" % (obs_values, obs_sds, obs_indices)
        )

        # Allow plugins to act before analysis begins
        if plugins:
            for plugin in plugins:
                plugin.before_analysis(
                    time, f.model_states, obs_indices, obs_values, obs_sds, f,
                )

        # If we have observations, then perform assimilation.
        # This updates f.model_states
        obs_indices += 1  # convert to 1-based indices for Fortran/PDAF
        f.assimilate(obs_indices, obs_values, obs_sds)

        # Allow plugins to act before analysis state is sent back to models
        for plugin in reversed(plugins):
            plugin.after_analysis(f.model_states)

        # Select the parts of the model state that the filter has acted upon
        for model_state, filter_state in model2filter_state_map:
            model_state[...] = filter_state

        # Return state to models
        MPI.Request.Waitall(
            [
                comm_model.Isend(
                    model_states[imodel, :], dest=imodel + 1, tag=shared.TAG_ANALYSIS
                )
                for imodel in range(nmodel)
            ]
        )

    # Tell ensemble members to integrate to the end of the configured simulation period
    for imodel in range(nmodel):
        comm_model.Send(
            [b"0000-00-00 00:00:00", MPI.CHARACTER],
            dest=imodel + 1,
            tag=shared.TAG_TIMESTR,
        )

    logger.info(
        "No more observations to assimilate; waiting for ensemble members to finish"
        " forecast-only remainder of the simulation..."
    )

    # Allow plugins and filter to clean-up
    for plugin in plugins:
        plugin.finalize()
    f.finalize()


if __name__ == "__main__":
    main()
