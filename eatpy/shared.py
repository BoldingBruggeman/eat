from typing import Mapping, Any, Tuple, Sequence, List
import datetime
import collections
import logging

import numpy as np
from mpi4py import MPI

# FOR MPI_SPLIT - MUST MATCH COLORS DEFINED IN EAT_CONFIG.F90
COLOR_OBS = 1
COLOR_MODEL = 2
COLOR_FILTER = 4

# TAGS - MUST MATCH TAGS DEFINED IN EAT_CONFIG.F90
TAG_TIMESTR = 1
TAG_NOBS = 1
TAG_IOBS = 2
TAG_OBS = 3
TAG_ANALYSIS = 1
TAG_FORECAST = 2


def setup_mpi(color: int) -> Tuple[MPI.Comm, MPI.Comm, MPI.Comm, MPI.Comm]:
    """Set up MPI communicators."""
    comm_model = MPI.COMM_WORLD.Split(
        color=color if color == COLOR_MODEL else MPI.UNDEFINED
    )
    comm_filter = MPI.COMM_WORLD.Split(
        color=color if color == COLOR_FILTER else MPI.UNDEFINED
    )

    comm_model_filter = MPI.COMM_WORLD.Split(
        color=MPI.UNDEFINED if color == COLOR_OBS else COLOR_MODEL + COLOR_FILTER,
        key=-1 if color == COLOR_FILTER else 1,
    )

    MPI.COMM_WORLD.Barrier()

    return (
        {COLOR_MODEL: comm_model, COLOR_FILTER: comm_filter}[color],
        comm_model_filter,
    )


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


class Filter:
    """Base class for filters.
    Derived classes must implement initialize and assimilate.
    """

    model_states: np.ndarray

    def initialize(
        self,
        comm: MPI.Comm,
        state_size: int,
        ensemble_size: int,
        plugins: Sequence["Plugin"],
    ):
        raise NotImplementedError

    def assimilate(self, iobs: np.ndarray, obs: np.ndarray, cvt_handler=None):
        raise NotImplementedError

    def finalize(self):
        pass


class Plugin:
    """Base class for plugins.
    initialize: receives a description of the memory layout of the state,
        including variable metadata.
    update: receives the current model forecast and analysis states.
    """

    def __init__(self, name=None):
        self.logger = logging.getLogger(
            name or "filter.plugin.%s" % self.__class__.__name__
        )

    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        pass

    def before_analysis(
        self,
        time: datetime.datetime,
        state: np.ndarray,
        iobs: np.ndarray,
        obs: np.ndarray,
        obs_sds: np.ndarray,
        filter: Filter,
    ):
        pass

    def after_analysis(self, state: np.ndarray):
        pass

    def finalize(self):
        pass


class TestPlugin(Plugin):
    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        print("TestPlugin.initialize")

    def before_analysis(self, *args, **kwargs):
        print("TestPlugin.before_analysis")

    def after_analysis(self, *args, **kwargs):
        print("TestPlugin.after_analysis")

    def finalize(self):
        print("TestPlugin.finalize")


class Controller:
    state_size: int

    def __init__(self):
        self.plugins: List[Plugin] = []

        # Ensure logging goes to console
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger()

        # Set up parallel communication
        self.comm, self.comm_model = setup_mpi(COLOR_FILTER)
        self.nmodel = self.comm_model.size - self.comm.size
        if not self.nmodel:
            self.logger.info("Running without models")

        self.configure_model()

        # Receive size of state vector from lowest-ranking model
        state_size = np.array(0, dtype="i4")
        if self.nmodel:
            self.comm_model.Recv(state_size, source=1, tag=MPI.ANY_TAG)
        self.model_states = np.empty((self.nmodel, state_size))

        self.variables = self.get_model_variables()

        # Provide all variables with pointers to the model data
        # and rename "start" key to "model_start"
        for info in self.variables.values():
            start = info["start"]
            stop = start + info["length"]
            info["data"] = info["model_data"] = self.model_states[:, start:stop]
            info["model_start"] = start
            del info["start"]

        self._plugins_initialized = False

    def add_plugin(self, plugin: Plugin):
        if self._plugins_initialized:
            raise Exception(
                "You cannot add any more plugins after observations have been added."
            )
        self.plugins.append(plugin)

    def get_model_variables(self) -> Mapping[str, Any]:
        """Get dictionary describing the variables that are part of the model state"""
        raise NotImplementedError

    def configure_model(self):
        """Custom communication with the model"""
        pass

    def initialize_plugins(self):
        """Initialize plugins.
        """
        if self._plugins_initialized:
            return

        # Plugins can manipulate the state layout (add or remove variables);
        # keep a copy of all possible variables, needed by the observation handler
        self.all_variables = self.variables.copy()
        for plugin in self.plugins:
            plugin.initialize(self.variables, self.nmodel)

        # Add any new variables created by plugins
        self.all_variables.update(self.variables)

        # Determine final state size (plugins may have added or removed variables)
        if not self.variables:
            raise Exception(
                "No variables left for the DA filter to act upon"
                " (possibly due to plugins filtering out all)"
            )
        self.logger.info(
            "The filter will affect the following %i state variables:"
            % len(self.variables)
        )
        self.state_size = 0
        for name, info in self.variables.items():
            model_loc = "not in model"
            if "model_start" in info:
                start = info["model_start"]
                stop = start + info["length"]
                model_loc = "%i:%i in model" % (start, stop,)
            info["start"] = self.state_size
            self.state_size += info["length"]
            info["stop"] = self.state_size
            self.logger.info(
                "  %s: %s (%s) [%i values, %s]"
                % (name, info["long_name"], info["units"], info["length"], model_loc)
            )

        self._plugins_initialized = True

    def run(self, filter: Filter):
        self.initialize_plugins()

        # Initialize filter.
        # This will allocate filter.model_states, which holds the state of all members
        filter.initialize(self.comm, self.state_size, self.nmodel, self.plugins)

        # Create a mapping between state of the model and state as seen by the filter
        model2filter_state_map = []
        for info in self.variables.values():
            info["filter_data"] = filter.model_states[:, info["start"] : info["stop"]]
            info["data"] = info["filter_data"]
            if "model_data" in info:
                model2filter_state_map.append((info["model_data"], info["filter_data"]))

        for time in self.observations():
            strtime = time.strftime("%Y-%m-%d %H:%M:%S").encode("ascii")
            self.logger.info("(-> model) {}".format(strtime))

            # Send time of next observation to all ensemble members,
            # so they can integrate to that point
            MPI.Request.Waitall(
                [
                    self.comm_model.Issend(
                        [strtime, MPI.CHARACTER], dest=imodel + 1, tag=TAG_TIMESTR
                    )
                    for imodel in range(self.nmodel)
                ]
            )

            # Start receiving state from models
            reqs = [
                self.comm_model.Irecv(
                    self.model_states[imodel, :], source=imodel + 1, tag=TAG_FORECAST
                )
                for imodel in range(self.nmodel)
            ]

            # While we are waiting for the model, collect all observations
            # for the current time
            self.collect_observations()

            # Ensure the model state has been received
            MPI.Request.Waitall(reqs)

            self.logger.info("Ensemble spread (root-mean-square differences):")
            for name, info in self.all_variables.items():
                d = info.get("model_data")
                if d is not None:
                    rms = np.sqrt(np.mean(np.var(d, axis=0)))
                    self.logger.info("  %s: %.3g %s" % (name, rms, info["units"]))

            # Select the parts of the model state that the filter will act upon
            for model_state, filter_state in model2filter_state_map:
                filter_state[...] = model_state

            # Get observations mapped to to state indices
            obs_indices, obs_values, obs_sds = self.get_observations(self.all_variables)
            self.logger.debug(
                "Observations: %s (sd %s) @ %s" % (obs_values, obs_sds, obs_indices)
            )

            # Allow plugins to act before analysis begins
            for plugin in self.plugins:
                plugin.before_analysis(
                    time, filter.model_states, obs_indices, obs_values, obs_sds, filter
                )

            if not np.isfinite(filter.model_states).all():
                self.logger.error(
                    "Non-finite values in ensemble state state sent to PDAF"
                    " (after plugin.before_analysis)"
                )
                raise Exception("non-finite ensemble state")
            if not np.isfinite(obs_values).all():
                self.logger.error(
                    "Non-finite values in observations sent to PDAF"
                    " (after plugin.before_analysis)"
                )
                raise Exception("non-finite observations")
            if not np.isfinite(obs_sds).all():
                self.logger.error(
                    "Non-finite values in observation s.d. sent to PDAF"
                    " (after plugin.before_analysis)"
                )
                raise Exception("non-finite observation errors")

            # Perform assimilation. This updates f.model_states
            obs_indices += 1  # convert to 1-based indices for Fortran/PDAF
            filter.assimilate(obs_indices, obs_values, obs_sds)

            # Allow plugins to act before analysis state is sent back to models
            for plugin in reversed(self.plugins):
                plugin.after_analysis(filter.model_states)

            # Select the parts of the model state that the filter has acted upon
            for model_state, filter_state in model2filter_state_map:
                model_state[...] = filter_state

            # Return state to models
            MPI.Request.Waitall(
                [
                    self.comm_model.Isend(
                        self.model_states[imodel, :], dest=imodel + 1, tag=TAG_ANALYSIS
                    )
                    for imodel in range(self.nmodel)
                ]
            )

        # Tell ensemble members to integrate to the end of the configured
        # simulation period
        for imodel in range(self.nmodel):
            self.comm_model.Send(
                [b"0000-00-00 00:00:00", MPI.CHARACTER],
                dest=imodel + 1,
                tag=TAG_TIMESTR,
            )

        self.logger.info("No more observations to assimilate")
        self.logger.info(
            "Waiting for ensemble members to finish forecast-only remainder"
            " of the simulation..."
        )

        # Allow plugins and filter to clean-up
        for plugin in self.plugins:
            plugin.finalize()
        filter.finalize()
