from typing import (
    Mapping,
    Any,
    Tuple,
    Sequence,
    List,
    Iterator,
    Optional,
    MutableMapping,
)
import datetime
import collections
import logging
import atexit

import numpy as np

import mpi4py.rc

mpi4py.rc.initialize = False
mpi4py.rc.finalize = False

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

MPI_initialized = False


def setup_mpi(color: int) -> Tuple[MPI.Comm, MPI.Comm]:
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

    def assimilate(self, iobs: np.ndarray, obs: np.ndarray, sds: np.ndarray):
        raise NotImplementedError

    def finalize(self):
        pass


class Plugin:
    """Base class for plugins.
    initialize: receives a description of the memory layout of the state,
        including variable metadata.
    update: receives the current model forecast and analysis states.
    """

    logger: logging.Logger

    def __init__(self, name: Optional[str] = None):
        self.name = name or self.__class__.__name__

    def initialize(self, variables: MutableMapping[str, Any], ensemble_size: int):
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


class Experiment:
    state_size: int

    def __init__(self, log_level=logging.INFO):
        self.plugins: List[Plugin] = []

        if not MPI.Is_initialized():
            global MPI_initialized
            MPI.Init()
            MPI_initialized = True

        # Ensure logging goes to console
        logging.basicConfig(level=log_level)
        self.logger = logging.getLogger()

        # Set up parallel communication
        self.comm, self.comm_model = setup_mpi(COLOR_FILTER)
        self.nmodel = self.comm_model.size - self.comm.size
        state_size = np.array(0, dtype=np.intc)
        if not self.nmodel:
            self.logger.info("Running without models")
            self.variables = collections.OrderedDict()
        else:
            self.configure_model()
            self.comm_model.Recv(state_size, source=1, tag=MPI.ANY_TAG)
            self.variables = collections.OrderedDict(self.get_model_variables())

        self.model_states = np.empty((self.nmodel, state_size))

        # Provide all variables with pointers to the model data
        # and rename "start" key to "model_start"
        for info in self.variables.values():
            start = info["start"]
            stop = start + info["length"]
            info["data"] = info["model_data"] = self.model_states[:, start:stop]
            info["model_start"] = start
            del info["start"]

        self._plugins_initialized = False

    def add_plugin(self, plugin: Plugin, name: Optional[str] = None):
        if not isinstance(plugin, Plugin):
            raise Exception("Plugins must be a subclass of eatpy.Plugin")
        if self._plugins_initialized:
            raise Exception(
                "You cannot add any more plugins after observations have been added."
            )
        if name is not None:
            plugin.name = name
        self.plugins.append(plugin)

    def get_model_variables(self) -> Mapping[str, Any]:
        """Get dictionary describing the variables that are part of the model state"""
        raise NotImplementedError

    def configure_model(self):
        """Custom communication with the model"""
        pass

    def initialize_plugins(self):
        """Initialize plugins."""
        if self._plugins_initialized:
            return

        # Plugins can manipulate the state layout (add or remove variables);
        # keep a copy of all possible variables, needed by the observation handler
        self.all_variables = collections.OrderedDict(self.variables)
        for plugin in self.plugins:
            if getattr(plugin, "name", None) is None:
                plugin.name = plugin.__class__.__name__
            plugin.logger = self.logger.getChild("plugins.%s" % plugin.name)
            plugin.logger.info("Initializing...")
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
                model_loc = "%i:%i in model" % (start, stop)
            info["start"] = self.state_size
            self.state_size += info["length"]
            info["stop"] = self.state_size
            self.logger.info(
                "  %s: %s (%s) [%i values, %s]"
                % (name, info["long_name"], info["units"], info["length"], model_loc)
            )

        self._plugins_initialized = True

    def observations(self) -> Iterator[datetime.datetime]:
        raise NotImplementedError

    def collect_observations(self):
        """Collect observations for the upcoming time while the models are running."""
        pass

    def get_observations(
        self, variables: Mapping[str, Any]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError

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

        timebuf = bytearray(19)

        # Create persistent MPI requests
        comm = self.comm_model
        send_time_reqs = []
        recv_state_reqs = []
        send_state_reqs = []
        for imodel in range(self.nmodel):
            rank = imodel + 1
            state = self.model_states[imodel, :]
            send_time_reqs.append(comm.Send_init(timebuf, dest=rank, tag=TAG_TIMESTR))
            recv_state_reqs.append(comm.Recv_init(state, source=rank, tag=TAG_FORECAST))
            send_state_reqs.append(comm.Send_init(state, dest=rank, tag=TAG_ANALYSIS))
        all_forecast_wait_reqs = send_state_reqs + send_time_reqs + recv_state_reqs
        all_final_wait_reqs = send_state_reqs + send_time_reqs

        # On the first iteration, we have no updated [analyzed] state yet,
        # so we skip sending that to the ensemble members
        forecast_reqs = forecast_wait_reqs = send_time_reqs + recv_state_reqs
        final_wait_reqs = send_time_reqs

        for time in self.observations():
            strtime = time.strftime("%Y-%m-%d %H:%M:%S")
            timebuf[:] = strtime.encode("ascii")
            self.logger.info(f"Forecasting until {strtime}")

            # Start forecast by ensemble members:
            # * send time of next observation, which members should integrate up to
            # * receive updated [forecasted] model state
            MPI.Prequest.Startall(forecast_reqs)

            # While we are waiting for the model, collect all observations
            # for the current time
            self.collect_observations()

            # Ensure the model-forecasted state has been received
            # On all but the first iteration, this also first waits
            # for the ensemble members to received the updated [analyzed] state
            MPI.Request.Waitall(forecast_wait_reqs)

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

            # Perform assimilation. This updates filter.model_states
            if obs_values.size > 0:
                filter.assimilate(obs_indices, obs_values, obs_sds)

            # Allow plugins to act before analysis state is sent back to models
            for plugin in reversed(self.plugins):
                plugin.after_analysis(filter.model_states)

            # Update the parts of the model state that the filter has acted upon
            for model_state, filter_state in model2filter_state_map:
                model_state[...] = filter_state

            # Start sending the the updated state
            # Also update the set of requests that the next Waitall will process,
            # so that it includes the send_state_reqs that we are now starting
            MPI.Prequest.Startall(send_state_reqs)
            forecast_wait_reqs = all_forecast_wait_reqs
            final_wait_reqs = all_final_wait_reqs

        # Tell ensemble members to integrate to the end of the configured
        # simulation period
        self.logger.info("No more observations to assimilate")
        timebuf[:] = b"0000-00-00 00:00:00"
        MPI.Prequest.Startall(send_time_reqs)

        # Allow plugins and filter to clean-up
        self.logger.info("Finalizing plugins...")
        for plugin in self.plugins:
            plugin.finalize()
        filter.finalize()

        self.logger.info(
            "Waiting for ensemble members to finish forecast-only remainder"
            " of the simulation..."
        )

        MPI.Request.Waitall(final_wait_reqs)


def finalize_MPI():
    if MPI_initialized and not MPI.Is_finalized():
        MPI.Finalize()


atexit.register(finalize_MPI)
