
from typing import Iterable
import argparse
import sys
import importlib
import datetime

import numpy
from mpi4py import MPI

from . import shared
from . import output
from . import _eat_filter_pdaf

class Filter:
    """Base class for filters. Derived classes must implement assimilate."""
    def __init__(self, model_states):
        self.model_states = model_states

    def assimilate(self, iobs: numpy.ndarray, obs: numpy.ndarray):
        raise NotImplementedError

    def finalize(self):
        pass

class PDAF(Filter):
    """Filter class that wraps PDAF."""
    def __init__(self, comm: MPI.Comm, state_size: int, ensemble_size: int):
        model_states = _eat_filter_pdaf.initialize(comm, state_size, ensemble_size)
        super().__init__(model_states)

    def assimilate(self, iobs: numpy.ndarray, obs: numpy.ndarray):
        _eat_filter_pdaf.assimilate(iobs, obs)

    def finalize(self):
        _eat_filter_pdaf.finalize()
        
def main(parse_args: bool=True, plugins: Iterable[shared.Plugin]=()):
    # Enable appending to plugins
    plugins = list(plugins)

    # Parse command line arguments if requested. This may append plugins, e.g., for output.
    if parse_args:
        parser = argparse.ArgumentParser()
        parser.add_argument('-o', '--output', help='NetCDF file to write forrcats/analysis state to')
        parser.add_argument('-p', '--plugin', help='Name of Python plug-in class. This must inherit from eatpy.shared.Plugin.', action='append', default=[])
        args = parser.parse_args()

        # Load user-specified plugins
        sys.path.insert(0, '.')
        for plugin in args.plugin:
            assert '.' in plugin, 'Plugins must be be prefixed with a module name, followed by a dot, for instance, eatpy.output.NetCDF.'
            smod = plugin.split('(', 1)[0].rsplit('.', 1)[0]
            mod = importlib.import_module(smod)
            p = eval(plugin[len(smod) + 1:], mod.__dict__)
            if not isinstance(p, shared.Plugin) and not issubclass(p, shared.Plugin):
                raise Exception('%s is not an instance or subclass of shared.Plugin' % (plugin,))
            if not isinstance(p, shared.Plugin):
                p = p()
            plugins.append(p)
        sys.path.pop(0)

        # Add a plugin for NetCDF output if --output was specified
        if args.output:
            plugins.insert(0, output.NetCDF(args.output))

    # Set up parallel communication
    comm, _, comm_obs, comm_model = shared.setup_mpi(shared.COLOR_FILTER)
    have_obs = comm_obs.size - comm.size > 0
    nmodel = comm_model.size - comm.size
    if not have_obs:
        print('Running without observation handler')
    if not nmodel:
        print('Running without models')

    # Receive size of state vector from lowest-ranking model
    state_size = numpy.array(0, dtype='i4')
    if nmodel:
        comm_model.Recv(state_size, source=1, tag=MPI.ANY_TAG)

    # Create filter
    f = PDAF(comm, state_size, nmodel)

    # Initialize plugins
    variables = dict(shared.parse_memory_map('da_variables.dat'))
    for plugin in plugins:
        plugin.initialize(variables, nmodel)

    nobs = numpy.array(-1, dtype='i4')
    while True:
        reqs = []
        if have_obs:
            # Receive current time and number of observations from observation handler
            timestr = numpy.empty((19,), dtype='|S1')
            comm_obs.Recv(timestr, source=0, tag=shared.TAG_TIMESTR)
            comm_obs.Recv(nobs, source=0, tag=shared.TAG_NOBS)

            # A negative number of observations indicates shutdown
            if nobs < 0:
                break

            # Set up arrays for observation indices and values.
            # If > 0, receive those from observation handler.
            iobs = numpy.empty((nobs,), dtype='i4')
            obs = numpy.empty((nobs,))
            if nobs > 0:
                reqs.append(comm_obs.Irecv(iobs, source=0, tag=shared.TAG_IOBS))
                reqs.append(comm_obs.Irecv(obs, source=0, tag=shared.TAG_OBS))

        # Receive state from models
        for imodel in range(nmodel):
            reqs.append(comm_model.Irecv(f.model_states[imodel, :], source=imodel + 1, tag=shared.TAG_FORECAST))
        MPI.Request.Waitall(reqs)

        # Allow plugins to act before analysis begins
        time = datetime.datetime.strptime(timestr.tobytes().decode('ascii'), '%Y-%m-%d %H:%M:%S')
        for plugin in plugins:
            plugin.before_analysis(time, f.model_states, iobs, obs)

        # If we have observations, then perform assimilation. This updates f.model_states
        if nobs > 0:
            f.assimilate(iobs, obs)

        # Allow plugins to act before analysis state is sent back to models
        for plugin in reversed(plugins):
            plugin.after_analysis(f.model_states)

        # Send state to models
        reqs = []
        for imodel in range(nmodel):
            reqs.append(comm_model.Isend(f.model_states[imodel, :], dest=imodel + 1, tag=shared.TAG_ANALYSIS))
        MPI.Request.Waitall(reqs)

    # Allow plugins and filter to clean-up
    for plugin in plugins:
        plugin.finalize()
    f.finalize()

if __name__ == '__main__':
    main()