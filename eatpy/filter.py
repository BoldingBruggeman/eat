import collections
from typing import Iterable, List, Optional
import argparse
import sys
import importlib
import copy
import logging

import numpy
from mpi4py import MPI

from . import shared
from . import output
from .gotm import obs
from . import _eat_filter_pdaf

class CvtHandler(shared.Plugin):
    def cvt(self, iter: int, state: numpy.ndarray, v_p: numpy.ndarray) -> numpy.ndarray:
        raise Exception('cvt called but not implemented; state shape = %s, v_p shape = %s' % (state.shape, v_p.shape,))
    def cvt_adj(self, iter: int, state: numpy.ndarray, Vv_p: numpy.ndarray) -> numpy.ndarray:
        raise Exception('cvt_adj called but not implemented; state shape = %s, Vv_p shape = %s' % (state.shape, Vv_p.shape,))
    def cvt_ens(self, iter: int, state: numpy.ndarray, v_p: numpy.ndarray) -> numpy.ndarray:
        raise Exception('cvt_ens called but not implemented; state shape = %s, v_p shape = %s, Vv_p = %s' % (state.shape, v_p.shape))
    def cvt_adj_ens(self, iter: int, state: numpy.ndarray, Vv_p: numpy.ndarray) -> numpy.ndarray:
        raise Exception('cvt_adj_ens called but not implemented; state shape = %s, Vv_p shape = %s' % (state.shape, Vv_p.shape))

class PDAF(shared.Filter):
    """Filter class that wraps PDAF."""
    def __init__(self, comm: MPI.Comm, state_size: int, ensemble_size: int, cvt_handler: Optional[CvtHandler]=None):
        model_states = _eat_filter_pdaf.initialize(comm, state_size, ensemble_size, cvt_handler or CvtHandler())
        super().__init__(model_states)

    def assimilate(self, iobs: numpy.ndarray, obs: numpy.ndarray, sds: numpy.ndarray):
        _eat_filter_pdaf.assimilate(iobs, obs, sds)

    def finalize(self):
        _eat_filter_pdaf.finalize()
        
def main(parse_args: bool=True, plugins: Iterable[shared.Plugin]=()):
    # Enable appending to plugins
    plugins = list(plugins)

    # Ensure logging goes to console
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('filter')

    obs_handler = obs.ObservationHandler(verbose=True)

    # Parse command line arguments if requested. This may append plugins, e.g., for output.
    if parse_args:
        parser = argparse.ArgumentParser()
        parser.add_argument('--output', help='NetCDF file to write forecast/analysis state to')
        parser.add_argument('-p', '--plugin', help='Name of Python plug-in class. This must inherit from eatpy.shared.Plugin.', action='append', default=[])
        obs_handler.add_arguments(parser)
        args = parser.parse_args()
        obs_handler.process_arguments(args)

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
    comm, comm_model = shared.setup_mpi(shared.COLOR_FILTER)
    nmodel = comm_model.size - comm.size
    if not nmodel:
        logger.info('Running without models')

    # Receive size of state vector from lowest-ranking model
    state_size = numpy.array(0, dtype='i4')
    if nmodel:
        comm_model.Recv(state_size, source=1, tag=MPI.ANY_TAG)
    model_states = numpy.empty((nmodel, state_size))

    # Get state layout from the model
    model_variables = collections.OrderedDict(shared.parse_memory_map('da_variables.dat'))

    # Initialize plugins. These may manipulate the state layout (add or remove variables)
    variables = copy.deepcopy(model_variables)
    for plugin in plugins:
        plugin.initialize(variables, nmodel)

    # Determine final state size (plugins may have added or removed variables)
    if not variables:
        raise Exception('No variables left to act upon (possibly due to plugins filtering out all)')
    logger.info('The filter will affect the following %i state variables:' % len(variables))
    state_size = 0
    for name, info in variables.items():
        model_loc = 'not in model'
        if 'start' in info:
            info['_data'] = model_states[:, info['start']:info['start'] + info['length']]
            model_loc = '%i:%i in model' % (info['start'], info['start'] + info['length'])
        info['start'] = state_size
        state_size += info['length']
        logger.info('- %s = %s (%s) [%i values, %s]' % (name, info['long_name'], info['units'], info['length'], model_loc))

    # Let the obsveration handler verify whether the model state and the set of observed variables are compatible
    obs_handler.initialize(variables, nmodel)

    # Create filter. This will allocate an array to hold the state of all members: f.model_states
    cvt_handler = None
    for plugin in plugins:
        if isinstance(plugin, CvtHandler):
            cvt_handler = plugin
    f = PDAF(comm, state_size, nmodel, cvt_handler)

    # Create a mapping between state fo the model and state as seen by the filter
    model2filter_state_map = []
    for info in variables.values():
        if '_data' in info:
            model2filter_state_map.append((info['_data'], f.model_states[:, info['start']:info['start'] + info['length']]))

    for time in obs_handler.observations():
        strtime = time.strftime('%Y-%m-%d %H:%M:%S').encode('ascii')
        logger.info('(-> model) {}'.format(strtime))
        MPI.Request.Waitall([comm_model.Issend([strtime, MPI.CHARACTER], dest=imodel + 1, tag=shared.TAG_TIMESTR) for imodel in range(nmodel)])

        # Start receiving state from models
        reqs = [comm_model.Irecv(model_states[imodel, :], source=imodel + 1, tag=shared.TAG_FORECAST) for imodel in range(nmodel)]

        # While we are waiting for the model, collect all observations for the current time
        obs_handler.collect_observations()

        # Ensure the model state has been received
        MPI.Request.Waitall(reqs)

        # Select the parts of the model state that the filter will act upon
        for model_state, filter_state in model2filter_state_map:
            filter_state[...] = model_state

        # Get observations mapped to to state indices
        all_obs_indices, all_obs_values, all_obs_sds = obs_handler.get_observations(model_variables, model_states, variables)

        # Allow plugins to act before analysis begins
        if plugins:
            for plugin in plugins:
                plugin.before_analysis(time, f.model_states, all_obs_indices, all_obs_values, all_obs_sds, f)

        # If we have observations, then perform assimilation. This updates f.model_states
        f.assimilate(all_obs_indices, all_obs_values, all_obs_sds)

        # Allow plugins to act before analysis state is sent back to models
        for plugin in reversed(plugins):
            plugin.after_analysis(f.model_states)

        # Select the parts of the model state that the filter has acted upon
        for model_state, filter_state in model2filter_state_map:
            model_state[...] = filter_state

        # Return state to models
        MPI.Request.Waitall([comm_model.Isend(model_states[imodel, :], dest=imodel + 1, tag=shared.TAG_ANALYSIS) for imodel in range(nmodel)])

    for imodel in range(nmodel):
        comm_model.Send([b'0000-00-00 00:00:00', MPI.CHARACTER], dest=imodel + 1, tag=shared.TAG_TIMESTR)
    logger.info('No more observations to assimilate; waiting for model to finish forecast-only remainder of the simulation...')

    # Allow plugins and filter to clean-up
    for plugin in plugins:
        plugin.finalize()
    f.finalize()

if __name__ == '__main__':
    main()
