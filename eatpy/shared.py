from typing import Mapping, Any, Tuple
import datetime
import collections
import logging

import numpy
from mpi4py import MPI

# FOR MPI_SPLIT - MUST MATCH COLORS DEFINED IN EAT_CONFIG.F90
COLOR_OBS=1
COLOR_MODEL=2
COLOR_FILTER=4

# TAGS - MUST MATCH TAGS DEFINED IN EAT_CONFIG.F90
TAG_TIMESTR=1
TAG_NOBS=1
TAG_IOBS=2
TAG_OBS=3
TAG_ANALYSIS=1
TAG_FORECAST=2

def setup_mpi(color: int) -> Tuple[MPI.Comm, MPI.Comm, MPI.Comm, MPI.Comm]:
    """Set up MPI communicators."""
    comm_model = MPI.COMM_WORLD.Split(color=color if color == COLOR_MODEL else MPI.UNDEFINED)
    comm_filter = MPI.COMM_WORLD.Split(color=color if color == COLOR_FILTER else MPI.UNDEFINED)

    comm_model_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED if color == COLOR_OBS else COLOR_MODEL + COLOR_FILTER, key=-1 if color == COLOR_FILTER else 1)

    MPI.COMM_WORLD.Barrier()

    return {COLOR_MODEL: comm_model, COLOR_FILTER: comm_filter}[color], comm_model_filter

def parse_memory_map(path: str):
    """Parse file describing layout of state in memory. This file is written by GOTM."""
    with open(path, 'r') as f:
        for l in f:
            name, long_name, units, category, dims, start, length = l.split('\t')
            dim2length = collections.OrderedDict()
            for dim in dims.split(','):
                n, l = dim.split('=')
                l = int(l)
                dim2length[n] = None if l == -1 else l
            yield name, {'long_name': long_name, 'units': units, 'dimensions': dim2length, 'start': int(start) - 1, 'length': int(length)}


class Filter:
    """Base class for filters. Derived classes must implement assimilate."""
    def __init__(self, model_states):
        self.model_states = model_states

    def assimilate(self, iobs: numpy.ndarray, obs: numpy.ndarray):
        raise NotImplementedError

    def finalize(self):
        pass

class Plugin:
    """Base class for plugins.
    initialize: receives a description of the memory layout of the state, including variable metadata.
    update: receives the current model forecast and analysis states.
    """
    def __init__(self, name=None):
        self.logger = logging.getLogger(name or 'filter.plugin.%s' % self.__class__.__name__)

    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        pass

    def before_analysis(self, time: datetime.datetime, state: numpy.ndarray, iobs: numpy.ndarray, obs: numpy.ndarray, obs_sds: numpy.ndarray, filter: Filter):
        pass

    def after_analysis(self, state: numpy.ndarray):
        pass

    def finalize(self):
        pass

class TestPlugin(Plugin):
    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        print('TestPlugin.initialize')

    def before_analysis(self, *args, **kwargs):
        print('TestPlugin.before_analysis')

    def after_analysis(self, *args, **kwargs):
        print('TestPlugin.after_analysis')

    def finalize(self):
        print('TestPlugin.finalize')

