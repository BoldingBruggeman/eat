from typing import Mapping, Any
import datetime
import collections

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

def setup_mpi(color: int):
    comm_obs = MPI.COMM_WORLD.Split(color=color if color == COLOR_OBS else MPI.UNDEFINED)
    comm_model = MPI.COMM_WORLD.Split(color=color if color == COLOR_MODEL else MPI.UNDEFINED)
    comm_filter = MPI.COMM_WORLD.Split(color=color if color == COLOR_FILTER else MPI.UNDEFINED)

    comm_obs_model = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED if color == COLOR_FILTER else COLOR_OBS + COLOR_MODEL, key=-1 if color == COLOR_OBS else MPI.COMM_WORLD.rank)
    comm_obs_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED if color == COLOR_MODEL else COLOR_OBS + COLOR_FILTER, key=-1 if color == COLOR_OBS else 1)
    comm_model_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED if color == COLOR_OBS else COLOR_MODEL + COLOR_FILTER, key=-1 if color == COLOR_FILTER else 1)

    MPI.COMM_WORLD.Barrier()

    return {COLOR_OBS: comm_obs, COLOR_MODEL: comm_model, COLOR_FILTER: comm_filter}[color], comm_obs_model, comm_obs_filter, comm_model_filter

def parse_memory_map(path: str):
      with open(path, 'r') as f:
         for l in f:
            name, long_name, units, category, dims, start, length = l.split('\t')
            dim2length = collections.OrderedDict()
            for dim in dims.split(','):
                n, l = dim.split('=')
                l = int(l)
                dim2length[n] = None if l == -1 else l
            yield name, {'long_name': long_name, 'units': units, 'dimensions': dim2length, 'start': int(start), 'length': int(length)}

class Plugin:
    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        pass

    def update(self, time: datetime.datetime, forecast: numpy.ndarray, analysis: numpy.ndarray):
        pass

    def finalize(self):        
        pass
