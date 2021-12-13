from typing import Mapping, Any, List, Tuple
import datetime

import numpy

from . import shared

class Log(shared.Plugin):
    def __init__(self, *variables):
        self.variables = frozenset(variables)
        self.slices: List[Tuple[int, int]] = []

    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        for name in self.variables:
            info = variables[name]
            self.slices.append((info['start'], info['start'] + info['length']))

    def before_analysis(self, time: datetime.datetime, state: numpy.ndarray, iobs: numpy.ndarray, obs: numpy.ndarray):
        for start, stop in self.slices:
            affected_obs = numpy.logical_and(iobs > start, iobs <= stop)
            if affected_obs.any():
                obs[affected_obs] = numpy.log10(obs[affected_obs])
            state[:, start:stop] = numpy.log10(state[:, start:stop])

    def after_analysis(self, state: numpy.ndarray):
        for start, stop in self.slices:
            state[:, start:stop] = 10. ** state[:, start:stop]
