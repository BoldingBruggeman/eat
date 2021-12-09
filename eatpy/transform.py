from typing import Mapping, Any
import datetime

import numpy

from . import shared

class LogTransform(shared.Plugin):
    def __init__(self, *variables):
        self.variables = frozenset(variables)
        self.slices = []

    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        for name in self.variables:
            info = variables[name]
            self.slices.append(slice(info['start'], info['start'] + info['length']))

    def before_analysis(self, time: datetime.datetime, state: numpy.ndarray):
        for s in self.slices:
            state[s] = numpy.log10(state[s])

    def after_analysis(self, time: datetime.datetime, state: numpy.ndarray):
        for s in self.slices:
            state[s] = 10. ** state[s]
