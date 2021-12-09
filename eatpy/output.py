from typing import List, Tuple, Mapping, Any
import datetime

import numpy
import netCDF4

from . import shared

class NetCDF(shared.Plugin):
    def __init__(self, path: str, sync_interval: int=1):
        self.nc = netCDF4.Dataset(path, 'w')
        self.variables: List[Tuple[netCDF4.Variable, int, int]] = []
        self.itime = 0
        self.sync_interval = sync_interval

    def initialize(self, variables: Mapping[str, Any], ensemble_size: int):
        self.nc.createDimension('da_step', 2)
        self.nc.createDimension('member', ensemble_size)
        for name, metadata in variables.items():
            istart, length = metadata['start'], metadata['length']
            for n, l in metadata['dimensions'].items():
                if n not in self.nc.dimensions:
                    self.nc.createDimension(n, l)
            dimnames = list(metadata['dimensions'].keys())[::-1]
            time_dependent = 'time' in metadata['dimensions']
            if time_dependent:
                dimnames.insert(1, 'da_step')
                dimnames.insert(2, 'member')
            ncvar = self.nc.createVariable(name, float, dimnames)
            ncvar.units = metadata['units']
            ncvar.long_name = metadata['long_name']
            self.variables.append((ncvar, istart, istart + length - 1, time_dependent))

    def update(self, time: datetime.datetime, forecast: numpy.ndarray, analysis: numpy.ndarray):
        for ncvar, istart, istop, time_dependent in self.variables:
            if time_dependent:
                ncvar[self.itime, 0, ...] = forecast[:, istart - 1:istop]
                ncvar[self.itime, 1, ...] = analysis[:, istart - 1:istop]
            elif self.itime == 0:
                ncvar[...] = forecast[0, istart - 1:istop]
        self.itime += 1
        if self.itime % self.sync_interval == 0:
            self.nc.sync()

    def finalize(self):
        self.nc.close()
