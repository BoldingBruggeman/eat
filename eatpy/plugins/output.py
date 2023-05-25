from typing import List, Tuple, MutableMapping, Any, Optional
import datetime

import netCDF4
import numpy as np

from .. import shared


class NetCDF(shared.Plugin):
    """Plugin that saves the forecast and analysis states to NetCDF."""

    def __init__(self, path: str, sync_interval: int = 1):
        self.nc = netCDF4.Dataset(path, "w")
        self.variables: List[Tuple[netCDF4.Variable, Any, bool]] = []
        self.itime = 0
        self.sync_interval = sync_interval
        self.reftime: Optional[datetime.datetime] = None
        self.nctime: Optional[netCDF4.Variable] = None

    def initialize(self, variables: MutableMapping[str, Any], ensemble_size: int):
        self.nc.createDimension("da_step", 2)
        self.nc.createDimension("member", ensemble_size)
        for name, metadata in variables.items():
            for n, l in metadata["dimensions"].items():
                if n not in self.nc.dimensions:
                    self.nc.createDimension(n, l)
            dimnames = list(metadata["dimensions"].keys())[::-1]
            time_dependent = "time" in metadata["dimensions"]
            if time_dependent:
                dimnames.insert(1, "da_step")
                dimnames.insert(2, "member")
            ncvar = self.nc.createVariable(name, float, dimnames)
            ncvar.units = metadata["units"]
            ncvar.long_name = metadata["long_name"]
            self.variables.append((ncvar, metadata, time_dependent))

    def before_analysis(self, time: datetime.datetime, *args, **kwargs):
        if self.nctime is None:
            self.reftime = time
            self.nctime = self.nc.createVariable("time", float, ("time",))
            self.nctime.units = "seconds since %s" % self.reftime.strftime(
                "%Y-%m-%d %H:%M:%S"
            )
        self.nctime[self.itime] = (time - self.reftime).total_seconds()
        for ncvar, metadata, time_dependent in self.variables:
            if time_dependent:
                ncvar[self.itime, 0, ...] = metadata["data"]
            elif self.itime == 0:
                ncvar[...] = metadata["data"][0, ...]

    def after_analysis(self, *args, **kwargs):
        for ncvar, metadata, time_dependent in self.variables:
            if time_dependent:
                ncvar[self.itime, 1, ...] = metadata["data"]
        self.itime += 1
        if self.itime % self.sync_interval == 0:
            self.nc.sync()

    def finalize(self):
        self.nc.close()


class NetCDFStats(shared.Plugin):
    """Plugin that saves the ensemble mean and sd of forecast and analysis
     states to NetCDF."""

    def __init__(self, path: str, sync_interval: int = 1):
        self.nc = netCDF4.Dataset(path, "w")
        self.variables: List[Tuple[netCDF4.Variable, Optional[netCDF4.Variable], Any]] = []
        self.itime = 0
        self.sync_interval = sync_interval
        self.reftime: Optional[datetime.datetime] = None
        self.nctime: Optional[netCDF4.Variable] = None

    def initialize(self, variables: MutableMapping[str, Any], ensemble_size: int):
        self.nc.createDimension("da_step", 2)
        for name, metadata in variables.items():
            for n, l in metadata["dimensions"].items():
                if n not in self.nc.dimensions:
                    self.nc.createDimension(n, l)
            dimnames = list(metadata["dimensions"].keys())[::-1]
            if "time" in metadata["dimensions"]:
                dimnames.insert(1, "da_step")
                ncmean = self.nc.createVariable(name + '_mean', float, dimnames)
                ncmean.units = metadata["units"]
                ncmean.long_name = 'mean ' + metadata["long_name"]
                ncsd = self.nc.createVariable(name + '_sd', float, dimnames)
                ncsd.units = metadata["units"]
                ncsd.long_name = 'sd ' + metadata["long_name"]
            else:
                ncmean = self.nc.createVariable(name, float, dimnames)
                ncmean.units = metadata["units"]
                ncmean.long_name = metadata["long_name"]
                ncsd = None
            self.variables.append((ncmean, ncsd, metadata))

    def before_analysis(self, time: datetime.datetime, *args, **kwargs):
        if self.nctime is None:
            self.reftime = time
            self.nctime = self.nc.createVariable("time", float, ("time",))
            self.nctime.units = "seconds since %s" % self.reftime.strftime(
                "%Y-%m-%d %H:%M:%S"
            )
        self.nctime[self.itime] = (time - self.reftime).total_seconds()
        for ncmean, ncsd, metadata in self.variables:
            if ncsd is not None:
                ncmean[self.itime, 0, ...] = metadata["data"].mean(axis=0)
                ncsd[self.itime, 0, ...] = np.std(metadata["data"], axis=0)
            elif self.itime == 0:
                ncmean[...] = metadata["data"][0, ...]

    def after_analysis(self, *args, **kwargs):
        for ncmean, ncsd, metadata in self.variables:
            if ncsd is not None:
                ncmean[self.itime, 1, ...] = metadata["data"].mean(axis=0)
                ncsd[self.itime, 1, ...] = np.std(metadata["data"], axis=0)
        self.itime += 1
        if self.itime % self.sync_interval == 0:
            self.nc.sync()

    def finalize(self):
        self.nc.close()
