import eatpy

gotm = eatpy.models.GOTM()

filter = eatpy.PDAF(eatpy.pdaf.FilterType.ESTKF, 0)

#gotm.add_plugin(eatpy.plugins.output.NetCDF('da.nc'))

gotm.add_observations("temp[-1]", "cci_sst.dat")

gotm.run(filter)