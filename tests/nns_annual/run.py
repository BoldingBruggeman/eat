import eatpy

gotm = eatpy.models.GOTM(diagnostics_in_state=('rho',))

filter = eatpy.PDAF(eatpy.pdaf.FilterType.ESTKF, 0)

#cvt = eatpy.pdaf.CvtHandler()
#cvt.dim_cvec=2
#gotm.add_plugin(cvt)
gotm.add_observations('temp[-1]', 'cci_sst.dat')

gotm.run(filter)