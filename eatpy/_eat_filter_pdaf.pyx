# cython: language_level=3

cimport cython
cimport mpi4py.MPI
cimport mpi4py.libmpi

cimport numpy
import numpy

cdef extern void ceat_pdaf_init(mpi4py.libmpi.MPI_Fint comm, int state_size, int ensemble_sizes, void* cvt_callback, double** p, int* stat) nogil
cdef extern void assimilation_pdaf() nogil
cdef extern void finish_pdaf() nogil
cdef extern void ceat_pdaf_set_observations(int nobs, int* iobs, double* obs, double* rms_obs) nogil

cvt_handler_ = None
state_ = None

def initialize(mpi4py.MPI.Comm comm, int state_size, int ensemble_size, cvt_handler=None):
    global cvt_handler_, state_
    cdef int stat
    cdef double* p
    cvt_handler_ = cvt_handler
    ceat_pdaf_init(mpi4py.libmpi.MPI_Comm_c2f(comm.ob_mpi), state_size, ensemble_size, &ccvt_callback, &p, &stat)
    assert stat == 0, 'init_pdaf failed'
    state_ = numpy.asarray(<double[:ensemble_size, :state_size:1]> p)
    return state_

def assimilate(int[::1] iobs not None, double[::1] obs not None, double[::1] rms_obs not None):
    ceat_pdaf_set_observations(iobs.shape[0], &iobs[0], &obs[0], &rms_obs[0])
    assimilation_pdaf()

def finalize():
    finish_pdaf()

cdef void ccvt_callback(int cb_type, int iter, int dim_p, int dim_ens, int dim_cvec_ens, double* ens_p, double* v_p, double* Vv_p) with gil:
    v_p_ = numpy.asarray(<double[:dim_cvec_ens:1]> v_p)
    Vv_p_ = numpy.asarray(<double[:dim_p:1]> Vv_p)
    ens_p_ = None
    if ens_p != NULL:
        ens_p_ = numpy.asarray(<double[:dim_ens, :dim_p:1]> ens_p)
        ens_p_.flags.writeable = False
    if cvt_handler_ is not None:
        if cb_type == 1:
            v_p_.flags.writeable = False
            result = cvt_handler_.cvt(iter, state_.mean(axis=0), v_p_)
            assert numpy.shape(result) == Vv_p_.shape, 'Array returned by cvt should have shape %s, but it has shape %s.' % (Vv_p_.shape, numpy.shape(result))
            Vv_p_[:] = result
        elif cb_type == 2:
            Vv_p_.flags.writeable = False
            result = cvt_handler_.cvt_adj(iter, state_.mean(axis=0), Vv_p_)
            assert numpy.shape(result) == v_p_.shape, 'Array returned by cvt_adj should have shape %s, but it has shape %s.' % (v_p_.shape, numpy.shape(result))
            v_p_[:] = result
        elif cb_type == 3:
            assert ens_p_ is not None
            v_p_.flags.writeable = False
            result = cvt_handler_.cvt_ens(iter, ens_p_, v_p_)
            assert numpy.shape(result) == Vv_p_.shape, 'Array returned by cvt_ens should have shape %s, but it has shape %s.' % (Vv_p_.shape, numpy.shape(result))
            Vv_p_[:] = result
        elif cb_type == 4:
            assert ens_p_ is not None
            Vv_p_.flags.writeable = False
            result = cvt_handler_.cvt_adj_ens(iter, ens_p_, Vv_p_)
            assert numpy.shape(result) == v_p_.shape, 'Array returned by cvt_adj_ens should have shape %s, but it has shape %s.' % (v_p_.shape, numpy.shape(result))
            v_p_[:] = result

