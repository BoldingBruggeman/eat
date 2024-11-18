# cython: language_level=3

cimport cython

# https://github.com/Unidata/netcdf4-python/blob/a371621e919b665460a71c25a70cd68491f455f1/include/mpi-compat.h#L19
cdef extern from *:
    """
    #include <mpi.h>

    #ifdef MSMPI_VER
    #define PyMPI_HAVE_MPI_Message 1
    #endif

    #if (MPI_VERSION < 3) && !defined(PyMPI_HAVE_MPI_Message)
    typedef void *PyMPI_MPI_Message;
    #define MPI_Message PyMPI_MPI_Message
    #endif
    
    #if (MPI_VERSION < 4) && !defined(PyMPI_HAVE_MPI_Session)
    typedef void *PyMPI_MPI_Session;
    #define MPI_Session PyMPI_MPI_Session
    #endif
    """

cimport mpi4py.MPI
cimport mpi4py.libmpi

#cimport numpy as np
import numpy as np

cdef extern void init_pdaf(mpi4py.libmpi.MPI_Fint comm, int state_size, int ensemble_sizes, int* stat) nogil
cdef extern void init_pdaf_with_param(int filtertype, int subtype, int* filter_param_i, int length_filter_param_i, double* filter_param_r, int length_filter_param_r, int comm, int screen, int* status_pdaf)
cdef extern void ceat_pdaf_get_ensemble_state_pointer(double** p)
cdef extern void ceat_set_3dvar_callback(void* cvt_callback) nogil
cdef extern void assimilation_pdaf() nogil
cdef extern void finish_pdaf() nogil
cdef extern void ceat_pdaf_set_observations(int nobs, int* iobs, double* obs, double* rms_obs) nogil

cvt_handler = None
state_ = None

ceat_set_3dvar_callback(&ccvt_callback)

def initialize(mpi4py.MPI.Comm comm, int state_size, int ensemble_size):
    global state_
    cdef int stat
    cdef double* p
    init_pdaf(mpi4py.libmpi.MPI_Comm_c2f(comm.ob_mpi), state_size, ensemble_size, &stat)
    assert stat == 0, 'init_pdaf failed'
    ceat_pdaf_get_ensemble_state_pointer(&p)
    assert stat == 0, 'ceat_pdaf_get_ensemble_state_pointer failed'
    state_ = np.asarray(<double[:ensemble_size, :state_size:1]> p)
    return state_

def initialize_with_params(int filtertype, int subtype, int[::1] filter_param_i, double[::1] filter_param_r, mpi4py.MPI.Comm comm, int screen):
    global state_
    cdef int status_pdaf
    cdef double* p
    init_pdaf_with_param(filtertype, subtype, &filter_param_i[0], filter_param_i.size, &filter_param_r[0], filter_param_r.size, mpi4py.libmpi.MPI_Comm_c2f(comm.ob_mpi), screen, &status_pdaf)
    assert status_pdaf == 0, 'init_pdaf_with_param failed'
    ceat_pdaf_get_ensemble_state_pointer(&p)
    assert status_pdaf == 0, 'ceat_pdaf_get_ensemble_state_pointer failed'
    state_ = np.asarray(<double[:filter_param_i[1], :filter_param_i[0]:1]> p)
    return state_

def assimilate(int[::1] iobs not None, double[::1] obs not None, double[::1] rms_obs not None):
    ceat_pdaf_set_observations(iobs.shape[0], &iobs[0], &obs[0], &rms_obs[0])
    assimilation_pdaf()

def finalize():
    finish_pdaf()

cdef void ccvt_callback(int cb_type, int iter, int dim_p, int dim_ens, int dim_cvec_ens, double* ens_p, double* v_p, double* Vv_p) with gil:
    v_p_ = np.asarray(<double[:dim_cvec_ens:1]> v_p)
    Vv_p_ = np.asarray(<double[:dim_p:1]> Vv_p)
    ens_p_ = None
    if ens_p != NULL:
        ens_p_ = np.asarray(<double[:dim_ens, :dim_p:1]> ens_p)
        ens_p_.flags.writeable = False
    if cvt_handler is not None:
        if cb_type == 1:
            v_p_.flags.writeable = False
            result = cvt_handler.cvt(iter, state_.mean(axis=0), v_p_)
            assert np.shape(result) == Vv_p_.shape, 'Array returned by cvt should have shape %s, but it has shape %s.' % (Vv_p_.shape, np.shape(result))
            Vv_p_[:] = result
        elif cb_type == 2:
            Vv_p_.flags.writeable = False
            result = cvt_handler.cvt_adj(iter, state_.mean(axis=0), Vv_p_)
            assert np.shape(result) == v_p_.shape, 'Array returned by cvt_adj should have shape %s, but it has shape %s.' % (v_p_.shape, np.shape(result))
            v_p_[:] = result
        elif cb_type == 3:
            assert ens_p_ is not None
            v_p_.flags.writeable = False
            result = cvt_handler.cvt_ens(iter, ens_p_, v_p_)
            assert np.shape(result) == Vv_p_.shape, 'Array returned by cvt_ens should have shape %s, but it has shape %s.' % (Vv_p_.shape, np.shape(result))
            Vv_p_[:] = result
        elif cb_type == 4:
            assert ens_p_ is not None
            Vv_p_.flags.writeable = False
            result = cvt_handler.cvt_adj_ens(iter, ens_p_, Vv_p_)
            assert np.shape(result) == v_p_.shape, 'Array returned by cvt_adj_ens should have shape %s, but it has shape %s.' % (v_p_.shape, np.shape(result))
            v_p_[:] = result

