# cython: language_level=3

cimport cython
cimport mpi4py.MPI

cimport numpy
import numpy

cdef extern void init_pdaf(int comm, int state_size, int ensemble_size, int* stat) nogil
cdef extern void assimilation_pdaf() nogil
cdef extern void eat_finish_pdaf() nogil
cdef extern void set_observations(int nobs, int* iobs, double* obs) nogil
cdef extern void get_model_states(double** p, int* stat) nogil

def initialize(mpi4py.MPI.Comm comm, int state_size, int ensemble_size):
    cdef int stat
    cdef double* p
    cdef int icomm
    icomm = <int>comm.ob_mpi
    init_pdaf(icomm, state_size, ensemble_size, &stat)
    assert stat == 0, 'init_pdaf failed'
    get_model_states(&p, &stat)
    assert stat == 0, 'get_model_states failed'
    return numpy.asarray(<double[:ensemble_size, :state_size:1]> p)

def assimilate(int[::1] iobs not None, double[::1] obs not None):
    set_observations(iobs.shape[0], &iobs[0], &obs[0])
    assimilation_pdaf()

def finalize():
    eat_finish_pdaf()