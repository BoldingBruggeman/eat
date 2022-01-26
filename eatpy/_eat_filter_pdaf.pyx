# cython: language_level=3

cimport cython
cimport mpi4py.MPI
cimport mpi4py.libmpi

cimport numpy
import numpy

cdef extern void ceat_pdaf_init(mpi4py.libmpi.MPI_Fint comm, int state_size, int ensemble_sizes, double** p, int* stat) nogil
cdef extern void assimilation_pdaf() nogil
cdef extern void finish_pdaf() nogil
cdef extern void ceat_pdaf_set_observations(int nobs, int* iobs, double* obs) nogil

def initialize(mpi4py.MPI.Comm comm, int state_size, int ensemble_size):
    cdef int stat
    cdef double* p
    ceat_pdaf_init(mpi4py.libmpi.MPI_Comm_c2f(comm.ob_mpi), state_size, ensemble_size, &p, &stat)
    assert stat == 0, 'init_pdaf failed'
    return numpy.asarray(<double[:ensemble_size, :state_size:1]> p)

def assimilate(int[::1] iobs not None, double[::1] obs not None):
    ceat_pdaf_set_observations(iobs.shape[0], &iobs[0], &obs[0])
    assimilation_pdaf()

def finalize():
    finish_pdaf()