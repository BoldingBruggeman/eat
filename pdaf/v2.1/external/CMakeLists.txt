add_library(optim)
set(SRC_OPTIM
    CG+_mpi/cgfam.f
    CG+_mpi/cgsearch.f
    CG+/cgfam.f
    CG+/cgsearch.f
    LBFGS/lbfgsb.f
    LBFGS/linpack.f
    LBFGS/timer.f
   )
list(TRANSFORM SRC_OPTIM PREPEND ${PDAF_BASE}/external/ )
target_sources(optim PRIVATE
               ${SRC_OPTIM}
              )
