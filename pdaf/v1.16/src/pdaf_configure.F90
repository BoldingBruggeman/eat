#include "pdaf_configure.h"
!-----------------------------------------------------------------------
!BOP
!
   program main
!
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   write(*,*) ''
   write(*,*) 'PDAF version: ', _VERSION_
   write(*,*) 'PDAF configuration options:'
#ifdef _PDAF_USE_MPI_
   write(*,*) 'Using MPI:',_MPI_VERSION_
   write(*,*) '  include:   ',_MPI_Fortran_INCLUDE_PATH_
   write(*,*) '  link:      ',_MPI_Fortran_LINK_LIBRARIES_
   write(*,*) '  libraries: ',_MPI_Fortran_LIBRARIES_
#else
   write(*,*) 'Not using MPI'
#endif
   write(*,*) 'Using BLAS:'
   write(*,*) '  libraries:    ',_BLAS_LIBRARIES_
   write(*,*) 'Using Lapack:'
   write(*,*) '  libraries:    ',_LAPACK_LIBRARIES_
   write(*,*) ''

end program

!EOC

!-----------------------------------------------------------------------
