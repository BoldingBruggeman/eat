! Copyright (C) 2021 Bolding & Bruggeman

!KB#define _F08_

module mpi_config

   !! General MPI routines shared by server, worker and obs_handler

   USE, INTRINSIC :: ISO_FORTRAN_ENV
#ifdef _F08_
   use mpi_f08
#else
   use mpi
#endif
   IMPLICIT NONE

   private

   public :: init_mpi_config, version_mpi_config
   integer, parameter, public :: color_server=1
   integer, parameter, public :: color_obs=2
   integer, parameter, public :: color_model=4

   integer, parameter, public :: signal_initialize=1
   integer, parameter, public :: signal_integrate=2
   integer, parameter, public :: signal_finalize=4
   integer, parameter, public :: signal_send_state=8
#ifdef _F08_
   TYPE(mpi_comm), public :: MPI_COMM_obs,MPI_COMM_model
#else
   integer, public :: MPI_COMM_obs,MPI_COMM_model
#endif
   integer, public :: rank, nprocs

   integer :: ierr

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_mpi_config(color)

   !! Initialize MPI, define communicators and set variables

   IMPLICIT NONE

   ! Subroutine arguments
   integer, intent(in) :: color

   ! Local variables
   CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: pname
   integer :: len
   integer :: n
   logical :: flag
!-------------------------------------------------------------------------
   call MPI_init(ierr)
   if(ierr /= MPI_SUCCESS) then
      write(0,*) 'Fatal error: unable to initialize MPI.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   ! Get number of processes
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   if (ierr /= MPI_SUCCESS) THEN
      write(error_unit,*) 'Fatal error: unable to get number of processes.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   !  Get rank of current in MPI_COMM_WORLD
   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
   if (ierr /= MPI_SUCCESS) THEN
      write(error_unit,*) 'Fatal error: unable to get RANK.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   ! Get the processor names
   call MPI_GET_PROCESSOR_NAME(pname,len,ierr)
   if(ierr /= MPI_SUCCESS) THEN
      write(error_unit,*) 'Fatal error: unable to get processor name.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   ! Setup inter/intra communicators
   if (iand(color,color_obs) == color_obs) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_obs,rank,MPI_COMM_obs,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,MPI_COMM_obs,ierr)
   end if

   if (iand(color,color_model) == color_model) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_model,rank,MPI_COMM_model,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,MPI_COMM_model,ierr)
   end if

   write(error_unit,*) 'Process ',rank,' of ',nprocs,' is alive on ',pname(1:len)

   call MPI_Barrier(MPI_COMM_WORLD,ierr)
end subroutine init_mpi_config

!-----------------------------------------------------------------------

subroutine version_mpi_config(mpi_version, library_version)

   !! Collect MPI and MPI Library version info

   IMPLICIT NONE

   ! Subroutine arguments
   character(len=*), intent(inout) :: mpi_version, library_version

   ! Local variables
   integer :: v,sv,len
!-------------------------------------------------------------------------
   call MPI_get_version(v,sv,ierr)
   write(mpi_version,'(i1,a1,i1)') v,'.',sv
   call MPI_Get_library_version(library_version,len,ierr)
end subroutine version_mpi_config

!-----------------------------------------------------------------------

end module mpi_config
