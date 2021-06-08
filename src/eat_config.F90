! Copyright (C) 2021 Bolding & Bruggeman

!KB#define _F08_

module eat_config

   !! General MPI routines shared by server, worker and obs_handler

   USE, INTRINSIC :: ISO_FORTRAN_ENV
#ifdef _F08_
   use mpi_f08
#else
   use mpi
#endif
   IMPLICIT NONE

   private

   public :: init_eat_config, version_mpi_config
   integer, parameter, public :: color_obs=1
   integer, parameter, public :: color_model=2
   integer, parameter, public :: color_filter=4
   integer, parameter, public :: debug=2048
   integer, parameter, public :: info=1024
   integer, parameter, public :: warn=512
   integer, parameter, public :: error=256
   integer, parameter, public :: fatal=128
   integer, parameter, public :: silent=64

#ifdef _F08_
   TYPE(mpi_comm), public :: MPI_COMM_obs,MPI_COMM_model
#else
   integer, public :: EAT_COMM_obs=MPI_COMM_NULL
   integer, public :: EAT_COMM_model=MPI_COMM_NULL
   integer, public :: EAT_COMM_filter=MPI_COMM_NULL
   integer, public :: EAT_COMM_obs_model=MPI_COMM_NULL
   integer, public :: EAT_COMM_obs_filter=MPI_COMM_NULL
   integer, public :: EAT_COMM_model_filter=MPI_COMM_NULL
#endif
   integer, public :: size_obs_comm=-1
   integer, public :: size_model_comm=-1
   integer, public :: size_filter_comm=-1
   integer, public :: size_obs_model_comm=-1
   integer, public :: size_obs_filter_comm=-1
   integer, public :: size_model_filter_comm=-1

   integer, public :: rank_obs_comm=-1
   integer, public :: rank_model_comm=-1
   integer, public :: rank_filter_comm=-1
   integer, public :: rank_obs_model_comm=-1
   integer, public :: rank_obs_filter_comm=-1
   integer, public :: rank_model_filter_comm=-1

   integer, public :: rank, nprocs

   integer :: ierr
   integer :: stderr=error_unit,stdout=output_unit

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_eat_config(color)

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
      write(stderr,*) 'Fatal error: unable to get number of processes.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   !  Get rank of current in MPI_COMM_WORLD
   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
   if (ierr /= MPI_SUCCESS) THEN
      write(stderr,*) 'Fatal error: unable to get RANK.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   ! Get the processor names
   call MPI_GET_PROCESSOR_NAME(pname,len,ierr)
   if(ierr /= MPI_SUCCESS) THEN
      write(stderr,*) 'Fatal error: unable to get processor name.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   if (color >= debug) then
      write(stderr,'(A,I4,A,I4,2A)') 'MPI_COMM_WORLD(process) ',rank,' of ',nprocs,' is alive on ',pname(1:len)
   end if

   ! Setup inter/intra communicators
   ! Observations only
   if (iand(color,color_obs) == color_obs) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_obs,rank,EAT_COMM_obs,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,EAT_COMM_obs,ierr)
   end if
   if (EAT_COMM_obs /= MPI_COMM_NULL) then
      call MPI_COMM_SIZE(EAT_COMM_obs,size_obs_comm,ierr)
      call MPI_COMM_RANK(EAT_COMM_obs,rank_obs_comm,ierr)
   else
      size_obs_comm=0
   end if

   ! Model only
   if (iand(color,color_model) == color_model) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_model,rank,EAT_COMM_model,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,EAT_COMM_model,ierr)
   end if
   if (EAT_COMM_model /= MPI_COMM_NULL) then
      call MPI_COMM_SIZE(EAT_COMM_model,size_model_comm,ierr)
      call MPI_COMM_RANK(EAT_COMM_model,rank_model_comm,ierr)
   else
      size_model_comm=0
   end if

   ! Filter only
   if (iand(color,color_filter) == color_filter) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_filter,rank,EAT_COMM_filter,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,EAT_COMM_filter,ierr)
   end if
   if (EAT_COMM_filter /= MPI_COMM_NULL) then
      call MPI_COMM_SIZE(EAT_COMM_filter,size_filter_comm,ierr)
      call MPI_COMM_RANK(EAT_COMM_filter,rank_filter_comm,ierr)
   else
      size_filter_comm=0
   end if

   ! Observations and model
   if (iand(color,color_obs) == color_obs .or. iand(color,color_model) == color_model) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_obs+color_model,rank,EAT_COMM_obs_model,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,EAT_COMM_obs_model,ierr)
   end if
   if (EAT_COMM_obs_model /= MPI_COMM_NULL) then
      call MPI_COMM_SIZE(EAT_COMM_obs_model,size_obs_model_comm,ierr)
   end if
   if (iand(color,color_obs) == color_obs .and. size_obs_comm == size_obs_model_comm) then
      EAT_COMM_obs_model = MPI_COMM_NULL
   end if
   if (iand(color,color_model) == color_model .and. size_model_comm == size_obs_model_comm) then
      EAT_COMM_obs_model = MPI_COMM_NULL
   end if
   if (EAT_COMM_obs_model /= MPI_COMM_NULL) then
      call MPI_COMM_RANK(EAT_COMM_obs_model,rank_obs_model_comm,ierr)
   else
      size_obs_model_comm=-1
   end if

   ! Observations and filter
   if (iand(color,color_obs) == color_obs .or. iand(color,color_filter) == color_filter) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_obs+color_filter,rank,EAT_COMM_obs_filter,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,EAT_COMM_obs_filter,ierr)
   end if
   if (EAT_COMM_obs_filter /= MPI_COMM_NULL) then
      call MPI_COMM_SIZE(EAT_COMM_obs_filter,size_obs_filter_comm,ierr)
   end if
   if (iand(color,color_obs) == color_obs .and. size_obs_comm == size_obs_filter_comm) then
      EAT_COMM_obs_filter = MPI_COMM_NULL
   end if
   if (iand(color,color_filter) == color_filter .and. size_filter_comm == size_obs_filter_comm) then
      EAT_COMM_obs_filter = MPI_COMM_NULL
   end if
   if (EAT_COMM_obs_filter /= MPI_COMM_NULL) then
      call MPI_COMM_RANK(EAT_COMM_obs_filter,rank_obs_filter_comm,ierr)
   else
      size_obs_filter_comm=-1
   end if

   ! Model and filter
   if (iand(color,color_model) == color_model .or. iand(color,color_filter) == color_filter) then
      call MPI_Comm_split(MPI_COMM_WORLD,color_model+color_filter,rank,EAT_COMM_model_filter,ierr)
   else
      call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,EAT_COMM_model_filter,ierr)
   end if
   if (EAT_COMM_model_filter /= MPI_COMM_NULL) then
      call MPI_COMM_SIZE(EAT_COMM_model_filter,size_model_filter_comm,ierr)
   end if
   if (iand(color,color_model) == color_model .and. size_model_comm == size_model_filter_comm) then
      EAT_COMM_model_filter = MPI_COMM_NULL
   end if
   if (iand(color,color_filter) == color_filter .and. size_filter_comm == size_model_filter_comm) then
      EAT_COMM_model_filter = MPI_COMM_NULL
   end if
   if (EAT_COMM_model_filter /= MPI_COMM_NULL) then
      call MPI_COMM_RANK(EAT_COMM_model_filter,rank_model_filter_comm,ierr)
   else
      size_model_filter_comm=-1
   end if

   if (iand(color,color_obs) == color_obs .and. color >= info) then
#if 0
      write(stderr,'(A,*(I5))') ' obs(ranks: O-OM-OF):    ',rank_obs_comm,rank_obs_model_comm,rank_obs_filter_comm
      write(stderr,'(A,*(I5))') ' obs(sizes: O-OM-OF):    ',size_obs_comm,size_obs_model_comm,size_obs_filter_comm
#else
      write(stderr,'(A,*(I5))') ' obs(ranks and sizes: O-OM-OF):    ',rank_obs_comm,rank_obs_model_comm,rank_obs_filter_comm, &
                                                                    size_obs_comm,size_obs_model_comm,size_obs_filter_comm
#endif
   end if

   if (iand(color,color_model) == color_model .and. color >= info) then
#if 0
      write(stderr,'(A,*(I5))') ' model(ranks: M-OM-MF):  ',rank_model_comm,rank_obs_model_comm,rank_model_filter_comm
      write(stderr,'(A,*(I5))') ' model(sizes: M-OM-OF):  ',size_model_comm,size_obs_model_comm,size_model_filter_comm
#else
      write(stderr,'(A,*(I5))') ' model(ranks and sizes: M-OM-MF):  ',rank_model_comm,rank_obs_model_comm,rank_model_filter_comm, &
                                                                    size_model_comm,size_obs_model_comm,size_model_filter_comm
#endif
   end if

   if (iand(color,color_filter) == color_filter .and. color >= info) then
#if 0
      write(stderr,'(A,*(I5))') ' filter(ranks: F-OF-MF): ',rank_filter_comm,    rank_obs_filter_comm,rank_model_filter_comm
      write(stderr,'(A,*(I5))') ' filter(sizes: F-OF-MF): ',size_filter_comm,    size_obs_filter_comm,size_model_filter_comm
#else
      write(stderr,'(A,*(I5))') ' filter(ranks and sizes: F-OF-MF): ',rank_filter_comm,rank_obs_filter_comm,rank_model_filter_comm,&
                                                                    size_filter_comm,size_obs_filter_comm,size_model_filter_comm
#endif
   end if

   call MPI_Barrier(MPI_COMM_WORLD,ierr)
end subroutine init_eat_config

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

end module eat_config
