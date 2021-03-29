! Copyright (C) 2021 Bolding & Bruggeman

program eat_observations

   !! A observation handler example program in Fortran.
   !! No real observation handling but illustrates communication
   !! with the server.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   IMPLICIT NONE

   integer :: ierr
   character(32), allocatable  :: obs_times(:)
   logical :: have_model=.true.
   logical :: have_filter=.true.
   integer :: nmodel=-1
!KB   integer :: nfilter=-1
!-----------------------------------------------------------------------

   call init_observations()
   call do_observations()
   call finish_observations()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_observations()

   !! Observation times are obtained from an external file - prepare to
   !! send to the server

   ! Local variables
   integer :: unit,ios
   character(len=32) :: buf
   integer :: i,n
!-----------------------------------------------------------------------
   call init_eat_config(color_obs)

   if (size_obs_comm == size_obs_model_comm) then
      write(error_unit,*) "obs(no model executable present)"
      have_model=.false.
      EAT_COMM_obs_model=MPI_COMM_NULL
   else
      call MPI_COMM_RANK(EAT_COMM_obs_model,rank_obs_model_comm,ierr)
      call MPI_COMM_SIZE(EAT_COMM_obs_model,nmodel,ierr)
      nmodel=nmodel-1 !KB - need size of obs comm
   end if

   if (size_obs_comm == size_obs_filter_comm) then
      write(error_unit,*) "obs(no filter executable present)"
      have_filter=.false.
      EAT_COMM_obs_filter=MPI_COMM_NULL
   else
      call MPI_COMM_RANK(EAT_COMM_obs_filter,rank_obs_filter_comm,ierr)
   end if

   write(error_unit,'(A,3I5)') ' obs(ranks: O-OM-OF):    ',rank_obs_comm,rank_obs_model_comm,rank_obs_filter_comm
   write(error_unit,'(A,3I5)') ' obs(sizes: O-OM-OF):    ',size_obs_comm,size_obs_model_comm,size_obs_filter_comm

   open(newunit=unit,file="obs_times.dat",status='old',action='read',iostat=ios)
   if (ios /= 0) stop 'init_obs(): unable to open obs_times.dat for reading'
   n=0
   do
      read(unit,*,iostat=ios) buf
      if (ios < 0) exit
      if(len(buf) > 0) n=n+1
   end do
   allocate(obs_times(n))
   rewind unit
   do i=1,n
      read(unit,*,iostat=ios) buf
      if (ios < 0) exit
      if(len(buf) > 0) obs_times(i) = trim(buf)
   end do
end subroutine init_observations

!-----------------------------------------------------------------------

subroutine do_observations()

   !! Handle observations and send size and observations to the server

   ! Local variables
   integer :: stat(MPI_STATUS_SIZE)
   integer :: m,n,nobs
   real(real64), allocatable :: obs(:)
   character(len=32) :: timestr,halt="0000-00-00 00:00:00"
   integer :: request
!-----------------------------------------------------------------------
   do n=1,size(obs_times)
      if (have_model) then
         do m=1,nmodel
            call MPI_SEND(obs_times(n),19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr)
         end do
      else
         write(error_unit,*) 'obs(-> time)  ',trim(obs_times(n))
      end if
      if (have_filter) then
         nobs=10000*n
         write(error_unit,'(A,I6)') ' obs(-> nobs)    ',nobs
         call MPI_SEND(nobs,1,MPI_INTEGER,1,1,EAT_COMM_obs_filter,ierr)
         if (nobs > 0) then
            if (.not. allocated(obs)) allocate(obs(nobs))
            if (allocated(obs) .and. nobs > size(obs)) then
                deallocate(obs)
                allocate(obs(nobs))
            end if
            CALL RANDOM_NUMBER(obs)
            call MPI_ISEND(obs,nobs,MPI_DOUBLE,1,1,EAT_COMM_obs_filter,request,ierr)
            call sleep(1)
            call MPI_WAIT(request,stat,ierr)
         end if
      end if
   end do
   if (have_model) then
      do m=1,nmodel
         call MPI_SEND(halt,19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr) !!!! nmodel
      end do
   else
      write(error_unit,*) 'obs(no more obs) -> halting'
   end if
   if (have_filter) then
      nobs=-1
      write(error_unit,'(A,I6)') ' obs(-> nobs)    ',nobs
      call MPI_SEND(nobs,1,MPI_INTEGER,1,1,EAT_COMM_obs_filter,ierr)
   end if
end subroutine do_observations

!-----------------------------------------------------------------------

subroutine finish_observations()

   !! Finish the eat_observation program by calling MPI_Finalize()

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_observations

!-----------------------------------------------------------------------

end program eat_observations
