! Copyright (C) 2021 Bolding & Bruggeman

program eat_obs_gotm

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
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info
   integer :: filter=1
!-----------------------------------------------------------------------

   call init_obs()
   call do_obs()
   call finish_obs()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_obs()

   !! Observation times are obtained from an external file - prepare to
   !! send to the server

   ! Local variables
   integer :: unit,ios
   character(len=32) :: buf
   integer :: i,n
   character(len=*), parameter :: nmlfile='eat_gotm.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   character(len=128)  :: obs_times_file='obs_times.dat'
   namelist /nml_eat_obs/ verbosity,obs_times_file
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_eat_obs)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'obs(read namelist)'
   end if

   call init_eat_config(color_obs+verbosity)

   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= info) write(stderr,*) "obs(no model executable present)"
      have_model=.false.
   else
      nmodel=size_obs_model_comm-size_obs_comm
   end if

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity <= info) write(stderr,*) "obs(no filter executable present)"
      have_filter=.false.
   end if

   open(newunit=unit,file=trim(obs_times_file),status='old',action='read',iostat=ios)
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
end subroutine init_obs

!-----------------------------------------------------------------------

subroutine do_obs()

   !! Handle observations and send size and observations to the server

   ! Local variables
   integer :: m ! ensemble counter
   integer :: n,nobs
   integer, allocatable :: iobs(:)
   real(real64), allocatable :: obs(:)
   integer :: requests(2)
   integer :: stats(MPI_STATUS_SIZE,2)
   character(len=32) :: timestr,halt="0000-00-00 00:00:00"
!-----------------------------------------------------------------------
   do n=1,size(obs_times)
      if (have_model) then
         if (verbosity >= info) write(stderr,*) 'obs(-> model)  ',trim(obs_times(n))
         do m=1,nmodel
            call MPI_SSEND(obs_times(n),19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr)
            if (ierr /= MPI_SUCCESS) then
               if (verbosity >= error) then
                  write(stderr,*) 'obs: failing to send to model process#: ',m
               end if
            end if
         end do
      end if

      if (have_filter) then
         nobs=10*n
         if (verbosity >= info) write(stderr,'(A,I6)') ' obs(-> filter) ',nobs
         call MPI_SSEND(nobs,1,MPI_INTEGER,filter,tag_nobs,EAT_COMM_obs_filter,ierr)
         if (ierr /= MPI_SUCCESS) then
            if (verbosity >= error) then
               write(stderr,*) 'obs: failing to send to filter process'
            end if
         end if
      end if
      if (nobs > 0) then
         if (.not. allocated(iobs)) allocate(iobs(nobs))
         if (allocated(iobs) .and. nobs > size(iobs)) then
             deallocate(iobs)
             allocate(iobs(nobs))
         end if
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (allocated(obs) .and. nobs > size(obs)) then
             deallocate(obs)
             allocate(obs(nobs))
         end if
         CALL RANDOM_NUMBER(obs)
         call MPI_ISEND(iobs(1:nobs),nobs,MPI_INTEGER,1,tag_iobs,EAT_COMM_obs_filter,requests(1),ierr)
         call MPI_ISEND(obs(1:nobs),nobs,MPI_DOUBLE,1,tag_obs,EAT_COMM_obs_filter,requests(2),ierr)
         call MPI_WAITALL(2,requests,stats,ierr)
         if (ierr /= MPI_SUCCESS) then
            if (verbosity >= error) write(stderr,*) 'obs: failing to wait for requests'
         end if
      end if
   end do

   ! Here we must NOT use MPI_SSEND()
   if (have_filter) then
      nobs=-1
      call MPI_SEND(nobs,1,MPI_INTEGER,filter,tag_nobs,EAT_COMM_obs_filter,ierr)
      if (ierr /= MPI_SUCCESS .and. verbosity >= error) write(stderr,*) 'obs: failing to send nobs=-1'
      if (verbosity >= info) write(stderr,'(A,I6)') ' obs(--> filter(exit))'
   end if
   if (have_model) then
      do m=1,nmodel
         call MPI_SEND(halt,19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr)
         if (ierr /= MPI_SUCCESS .and. verbosity >= error) write(stderr,*) 'obs: failing to send - halt - to member#',m
      end do
      if (verbosity >= info) write(stderr,*) 'obs(--> model(exit))'
   end if
   if (verbosity >= info) write(stderr,*) 'obs(exit)'
end subroutine do_obs

!-----------------------------------------------------------------------

subroutine finish_obs()

   !! Finish the eat_observation program by calling MPI_Finalize()

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_obs

!-----------------------------------------------------------------------

end program eat_obs_gotm
