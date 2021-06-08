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
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info
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
   character(len=*), parameter :: nmlfile='eat_seamless.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   logical :: all_verbose=.true.
   character(len=64) :: obs_times_file="obs_times.dat"
   namelist /nml_seamless_obs/ verbosity,all_verbose,obs_times_file
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_seamless_obs)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'obs(read namelist)'
   end if

   call init_eat_config(color_obs+verbosity)

   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "obs(no model program present)"
      have_model=.false.
   else
      nmodel=size_obs_model_comm-size_obs_comm
   end if

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "obs(no filter program present)"
      have_filter=.false.
   end if

!KB   if (.not. all_verbose .and. member /= 0) then
!KB      verbosity=silent
!KB   end if

   open(newunit=unit,file=trim(obs_times_file),status='old',action='read',iostat=ios)
   if (verbosity >= fatal) then
      if (ios /= 0) stop 'obs(unable to open obs_times.dat for reading)'
   else
      if (ios /= 0) stop 1
   end if
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
      if (verbosity >= info) write(stderr,*) 'obs(-> time)  ',trim(obs_times(n))
      if (have_model) then
         do m=1,nmodel
            call MPI_SEND(obs_times(n),19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr)
         end do
      end if
      if (have_filter) then
         nobs=10000*n
         if (verbosity >= info) write(stderr,'(A,I6)') ' obs(-> nobs)    ',nobs
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
      if (verbosity >= info) write(stderr,*) 'obs(no more obs) -> halting'
   end if
   if (have_filter) then
      nobs=-1
      if (verbosity >= info) write(stderr,'(A,I6)') ' obs(-> nobs)    ',nobs
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
