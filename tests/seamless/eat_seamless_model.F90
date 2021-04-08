! Copyright (C) 2021 Bolding & Bruggeman

program eat_seamless_model

   !! An implementation of a test in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use datetime_module, only: datetime, timedelta, strptime
   IMPLICIT NONE

   integer, parameter :: signal_initialize=1
   integer, parameter :: signal_integrate=2
   integer, parameter :: signal_finalize=4
   integer, parameter :: signal_send_state=8

   integer :: ierr
   integer :: member
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request
   integer :: state_size=1234
   character(len=19) :: start,stop,timestr
   real(real64) :: timestep
   real(real64), allocatable :: state(:)
   logical :: ensemble_only=.false.
   integer :: signal

   ! Most of this must go to a model specific file
   character(len=256), parameter :: time_format='%Y-%m-%d %H:%M:%S'
   TYPE(datetime), save :: sim_start, sim_stop
   TYPE(datetime) :: start_time,stop_time
   logical :: have_obs=.true.
   logical :: have_filter=.true.
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info
   character(len=*), parameter :: nmlfile='eat_seamless.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   logical :: all_verbose=.true.
   namelist /nml_seamless_model/ verbosity,all_verbose
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_seamless_model)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'model(read namelist)'
   end if

   call init_eat_config(color_model+verbosity)

   member=rank_model_comm

   if (.not. all_verbose .and. member /= 0) then
      verbosity=silent
   end if

   call pre_model_initialize()

   do
      call signal_setup()

      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_seamless()
         call post_model_initialize()
!KB         if (iand(signal,signal_send_state) == signal_send_state) then
!KB            allocate(state(state_size))
!KB         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_seamless()
         call post_model_integrate()
      end if

!KB      if (iand(signal,signal_send_state) == signal_send_state) then
!KB         CALL RANDOM_NUMBER(state)
!KB         state=state+rank-1
!KBwrite(stderr,*) 'CCC3 ',member,state_size
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model_filter,request,ierr)
         call MPI_WAIT(request,stat,ierr)
!KBwrite(stderr,*) 'CCC4 ',member
!KB      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_seamless()
         exit
      end if
   end do
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine signal_setup() ! setup signal
   logical :: first=.true.

   if (ensemble_only) then
      signal=signal_initialize+signal_integrate+signal_finalize
   else
      if (first) then
         signal=signal+signal_integrate
         first=.false.
!KBwrite(stderr,*) 'CCC1 ',rank_model_comm,rank_model_filter_comm,state_size
! go to post_model_initialize
         if (rank_model_comm == 0) then
            call MPI_SEND(state_size,1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
         end if
!KBwrite(stderr,*) 'CCC2 ',rank_model_comm
      else
         signal=signal_integrate
      end if

      call MPI_RECV(timestr,19,MPI_CHARACTER,0,MPI_ANY_TAG,EAT_COMM_obs_model,stat,ierr)
      if (ierr /= MPI_SUCCESS) then
         call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      if(trim(timestr) == "0000-00-00 00:00:00") then
         signal=signal+signal_finalize
      end if
      member=stat(MPI_TAG)
   end if
end subroutine signal_setup

!-----------------------------------------------------------------------

subroutine pre_model_initialize()
   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no observation program present)"
      have_obs=.false.
      ensemble_only=.true.
   else
      signal=signal_initialize+signal_send_state
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no filter program present)"
      have_filter=.false.
   end if
end subroutine pre_model_initialize

!-----------------------------------------------------------------------

subroutine post_model_initialize()
   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   start_time=sim_start
   stop_time=sim_stop
   if (verbosity >= debug) then
      write(stderr,*) 'model(sim_start) ',sim_start%isoformat()
      write(stderr,*) 'model(sim_stop)  ',sim_stop%isoformat()
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
      signal=signal_initialize+signal_integrate
   end if
end subroutine post_model_initialize

!-----------------------------------------------------------------------

subroutine pre_model_integrate()
   logical :: first=.true.
   if (.not. ensemble_only) then
      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
      else
         stop_time=strptime(trim(timestr),time_format)
      end if
   end if
end subroutine pre_model_integrate

!-----------------------------------------------------------------------

subroutine post_model_integrate()
   start_time=stop_time
end subroutine post_model_integrate

!-----------------------------------------------------------------------

subroutine initialize_seamless()
   start="1998-01-01 00:00:00"
   stop="1999-01-01 00:00:00"
!KB   state_size=1234
   allocate(state(state_size))
   CALL RANDOM_NUMBER(state)
   if (verbosity >= info) write(stderr,*) 'model(initialize): '
end subroutine initialize_seamless

subroutine integrate_seamless()
   integer :: n=0
   if (verbosity >= info) write(stderr,'(A,A,A,A)') ' model(integrate): ',start_time%isoformat(),' -> ',stop_time%isoformat()
   CALL RANDOM_NUMBER(state)
   state=state+n*(rank_model_comm+1)
   n=n+1
   call sleep(1)
end subroutine integrate_seamless

subroutine finalize_seamless()
   if (verbosity >= info) write(stderr,*) 'model(finalize)'
end subroutine finalize_seamless

!-----------------------------------------------------------------------

end program eat_seamless_model
