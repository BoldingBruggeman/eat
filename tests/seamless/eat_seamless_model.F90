! Copyright (C) 2021 Bolding & Bruggeman

module seamless_model_m

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use eat_model, only: type_eat_model
   use eat_model, only: signal_initialize,signal_integrate,signal_finalize,signal_send_state

   IMPLICIT NONE

   type, extends(type_eat_model) :: type_seamless_model
      character(len=19) :: start,stop
   contains
      procedure :: pre_initialize => pre_initialize_seamless
      procedure :: initialize => initialize_seamless
      procedure :: post_initialize => post_initialize_seamless
      procedure :: pre_integrate => pre_integrate_seamless
      procedure :: integrate => integrate_seamless
      procedure :: post_integrate => post_integrate_seamless
      procedure :: pre_finalize => pre_finalize_seamless
      procedure :: finalize => finalize_seamless
      procedure :: post_finalize => post_finalize_seamless
   end type type_seamless_model

contains

module subroutine pre_initialize_seamless(self)
   class(type_seamless_model), intent(inout) :: self
   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (self%verbosity >= warn) write(self%stderr,*) "model(no observation program present)"
      self%have_obs=.false.
!KB      ensemble_only=.true.
   else
!KB      signal=signal_initialize+signal_send_state
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (self%verbosity >= warn) write(self%stderr,*) "model(no filter program present)"
      self%have_filter=.false.
   end if
end subroutine pre_initialize_seamless

!-----------------------------------------------------------------------
module subroutine initialize_seamless(self)
   class(type_seamless_model), intent(inout) :: self
   if (self%verbosity >= info) write(self%stderr,*) 'model(initialize): '
   self%start="1998-01-01 00:00:00"
   self%stop="1999-01-01 00:00:00"
#if 0
   allocate(state(state_size))
   CALL RANDOM_NUMBER(state)
#endif
end subroutine initialize_seamless

!-----------------------------------------------------------------------
module subroutine post_initialize_seamless(self)
   class(type_seamless_model), intent(inout) :: self
#if 0
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
#endif
end subroutine post_initialize_seamless

!-----------------------------------------------------------------------
module subroutine pre_integrate_seamless(self)
   class(type_seamless_model), intent(inout) :: self
#if 0
   if (.not. ensemble_only) then
      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
      else
         stop_time=strptime(trim(timestr),time_format)
      end if
   end if
#endif
end subroutine pre_integrate_seamless

!-----------------------------------------------------------------------
module subroutine integrate_seamless(self)
   class(type_seamless_model), intent(inout) :: self
   integer :: n=0
#if 0
   if (verbosity >= info) write(self%stderr,'(A,A,A,A)') ' model(integrate): ',start_time%isoformat(),' -> ',stop_time%isoformat()
   CALL RANDOM_NUMBER(state)
   state=state+n*(rank_model_comm+1)
   n=n+1
   call sleep(1)
#endif
end subroutine integrate_seamless

!-----------------------------------------------------------------------
module subroutine post_integrate_seamless(self)
   class(type_seamless_model), intent(inout) :: self
#if 0
   if (.not. ensemble_only) then
      start_time=stop_time
   end if
#endif
end subroutine post_integrate_seamless

!-----------------------------------------------------------------------
module subroutine pre_finalize_seamless(self)
   class(type_seamless_model), intent(inout) :: self
#if 0
   if (.not. ensemble_only) then
      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
      else
         stop_time=strptime(trim(timestr),time_format)
      end if
   end if
#endif
end subroutine pre_finalize_seamless

!-----------------------------------------------------------------------
module subroutine finalize_seamless(self)
   class(type_seamless_model), intent(inout) :: self
   if (self%verbosity >= info) write(self%stderr,*) 'model(finalize)'
end subroutine finalize_seamless

!-----------------------------------------------------------------------
module subroutine post_finalize_seamless(self)
   class(type_seamless_model), intent(inout) :: self
#if 0
   if (.not. ensemble_only) then
      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
      else
         stop_time=strptime(trim(timestr),time_format)
      end if
   end if
#endif
end subroutine post_finalize_seamless

end module seamless_model_m

!-----------------------------------------------------------------------

program eat_seamless_model

   !! An implementation of a test in an ensemble context
   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use eat_model, only: signal_initialize,signal_integrate,signal_finalize,signal_send_state
   use seamless_model_m, only: type_seamless_model
   use datetime_module, only: datetime, timedelta, strptime
   IMPLICIT NONE

   integer :: ierr
   integer :: member
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request
   integer :: state_size=1234
!KB   character(len=19) :: start,stop
   character(len=19) :: timestr
   real(real64) :: timestep
   real(real64), allocatable :: state(:)
   logical :: ensemble_only=.true. !KB
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
   class(type_seamless_model), allocatable :: model
!-----------------------------------------------------------------------
   allocate (model)

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

   call model%pre_initialize()

   do
!KB      call signal_setup()
      if (verbosity >= debug) write(stderr,*) 'model(signal) ',signal

      if (iand(signal,signal_initialize) == signal_initialize) then
         call model%initialize()
         call model%post_initialize()
!KB         if (iand(signal,signal_send_state) == signal_send_state) then
            allocate(state(state_size))
!KB         end if
      else
         if (model%have_filter) then
            call MPI_IRECV(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model_filter,request,ierr)
            call MPI_WAIT(request,stat,ierr)
         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call model%pre_integrate()
         call model%integrate()
         call model%post_integrate()
      end if

      if (model%have_filter .and. iand(signal,signal_send_state) == signal_send_state) then
         CALL RANDOM_NUMBER(state)
         state=state+rank-1
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model_filter,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call model%finalize()
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
         if (rank_model_comm == 0) then
            call MPI_SEND(state_size,1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
         end if
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

end program eat_seamless_model
