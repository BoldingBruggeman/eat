! Copyright (C) 2021 Bolding & Bruggeman

program eat_model_seamless

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
   logical :: flag
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request
   integer :: state_size
   character(len=19) :: start,stop,timestr
   real(real64) :: timestep
   integer :: MinN,MaxN
   real(real64), allocatable :: state(:)
   logical :: ensemble_only=.false.
   integer :: signal

   ! Most of this must go to a model specific file
   character(len=256), parameter :: time_format='%Y-%m-%d %H:%M:%S'
   TYPE(datetime), save :: sim_start, sim_stop
   TYPE(datetime) :: start_time,stop_time
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=warn
!-----------------------------------------------------------------------
   call init_eat_config(color_model+verbosity)

   call MPI_COMM_RANK(EAT_COMM_model,member,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
   end if

   call pre_eat_model_initialize()

   do
      call signal_setup()

      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_seamless()
         call post_eat_model_initialize()
         state_size=1234 !!!!KB
         if (iand(signal,signal_send_state) == signal_send_state) then
            allocate(state(state_size))
         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_eat_model_integrate()
         call integrate_seamless()
         call post_eat_model_integrate()
      end if

      if (iand(signal,signal_send_state) == signal_send_state) then
         CALL RANDOM_NUMBER(state)
         state=state+rank-1
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_seamless()
         exit
      end if
   end do
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine pre_eat_model_initialize()
   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= info) write(stderr,*) "No observation program present"
      ensemble_only=.true.
   else
      signal=signal_initialize
   end if

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity >= info) then
         write(stderr,*) "No filter program present"
      end if
   end if
end subroutine pre_eat_model_initialize

!-----------------------------------------------------------------------

subroutine post_eat_model_initialize()
   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   if (verbosity >= info) then
      write(stderr,*) 'model(sim_start) ',sim_start%isoformat()
      write(stderr,*) 'model(sim_stop)  ',sim_stop%isoformat()
      write(stderr,*) 'model(timestep)  ',timestep
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
      signal=signal_initialize+signal_integrate
   end if
end subroutine post_eat_model_initialize

!-----------------------------------------------------------------------

subroutine signal_setup() ! setup signal
   logical :: first=.true.

   if (ensemble_only) then
      signal=signal_initialize+signal_integrate+signal_finalize
   else
      if (first) then
         signal=signal+signal_integrate
         first=.false.
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

subroutine pre_eat_model_integrate()
   logical :: first=.true.
   TYPE(timedelta) :: td
   if (.not. ensemble_only) then
      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
      else
         stop_time=strptime(trim(timestr),time_format)
      end if
      td = stop_time-sim_start
      MaxN=int(td%total_seconds()/timestep)
   end if
end subroutine pre_eat_model_integrate

!-----------------------------------------------------------------------

subroutine post_eat_model_integrate()
   if (.not. ensemble_only) then
      MinN=MaxN+1
   end if
end subroutine post_eat_model_integrate

!-----------------------------------------------------------------------

subroutine initialize_seamless()
   start="1998-01-01 00:00:00"
   stop="1999-01-01 00:00:00"
   timestep=3600
   if (verbosity >= info) write(stderr,*) 'model(initialize): '
   MinN=1
   MaxN=8760
end subroutine initialize_seamless
subroutine integrate_seamless()
   if (verbosity >= info) write(stderr,'(A,I5,A,I5)') ' model(integrate): ',MinN,' -> ',MaxN
   call sleep(1)
end subroutine integrate_seamless
subroutine finalize_seamless()
   if (verbosity >= info) write(stderr,*) 'model(finalize)'
end subroutine finalize_seamless

!-----------------------------------------------------------------------

end program eat_model_seamless
