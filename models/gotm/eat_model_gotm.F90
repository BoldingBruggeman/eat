! Copyright (C) 2021 Bolding & Bruggeman

program eat_model_gotm

   !! An implementation of GOTM in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
!KB For version 7
!KB   use gotm, only: initialize_gotm, integrate_gotm, finalize_gotm
   use gotm, only: initialize_gotm => init_gotm
   use gotm, only: integrate_gotm => time_loop
   use gotm, only: finalize_gotm => clean_up
   use time, only: start,stop,timestep
   use time, only: MinN,MaxN
   use datetime_module, only: datetime, timedelta, clock, strptime
   IMPLICIT NONE

   integer, parameter :: signal_initialize=1
   integer, parameter :: signal_integrate=2
   integer, parameter :: signal_finalize=4
   integer, parameter :: signal_send_state=8

   integer :: ierr
   integer :: member
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request
!KB   integer :: state_size
   integer :: state_size=1234 !!!!KB
   character(len=19) :: timestr
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
   character(len=*), parameter :: nmlfile='eat_gotm.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   logical :: all_verbose=.false.
   namelist /nml_model_gotm/ verbosity,all_verbose
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_model_gotm)
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
!KB      call model%signal_setup()
      call signal_setup()
      if (verbosity >= debug) write(stderr,*) 'model(signal) ',signal

      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_gotm()
         call post_model_initialize()
         if (iand(signal,signal_send_state) == signal_send_state) then
            allocate(state(state_size))
         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_gotm()
         call post_model_integrate()
      end if

      if (iand(signal,signal_send_state) == signal_send_state) then
         CALL RANDOM_NUMBER(state)
         state=state+rank-1
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_gotm()
         exit
      end if
   end do
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine signal_setup()
   logical :: first=.true.

   if (ensemble_only) then
      signal=signal_initialize+signal_integrate+signal_finalize
   else
      if (first) then
         signal=signal+signal_integrate
         MinN=1
         first=.false.
      else
         signal=signal_integrate
         MinN=MaxN+1
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
      signal=signal_initialize
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no filter program present)"
      have_filter=.false.
   end if

   output: block
      use gotm, only: yaml_file,output_id
      character(len=128) :: fname,strbuf
      write(output_id, "(A,I0.4)") '_', member+1
      write(strbuf, "(A,I0.4)") 'gotm_', member+1
      yaml_file = TRIM(strbuf) // '.yaml'
      fname = TRIM(strbuf) // '.stderr'
      open(stderr,file=fname)
      fname = TRIM(strbuf) // '.stdout'
      open(output_unit,file=fname)
   end block output
end subroutine pre_model_initialize

!-----------------------------------------------------------------------

subroutine post_model_initialize()
   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   if (verbosity >= debug) then
      write(stderr,*) 'model(sim_start) ',sim_start%isoformat()
      write(stderr,*) 'model(sim_stop)  ',sim_stop%isoformat()
      write(stderr,*) 'model(timestep)  ',timestep
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
      signal=signal_initialize+signal_integrate
   end if
end subroutine post_model_initialize

!-----------------------------------------------------------------------

subroutine pre_model_integrate()
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
end subroutine pre_model_integrate

!-----------------------------------------------------------------------

subroutine post_model_integrate()
   if (.not. ensemble_only) then
      MinN=MaxN+1
   end if
end subroutine post_model_integrate

!-----------------------------------------------------------------------

end program eat_model_gotm
