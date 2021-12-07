! Copyright (C) 2021 Bolding & Bruggeman

program eat_model_gotm

   !! An implementation of GOTM in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use gotm, only: initialize_gotm, integrate_gotm, finalize_gotm
   use time, only: start,stop,timestep,julianday,fsecondsofday
   use time, only: MinN,MaxN
   use datetime_module, only: datetime, timedelta, clock, strptime
   use field_manager
   use output_manager_core
   use output_manager
   use memory_output
   use register_all_variables, only: fm
   IMPLICIT NONE

   integer :: ierr
   integer :: member
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request

   class (type_memory_file), pointer :: memory_file

   character(len=19) :: timestr
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
   logical :: shared_gotm_yaml=.true.
   logical :: shared_restart_file=.true.
   namelist /nml_eat_model/ verbosity,all_verbose,shared_gotm_yaml,shared_restart_file

   type (type_field_set) :: field_set
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_eat_model)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'model(read namelist)'
   end if

   call init_eat_config(color_model+verbosity)

   member=rank_model_comm

   if (.not. all_verbose .and. member /= 0) verbosity=silent

   call pre_model_initialize()

   call initialize_gotm()
   call post_model_initialize()
   if (have_obs) call MPI_Barrier(EAT_COMM_obs_model,ierr)

   do
      call signal_setup()
!signal = signal-signal_recv_state
      if (verbosity >= debug) write(stderr,*) 'model(signal) ',signal

      if (iand(signal,signal_initialize) == signal_initialize) then
      else
         if (have_filter .and. iand(signal,signal_recv_state) == signal_recv_state) then
            call MPI_IRECV(memory_file%data,size(memory_file%data),MPI_DOUBLE,0,tag_analysis,EAT_COMM_model_filter,request,ierr)
            if(ierr /= MPI_SUCCESS) THEN
               write(stderr,*) 'Fatal error (MODEL): Unable to receive: ',member; call flush(stderr)
               call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
            end if
            call MPI_WAIT(request,stat,ierr)
            if(ierr /= MPI_SUCCESS) THEN
               write(stderr,*) 'Fatal error (MODEL): Unable to wait: ',member; call flush(stderr)
               call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
            end if
            call memory_file%restore()
         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_gotm()
         call post_model_integrate()
      end if

      if (have_filter .and. iand(signal,signal_send_state) == signal_send_state) then
         call memory_file%save(julianday,int(fsecondsofday),int(mod(fsecondsofday,1._real64)*1000000))
         call MPI_ISEND(memory_file%data,size(memory_file%data),MPI_DOUBLE,0,tag_forecast,EAT_COMM_model_filter,request,ierr)
         if(ierr /= MPI_SUCCESS) THEN
            write(stderr,*) 'Fatal error (MODEL): Unable to send: ',member; call flush(stderr)
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
         call MPI_WAIT(request,stat,ierr)
         if(ierr /= MPI_SUCCESS) THEN
            write(stderr,*) 'Fatal error (MODEL): Unable to send: ',member; call flush(stderr)
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_gotm()
         exit
      end if
      call flush(stderr)
   end do
   if (verbosity >= info) write(stderr,*) 'model(exit)'
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine signal_setup()
   logical, save :: first=.true.
   character(len=64) :: fn='da_variables.dat'

   if (ensemble_only) then
      signal=signal_initialize+signal_integrate+signal_finalize
   else
      if (first) then
         first=.false.
         signal=signal_initialize+signal_integrate
         MinN=1
         if (have_filter .and. rank_model_comm == 0) then
            call MPI_SSEND(size(memory_file%data),1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
         end if
         if (have_filter) signal=signal+signal_send_state
      else
         signal=signal_integrate
         if (have_filter) signal=signal+signal_send_state+signal_recv_state
         MinN=MaxN+1
      end if
      call MPI_RECV(timestr,19,MPI_CHARACTER,0,MPI_ANY_TAG,EAT_COMM_obs_model,stat,ierr)
      if (verbosity >= debug) write(stderr,*) 'model(<-- time)  ',trim(timestr)
      if (ierr /= MPI_SUCCESS) then
         call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      if(trim(timestr) == "0000-00-00 00:00:00") then
         signal=signal+signal_finalize
         if (have_filter) signal=signal-signal_send_state
      end if
   end if
end subroutine signal_setup

!-----------------------------------------------------------------------

subroutine pre_model_initialize()
   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no observation program present)"
      have_obs=.false.
      ensemble_only=.true.
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no filter program present)"
      have_filter=.false.
   end if

   output: block
      use gotm, only: yaml_file,restart_file,output_id
      character(len=128) :: fname,strbuf
      write(strbuf, "(A,I0.4)") 'gotm_', member+1
      fname = TRIM(strbuf) // '.stderr'
      open(stderr,file=fname)
      fname = TRIM(strbuf) // '.stdout'
      open(output_unit,file=fname)
      if ( .not. shared_gotm_yaml) then
         yaml_file = TRIM(strbuf) // '.yaml'
      end if
      write(output_id, "(A,I0.4)") '_', member+1
      if ( .not. shared_restart_file) then
         write(restart_file, "(A,I0.4)") 'restart_', member+1
!KB         restart_file = TRIM(strbuf) // '.nc'
      end if
   end block output
end subroutine pre_model_initialize

!-----------------------------------------------------------------------

subroutine post_model_initialize()
   type (type_output_item),  pointer :: item
   integer :: ios
   integer, parameter :: unit = 250
   character(len=64) :: fn='da_variables.dat'

   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   if (verbosity >= debug) then
      write(stderr,*) 'model(sim_start) ',sim_start%isoformat()
      write(stderr,*) 'model(sim_stop)  ',sim_stop%isoformat()
      write(stderr,*) 'model(timestep)  ',timestep
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
   end if
   if (have_obs) then
      allocate(memory_file)
      call output_manager_add_file(fm, memory_file)
      allocate(item)
      item%name = 'state'
      item%output_level = output_level_debug
      call memory_file%append_item(item)
      call output_manager_start(julianday,int(fsecondsofday),int(mod(fsecondsofday,1._real64)*1000000),0)
      if (member == 0) then
         open(unit, file=trim(fn), action='write', status='replace', iostat=ios)
         call memory_file%write_metadata(unit)
         close(unit=unit)
      end if
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

!-----------------------------------------------------------------------
