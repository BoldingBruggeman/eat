! Copyright (C) 2021 Bolding & Bruggeman

program eat_model_gotm

   !! An implementation of GOTM in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use gotm, only: initialize_gotm, integrate_gotm, finalize_gotm, generate_restart_file
   use cmdline
   use time, only: start,stop,timestep,julianday,fsecondsofday
   use time, only: MinN,MaxN
   use datetime_module, only: datetime, timedelta, clock, strptime
   use field_manager
   use output_manager_core
   use output_manager
   use memory_output
   use register_all_variables, only: fm
   use fabm_types, only: fabm_parameter_pointers

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
   character(len=128) :: extra_state(64)
   character(len=128) :: fabm_parameters_in_state(64)
   logical :: fileexists
   integer :: nmlunit,outunit
   logical :: all_verbose=.false.
   logical :: shared_gotm_yaml=.true.
   logical :: shared_restart_file=.true.
   namelist /nml_eat_model/ verbosity,all_verbose,shared_gotm_yaml,shared_restart_file,extra_state, fabm_parameters_in_state

   type (type_field_set) :: field_set
   integer :: i, n
   character(len=1024) :: arg
!-----------------------------------------------------------------------
   extra_state = ''
   fabm_parameters_in_state = ''
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_eat_model)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'model(read namelist)'
   end if

   n = command_argument_count()
   i = 1
   do while (i <= n)
      call get_command_argument(i, arg)
      select case (arg)
      case ('--separate_gotm_yaml')
         shared_gotm_yaml = .false.
      case ('--separate_restart_file')
         shared_restart_file = .false.
      case ('--generate_restart_file')
         generate_restart_file = .true.
      end select
      i = i+1
   end do

   call init_eat_config(color_model+verbosity)

   member=rank_model_comm

   if (.not. all_verbose .and. member /= 0) verbosity=silent

   call pre_model_initialize()

   extra_state = ''
   fabm_parameters_in_state = ''
   if (have_filter) then
      call MPI_RECV(n,1,MPI_INTEGER,0,MPI_ANY_TAG,EAT_COMM_model_filter,stat,ierr)
      if(ierr /= MPI_SUCCESS) call fatal_error('Unable to receive number of diagnostics to be added to state')
      write(stderr,*) 'Adding ', n, ' diagnostics to the model state as seen by filter'
      do i = 1, n
         call MPI_RECV(extra_state(i),len(extra_state(i)),MPI_CHARACTER,0,MPI_ANY_TAG,EAT_COMM_model_filter,stat,ierr)
         if(ierr /= MPI_SUCCESS) call fatal_error('Unable to receive diagnostic to be added to state')
         write(stderr,*) '- ', trim(extra_state(i))
      end do

      call MPI_RECV(n,1,MPI_INTEGER,0,MPI_ANY_TAG,EAT_COMM_model_filter,stat,ierr)
      if(ierr /= MPI_SUCCESS) call fatal_error('Unable to receive number of FABM parameter to be added to state')
      write(stderr,*) 'Adding ', n, ' FABM parameters to the model state as seen by filter'
      do i = 1, n
         call MPI_RECV(fabm_parameters_in_state(i),len(fabm_parameters_in_state(i)),MPI_CHARACTER,0,MPI_ANY_TAG, &
            EAT_COMM_model_filter,stat,ierr)
         if(ierr /= MPI_SUCCESS) call fatal_error('Unable to receive FABM parameter to be added to state')
         write(stderr,*) '- ', trim(fabm_parameters_in_state(i))
      end do
   end if

   fabm_parameter_pointers = any(fabm_parameters_in_state /= '')

   if (ensemble_only) then
      signal=signal_initialize+signal_integrate+signal_finalize
   else
      signal=signal_initialize
   end if
   do
      if (verbosity >= debug) write(stderr,*) 'model(signal) 1 ',signal
      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_gotm()
         call post_model_initialize()
         ! if (have_obs) call MPI_Barrier(EAT_COMM_obs_model,ierr)
      end if

      if (verbosity >= debug) write(stderr,*) 'model(signal) 2 ',signal
      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_gotm()
         call post_model_integrate()
      end if

      if (verbosity >= debug) write(stderr,*) 'model(signal) 3 ',signal
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

subroutine pre_model_initialize()
   ! if (EAT_COMM_obs_model == MPI_COMM_NULL) then
   !    if (verbosity >= warn) write(stderr,*) "model(no observation program present)"
   !    have_obs=.false.
   !    ensemble_only=.true.
   ! end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no filter program present)"
      have_filter=.false.
      have_obs=.false.
      ensemble_only=.true.
      if (nprocs == 1) then
         ! Behave like normal (non-EAT) GOTM running in serial
         if (.not. generate_restart_file) call parse_cmdline('eat-gotm')
         return
      end if
   end if

   output: block
      use gotm, only: yaml_file,restart_file,output_id,force_restart_offline
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
         force_restart_offline = .true.
      end if
   end block output
end subroutine pre_model_initialize

!-----------------------------------------------------------------------

subroutine post_model_initialize()
   use yaml_settings, only: type_key_value_pair, type_real_setting, format_real
   use gotm_fabm, only: fabm_model => model

   type (type_output_item), pointer :: item
   integer :: ios, i, j
   integer, parameter :: unit = 250
   character(len=64) :: fn='da_variables.dat'
   character(len=:), allocatable :: str_scale_factor
   class (type_key_value_pair), pointer :: pair

   if (have_filter .and. rank_model_comm == 0) then
      call MPI_SEND(start,19,MPI_CHARACTER,0,0,EAT_COMM_model_filter,ierr)
      if(ierr /= MPI_SUCCESS) call fatal_error('Unable to send time of simulation start')
      call MPI_SEND(stop,19,MPI_CHARACTER,0,0,EAT_COMM_model_filter,ierr)
      if(ierr /= MPI_SUCCESS) call fatal_error('Unable to send time of simulation stop')
   end if

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
      ! Create output object that will manage memory blob with (extended) model state
      allocate(memory_file)
      call output_manager_add_file(fm, memory_file)

      ! Add model state
      allocate(item)
      item%name = 'state'
      item%output_level = output_level_debug
      call memory_file%append_item(item)

      ! Optionally extend state with user-selected variables
      do i = 1, size(fabm_parameters_in_state)
         if (fabm_parameters_in_state(i) /= '') then
            if (.not. associated(fabm_model)) call fatal_error('Cannot add FABM parameters to state &
                  &because FABM is not used in this GOTM configuration')
            pair => fabm_model%settings%get_node(trim(fabm_parameters_in_state(i)))
            select type (value => pair%value)
            class is (type_real_setting)
               allocate(item)
               str_scale_factor = ''
               if (value%scale_factor /= 1.0) str_scale_factor = '*' // format_real(value%scale_factor)
               item%name = trim(fabm_parameters_in_state(i))
               do j=1,len_trim(item%name)
                  if (item%name(j:j) == '/') item%name(j:j) = '_'
               end do
               call fm%register(item%name, value%units//str_scale_factor, value%long_name, data0d=value%pvalue, field=item%field)
               call memory_file%append_item(item)
            class default
               call fatal_error(trim(fabm_parameters_in_state(i)) // ' is not a real-valued FABM parameter')
            end select
         end if
      end do

      ! Optionally extend state with user-selected variables
      do i = 1, size(extra_state)
         if (extra_state(i) /= '') then
            allocate(item)
            item%name = trim(extra_state(i))
            item%field => memory_file%field_manager%select_for_output(extra_state(i))
            call memory_file%append_item(item)
         end if
      end do

      ! Start output manager - this initializes all output "files" and thereby ensures the memory blob in memory_file is allocated
      call output_manager_start(julianday,int(fsecondsofday),int(mod(fsecondsofday,1._real64)*1000000),0)

      ! Write file describing the structure of the memory blob (rank 0 only)
      if (member == 0) then
         open(unit, file=trim(fn), action='write', status='replace', iostat=ios)
         call memory_file%write_metadata(unit)
         close(unit=unit)
      end if

      if (have_filter .and. rank_model_comm == 0) then
         call MPI_SSEND(size(memory_file%data),1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
         if(ierr /= MPI_SUCCESS) call fatal_error('Unable to send size of model state')
      end if

      signal=signal_integrate
   end if
end subroutine post_model_initialize

!-----------------------------------------------------------------------

subroutine pre_model_integrate()
   TYPE(timedelta) :: td

   if (ensemble_only) return

   if (have_filter .and. iand(signal,signal_recv_state) == signal_recv_state) then
      call MPI_RECV(memory_file%data,size(memory_file%data),MPI_DOUBLE,0,tag_analysis,EAT_COMM_model_filter,stat,ierr)
      if(ierr /= MPI_SUCCESS) call fatal_error('Unable to receive model state')
      call memory_file%restore()
      signal=signal-signal_recv_state
   end if

   call MPI_RECV(timestr,19,MPI_CHARACTER,0,MPI_ANY_TAG,EAT_COMM_model_filter,stat,ierr)
   if(ierr /= MPI_SUCCESS) call fatal_error('Unable to receive time of next forecast')
   if (verbosity >= debug) write(stderr,*) 'model(<-- time)  ',trim(timestr)

   if(trim(timestr) == "0000-00-00 00:00:00") then
      stop_time=sim_stop
      signal=signal+signal_finalize
   else
      stop_time=strptime(trim(timestr),time_format)
      if (have_filter) signal=signal+signal_send_state
   end if
   td = stop_time-sim_start
   MaxN=int(td%total_seconds()/timestep)
   if (MaxN < 0) call fatal_error('New stop time ' // timestr // ' precedes start time &
      &(overall simulation start or end of previous time slice)')
end subroutine pre_model_integrate

!-----------------------------------------------------------------------

subroutine post_model_integrate()
   if (ensemble_only) return

   if (have_filter .and. iand(signal,signal_send_state) == signal_send_state) then
      call memory_file%save(julianday,int(fsecondsofday),int(mod(fsecondsofday,1._real64)*1000000))
      call MPI_SEND(memory_file%data,size(memory_file%data),MPI_DOUBLE,0,tag_forecast,EAT_COMM_model_filter,ierr)
      if(ierr /= MPI_SUCCESS) call fatal_error('Unable to send model state')
      signal=signal-signal_send_state
   end if
   if (have_filter) then
      signal=signal+signal_recv_state
   end if
   MinN=MaxN+1
end subroutine post_model_integrate

subroutine fatal_error(msg)
   character(len=*), intent(in) :: msg

   write(stderr,*) 'Fatal error: ', msg
   call flush(stderr)
   call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
end subroutine fatal_error

!-----------------------------------------------------------------------

end program eat_model_gotm

!-----------------------------------------------------------------------
