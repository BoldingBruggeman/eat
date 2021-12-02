! Copyright (C) 2021 Bolding & Bruggeman
!KB#define _WORKSHOP_

program eat_model_gotm

   !! An implementation of GOTM in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use gotm, only: initialize_gotm, integrate_gotm, finalize_gotm
   use time, only: start,stop,timestep
   use time, only: MinN,MaxN
#ifdef _WORKSHOP_
   use meanflow, only: T
#endif
   use datetime_module, only: datetime, timedelta, clock, strptime
   use field_manager
   use output_manager
   use register_all_variables, only: fm
   IMPLICIT NONE

   integer :: ierr
   integer :: member
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request
   integer :: state_size=1657 !!!!KB
   real(real64), allocatable :: state(:)
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
   logical :: all_verbose=.true.
   namelist /nml_eat_model/ verbosity,all_verbose

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

   do
      call signal_setup()
      if (verbosity >= debug) write(stderr,*) 'model(signal) ',signal

      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_gotm()
         call post_model_initialize()
         if (iand(signal,signal_send_state) == signal_send_state) then
            call state_vector_size(state_size)
            allocate(state(state_size))
         end if
      else
         if (have_filter .and. iand(signal,signal_recv_state) == signal_recv_state) then
            call MPI_IRECV(state,state_size,MPI_DOUBLE,0,tag_analysis,EAT_COMM_model_filter,request,ierr)
            call MPI_WAIT(request,stat,ierr)
            call state_to_fields()
         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_gotm()
         call post_model_integrate()
      end if

      if (have_filter .and. iand(signal,signal_send_state) == signal_send_state) then
         call fields_to_state()
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,tag_forecast,EAT_COMM_model_filter,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_gotm()
         exit
      end if
   end do
   if (verbosity >= info) write(stderr,*) 'model(exit)'
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine da_layout(fn)
   use memory_output
   use output_manager_core

   IMPLICIT NONE
   character(len=*), intent(in) :: fn

   class (type_memory_file), pointer :: file
   type (type_output_item),  pointer :: item
   integer :: ios
   integer, parameter :: unit = 250
   integer :: jul,secs
   logical :: success
   real(real64) :: fsecs=0.

   allocate(file)
   call output_manager_add_file(fm,file)

   allocate(item)
   item%name = 'state'
   item%output_level = output_level_debug
   call file%append_item(item)

   call read_time_string(start,jul,secs,success)
   call output_manager_start(jul,secs,int(mod(1._real64*secs,1._real64)*1000000),0)

   open(unit, file=trim(fn), action='write', status='replace', iostat=ios)
   call file%write_metadata(unit)
   close(unit=unit)
end subroutine da_layout

!-----------------------------------------------------------------------

subroutine signal_setup()
   logical :: first=.true.
   character(len=64) :: fn='da_variables.dat'

   if (ensemble_only) then
      signal=signal_initialize+signal_integrate+signal_finalize
   else
      if (first) then
         first=.false.
         signal=signal_initialize+signal_integrate
         MinN=1
         if (have_filter .and. rank_model_comm == 0) then
            call MPI_SSEND(state_size,1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
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
!KB      member=stat(MPI_TAG)
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
#ifdef _WORKSHOP_
T = 0.01*member+T
#endif
!KB   call da_layout(fn)
#if 0
write(0,*) 'member= ',member
   if (have_obs .and. member == 0) then
write(0,*) 'model: ready to write layout'
      call da_layout(fn)
write(0,*) 'model: layout is written'
write(0,*) 'model: ready to send layout'
      call MPI_SEND(fn,64,MPI_CHARACTER,0,member,EAT_COMM_obs_model,ierr)
write(0,*) 'model: layout is sent'
   end if
#endif
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
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine state_vector_size(state_size)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: state_size
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   class (type_field_set_member),pointer    :: member
!EOP
!-------------------------------------------------------------------------
!BOC
   field_set = fm%get_state()
   member => field_set%first
   state_size = 0
   do while (associated(member))
      if (verbosity >= info)  write(stderr,*) '   state vector:  ',trim(member%field%name)
      if (associated(member%field%data%p0d)) then
         state_size = state_size+1
      elseif (associated(member%field%data%p1d)) then
         ! Depth-dependent variable with data pointed to by member%field%data%p1d
         state_size = state_size+size(member%field%data%p1d)
      else
         stop 'no data assigned to state field'
      end if
      member => member%next
   end do
   if (verbosity >= info)  write(stderr,*) '   state vector size:  ',state_size
   end subroutine state_vector_size
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save the model results in file
!
! !INTERFACE:
   subroutine state_to_fields()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                  :: n,len
   class (type_field_set_member),pointer    :: member
!EOP
!-------------------------------------------------------------------------
!BOC
   member => field_set%first
   n = 1
   len = 0
   do while (associated(member))
      ! This field is part of the model state. Its name is member%field%name.
      if (associated(member%field%data%p0d)) then
         ! Depth-independent variable with data pointed to by child%field%data%p0d
         member%field%data%p0d = state(n)
         n = n + 1
      elseif (associated(member%field%data%p1d)) then
         ! Depth-dependent variable with data pointed to by member%field%data%p1d
         len = size(member%field%data%p1d)
         member%field%data%p1d = state(n:n+len-1)
         n = n + len
      else
         stop 'no data assigned to state field'
      end if
      member => member%next
   end do
   end subroutine state_to_fields
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Data from fields to state vector
!
! !INTERFACE:
   subroutine fields_to_state()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                  :: n,len
   class (type_field_set_member),pointer    :: member
!EOP
!-------------------------------------------------------------------------
!BOC
   member => field_set%first
   n = 1
   len = 0
   do while (associated(member))
      ! This field is part of the model state. Its name is member%field%name.
      if (associated(member%field%data%p0d)) then
         ! Depth-independent variable with data pointed to by member%field%data%data%p0d
         state(n) = member%field%data%p0d
         n = n+1
      elseif (associated(member%field%data%p1d)) then
         ! Depth-dependent variable with data pointed to by member%field%data%p1d
         len = size(member%field%data%p1d)
         state(n:n+len-1) = member%field%data%p1d
         n = n+len
      else
         stop 'no data assigned to state field'
      end if
      member => member%next
   end do
   end subroutine fields_to_state
!EOC

!-----------------------------------------------------------------------

end program eat_model_gotm

!-----------------------------------------------------------------------
