! Copyright (C) 2021 Bolding & Bruggeman

program eat_2d_model

   !! An re-implementation of the 'model' from here:
   !! http://pdaf.awi.de/files/pdaf_tutorial_onlineserial.pdf
   !! Used for testing the analysis step

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use eat_2d_data
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
   integer :: state_size=nx*ny
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
   character(len=*), parameter :: nmlfile='eat_2d.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   logical :: all_verbose=.true.
   namelist /nml_2d_model/ verbosity,all_verbose
   character(len=64) :: fn
   integer :: total_steps
   integer :: N=0
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_2d_model)
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
      if (verbosity >= debug) write(stderr,*) 'model(signal) ',signal

      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_2d()
         call post_model_initialize()
!KB         if (iand(signal,signal_send_state) == signal_send_state) then
!KB            allocate(state(state_size))
!KB         end if
#if 1
      else
         if (have_filter) then
            call MPI_IRECV(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model_filter,request,ierr)
            call MPI_WAIT(request,stat,ierr)
         end if
#endif
      end if
!KB      if (verbosity >= info) write(stderr,*) 'model1(signal) ',signal

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_2d()
         call post_model_integrate()
      end if
!KB      if (verbosity >= info) write(stderr,*) 'model2(signal) ',signal

      if (have_filter .and. iand(signal,signal_send_state) == signal_send_state) then
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model_filter,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if
!KB      if (verbosity >= info) write(stderr,*) 'model3(signal) ',signal

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_2d()
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
         signal=signal_initialize+signal_integrate+signal_send_state
         first=.false.
         if (rank_model_comm == 0 .and. have_filter) then
            call MPI_SSEND(state_size,1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
         end if
      else
         signal=signal_integrate+signal_send_state
      end if

      call MPI_RECV(timestr,19,MPI_CHARACTER,0,MPI_ANY_TAG,EAT_COMM_obs_model,stat,ierr)
      if (verbosity >= debug) write(stderr,*) 'model(<- time)  ',trim(timestr)
      if (ierr /= MPI_SUCCESS) then
         call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      if(trim(timestr) == "0000-00-00 00:00:00") then
         signal=signal+signal_finalize
         if (have_filter) then
            signal=signal-signal_send_state
         end if
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
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no filter program present)"
      have_filter=.false.
   end if
end subroutine pre_model_initialize

!-----------------------------------------------------------------------

subroutine initialize_2d()
   if (verbosity >= info) write(stderr,*) 'model(initialize): '
   start="1998-01-01 00:00:00"
   stop="1999-01-01 00:00:00"
   ! reproduce the fields read in the PDAF tutorial
   if (member == 1) then
      call true_field(member)
      write(fn,'(A,I0.4,A)') 'true_',0,'.dat'
      call write_field(fn,field)
   end if

   call true_field(member,nmember=size_model_comm)
   write(fn,'(A,*(I0.2,A))') 'ens_',member,'_step',0,'_ini.dat'
   call write_field(fn,field)
end subroutine initialize_2d

!-----------------------------------------------------------------------

subroutine post_model_initialize()
   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   if (verbosity >= info) then
      write(stderr,'(x,4A)') 'model(sim_start->sim_stop) ',sim_start%isoformat(),' -> ',sim_stop%isoformat()
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
      signal=signal-signal_initialize
   end if
   if (have_filter) then
      write(fn,'(A,*(I0.2,A))') 'ens_',member,'_step',N,'_ini.dat'
   end if
   state=reshape(field,(/state_size/))
end subroutine post_model_initialize

!-----------------------------------------------------------------------

subroutine pre_model_integrate()
   if (.not. ensemble_only) then
      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
      else
         stop_time=strptime(trim(timestr),time_format)
      end if
      if (have_filter) then
#if 0
         call MPI_IRECV(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model_filter,request,ierr)
         call MPI_WAIT(request,stat,ierr)
#endif
!KB         field=reshape(state,(/nx,ny/))
         write(fn,'(A,*(I0.2,A))') 'ens_',member,'_step',N,'_ana.dat'
         call write_field(fn,field)
      end if
   end if
end subroutine pre_model_integrate

!-----------------------------------------------------------------------

subroutine integrate_2d()
   integer :: nsteps=2
   if (verbosity >= info) write(stderr,'(A,I4)') ' model(integrate): ',N
      write(stderr,'(x,4A)') 'model(start_time->stop_time) ',start_time%isoformat(),' -> ',stop_time%isoformat()
   if (.not. ensemble_only) then
      N=N+nsteps
   end if
   call update_field(nsteps=nsteps)
!KB   call sleep(1)
end subroutine integrate_2d

!-----------------------------------------------------------------------

subroutine post_model_integrate()
   start_time=stop_time
   if (member == 1 .and. size_model_comm == 1) then
      write(fn,'(A,I0.4,A)') 'true_',N,'.dat'
   else
      write(fn,'(A,*(I0.2,A))') 'ens_',member,'_step',N,'_for.dat'
   end if
   call write_field(fn,field)
end subroutine post_model_integrate

!-----------------------------------------------------------------------

subroutine finalize_2d()
   if (verbosity >= info) write(stderr,*) 'model(finalize)'
end subroutine finalize_2d

!-----------------------------------------------------------------------

end program eat_2d_model
