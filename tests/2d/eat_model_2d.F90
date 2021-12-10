! Copyright (C) 2021 Bolding & Bruggeman

program eat_model_2d

   !! An re-implementation of the 'model' from here:
   !! http://pdaf.awi.de/files/pdaf_tutorial_onlineserial.pdf
   !! Used for testing the analysis step

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use fields_2d
   use datetime_module, only: datetime, timedelta, strptime
   IMPLICIT NONE

   integer :: ierr
   integer :: member
   integer :: stat(MPI_STATUS_SIZE)
   integer :: request
   integer :: state_size=nx*ny
   character(len=19) :: start="1998-01-01 00:00:00"
   character(len=19) :: stop="1998-01-02 12:00:00"
   character(len=19) :: timestr
   character(len=14) :: timestamp
   real(real64) :: timestep=3600._real64
   real(real64), allocatable :: state(:)
   logical :: ensemble_only=.false.
   integer :: signal

   ! Most of this must go to a model specific file
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
   namelist /nml_eat_model/ verbosity,all_verbose,start,stop
   character(len=64) :: fn
   integer :: total_steps
   integer :: N=0
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

   if (.not. all_verbose .and. member /= 0) then
      verbosity=silent
   end if

   call pre_model_initialize()

   signal=signal_initialize
   do
      if (verbosity >= debug) write(stderr,*) 'model(signal) 1 ',signal
      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_2d()
         call post_model_initialize()
      end if

      if (verbosity >= debug) write(stderr,*) 'model(signal) 2 ',signal
      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_model_integrate()
         call integrate_2d()
         call post_model_integrate()
      end if

      if (verbosity >= debug) write(stderr,*) 'model(signal) 3 ',signal
      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_2d()
         exit
      end if
   end do
   if (verbosity >= info) write(stderr,*) 'model(exit)'
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

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
   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   timestamp = sim_start%strftime("%Y%m%d%H%M%S")
   ! reproduce the fields read in the PDAF tutorial
   if (member == 1 .and. size_model_comm == 1) then
      call true_field(member)
      write(fn,'(3A)') 'true_',timestamp,'.dat'
      call write_field(fn,field)
   end if

   call true_field(member,nmember=size_model_comm)
   write(fn,'(A,I0.2,A,A,A)') 'ens_',member,'_ini_',timestamp,'.dat'
   call write_field(fn,field)
end subroutine initialize_2d

!-----------------------------------------------------------------------

subroutine post_model_initialize()
   if (have_filter .and. rank_model_comm == 0) then
      call MPI_SSEND(state_size,1,MPI_INTEGER,0,10,EAT_COMM_model_filter,ierr)
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
      signal=signal_integrate
      if (have_filter) then
         if (.not. allocated(state)) allocate(state(state_size))
      end if
   end if
!KB   if (have_filter) then
!KB      write(fn,'(A,*(I0.2,A))') 'ens_',member,'_step',N,'_ini.dat'
!KB   end if
end subroutine post_model_initialize

!-----------------------------------------------------------------------

subroutine pre_model_integrate()
   if (.not. ensemble_only) then

      call MPI_RECV(timestr,19,MPI_CHARACTER,0,MPI_ANY_TAG,EAT_COMM_obs_model,stat,ierr)
      if (ierr /= MPI_SUCCESS) then
         call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      if (verbosity >= info) write(stderr,*) 'model(<-- time)  ',trim(timestr)

      if (have_filter .and. iand(signal,signal_recv_state) == signal_recv_state) then
         if (verbosity >= debug) write(stderr,*) 'model: before receiving state'
         call MPI_IRECV(state,state_size,MPI_DOUBLE,0,tag_analysis,EAT_COMM_model_filter,request,ierr)
         if (verbosity >= debug) write(stderr,*) 'model: waiting receiving state'
         call MPI_WAIT(request,stat,ierr)
         if (verbosity >= debug) write(stderr,*) 'model: after receiving state'
         field=reshape(state,(/nx,ny/))
         signal=signal-signal_recv_state
      end if

      if(trim(timestr) == "0000-00-00 00:00:00") then
         stop_time=sim_stop
         signal=signal+signal_finalize
      else
         stop_time=strptime(trim(timestr),time_format)
         if (have_filter) signal=signal+signal_send_state
      end if

      timestamp = start_time%strftime("%Y%m%d%H%M%S")
      if (have_filter) then
         write(fn,'(A,I0.2,A,A,A)') 'ens_',member,'_ana_',timestamp,'.dat'
      else
         if (member == 1 .and. size_model_comm == 1) then
            write(fn,'(3A)') 'true_',timestamp,'.dat'
         end if
      end if
      call write_field(fn,field)
   end if
end subroutine pre_model_integrate

!-----------------------------------------------------------------------

subroutine integrate_2d()
   integer :: nsteps
   type(timedelta) :: td
   td = stop_time-start_time
   nsteps=nint(td%total_seconds()/timestep)
   if (verbosity >= info) write(stderr,'(x,2(A,I4))') 'model(integrate): ',N,' --> ',N+nsteps
   if (verbosity >= debug) write(stderr,'(x,4A)') 'model(start_time->stop_time) ', &
                                        start_time%isoformat(),' --> ',stop_time%isoformat()
   if (.not. ensemble_only) then
      N=N+nsteps
   end if
   call update_field(nsteps=nsteps)
!KB   call sleep(1)
end subroutine integrate_2d

!-----------------------------------------------------------------------

subroutine post_model_integrate()

   if (have_filter .and. iand(signal,signal_send_state) == signal_send_state) then
      state=reshape(field,(/state_size/))
      call MPI_ISEND(state,state_size,MPI_DOUBLE,0,tag_forecast,EAT_COMM_model_filter,request,ierr)
      if(ierr /= MPI_SUCCESS) THEN
         write(stderr,*) 'Fatal error (MODEL): Unable to send: ',member
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      call MPI_WAIT(request,stat,ierr)
      if(ierr /= MPI_SUCCESS) THEN
         write(stderr,*) 'Fatal error (MODEL): Unable to send: ',member
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      signal=signal-signal_send_state
   end if
   if (have_filter) then
      signal=signal+signal_recv_state
   end if
   timestamp = stop_time%strftime("%Y%m%d%H%M%S")
   start_time=stop_time
   write(fn,'(A,I0.2,A,A,A)') 'ens_',member,'_for_',timestamp,'.dat'
   call write_field(fn,field)
end subroutine post_model_integrate

!-----------------------------------------------------------------------

subroutine finalize_2d()
   if (verbosity >= info) write(stderr,*) 'model(finalize)'
end subroutine finalize_2d

!-----------------------------------------------------------------------

end program eat_model_2d
