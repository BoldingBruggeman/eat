! Copyright (C) 2021 Bolding & Bruggeman

program eat_model_2d

   !! An re-implementation of the 'model' from here:
   !! http://pdaf.awi.de/files/pdaf_tutorial_onlineserial.pdf
   !! Used for testing the analysis step

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
   real(real64), allocatable :: state(:)
   logical :: ensemble_only=.false.
   integer :: signal

   ! Most of this must go to a model specific file
   character(len=256), parameter :: time_format='%Y-%m-%d %H:%M:%S'
   TYPE(datetime), save :: sim_start, sim_stop
   TYPE(datetime) :: start_time,stop_time
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info

   integer, parameter ::nx=36,ny=18
   integer :: total_steps=5
   integer :: N=0
   real(real64) :: field(ny,nx)
   character(len=*), parameter :: nmlfile='eat_2d.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   integer :: i
   logical :: all_verbose=.true.
   namelist /nml_2d_model/ verbosity,all_verbose
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_2d_model)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'model(read namelist)'
   end if

   call init_eat_config(color_model+verbosity)

   call MPI_COMM_RANK(EAT_COMM_model,member,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
   end if

   if (.not. all_verbose .and. member /= 0) then
      verbosity=silent
   end if

   call pre_eat_model_initialize()

   do
      call signal_setup()
      if (verbosity >= debug) write(stderr,*) 'model(signal) ',signal

      if (iand(signal,signal_initialize) == signal_initialize) then
         call initialize_2d()
         call post_eat_model_initialize()
         state_size=nx*ny
         if (iand(signal,signal_send_state) == signal_send_state) then
            allocate(state(state_size))
         end if
      end if

      if (iand(signal,signal_integrate) == signal_integrate) then
         call pre_eat_model_integrate()
         call integrate_2d()
         call post_eat_model_integrate()
      end if

      if (iand(signal,signal_send_state) == signal_send_state) then
         CALL RANDOM_NUMBER(state)
         state=state+rank-1
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,EAT_COMM_model,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if

      if (iand(signal,signal_finalize) == signal_finalize) then
         call finalize_2d()
         exit
      end if
   end do
   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine pre_eat_model_initialize()
   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "model(no observation program present)"
      ensemble_only=.true.
   else
      signal=signal_initialize
   end if
   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) then
         write(stderr,*) "model(no filter program present)"
      end if
   end if
end subroutine pre_eat_model_initialize

!-----------------------------------------------------------------------

subroutine post_eat_model_initialize()
   integer :: j
   character(len=32) :: fn
   sim_start = strptime(trim(start), time_format)
   sim_stop  = strptime(trim(stop), time_format)
   if (verbosity >= info) then
      write(stderr,'(x,2A,I4)') 'model(sim_start,total_steps) ',sim_start%isoformat(),total_steps
   end if
   if (.not. ensemble_only) then
      start_time=sim_start
      signal=signal_initialize+signal_integrate
   else
      write(fn,'(A,I0.4,A)') 'model_2d_',member,'.dat'
      open(newunit=outunit,file=fn)
      write(outunit,'(A,I4)') 'N= ',N
      do j=1,ny
         write(outunit,'(36F8.4)') field(j, :)
      end do
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
   end if
   N=N+1
end subroutine pre_eat_model_integrate

!-----------------------------------------------------------------------

subroutine post_eat_model_integrate()
   integer :: j
   if (.not. ensemble_only) then
      N=N+1
   else
      write(outunit,'(A,I4)') 'N= ',N
      do j=1,ny
         write(outunit,'(36F8.4)') field(j, :)
      end do
   end if
end subroutine post_eat_model_integrate

!-----------------------------------------------------------------------

subroutine initialize_2d()
   integer :: i,j
   real(real64),parameter :: pi=acos(-1._real64)
   real(real64),parameter :: phi=asin(0.5_real64)
   if (verbosity >= info) write(stderr,*) 'model(initialize): '
   start='1998-01-01 00:00:00'
   ! reproduce the fields read in the PDAF tutorial
   DO j = 1, nx
      DO i = ny,1,-1
         field(i, j) = sin(2._real64*pi*(j+2*(i-1)-1)/nx+phi)
      END DO
   END DO
end subroutine initialize_2d

subroutine integrate_2d()
   integer :: i,j
   real(real64) :: store
   if (verbosity >= info) write(stderr,'(A,I5,A,I5)') ' model(integrate): ',N,' of ',total_steps
   ! *** Time step: Shift field vertically ***
   DO j = 1, nx
      store = field(ny, j)
      DO i = ny-1,1,-1
         field(i+1, j) = field(i, j)
      END DO
      field(1, j) = store
   END DO
   call sleep(1)
end subroutine integrate_2d

subroutine finalize_2d()
   if (verbosity >= info) write(stderr,*) 'model(finalize)'
end subroutine finalize_2d

!-----------------------------------------------------------------------

end program eat_model_2d
