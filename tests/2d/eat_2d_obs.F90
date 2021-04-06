! Copyright (C) 2021 Bolding & Bruggeman

program eat_2d_obs

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
   integer :: verbosity=debug

   integer, parameter ::nx=36,ny=18
   integer :: total_steps=5
   integer :: N=0
   real(real64) :: field(ny,nx)
   integer :: outunit
   integer :: i
   character(len=*), parameter :: nmlfile='eat_2d.nml'
   logical :: fileexists
   integer :: nmlunit
   logical :: all_verbose=.true.
   namelist /nml_2d_obs/ verbosity,all_verbose
!-----------------------------------------------------------------------
   inquire(FILE="eat_2d.nml", EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_2d_obs)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'obs(read namelist)'
   end if
   call init_eat_config(color_model+verbosity)

!-----------------------------------------------------------------------

end program eat_2d_obs
