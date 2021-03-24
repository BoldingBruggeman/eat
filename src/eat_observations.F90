! Copyright (C) 2021 Bolding & Bruggeman

program eat_observations

   !! A observation handler example program in Fortran.
   !! No real observation handling but illustrates communication 
   !! with the server.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use mpi_config
   IMPLICIT NONE

   integer :: ierr
   character(32), allocatable  :: obs_times(:)
!-----------------------------------------------------------------------

   call init_observations()
   call do_observations()
   call finish_observations()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_observations()

   !! Observations are obtained from an external file - prepare to send to the server

   ! Local variables 
   integer :: unit,ios
   character(len=32) :: buf
   integer :: i,n
!-----------------------------------------------------------------------
   call init_mpi_config(obs_color)

   open(newunit=unit,file="obs_times.dat",status='old',action='read',iostat=ios)
   if (ios /= 0) stop 'init_obs(): unable to open obs_times.dat for reading'
   n=0
   do
      read(unit,*,iostat=ios) buf
      if (ios < 0) exit
      if(len(buf) > 0) n=n+1
   end do
   allocate(obs_times(n))
   rewind unit
   do i=1,n
      read(unit,*,iostat=ios) buf
      if (ios < 0) exit
      if(len(buf) > 0) obs_times(i) = trim(buf)
   end do
end subroutine init_observations

!-----------------------------------------------------------------------

subroutine do_observations()

   !! Handle observations and send size and observations to the server

   ! Local variables 
   integer :: stat(MPI_STATUS_SIZE)
   integer :: n,nobs
   real(real64), allocatable :: obs(:)
   character(len=32) :: timestr,halt="0000-00-00 00:00:00"
!-----------------------------------------------------------------------
   do n=1,size(obs_times)
      call MPI_SEND(obs_times(n),19,MPI_CHARACTER,0,1,MPI_COMM_obs,ierr)
      nobs=10000*n
      call MPI_SEND(nobs,1,MPI_INTEGER,0,1,MPI_COMM_obs,ierr)
      if (nobs > 0) then
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (allocated(obs) .and. nobs > size(obs)) then
             deallocate(obs)
             allocate(obs(nobs))
         end if
         CALL RANDOM_NUMBER(obs)
         call MPI_SEND(obs,nobs,MPI_DOUBLE,0,1,MPI_COMM_obs,ierr)
      end if
   end do
   call MPI_SEND(halt,16,MPI_CHARACTER,0,1,MPI_COMM_obs,ierr)
end subroutine do_observations

!-----------------------------------------------------------------------

subroutine finish_observations()

   !! Finish the eat_observation program by calling MPI_Finalize()

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_observations

!-----------------------------------------------------------------------

end program eat_observations
