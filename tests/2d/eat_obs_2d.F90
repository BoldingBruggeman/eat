! Copyright (C) 2021 Bolding & Bruggeman

program eat_obs_2d

   !! The observation handler example program for the 2d test case

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use fields_2d
   use datetime_module, only: datetime, strptime
   IMPLICIT NONE

   integer :: ierr
   character(len=128)  :: obs_times_file='obs_2d_times.dat'
   character(32), allocatable  :: obs_times(:)
   integer, allocatable :: iobs(:)
   real(real64), allocatable :: obs(:)
   logical :: have_model=.true.
   logical :: have_filter=.true.
   integer :: nmodel=-1
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info
   integer :: filter=1
!-----------------------------------------------------------------------

   call init_observations()
   call do_observations()
   call finish_observations()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_observations()

   !! Observation times are obtained from an external file - prepare to
   !! send to the server

   ! Local variables
   integer :: unit,ios
   character(len=32) :: buf
   integer :: i,n
   character(len=*), parameter :: nmlfile='eat_2d.nml'
   logical :: fileexists
   integer :: nmlunit,outunit
   logical :: all_verbose=.true.
   namelist /nml_obs_2d/ verbosity,all_verbose,obs_times_file,nobsmin,nobsmax,obsstd
!-----------------------------------------------------------------------
   inquire(FILE=nmlfile,EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile)
      read(unit=nmlunit,nml=nml_obs_2d)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'obs(read namelist)'
   end if

   if (nobsmin > nobsmax) then
      stop "init_observations(): nobsmin must be <= nobsmax"
   end if

   call init_eat_config(color_obs+verbosity)

   if (EAT_COMM_obs_model == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "obs(no model program present)"
      have_model=.false.
   else
      nmodel=size_obs_model_comm-size_obs_comm
   end if

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity >= warn) write(stderr,*) "obs(no filter program present)"
      have_filter=.false.
   end if

   open(newunit=unit,file=obs_times_file,status='old',action='read',iostat=ios)
   if (verbosity >= fatal) then
      if (ios /= 0) stop 'obs(unable to open obs_times.dat for reading)'
   else
      if (ios /= 0) stop 1
   end if
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
   integer :: stats(MPI_STATUS_SIZE,2)
   integer :: m ! model counter
   integer :: n,nobs
   real(real64) :: x
   character(len=32) :: timestr,halt="0000-00-00 00:00:00"
   integer :: requests(2)
!-----------------------------------------------------------------------
   do n=1,size(obs_times)
      if (verbosity >= info) write(stderr,*) 'obs(--> time)  ',trim(obs_times(n))
      if (have_model) then
         do m=1,nmodel
            call MPI_SSEND(obs_times(n),19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr)
            if (ierr /= MPI_SUCCESS) then
               if (verbosity >= error) then
                  write(stderr,*) 'obs: failing to send to model process#: ',m
               end if
            end if
         end do
      end if

      call random_number(x)
      nobs=nobsmin+FLOOR((nobsmax+1-nobsmin)*x)
      if (verbosity >= info) write(stderr,'(A,I6)') ' obs(--> nobs)    ',nobs
      if (have_filter) then
         call MPI_SSEND(nobs,1,MPI_INTEGER,filter,1,EAT_COMM_obs_filter,ierr)
         if (ierr /= MPI_SUCCESS) then
            if (verbosity >= error) then
               write(stderr,*) 'obs: failing to send to filter process'
            end if
         end if
      end if
      if (nobs > 0) then
         if (.not. allocated(iobs)) allocate(iobs(nobs))
         if (allocated(iobs) .and. nobs > size(iobs)) then
             deallocate(iobs)
             allocate(iobs(nobs))
         end if
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (allocated(obs) .and. nobs > size(obs)) then
             deallocate(obs)
             allocate(obs(nobs))
         end if
         call get_obs(n,obsstd,iobs(1:nobs),obs(1:nobs))
         block
         character(len=64) :: fn
         TYPE(datetime) :: t
         character(len=14) :: timestamp
         t = strptime(trim(obs_times(n)),time_format)
         timestamp = t%strftime("%Y%m%d%H%M%S")
         write(fn,'(3A)') 'obs_',timestamp,'.dat'
         call write_field(trim(fn),obsfield)
         end block
         if (verbosity >= info)  write(stderr,'(A,F10.6)') ' obs(--> obs)    ',sum(obs)/nobs
      end if
      if (have_filter .and. nobs > 0) then
         call MPI_ISEND(iobs(1:nobs),nobs,MPI_INTEGER,filter,1,EAT_COMM_obs_filter,requests(1),ierr)
         call MPI_ISEND(obs(1:nobs),nobs,MPI_DOUBLE,filter,1,EAT_COMM_obs_filter,requests(2),ierr)
!KB         call sleep(1)
         call MPI_WAITALL(2,requests,stats,ierr)
         if (ierr /= MPI_SUCCESS) then
            if (verbosity >= error) write(stderr,*) 'obs: failing to wait for requests'
         end if
      end if
   end do

   ! Here we must NOT use MPI_SSEND()
   if (have_filter) then
      nobs=-1
      if (verbosity >= info) write(stderr,'(A,I6)') ' obs(--> nobs)    ',nobs
      call MPI_SEND(nobs,1,MPI_INTEGER,filter,1,EAT_COMM_obs_filter,ierr)
   end if
   if (have_model) then
      if (verbosity >= info) write(stderr,*) 'obs(--> exit)'
      do m=1,nmodel
         call MPI_SEND(halt,19,MPI_CHARACTER,m,m,EAT_COMM_obs_model,ierr)
      end do
   end if
end subroutine do_observations

!-----------------------------------------------------------------------

subroutine finish_observations()

   !! Finish the eat_obs_2d program by calling MPI_Finalize()

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_observations

!-----------------------------------------------------------------------

end program eat_obs_2d
