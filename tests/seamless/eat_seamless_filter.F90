! Copyright (C) 2021 Bolding & Bruggeman

program eat_seamless_filter

   !! A simple filter that does exchanges with model and observations -
   !! if present

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   IMPLICIT NONE

   integer :: ierr
   logical :: have_obs=.true.
   logical :: have_model=.true.
   integer :: state_size,ensemble_size
   integer, allocatable :: model_reqs(:)
   integer, allocatable :: model_stats(:,:)
   real(real64), allocatable :: model_states(:,:)
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info
!-----------------------------------------------------------------------

   call init_filter()
   call do_filter()
   call finish_filter()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_filter()

   !! Initialize EAT/PDAF component

   ! Local variables
   integer :: stat(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------
   call init_eat_config(color_filter+verbosity)

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      write(stderr,*) "filter(no model executable present)"
      have_model=.false.
   else
!KBwrite(stderr,*) 'AAA1 ',state_size,ensemble_size
!KBstate_size=1234
      call MPI_RECV(state_size,1,MPI_INTEGER,1,MPI_ANY_TAG,EAT_COMM_model_filter,stat,ierr)
      ensemble_size=size_model_filter_comm-size_filter_comm
!KBwrite(stderr,*) 'AAA2 ',state_size,ensemble_size
      allocate(model_states(state_size,ensemble_size))
      allocate(model_reqs(ensemble_size))
      allocate(model_stats(MPI_STATUS_SIZE,ensemble_size))
   end if

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      write(stderr,*) "filter(no observation executable present)"
      have_obs=.false.
   end if

end subroutine init_filter

!-----------------------------------------------------------------------

subroutine do_filter()

   !! Receive observation and states and do the filter calculation

   ! Local variables
   integer :: recv_signal(4)
   integer :: nobs
   real(real64), allocatable :: states(:,:)
   real(real64), allocatable :: obs(:)
   integer :: stat(MPI_STATUS_SIZE)
   integer :: obs_request
   integer :: m
!-----------------------------------------------------------------------
   do
      if (have_model) then
         do m=1,ensemble_size
!KBwrite(stderr,*) 'BBB1 ',m
            call MPI_IRECV(model_states(:,m),state_size,MPI_DOUBLE,m,MPI_ANY_TAG,EAT_COMM_model_filter,model_reqs(m),ierr)
!KBwrite(stderr,*) 'BBB2 ',m,state_size
         end do
      end if

      if (have_obs) then
         call MPI_RECV(nobs,1,MPI_INTEGER,0,1,EAT_COMM_obs_filter,stat,ierr)
         if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<- nobs) ',nobs
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (nobs > 0) then
            if (nobs > size(obs)) then
               deallocate(obs)
               allocate(obs(nobs))
            end if
            call MPI_IRECV(obs,nobs,MPI_DOUBLE,0,1,EAT_COMM_obs_filter,obs_request,ierr)
         else
            exit
         end if
      end if

      if (have_model) then
         call MPI_WAITALL(ensemble_size,model_reqs,model_stats(:,:),ierr)
         do m=1,ensemble_size
            if (verbosity >= info) write(stderr,'(x,A,I4,F8.5)') 'filter(<- state)',m,sum(model_states(:,m))/state_size
         end do
      end if
      if (have_obs) then
         call MPI_WAIT(obs_request,stat,ierr)
         if (verbosity >= info) write(stderr,*) 'filter(<- obs) ',sum(obs)/nobs
      end if
   end do
end subroutine do_filter

!-----------------------------------------------------------------------

subroutine finish_filter()

   !! Cleanup and finalize the EAT/PDAF component

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_filter

!-----------------------------------------------------------------------

end program eat_seamless_filter
