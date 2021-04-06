! Copyright (C) 2021 Bolding & Bruggeman

program eat_pdaf

   !! A wrapper around the 'off_line PDAF' implmentation to keep it alive during
   !! ensemble simulations

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use pdaf_mod_filter
   use eat_config
   IMPLICIT NONE

   integer :: ierr
   logical :: have_obs=.true.
   logical :: have_model=.true.
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=warn
!-----------------------------------------------------------------------

   call eat_init_pdaf()
   call eat_do_pdaf()
   call eat_finish_pdaf()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine eat_init_pdaf()

   !! Initialize EAT/PDAF component

   ! Local variables
!-----------------------------------------------------------------------
   call init_eat_config(color_filter+verbosity)

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity >= info) write(stderr,*) "filter(no observation executable present)"
      have_obs=.false.
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= info) write(stderr,*) "filter(no model executable present)"
      have_model=.false.
   end if
end subroutine eat_init_pdaf

!-----------------------------------------------------------------------

subroutine eat_do_pdaf()

   !! Handle observations and send size and observations to the server

   ! Local variables
   integer :: recv_signal(4)
   integer :: nensemble
   integer :: state_size
   integer :: nobs
   real(real64), allocatable :: states(:,:)
   real(real64), allocatable :: obs(:)
   integer :: stat(MPI_STATUS_SIZE)
   integer :: filter_reqs(2)
   integer :: filter_stats(MPI_STATUS_SIZE,2)
!-----------------------------------------------------------------------
   do
#if 0
      call MPI_RECV(recv_signal,4,MPI_INTEGER,0,1,EAT_COMM_filter,stat,ierr)
      if (ierr /= MPI_SUCCESS) then
         call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      if (verbosity >= debug) write(stderr,*) 'filter - signal ',recv_signal
      if (recv_signal(1) == -1) then
         exit
      else
         state_size=recv_signal(2)
         nensemble=recv_signal(3)
         nobs=recv_signal(4)
         if (.not. allocated(states)) allocate(states(state_size,nensemble))
         if (state_size > size(states,1) .or. nensemble > size(states,2)) then
            deallocate(states)
            allocate(states(state_size,nensemble))
         end if
         call MPI_IRECV(states,nensemble*state_size,MPI_DOUBLE,0,1,EAT_COMM_filter,filter_reqs(1),ierr)
#endif
      if (have_obs) then
         call MPI_RECV(nobs,1,MPI_INTEGER,0,1,EAT_COMM_obs_filter,stat,ierr)
         if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<- nobs) ',nobs
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (nobs > 0) then
            if (nobs > size(obs)) then
               deallocate(obs)
               allocate(obs(nobs))
            end if
            call MPI_IRECV(obs,nobs,MPI_DOUBLE,0,1,EAT_COMM_obs_filter,filter_reqs(1),ierr)
            call MPI_WAITALL(1,filter_reqs,filter_stats,ierr)
            if (verbosity >= info) write(stderr,*) 'filter(<- obs) ',sum(obs)/nobs
         else
            exit
         end if
      end if
   end do
end subroutine eat_do_pdaf

!-----------------------------------------------------------------------

subroutine eat_finish_pdaf()

   !! Cleanup and finalize the EAT/PDAF component

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine eat_finish_pdaf

!-----------------------------------------------------------------------

end program eat_pdaf
