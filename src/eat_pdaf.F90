! Copyright (C) 2021 Bolding & Bruggeman

program eat_pdaf

   !! A wrapper around the 'off_line PDAF' implmentation to keep it alive during
   !! ensemble simulations

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use mpi_config
   IMPLICIT NONE

   integer :: ierr
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
   call init_mpi_config(color_filter)
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
      call MPI_RECV(recv_signal,4,MPI_INTEGER,0,1,MPI_COMM_filter,stat,ierr)
      if (ierr /= MPI_SUCCESS) then
         call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      write(error_unit,*) 'filter - signal ',recv_signal
      if (recv_signal(1) == -1) then
         exit
      else
         state_size=recv_signal(2)
         nensemble=recv_signal(3)
         nobs=recv_signal(4)
#if 1
         if (.not. allocated(states)) allocate(states(state_size,nensemble))
         if (state_size > size(states,1) .or. nensemble > size(states,2)) then
            deallocate(states)
            allocate(states(state_size,nensemble))
         end if
         call MPI_IRECV(states,nensemble*state_size,MPI_DOUBLE,0,1,MPI_COMM_filter,filter_reqs(1),ierr)
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (nobs > size(obs)) then
            deallocate(obs)
            allocate(obs(nobs))
         end if
         call MPI_IRECV(obs,nobs,MPI_DOUBLE,0,1,MPI_COMM_filter,filter_reqs(2),ierr)
         call MPI_WAITALL(2,filter_reqs,filter_stats,ierr)
         write(error_unit,*) 'filter - recv obs ',sum(obs)/nobs
#endif
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
