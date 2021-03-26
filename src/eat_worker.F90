! Copyright (C) 2021 Bolding & Bruggeman

#define NO_GOTM

program eat_worker

   !! An implementation of GOTM in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use mpi_config
#ifndef NO_GOTM
   use gotm, only: initialize_gotm, integrate_gotm, finalize_gotm
#endif
   use time, only: MinN,MaxN
   IMPLICIT NONE

   integer :: ierr
   integer :: member
   logical :: flag
   integer :: stat(MPI_STATUS_SIZE)
   integer :: recv_signal(4)
   integer :: request
   integer :: state_size
   real(real64), allocatable :: state(:)
!-----------------------------------------------------------------------
   call init_mpi_config(color_model)

   do
      call MPI_RECV(recv_signal,4,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_model,stat,ierr)
      member=stat(MPI_TAG)

      if (iand(recv_signal(1),signal_initialize) == signal_initialize) then
         output: block
         use gotm, only: yaml_file,output_id
         character(len=128) :: fname,strbuf
         write(output_id, "(A,I0.4)") '_', member
         write(strbuf, "(A,I0.4)") 'gotm_', member
         yaml_file = TRIM(strbuf) // '.yaml'
         fname = TRIM(strbuf) // '.stderr'
         open(error_unit,file=fname)
         fname = TRIM(strbuf) // '.stdout'
         open(output_unit,file=fname)
         end block output

         call initialize_gotm()
         state_size=1234 !!!!KB
         if (iand(recv_signal(1),signal_send_state) == signal_send_state) then
            allocate(state(state_size))
         end if
      end if

      if (iand(recv_signal(1),signal_integrate) == signal_integrate) then
         if (recv_signal(2) /= -1) MaxN=recv_signal(2)
         call integrate_gotm()
         MinN=MaxN+1
      end if

      if (iand(recv_signal(1),signal_send_state) == signal_send_state) then
         CALL RANDOM_NUMBER(state)
         state=state+rank-1
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,MPI_COMM_model,request,ierr)
         call MPI_WAIT(request,stat,ierr)
      end if

      if (iand(recv_signal(1),signal_finalize) == signal_finalize) then
         call finalize_gotm()
         exit
      end if
   end do

   call MPI_Finalize(ierr)

!-----------------------------------------------------------------------

#ifdef NO_GOTM
contains

subroutine initialize_gotm()
   use gotm, only: yaml_file
   use time, only: MinN,MaxN
   write(error_unit,*) 'initialize_gotm()'
   MinN=1
   write(error_unit,*) trim(yaml_file)
end subroutine initialize_gotm
subroutine integrate_gotm()
   use time, only: MinN,MaxN
   write(error_unit,*) 'integrate_gotm()'
   write(error_unit,*) MinN,MaxN,MaxN-MinN+1
   call sleep(1)
end subroutine integrate_gotm
subroutine finalize_gotm()
   write(error_unit,*) 'finalize_gotm()'
end subroutine finalize_gotm
#endif

!-----------------------------------------------------------------------

end program eat_worker
