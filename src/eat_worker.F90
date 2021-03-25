! Copyright (C) 2021 Bolding & Bruggeman

#define NO_GOTM
#define _ASYNC_

program eat_worker

   !! An implementation of GOTM in an ensemble context

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use mpi_config
#ifndef NO_GOTM
   use gotm, only: init_gotm, time_loop, clean_up
#endif
   IMPLICIT NONE

   integer :: ierr
!-----------------------------------------------------------------------

   call init_worker()
   call do_worker()
   call finish_worker()

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_worker()

   !! An implementation of GOTM in an ensemble context

!-----------------------------------------------------------------------
   call init_mpi_config(color_model)
end subroutine init_worker

!-----------------------------------------------------------------------

subroutine do_worker()

   !! An implementation of GOTM in an ensemble context

   ! Subroutine arguments

   use gotm, only: yaml_file,output_id
   use time, only: MinN,MaxN

   ! Local variables
   integer :: state_size
   real(real64), allocatable :: state(:)
   integer :: member
   logical :: flag
   integer :: stat(MPI_STATUS_SIZE)
   character(len=128) :: fname,strbuf
   integer :: recv_signal(5)
   integer :: request
!-----------------------------------------------------------------------
   do
      call MPI_RECV(recv_signal,5,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_model,stat,ierr)
      member=stat(MPI_TAG)

      if (iand(recv_signal(1),signal_initialize) == signal_initialize) then
         write(output_id, "(A,I0.4)") '_', member
         write(strbuf, "(A,I0.4)") 'gotm_', member
         yaml_file = TRIM(strbuf) // '.yaml'
         fname = TRIM(strbuf) // '.stderr'
         open(error_unit,file=fname)
         fname = TRIM(strbuf) // '.stdout'
         open(output_unit,file=fname)
         call init_gotm()
         state_size=1234 !!!!KB
         if (iand(recv_signal(1),signal_send_state) == signal_send_state) then
            allocate(state(state_size))
         end if
!         if (MinN /= recv_signal(4)) stop 'MinN diffrent between server and worker'
!         if (MaxN /= recv_signal(5)) stop 'MaxN diffrent between server and worker'
      end if

      if (iand(recv_signal(1),signal_integrate) == signal_integrate) then
         if (recv_signal(2) /= -1) MinN=recv_signal(2)
         if (recv_signal(3) /= -1) MaxN=recv_signal(3)
         call time_loop()
      end if

      if (iand(recv_signal(1),signal_send_state) == signal_send_state) then
         CALL RANDOM_NUMBER(state)
         state=state+rank-1
#ifdef _ASYNC_
         call MPI_ISEND(state,state_size,MPI_DOUBLE,0,member,MPI_COMM_model,request,ierr)
         call MPI_WAIT(request,stat,ierr)
#else
         call MPI_SEND(state,state_size,MPI_DOUBLE,0,member,MPI_COMM_model,ierr)
#endif
      end if

      if (iand(recv_signal(1),signal_finalize) == signal_finalize) then
         call clean_up()
         exit
      end if
   end do
end subroutine do_worker

!-----------------------------------------------------------------------

subroutine finish_worker()

   !! An implementation of GOTM in an ensemble context

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_worker

!-----------------------------------------------------------------------

#ifdef NO_GOTM
subroutine init_gotm()
   use gotm, only: yaml_file
   write(error_unit,*) trim(yaml_file)
end subroutine init_gotm
subroutine time_loop()
   use time, only: MinN,MaxN
   write(error_unit,*) MinN,MaxN,MaxN-MinN+1
   call sleep(1)
end subroutine time_loop
subroutine clean_up()
   write(error_unit,*) 'cleaning up'
end subroutine clean_up
#endif

!-----------------------------------------------------------------------

end program eat_worker
