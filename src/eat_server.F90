! Copyright (C) 2021 Bolding & Bruggeman

program eat_server

   !! An ensemble/assimilation server

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use mpi_config
   IMPLICIT NONE

   integer :: ierr
   integer :: nensemble
   integer :: state_size=1234 !!!!!KB
   real(real64), allocatable :: all_states(:,:)
   logical :: do_observations=.false.
   logical :: do_ensemble=.false.
   logical :: do_assimilation=.false.
   integer :: nreqs
   integer, allocatable :: all_reqs(:)
   integer, allocatable :: all_stats(:,:)
!-----------------------------------------------------------------------
   call init_server()
   call do_server()
   call finish_server()
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_server()

   !! A observation handler stub

   use gotm, only: initialize_gotm,yaml_file

   ! Local variables
   logical :: flag
   integer :: n
   character(len=128) :: fname,strbuf
!-----------------------------------------------------------------------
   call cmdline()

#ifdef _PDAF_
   call init_mpi_config(color_obs+color_model+color_filter)
#else
   call init_mpi_config(color_obs+color_model)
#endif

   write(strbuf, "(A,I0.4)") 'gotm_', 0
   yaml_file = TRIM(strbuf) // '.yaml'
   fname = TRIM(strbuf) // '.stderr'
   open(error_unit,file=fname)
   fname = TRIM(strbuf) // '.stdout'
   open(output_unit,file=fname)
   call initialize_gotm()
   close(error_unit)
   close(output_unit)

   open(error_unit,file='eat_server.log')

   call MPI_Comm_size(MPI_COMM_model,nensemble,ierr)
   nensemble=nensemble-1
   if (nensemble > 1) then
      do_ensemble=.true.
      nreqs=nensemble
      call MPI_COMM_TEST_INTER(MPI_COMM_model,flag,ierr)
      if (flag) then
         write(error_unit,*) 'Inter-communicator: server-worker created'
      else
         write(error_unit,*) 'Intra-communicator: server-worker created'
      end if
   end if

   call MPI_Comm_size(MPI_COMM_obs,n,ierr)
   !KBn=n-1
   if (n == 2) then
      do_observations=.true.
      nreqs=nreqs+1
      call MPI_COMM_TEST_INTER(MPI_COMM_obs,flag,ierr)
      if (flag) then
         write(error_unit,*) 'Inter-communicator: server-obs created'
      else
         write(error_unit,*) 'Intra-communicator: server-obs created'
      end if
   end if

   if (do_ensemble .and. do_observations) do_assimilation=.true.

   if (do_assimilation) then
      allocate(all_states(state_size,nensemble))
      allocate(all_reqs(nreqs))
      allocate(all_stats(MPI_STATUS_SIZE,nreqs))
   end if
   if (.not. do_ensemble .and. .not. do_assimilation) then
      write(error_unit,*) './server_handler -h'
      stop 0
   end if
end subroutine init_server

!-----------------------------------------------------------------------

subroutine do_server()

   !! A observation handler stub

   ! Local variables
   integer :: stat(MPI_STATUS_SIZE)
   integer :: nobs
   real(real64), allocatable :: obs(:)
!KB   integer :: i,n
   character(len=32) :: timestr,t1='',t2=''
   logical :: halt=.false.
   integer :: filter_reqs(2)
   integer :: filter_stats(MPI_STATUS_SIZE,2)
   integer :: send_signal(4)
!-----------------------------------------------------------------------
   if (do_ensemble .and. .not. do_assimilation) then
      write(error_unit,'(a)') 'Only doing ensemble:'
      call ensemble_integration()
   else
      write(error_unit,'(a)') 'Doing assimilation:'
      do
         ! receive next observation time
         call MPI_RECV(t2,19,MPI_CHARACTER,1,1,MPI_COMM_obs,stat,ierr)
         if (ierr /= MPI_SUCCESS) then
            call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
         end if
         if(trim(t2) == "0000-00-00 00:00:00") then
            t2=''
            halt=.true.
         end if
         call ensemble_integration(t1,t2,all_states)

         if (.not. halt) then
            ! recieve observations
            call MPI_RECV(nobs,1,MPI_INTEGER,1,1,MPI_COMM_obs,stat,ierr)
            if (ierr /= MPI_SUCCESS) then
               call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
            end if
            if (nobs > 0) then
               if (.not. allocated(obs)) allocate(obs(nobs))
               if (allocated(obs) .and. nobs > size(obs)) then
                   deallocate(obs)
                   allocate(obs(nobs))
               end if
               call MPI_IRECV(obs,nobs,MPI_DOUBLE,1,1,MPI_COMM_obs,all_reqs(nreqs),ierr)
               write(error_unit,'(a,i6)') 'observations(receive) --> ',nobs
               write(error_unit,*) 'server - recv obs ',sum(obs)/nobs
            end if
            ! waiting for ensembles and observations to arrive
            if (allocated(all_states)) then
               call MPI_WAITALL(nreqs,all_reqs,all_stats(:,:),ierr)
               if (ierr /= MPI_SUCCESS) then
                  call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
               end if
            end if
         end if
#ifdef _PDAF_
         if (.not. halt) then
            write(error_unit,'(a)') 'pdaf(filter)'
            send_signal(1)=1
            send_signal(2)=state_size
            send_signal(3)=nensemble
            send_signal(4)=nobs
            call MPI_SEND(send_signal,4,MPI_INTEGER,1,1,MPI_COMM_filter,ierr)
#if 1
            call MPI_ISEND(all_states,nensemble*state_size,MPI_DOUBLE,1,1,MPI_COMM_filter,filter_reqs(1),ierr)
            call MPI_ISEND(obs,nobs,MPI_DOUBLE,1,1,MPI_COMM_filter,filter_reqs(2),ierr)
            call MPI_WAITALL(2,filter_reqs,filter_stats,ierr)
#endif
         else
            send_signal(1)=-1
            call MPI_SEND(send_signal,4,MPI_INTEGER,1,1,MPI_COMM_filter,ierr)
         end if
#endif
         if(halt) exit
         t1=t2
      end do
   end if
end subroutine do_server

!-----------------------------------------------------------------------

subroutine finish_server()

   !! A observation handler stub

!-----------------------------------------------------------------------
   call MPI_Finalize(ierr)
end subroutine finish_server

!-----------------------------------------------------------------------

subroutine cmdline()
   character(len=1024) :: arg
   integer :: n,i
!-----------------------------------------------------------------------
   n = command_argument_count()
   i = 1
   do while (i <= n)
      call get_command_argument(i, arg)
      select case (arg)
      case ('-h', '--help')
         call print_help()
         stop 0
      end select
      i = i+1
   end do
end subroutine cmdline

!-----------------------------------------------------------------------

subroutine print_help()
   print '(a)', 'Usage: server_handler [OPTIONS]'
   print '(a)', '   or: server_handler -n 1 ./server_handler : -n <n> worker_handler'
   print '(a)', '       ensemble only execution'
   print '(a)', '   or: server_handler -n 1 ./server_handler : -n 1 observation_handler : -n <n> worker_handler'
   print '(a)', '       data assimilation execution'
   print '(a)', ''
   print '(a)', 'Options:'
   print '(a)', ''
   print '(a)', '  -h, --help                print usage information and exit'
   print '(a)', ''
end subroutine print_help

!-----------------------------------------------------------------------

subroutine ensemble_integration(t1,t2,all_states)

   !! A observation handler stub

   use time, only: timestep,jul1,secs1,jul2,secs2
   use time, only: MinN,MaxN
   use time, only: read_time_string,time_diff

   ! Subroutine arguments
   character(len=*), intent(in), optional :: t1,t2
   real(real64), intent(inout), optional :: all_states(:,:)

   ! Local variables
   integer :: stat(MPI_STATUS_SIZE)
   integer :: send_signal(4)
   integer :: member
   integer :: julday,secs
   integer :: nsecs
   integer, save :: n
!-----------------------------------------------------------------------
   send_signal(1)=signal_integrate
   send_signal(2)=-1   ! new MaxN
   send_signal(3)=MinN ! Global MinN
   send_signal(4)=MaxN ! Global MaxN
   if (present(t1)) then
      if (len(trim(t1)) == 0) then
         send_signal(1)=send_signal(1)+signal_initialize
      end if
   else
      send_signal(1)=send_signal(1)+signal_initialize
   end if
   if (present(t2)) then
      if (len(trim(t2)) .gt. 0) then
         call read_time_string(t2,julday,secs)
         nsecs = time_diff(julday,secs,jul1,secs1)
         send_signal(2) = nint(nsecs/timestep)
      else
         send_signal(1)=send_signal(1)+signal_finalize
      end if
   else
      send_signal(1)=send_signal(1)+signal_finalize
   end if

   if (present(all_states)) then
      send_signal(1) = send_signal(1)+signal_send_state
   end if

   if (iand(send_signal(1),signal_initialize) == signal_initialize) then
      write(error_unit,'(a,i4,a)') 'ensemble(initialize)  --> ',nensemble,' members'
   end if
   if (iand(send_signal(1),signal_integrate) == signal_integrate) then
      n=0
      write(error_unit,'(a,i4,a,i6)') 'ensemble(integrate)   --> ',nensemble,' members - n=',n
   end if
   if (iand(send_signal(1),signal_finalize) == signal_finalize) then
      write(error_unit,'(a,i4,a)') 'ensemble(finalize)    --> ',nensemble,' members'
      send_signal(2) = MaxN
   end if

   do member=1,nensemble
      call MPI_SEND(send_signal,4,MPI_INTEGER,member,member,MPI_COMM_model,ierr)
      if (present(all_states)) then
         call MPI_IRECV(all_states(:,member),state_size,MPI_DOUBLE,member,MPI_ANY_TAG,MPI_COMM_model,all_reqs(member),ierr)
         write(error_unit,*) 'server - recv state ',member,sum(all_states(:,member))/state_size
      end if
   end do
end subroutine ensemble_integration

!-----------------------------------------------------------------------

subroutine do_filter(framework,all_states,obs)

   ! Subroutine arguments
   integer, intent(in) :: framework
   real(real64), intent(inout) :: all_states(:,:)
   real(real64), intent(in) :: obs(:)
!KB   call do_filter(frame_work,all_states,obs)
end subroutine do_filter

!-----------------------------------------------------------------------

end program eat_server
