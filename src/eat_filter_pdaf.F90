! Copyright (C) 2021 Bolding & Bruggeman

#undef _USE_PDAF_
#define _USE_PDAF_

program eat_filter_pdaf

   !! A wrapper around the 'off_line PDAF' implmentation to keep it alive during
   !! ensemble simulations

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use pdaf_wrapper
   IMPLICIT NONE

   integer :: nobs=0
   logical :: have_obs=.true.
   logical :: have_model=.true.
   integer :: state_size,ensemble_size
   integer, allocatable :: reqs(:)
   integer, allocatable :: stats(:,:)
   real(real64), pointer, contiguous :: model_states(:,:)
   integer, parameter :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info

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
   integer :: stat(MPI_STATUS_SIZE)
   integer :: ierr

   logical :: fileexists
   character(len=*), parameter :: nmlfile='eat_pdaf.nml'
   integer :: nmlunit
   namelist /nml_filter_pdaf/ verbosity
!-----------------------------------------------------------------------
   INQUIRE(FILE=nmlfile, EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile,status='old',action='read')
      read(nmlunit,nml=nml_filter_pdaf)
      close(nmlunit)
      if (verbosity >= warn) write(stderr,*) 'filter(read namelist)'
   end if

   call init_eat_config(color_filter+verbosity)

   if (EAT_COMM_obs_filter == MPI_COMM_NULL) then
      if (verbosity >= info) write(stderr,*) "filter(no observation executable present)"
      have_obs=.false.
   end if

   if (EAT_COMM_model_filter == MPI_COMM_NULL) then
      if (verbosity >= info) write(stderr,*) "filter(no model executable present)"
      have_model=.false.
   else
      call MPI_RECV(state_size,1,MPI_INTEGER,1,MPI_ANY_TAG,EAT_COMM_model_filter,stat,ierr)
      if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<-- state_size) ',state_size
      ensemble_size=size_model_filter_comm-size_filter_comm
      allocate(reqs(2+ensemble_size))
      allocate(stats(MPI_STATUS_SIZE,2+ensemble_size))
   end if

#ifdef _USE_PDAF_
   CALL init_pdaf(EAT_COMM_filter, state_size, ensemble_size, model_states, ierr)
   if (ierr /= 0) then
      call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
   else
      write(error_unit,*) 'filter(PDAF is initialized): ',shape(model_states)
   end if
#else
   allocate(model_states(state_size,ensemble_size))
#endif
end subroutine eat_init_pdaf

!-----------------------------------------------------------------------

subroutine eat_do_pdaf()

   !! Get observations and states and do the PDAF/assimilation step

   ! Local variables
   real(real64), allocatable :: states(:,:)
   integer :: stat(MPI_STATUS_SIZE)
   integer :: obs_stats(MPI_STATUS_SIZE,2)
   integer :: obs_requests(2)
   integer :: m
   integer :: ierr
   character(len=19) :: timestr
!-----------------------------------------------------------------------
   do
      if (have_obs) then
         call MPI_RECV(timestr,19,MPI_CHARACTER,0,tag_timestr,EAT_COMM_obs_filter,stat,ierr)
         call MPI_RECV(nobs,1,MPI_INTEGER,0,tag_nobs,EAT_COMM_obs_filter,stat,ierr)
         if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<-- nobs) ',nobs
      end if

      if (have_obs .and. nobs > 0) then
         if (.not. associated(iobs)) allocate(iobs(nobs))
         if (nobs > size(iobs)) then
            deallocate(iobs)
            allocate(iobs(nobs))
         end if
         if (.not. associated(obs)) allocate(obs(nobs))
         if (nobs > size(obs)) then
            deallocate(obs)
            allocate(obs(nobs))
         end if
         call MPI_IRECV(iobs(1:nobs),nobs,MPI_INTEGER,0,tag_iobs,EAT_COMM_obs_filter,reqs(1),ierr)
         call MPI_IRECV(obs(1:nobs),nobs,MPI_DOUBLE,0,tag_obs,EAT_COMM_obs_filter,reqs(2),ierr)
         if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<-- iobs, obs)'
      end if

      if (have_model .and. nobs > 0) then
!KB      if (have_model) then
         do m=1,ensemble_size
            call MPI_IRECV(model_states(:,m),state_size,MPI_DOUBLE,m,tag_forecast,EAT_COMM_model_filter,reqs(2+m),ierr)
            if(ierr /= MPI_SUCCESS) THEN
               write(stderr,*) 'Fatal error (PDAF): Unable to receive: ',m
               call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
            end if
         end do
      end if

     if (verbosity >= info) write(stderr,'(x,A)') 'filter(<-- state)'
     call MPI_WAITALL(ensemble_size,reqs,stats(:,:),ierr)

      if (have_obs .and. have_model .and. nobs > 0) then
#ifdef _USE_PDAF_
         ! Begin PDAF specific part
         ! from .../tutorial/classical/offline_2D_serial/main_offline.F90
         if (verbosity >= info) write(stderr,'(x,A)') 'filter(--> PDAF)'
         CALL assimilation_pdaf()
         if (verbosity >= info) write(stderr,'(x,A)') 'filter(<-- PDAF)'
         ! End PDAF specific part
#endif

         do m=1,ensemble_size
            call MPI_ISEND(model_states(:,m),state_size,MPI_DOUBLE,m,tag_analysis,EAT_COMM_model_filter,reqs(2+m),ierr)
            if(ierr /= MPI_SUCCESS) THEN
               write(stderr,*) 'Fatal error (PDAF): Unable to send: ',m
               call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
            end if
         end do
         m=2+ensemble_size
         call MPI_WAITALL(ensemble_size,reqs(3:m),stats(:,3:m),ierr)
         if(ierr /= MPI_SUCCESS) THEN
            write(stderr,*) 'Fatal error (PDAF): Unable to wait'; call flush(stderr)
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
         if (verbosity >= info) write(stderr,'(x,A)') 'filter(--> state)'
      end if

      if (nobs < 0) then
         if (verbosity >= info) write(stderr,*) 'filter(exit)'
         exit
      end if
   end do
end subroutine eat_do_pdaf

!-----------------------------------------------------------------------

subroutine eat_finish_pdaf()

   !! Cleanup and finalize the EAT/PDAF component

!-----------------------------------------------------------------------
   integer :: ierr
   call finish_pdaf()
   call MPI_Finalize(ierr)
end subroutine eat_finish_pdaf

end program eat_filter_pdaf
