! Copyright (C) 2021 Bolding & Bruggeman

#define _USE_PDAF_
!KB#undef _USE_PDAF_

program eat_filter_pdaf

   !! A wrapper around the 'off_line PDAF' implmentation to keep it alive during
   !! ensemble simulations

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use mpi
   use eat_config
   use pdaf_mod_filter
   use PDAF_interfaces_module
   IMPLICIT NONE

   integer :: nobs
   integer, allocatable :: iobs(:)
   real(real64), allocatable :: obs(:)

   integer :: ierr
   logical :: have_obs=.true.
   logical :: have_model=.true.
   integer :: state_size,ensemble_size
   integer, allocatable :: model_reqs(:)
   integer, allocatable :: model_stats(:,:)
#ifdef _USE_PDAF_
   real(real64), pointer :: model_states(:,:) => null()
#else
   real(real64), allocatable :: model_states(:,:)
#endif
   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info

   logical :: fileexists
   character(len=*), parameter :: nmlfile='eat_pdaf.nml'
   integer :: nmlunit

   ! PDAF variables
   integer :: filtertype=6
!KB
REAL(real64) :: timenow
integer :: doexit,steps
   real(real64) :: rms_obs = 0.05    ! Observation error standard deviation

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
      allocate(model_reqs(ensemble_size))
      allocate(model_stats(MPI_STATUS_SIZE,ensemble_size))
   end if

#ifdef _USE_PDAF_
   CALL init_pdaf(ierr)
   call PDAF_set_ens_pointer(model_states,ierr)
   if (ierr /= 0) then
      call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
   else
      write(error_unit,*) 'filter(PDAF is initialized)'
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
!-----------------------------------------------------------------------
   do
      if (have_obs) then
         call MPI_RECV(nobs,1,MPI_INTEGER,0,tag_nobs,EAT_COMM_obs_filter,stat,ierr)
         if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<-- nobs) ',nobs
      end if

      if (have_model .and. nobs > 0) then
         do m=1,ensemble_size
            call MPI_IRECV(model_states(:,m),state_size,MPI_DOUBLE,m,tag_forecast,EAT_COMM_model_filter,model_reqs(m),ierr)
         end do
      end if

      if (have_obs .and. nobs > 0) then
         if (.not. allocated(iobs)) allocate(iobs(nobs))
         if (nobs > size(iobs)) then
            deallocate(iobs)
            allocate(iobs(nobs))
         end if
         if (.not. allocated(obs)) allocate(obs(nobs))
         if (nobs > size(obs)) then
            deallocate(obs)
            allocate(obs(nobs))
         end if
         if (verbosity >= info) write(stderr,'(A,I6)') ' filter(<-- iobs, obs)'
         call MPI_IRECV(iobs(1:nobs),nobs,MPI_INTEGER,0,tag_iobs,EAT_COMM_obs_filter,obs_requests(1),ierr)
         call MPI_IRECV(obs(1:nobs),nobs,MPI_DOUBLE,0,tag_obs,EAT_COMM_obs_filter,obs_requests(2),ierr)
      end if

      if (have_model .and. nobs > 0) then
         if (verbosity >= info) write(stderr,'(x,A)') 'filter(<-- state)'
         call MPI_WAITALL(ensemble_size,model_reqs,model_stats(:,:),ierr)
         if (verbosity >= debug) then
            do m=1,ensemble_size
               write(stderr,'(x,A,I4,*(F10.5))') 'filter(<-- state)',m,sum(model_states(:,m))/state_size
            end do
         end if
      end if

      if (have_obs .and. nobs > 0) then
         call MPI_WAITALL(2,obs_requests,obs_stats,ierr)
         if (verbosity >= debug) write(stderr,'(A,F10.6)') ' filter(<-- obs) ',sum(obs)/nobs
      end if

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
            call MPI_ISEND(model_states(:,m),state_size,MPI_DOUBLE,m,tag_analysis,EAT_COMM_model_filter,model_reqs(m),ierr)
         end do
         call MPI_WAITALL(ensemble_size,model_reqs,model_stats,ierr)
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
#ifdef _USE_PDAF__
   CALL finalize_pdaf(0) ! Basically CALL PDAF_deallocate()
#endif
   call MPI_Finalize(ierr)
end subroutine eat_finish_pdaf

!-----------------------------------------------------------------------

! Below are the routines implemented to link to the PDAF library - pdaf-d.


SUBROUTINE init_pdaf(stat)

#if 0
  USE mod_parallel, &     ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       rms_obs, incremental, covartype, type_forget, forget, &
       rank_analysis_enkf, locweight, local_range, srange, &
       filename, type_trans, type_sqrt
#endif

   IMPLICIT NONE
   integer, intent(out) :: stat

!KB
integer :: dim_state_p

  integer :: comm_couple, comm_filter, comm_model, n_modeltasks, task_id
  logical :: filterpe=.true.
!KB

! Local variables
   INTEGER :: filter_param_i(7) ! Integer parameter array for filter
   REAL(REAL64)    :: filter_param_r(2) ! Real parameter array for filter
   INTEGER :: status_pdaf       ! PDAF status flag

   ! External subroutines
!KB  EXTERNAL :: init_ens_pdaf  ! Ensemble initialization

   integer :: &
   screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
!KB  integer :: &
!KB  filtertype = 6    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
   integer :: &
   dim_ens = 9      ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
   integer :: &
   subtype = 1       ! (5) Offline mode !KB 5
   integer :: &
   type_trans = 0    ! Type of ensemble transformation
                     !   SEIK/LSEIK and ESTKF/LESTKF:
                     !     (0) use deterministic omega
                     !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                     !     (2) use product of (0) with random orthonormal matrix with
                     !         eigenvector (1,...,1)^T
                     !   ETKF/LETKF:
                     !     (0) use deterministic symmetric transformation
                     !     (2) use product of (0) with random orthonormal matrix with
                     !         eigenvector (1,...,1)^T
   integer :: &
   type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                     !   (0) fixed
                     !   (1) global adaptive
                     !   (2) local adaptive for LSEIK/LETKF/LESTKF
   real :: &
   forget  = 1.0     ! Forgetting factor
   real :: &
   type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
   integer :: &
   incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
   integer :: &
   covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                     !   (0) for dim_ens^-1 (old SEIK)
                     !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
                     !   This parameter has also to be set internally in PDAF_init.
   integer :: &
   rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                     ! in analysis of EnKF; (0) for analysis w/o eigendecomposition


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** specifications for observations ***
!KB  real :: &
!KB  rms_obs = 0.5    ! Observation error standard deviation
                   ! for the Gaussian distribution
! *** Localization settings
   integer :: &
   locweight = 0     ! Type of localizating weighting
                     !   (0) constant weight of 1
                     !   (1) exponentially decreasing with SRANGE
                     !   (2) use 5th-order polynomial
                     !   (3) regulated localization of R with mean error variance
                     !   (4) regulated localization of R with single-point error variance
   integer :: &
   local_range = 0  ! Range in grid points for observation domain in local filters
!KB  srange = local_range  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting

! *** File names
   character(len=128) :: &
   filename = 'output.dat'

   namelist /nml_config_pdaf/ screen,filtertype,subtype,dim_ens

   INQUIRE(FILE=nmlfile, EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile,status='old',action='read')
      read(nmlunit,nml=nml_config_pdaf)
      close(nmlunit)
   end if

   ! initialize variables for call to PDAF_init()
   dim_state_p=state_size
   comm_couple=EAT_COMM_filter ! suggested by Lars
   comm_filter=EAT_COMM_filter
   comm_model=EAT_COMM_filter ! suggested by Lars
   filterpe=.true.
   task_id=1
   n_modeltasks=1

   whichinit: IF (filtertype == 2) THEN
      ! *** EnKF with Monte Carlo init ***
      filter_param_i(1) = dim_state_p ! State dimension
      filter_param_i(2) = dim_ens     ! Size of ensemble
      filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
      filter_param_i(4) = incremental ! Whether to perform incremental analysis
      filter_param_i(5) = 0           ! Smoother lag (not implemented here)
      filter_param_r(1) = forget      ! Forgetting factor
      call PDAF_set_comm_pdaf(EAT_COMM_filter)
      CALL PDAF_init(filtertype, subtype, 0, &
           filter_param_i, 6,&
           filter_param_r, 2, &
           COMM_model, COMM_filter, COMM_couple, &
           task_id, n_modeltasks, filterpe, init_ens_pdaf, &
           screen, status_pdaf)
   ELSE
      ! *** All other filters                       ***
      ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
      filter_param_i(1) = dim_state_p ! State dimension
      filter_param_i(2) = dim_ens     ! Size of ensemble
      filter_param_i(3) = 0           ! Smoother lag (not implemented here)
      filter_param_i(4) = incremental ! Whether to perform incremental analysis
      filter_param_i(5) = type_forget ! Type of forgetting factor
      filter_param_i(6) = type_trans  ! Type of ensemble transformation
      filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
      filter_param_r(1) = forget      ! Forgetting factor
      call PDAF_set_comm_pdaf(EAT_COMM_filter)
      CALL PDAF_init(filtertype, subtype, 0, &
           filter_param_i, 7,&
           filter_param_r, 2, &
           COMM_model, COMM_filter, COMM_couple, &
           task_id, n_modeltasks, filterpe, init_ens_pdaf, &
           screen, status_pdaf)
   END IF whichinit

   CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
        distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

   stat=status_pdaf
END SUBROUTINE init_pdaf

!-----------------------------------------------------------------------

SUBROUTINE assimilation_pdaf()

#if 0
   USE mod_parallel, &    ! Parallelization
        ONLY: mype_world, abort_parallel
   USE mod_assimilation, & ! airables for assimilation
        ONLY: filtertype
#endif

   IMPLICIT NONE

!KB
!KB   integer :: mype_world, filtertype
!KB   integer :: mype_world=rank_filter_comm
!KB

!
#if 0
   EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector from model fields
        init_dim_obs_pdaf, &            ! Initialize dimension of observation vector
        obs_op_pdaf, &                  ! Implementation of the Observation operator
        init_obs_pdaf                   ! Routine to provide vector of measurements
! ! Subroutine used in SEIK and ETKF
   EXTERNAL :: prepoststep_ens_pdaf, &  ! User supplied pre/poststep routine
        init_obsvar_pdaf                ! Initialize mean observation error variance
! ! Subroutines used in EnKF
   EXTERNAL :: add_obs_error_pdaf, &    ! Add obs. error covariance R to HPH in EnKF
        init_obscovar_pdaf              ! Initialize obs error covar R in EnKF
! ! Subroutine used in SEIK and ETKF
   EXTERNAL :: prodRinvA_pdaf           ! Provide product R^-1 A for some matrix A
! ! Subroutines used in LSEIK and LETKF
   EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
        init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
        init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
        g2l_state_pdaf, &               ! Get state on local ana. domain from global state
        l2g_state_pdaf, &               ! Init global state from state on local analysis domain
        g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
        init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
        prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some matrix A (for LSEIK)
        init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
        init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
        obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
        init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain
#endif

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state (possible, but not required!)
! Calls: PDAF_put_state_seik
! Calls: PDAF_put_state_enkf
! Calls: PDAF_put_state_lseik
! Calls: PDAF_put_state_etkf
! Calls: PDAF_put_state_letkf
! Calls: MPI_barrier (MPI)
!EOP

! local variables
   INTEGER :: status    ! Status flag for filter routines


! ************************
! *** Perform analysis ***
! ************************

! *** Note on PDAF_get_state for offline implementation: ***
! *** For the offline mode of PDAF the call to           ***
! *** PDAF_get_state is not required as no forecasting   ***
! *** is performed in this mode. However, it is save     ***
! *** to call PDAF_get_state, even it is not necessary.  ***
! *** The functionality of PDAF_get_state is deactived   ***
! *** for the offline mode.                              ***

   call PDAF_force_analysis() ! Suggested by Lars

   IF (filtertype == 1) THEN
      CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
           init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
   ELSE IF (filtertype == 2) THEN
      CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
           init_obs_pdaf, prepoststep_ens_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
           status)
   ELSE IF (filtertype == 3) THEN
      CALL PDAF_put_state_lseik( &
           collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
           init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
           prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
           init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
           g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
   ELSE IF (filtertype == 4) THEN
      CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
           init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
   ELSE IF (filtertype == 5) THEN
      CALL PDAF_put_state_letkf( &
           collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
           init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
           prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
           init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
           g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
   ELSE IF (filtertype == 6) THEN
      CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
           init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
   ELSE IF (filtertype == 7) THEN
      CALL PDAF_put_state_lestkf( &
           collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
           init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
           prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
           init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
           g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
   END IF
   if (verbosity >= debug) write(stderr,*) 'PDAF PUT STATUS= ',status

   CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
        distribute_state_pdaf, prepoststep_ens_pdaf, status)
   if (verbosity >= debug) write(stderr,*) 'PDAF GET STATUS= ',status
END SUBROUTINE assimilation_pdaf

!-----------------------------------------------------------------------

! Below a number of callback routines implemented by the user.
! At present it is unclear how many - and which - are needed for the EAT implementation.

#define _CLASSICAL_OFFLINE_SERIAL_
#ifdef _CLASSICAL_OFFLINE_SERIAL_
SUBROUTINE init_ens_pdaf(filtertype,dim_p,dim_ens,state_p,Uinv,ens_p,flag)
! !ARGUMENTS:
   INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
   INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
   INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
   REAL(REAL64), INTENT(inout) :: state_p(dim_p)          ! PE-local model state
   ! It is not necessary to initialize the array 'state_p' for SEIK.
   ! It is available here only for convenience and can be used freely.
   REAL(REAL64), INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
   REAL(REAL64), INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
   INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

   if (verbosity >= debug) write(stderr,*) 'init_ens_pdaf() ',filtertype,dim_p,dim_ens
   ! copy all_states(:,:) to ens_p
!KB   ens_p=model_states
END SUBROUTINE init_ens_pdaf

! Routine to collect a state vector from model fields
SUBROUTINE collect_state_pdaf(dim_p, state_p)
   INTEGER, INTENT(in) :: dim_p
   INTEGER, INTENT(inout) :: state_p(dim_p)
   ! can be empty as all states have already been collected
   if (verbosity >= debug) write(stderr,*) 'collect_state_pdaf()'
END SUBROUTINE collect_state_pdaf

! Routine to distribute a state vector to model fields
SUBROUTINE distribute_state_pdaf(dim_p, state_p)
   INTEGER, INTENT(in) :: dim_p
   INTEGER, INTENT(inout) :: state_p(dim_p)
   ! can be empty as all states have already been distrubuted
   if (verbosity >= debug) write(stderr,*) 'distribute_state_pdaf()'
END SUBROUTINE distribute_state_pdaf

! Initialize dimension of observation vector
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)
   INTEGER, INTENT(in)  :: step       ! Current time step
   INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector
   dim_obs_p = nobs
   if (verbosity >= debug) write(stderr,*) 'init_dim_obs_pdaf() ',step,dim_obs_p
END SUBROUTINE init_dim_obs_pdaf

! Implementation of the Observation operator
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)
   INTEGER, INTENT(in) :: step               ! Currrent time step
   INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
   INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
   REAL(REAL64), INTENT(in)    :: state_p(dim_p)     ! PE-local model state
   REAL(REAL64), INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state

   integer :: i

   if (verbosity >= debug) write(stderr,*) 'obs_op_pdaf() ',dim_p, dim_obs_p
   DO i = 1, dim_obs_p
      m_state_p(i) = state_p(iobs(i))
   END DO
write(0,*) 'AAAA',m_state_p
END SUBROUTINE obs_op_pdaf

! Routine to provide vector of measurements
SUBROUTINE init_obs_pdaf(step, dim_obs_f, observation_f)
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
  REAL(REAL64), INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector
  if (verbosity >= debug) write(stderr,*) 'init_obs_pdaf() ',dim_obs_f,size(observation_f)
  observation_f = obs(1:dim_obs_f)
write(0,*) 'BBBB',observation_f
END SUBROUTINE init_obs_pdaf


! ! Subroutine used in SEIK and ETKF
! User supplied pre/poststep routine
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p,state_p, Uinv, ens_p, flag)
   INTEGER, INTENT(in) :: step        ! Current time step (not relevant for offline mode)
   INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
   INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
   INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
   INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
   REAL(REAL64), INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
   ! The array 'state_p' is not generally not initialized in the case of SEIK.
   ! It can be used freely here.
   REAL(REAL64), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
   REAL(REAL64), INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
   INTEGER, INTENT(in) :: flag        ! PDAF status flag

! *** local variables ***
  INTEGER :: i, j, member             ! Counters
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                ! estimated RMS error
  REAL, ALLOCATABLE :: variance(:)    ! model state variances
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)

   if (verbosity >= debug) write(stderr,*) 'prepoststep_ens_pdaf()',dim_p, dim_ens, dim_ens_p, dim_obs_p

! **********************
! *** INITIALIZATION ***
! **********************

  IF (firsttime) THEN
     WRITE (*, '(8x, a)') 'Analyze initial state ensemble'
     anastr = 'ini'
  ELSE
     IF (step<0) THEN
        WRITE (*, '(8x, a)') 'Analyze and write forecasted state ensemble'
        anastr = 'for'
     ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
        anastr = 'ana'
     END IF
  END IF

  ! Allocate fields
  ALLOCATE(variance(dim_p))

  ! Initialize numbers
  rmserror_est  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! ****************************************************************
! *** Perform prepoststep. The state and error information is  ***
! *** completely in the ensemble.                              ***
! ****************************************************************

  ! *** Compute mean state
  WRITE (*, '(8x, a)') '--- compute ensemble mean'

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

  write (*,*) anastr, 'state ', state_p(1:10)


  ! *** Compute sampled variances ***
  variance(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        variance(j) = variance(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  variance(:) = invdim_ensm1 * variance(:)

   block
   real(real64) :: s=0._real64,ss=0._real64
   real(real64) :: mean,var,std
   real(real64) :: x
   s=0._real64
   ss=0._real64
   DO i = 1, dim_obs_p
      x = obs(i)-state_p(iobs(i))
      s = s+x
      ss = ss+x*x
   END DO
   mean=s/dim_obs_p
   var=(ss-s*s/dim_obs_p)/(dim_obs_p-1)
   std=sqrt(var)
   if (verbosity >= info) then
      write(stderr,*) 'prepoststep_ens_pdaf() - mean, var, std',steps,mean,var,std
   end if
   end block


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
  DO i = 1, dim_p
     rmserror_est = rmserror_est + variance(i)
  ENDDO
  rmserror_est = SQRT(rmserror_est / dim_p)

  DEALLOCATE(variance)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(12x, a, es12.4)') &
       'RMS error according to sampled variance: ', rmserror_est

#if 0
   block
   integer :: u
   integer :: i=1,n
   character(len=32) :: fn
   if (step > 0) then
      write(fn,'(A,I0.4,A)') 'mean_',i,'.dat'
      open(newunit=u,file=trim(fn))
      write(u,*) 'step# ',i
      do n=1,dim_p
         write(u,'(F9.6)') sum(model_states(n,:))/dim_ens
      end do
      close(u)
      i=i+1
   end if
   end block
#endif

END SUBROUTINE prepoststep_ens_pdaf

! Initialize mean observation error variance
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)
   INTEGER, INTENT(in) :: step          ! Current time step
   INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
   REAL(REAL64), INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
   REAL(REAL64), INTENT(out)   :: meanvar       ! Mean observation error variance
   if (verbosity >= debug) write(stderr,*) 'init_obsvar_pdaf() ',rms_obs
   meanvar = rms_obs ** 2
END SUBROUTINE init_obsvar_pdaf

SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)
   INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
   INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
   INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
   REAL(REAL64), INTENT(out)    :: time     ! Current model (physical) time

   if (verbosity >= debug) write(stderr,*) 'next_observation_pdaf() '
   nsteps=1
   doexit=1
END SUBROUTINE next_observation_pdaf

! ! Subroutines used in EnKF
! Add obs. error covariance R to HPH in EnKF
SUBROUTINE add_obs_error_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-2,ierr)
END SUBROUTINE add_obs_error_pdaf

! Initialize obs error covar R in EnKF
SUBROUTINE init_obscovar_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_obscovar_pdaf


! ! Subroutine used in SEIK and ETKF
! Provide product R^-1 A for some matrix A
SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p)
   INTEGER, INTENT(in) :: step                ! Current time step
   INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
   INTEGER, INTENT(in) :: rank                ! Rank of initial covariance matrix
   REAL(REAL64), INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
   REAL(REAL64), INTENT(in)    :: A_p(dim_obs_p,rank) ! Input matrix from SEEK_ANALYSIS
   REAL(REAL64), INTENT(out)   :: C_p(dim_obs_p,rank) ! Output matrix

   integer :: i,j
   REAL(REAL64) :: ivariance_obs

   if (verbosity >= info) write(stderr,*) 'prodRinvA_pdaf() ',dim_obs_p,rank,rms_obs

   ivariance_obs = 1.0 / rms_obs ** 2
   DO j = 1, rank
      DO i = 1, dim_obs_p
         C_p(i, j) = ivariance_obs * A_p(i, j)
      END DO
   END DO

!   if (verbosity >= info) write(stderr,*) 'prodRinvA_pdaf() A_P',A_p

END SUBROUTINE prodRinvA_pdaf


! ! Subroutines used in LSEIK and LETKF
! Provide number of local analysis domains
SUBROUTINE init_n_domains_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_n_domains_pdaf

! Initialize state dimension for local ana. domain
SUBROUTINE init_dim_l_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_dim_l_pdaf

! Initialize dim. of obs. vector for local ana. domain
SUBROUTINE init_dim_obs_l_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_dim_obs_l_pdaf

! Get state on local ana. domain from global state
SUBROUTINE g2l_state_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE g2l_state_pdaf

! Init global state from state on local analysis domain
SUBROUTINE g2l_obs_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE g2l_obs_pdaf

! Restrict a global obs. vector to local analysis domain
SUBROUTINE l2g_state_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE l2g_state_pdaf

! Provide vector of measurements for local ana. domain
SUBROUTINE init_obs_l_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_obs_l_pdaf

! Provide product R^-1 A for some matrix A (for LSEIK)
SUBROUTINE prodRinvA_l_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE prodRinvA_l_pdaf

! Initialize local mean observation error variance
SUBROUTINE init_obsvar_l_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_obsvar_l_pdaf

! Provide full vector of measurements for PE-local domain
SUBROUTINE init_obs_f_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_obs_f_pdaf

! Obs. operator for full obs. vector for PE-local domain
SUBROUTINE obs_op_f_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE obs_op_f_pdaf

! Get dimension of full obs. vector for PE-local domain
SUBROUTINE init_dim_obs_f_pdaf()
call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE init_dim_obs_f_pdaf
#endif

end program eat_filter_pdaf
