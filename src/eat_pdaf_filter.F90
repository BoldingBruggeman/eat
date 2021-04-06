! Copyright (C) 2021 Bolding & Bruggeman

program eat_pdaf_filter

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

   ! PDAF variables
   integer :: filtertype=6
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
   integer :: ierr
!-----------------------------------------------------------------------
   call init_eat_config(color_filter+verbosity)
   CALL init_pdaf(ierr)
   if (ierr /= 0) then
      call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
   else
      write(error_unit,*) 'filter(PDAF is initialized)'
   end if

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

   !! Get observations and states and do the PDAF/assimilation step

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

      ! Begin PDAF specific part
      ! from .../tutorial/classical/offline_2D_serial/main_offline.F90
      CALL assimilation_pdaf()
      ! End PDAF specific part
   end do
end subroutine eat_do_pdaf

!-----------------------------------------------------------------------

subroutine eat_finish_pdaf()

   !! Cleanup and finalize the EAT/PDAF component

!-----------------------------------------------------------------------
#if 0
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

!KB Missing - comm_couple, comm_filter, comm_model, filterpe, n_modeltasks, task_id
!KB all known from eat_config - with different names
 integer :: comm_couple, comm_filter, comm_model, filterpe, n_modeltasks, task_id
!KB

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

  ! External subroutines
!KB  EXTERNAL :: init_ens_offline  ! Ensemble initialization

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
  dim_ens = 9       ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  integer :: &
  subtype = 5       ! (5) Offline mode
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
  real :: &
  rms_obs = 0.5    ! Observation error standard deviation
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

  logical :: fileexists
  character(len=*), parameter :: nmlfile='eat_pdaf_config.nml'
  integer :: nmlunit

!KB read some configuration from namelist
   namelist /eat_pdaf/ screen,filtertype

   INQUIRE(FILE=nmlfile, EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile,status='old',action='read')
      read(nmlunit,nml=eat_pdaf)
      close(nmlunit)
   end if
   write(*,*) screen,filtertype
!KB

  whichinit: IF (filtertype == 2) THEN
     ! *** EnKF with Monte Carlo init ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = 0           ! Smoother lag (not implemented here)
     filter_param_r(1) = forget      ! Forgetting factor
#if 1
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
#endif
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
#if 1
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
#endif
  END IF whichinit

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
  EXTERNAL :: prepoststep_ens_offline, & ! User supplied pre/poststep routine
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

  IF (filtertype == 1) THEN
     CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 2) THEN
     CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, add_obs_error_pdaf, init_obscovar_pdaf, &
          status)
  ELSE IF (filtertype == 3) THEN
     CALL PDAF_put_state_lseik( &
          collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_offline, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
  ELSE IF (filtertype == 4) THEN
     CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 5) THEN
     CALL PDAF_put_state_letkf( &
          collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_offline, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
  ELSE IF (filtertype == 6) THEN
     CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 7) THEN
     CALL PDAF_put_state_lestkf( &
          collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_offline, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
  END IF
END SUBROUTINE assimilation_pdaf

!-----------------------------------------------------------------------

! Below a number of callback routines implemented by the user.
! At present it is unclear how many - and which - are needed for the EAT implementation.


SUBROUTINE init_ens_offline(dim_p, state_p)
integer, intent(in) :: dim_p
integer, intent(inout) :: state_p(dim_p)
END SUBROUTINE init_ens_offline

! Routine to collect a state vector from model fields
SUBROUTINE collect_state_pdaf(dim_p, state_p)
integer, intent(in) :: dim_p
integer, intent(inout) :: state_p(dim_p)
END SUBROUTINE collect_state_pdaf

! Initialize dimension of observation vector
SUBROUTINE init_dim_obs_pdaf()
END SUBROUTINE init_dim_obs_pdaf

! Implementation of the Observation operator
SUBROUTINE obs_op_pdaf()
END SUBROUTINE obs_op_pdaf

! Routine to provide vector of measurements
SUBROUTINE init_obs_pdaf()
END SUBROUTINE init_obs_pdaf


! ! Subroutine used in SEIK and ETKF
! User supplied pre/poststep routine
SUBROUTINE prepoststep_ens_offline()
END SUBROUTINE prepoststep_ens_offline

! Initialize mean observation error variance
SUBROUTINE init_obsvar_pdaf()
END SUBROUTINE init_obsvar_pdaf


! ! Subroutines used in EnKF
! Add obs. error covariance R to HPH in EnKF
SUBROUTINE add_obs_error_pdaf()
END SUBROUTINE add_obs_error_pdaf

! Initialize obs error covar R in EnKF
SUBROUTINE init_obscovar_pdaf()
END SUBROUTINE init_obscovar_pdaf


! ! Subroutine used in SEIK and ETKF
! Provide product R^-1 A for some matrix A
SUBROUTINE prodRinvA_pdaf()
END SUBROUTINE prodRinvA_pdaf


! ! Subroutines used in LSEIK and LETKF
! Provide number of local analysis domains
SUBROUTINE init_n_domains_pdaf()
END SUBROUTINE init_n_domains_pdaf

! Initialize state dimension for local ana. domain
SUBROUTINE init_dim_l_pdaf()
END SUBROUTINE init_dim_l_pdaf

! Initialize dim. of obs. vector for local ana. domain
SUBROUTINE init_dim_obs_l_pdaf()
END SUBROUTINE init_dim_obs_l_pdaf

! Get state on local ana. domain from global state
SUBROUTINE g2l_state_pdaf()
END SUBROUTINE g2l_state_pdaf

! Init global state from state on local analysis domain
SUBROUTINE g2l_obs_pdaf()
END SUBROUTINE g2l_obs_pdaf

! Restrict a global obs. vector to local analysis domain
SUBROUTINE l2g_state_pdaf()
END SUBROUTINE l2g_state_pdaf

! Provide vector of measurements for local ana. domain
SUBROUTINE init_obs_l_pdaf()
END SUBROUTINE init_obs_l_pdaf

! Provide product R^-1 A for some matrix A (for LSEIK)
SUBROUTINE prodRinvA_l_pdaf()
END SUBROUTINE prodRinvA_l_pdaf

! Initialize local mean observation error variance
SUBROUTINE init_obsvar_l_pdaf()
END SUBROUTINE init_obsvar_l_pdaf

! Provide full vector of measurements for PE-local domain
SUBROUTINE init_obs_f_pdaf()
END SUBROUTINE init_obs_f_pdaf

! Obs. operator for full obs. vector for PE-local domain
SUBROUTINE obs_op_f_pdaf()
END SUBROUTINE obs_op_f_pdaf

! Get dimension of full obs. vector for PE-local domain
SUBROUTINE init_dim_obs_f_pdaf()
END SUBROUTINE init_dim_obs_f_pdaf

end program eat_pdaf_filter
