! Copyright (C) 2021 Bolding & Bruggeman

#undef _USE_PDAF_
#define _USE_PDAF_

module pdaf_wrapper

   !! A wrapper around 'PDAF'

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use eat_config, only: info, debug
   use PDAF_mod_filter
   use PDAF_interfaces_module
   use iso_c_binding, only: c_ptr, c_null_ptr, c_loc, c_int, c_double

   IMPLICIT NONE

   private

   public init_pdaf, assimilation_pdaf, finish_pdaf, iobs, obs, pcvt_callback

   integer :: stderr=error_unit,stdout=output_unit
   integer :: verbosity=info

   integer, pointer, contiguous :: iobs(:) => null()
      !! map observations to state
   real(real64), pointer, contiguous :: obs(:) => null()
      !! observation vector

   integer :: filtertype=6
      !! PDAF filter selection
      ! filtertype = 6 ! Type of filter
                       !  (1) SEIK
                       !  (2) EnKF
                       !  (3) LSEIK
                       !  (4) ETKF
                       !  (5) LETKF
                       !  (6) ESTKF
                       !  (7) LESTKF
                       ! (13) 3D-var
   integer :: subtype = 0
      !! PDAF subtype selection
      !! PDAF 3dvar subtype selection
                       ! (0) Parameterized 3D-Var
                       ! (1) Ensemble 3D-Var using LETKF for ensemble transformation
                       ! (4) Ensemble 3D-Var using global ETKF for ensemble transformation
                       ! (6) Hybrid 3D-Var using LETKF for ensemble transformation
                       ! (7) Hybrid 3D-Var using global ETKF for ensemble transformation
   REAL(real64) :: timenow
   integer :: doexit,steps
   real(real64) :: rms_obs = 0.05
     !! Observation error standard deviation !KB shall be replaced

   interface
      subroutine cvt_callback_interface(cb_type, iter, dim_p, dim_cvec, dim_cvec_ens, ens_p, v_p, Vv_p) bind(c)
         import C_PTR, C_INT
         INTEGER(C_INT), VALUE, INTENT(IN) :: cb_type
         INTEGER(C_INT), VALUE, INTENT(IN) :: iter
           !!< Iteration of optimization
         INTEGER(C_INT), VALUE, INTENT(IN) :: dim_p
           !!< PE-local observation dimension
         INTEGER(C_INT), VALUE, INTENT(IN) :: dim_cvec
           !!< Dimension of control vector
         INTEGER(C_INT), VALUE, INTENT(IN) :: dim_cvec_ens
           !!< Dimension of control vector
         TYPE(C_PTR), VALUE, INTENT(IN) :: ens_p
           !!< PE-local ensemble
         TYPE(C_PTR), VALUE, INTENT(IN) :: v_p
           !!< PE-local control vector
         TYPE(C_PTR), VALUE, INTENT(IN) :: Vv_p
           !!< PE-local result vector
      end subroutine
   end interface
   procedure (cvt_callback_interface), pointer, save :: pcvt_callback => null()
      !! user implemented routine (in Python)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

! Below are the routines implemented to link to the PDAF library - pdaf-d.

SUBROUTINE init_pdaf(EAT_COMM_filter, state_size, ensemble_size, model_states, stat)

      !! Initialize various variable and call PDAF_init()

   integer, intent(in), value :: EAT_COMM_filter, state_size, ensemble_size
   real(real64), pointer, contiguous :: model_states(:,:)
   integer, intent(out) :: stat

   INTEGER :: filter_param_i(7)
     !! Integer parameter array for filter
   REAL(REAL64) :: filter_param_r(2)
     !! Real parameter array for filter
   INTEGER :: status_pdaf
     !! PDAF status flag

   integer :: screen = 2
     !! Write screen output (1) for output, (2) add timings
   integer :: dim_state_p = -1
   integer :: comm_couple, comm_filter, comm_model, n_modeltasks, task_id
   logical :: filterpe=.true.

   ! Filter specific variables
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
   type_sqrt = 0    ! Type of transform matrix square-root
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

   ! File handling
   logical :: fileexists
   character(len=128) :: filename = 'output.dat'
   character(len=*), parameter :: nmlfile='eat_pdaf.nml'
   integer :: nmlunit

   namelist /nml_config_pdaf/ screen, filtertype, subtype, &
                              type_trans, type_forget, forget, &
                              dim_cvec, dim_cvec_ens, beta_3dvar, &
                              type_opt

   ! values must be provided via namelist
   dim_cvec=-1
   dim_cvec_ens=ensemble_size
   type_opt=0

   INQUIRE(FILE=nmlfile, EXIST=fileexists)
   if (fileexists) then
      open(newunit=nmlunit,file=nmlfile,status='old',action='read')
      read(nmlunit,nml=nml_config_pdaf)
      close(nmlunit)
   end if

   ! initialize variables for call to PDAF_init()
   dim_ens = ensemble_size
   dim_state_p=state_size
   comm_couple=EAT_COMM_filter ! suggested by Lars
   comm_filter=EAT_COMM_filter
   comm_model=EAT_COMM_filter ! suggested by Lars
   filterpe=.true.
   task_id=1
   n_modeltasks=1

   call PDAF_set_comm_pdaf(EAT_COMM_filter)
   select case (filtertype)
      case (2)
         ! *** EnKF with Monte Carlo init ***
         filter_param_i(1) = dim_state_p ! State dimension
         filter_param_i(2) = dim_ens     ! Size of ensemble
         filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
         filter_param_i(4) = incremental ! Whether to perform incremental analysis
         filter_param_i(5) = 0           ! Smoother lag (not implemented here)
         filter_param_r(1) = forget      ! Forgetting factor
         CALL PDAF_init(filtertype, subtype, 0, &
              filter_param_i, 6,&
              filter_param_r, 2, &
              COMM_model, COMM_filter, COMM_couple, &
              task_id, n_modeltasks, filterpe, init_ens_pdaf, &
              screen, status_pdaf)
      case (1,3,4,5,6,7)
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
         CALL PDAF_init(filtertype, subtype, 0, &
              filter_param_i, 7,&
              filter_param_r, 2, &
              COMM_model, COMM_filter, COMM_couple, &
              task_id, n_modeltasks, filterpe, init_ens_pdaf, &
              screen, status_pdaf)
      case (13)
            filter_param_i(1) = dim_state_p    ! State dimension
            filter_param_i(2) = dim_ens        ! Size of ensemble
            filter_param_i(3) = type_opt       ! Choose type of optimizer
            filter_param_i(4) = dim_cvec       ! Dimension of control vector (parameterized part)
            filter_param_i(5) = dim_cvec_ens   ! Dimension of control vector (ensemble part)
            filter_param_r(1) = forget         ! Forgetting factor
            filter_param_r(2) = beta_3dvar     ! Hybrid weight for hybrid 3D-Var
            select case (subtype)
               case (0) ! parameterized 3D-Var
                  if (dim_cvec < 0) then
                     call abort('init_pdaf(): dim_cvec < 0')
                  end if
                  CALL PDAF_init(filtertype, subtype, 0, &
                       filter_param_i, 5,&
                       filter_param_r, 1, &
                       COMM_model, COMM_filter, COMM_couple, &
                       task_id, n_modeltasks, filterpe, init_3dvar_pdaf, &
                       screen, status_pdaf)
               case (4,7) ! Ensemble or hybrid 3D-Var - local variations not supported for now
                  if (dim_cvec_ens < 0) then
                     call abort('init_pdaf(): dim_cvec_ens < 0')
                  end if
                  CALL PDAF_init(filtertype, subtype, 0, &
                       filter_param_i, 5,&
                       filter_param_r, 2, &
                       COMM_model, COMM_filter, COMM_couple, &
                       task_id, n_modeltasks, filterpe, init_ens_pdaf, &
                       screen, status_pdaf)
            end select
      case default
         stop 'init_pdaf(): Non-valid filtertype'
    end select
    if (status_pdaf /= 0) stop 'init_pdaf(): status_pdaf /= 0'

   CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
        distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

   stat=status_pdaf

   call PDAF_set_ens_pointer(model_states, stat)
END SUBROUTINE init_pdaf

!-----------------------------------------------------------------------

SUBROUTINE assimilation_pdaf() bind(c)

  !! local variables
  INTEGER :: status_pdaf    ! Status flag for filter routines

   call PDAF_force_analysis() ! Suggested by Lars

   select case (filtertype)
      case (1)
         CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
              init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status_pdaf)
      case (2)
         CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
              init_obs_pdaf, prepoststep_ens_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
              status_pdaf)
      case (3)
         CALL PDAF_put_state_lseik( &
              collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
              init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
              prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
              init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
              g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status_pdaf)
      case (4)
         CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
              init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status_pdaf)
      case (5)
         CALL PDAF_put_state_letkf( &
              collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
              init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
              prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
              init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
              g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status_pdaf)
      case (6)
         CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
              init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status_pdaf)
      case (7)
         CALL PDAF_put_state_lestkf( &
              collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
              init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
              prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
              init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
              g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status_pdaf)
      case (13)
         ! From .../PDAF_V1.16_var_devel/models/lorenz96/assimilation_pdaf.F90
! tkdiff ./tutorial/3dvar/online_2D_serialmodel/prepoststep_3dvar_pdaf.F90 ./tutorial/3dvar/online_2D_serialmodel/prepoststep_ens_pdaf.F90
         select case (subtype)
            case (0)
               CALL PDAF_put_state_3dvar(collect_state_pdaf, &
                    init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prodRinvA_pdaf, &
                    cvt_pdaf, cvt_adj_pdaf, obs_op_pdaf, obs_op_adj_pdaf, &
                    prepoststep_3dvar_pdaf, status_pdaf)
            case (1)
#if 0
               ! Ensemble 3D-Var with local ESTKF update of ensemble perturbations
               CALL PDAF_put_state_en3dvar_lestkf(collect_state_pdaf, &
                    init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prodRinvA_pdaf, &
                    cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
                    init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, &
                    prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                    init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                    g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, &
                    prepoststep_pdaf, status_pdaf)
#endif
            case (4)
               ! Ensemble 3D-Var with global ESTKF update of ensemble perturbations
               CALL PDAF_put_state_en3dvar_estkf(collect_state_pdaf, &
                    init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prodRinvA_pdaf, &
                    cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_pdaf, obs_op_adj_pdaf, &
                    init_obsvar_pdaf, &
                    prepoststep_ens_pdaf, status_pdaf)
            case (6)
#if 0
               ! Hybrid 3D-Var with local ESTKF update of ensmeble perturbations
               CALL PDAF_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
                    init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prodRinvA_pdaf, &
                    cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
                    obs_op_lin_pdaf, obs_op_adj_pdaf, &
                    init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, &
                    prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                    init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                    g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, &
                    prepoststep_pdaf, status_pdaf)
#endif
            case (7)
               ! Hybrid 3D-Var with global ESTKF update of ensemble perturbations
               CALL PDAF_put_state_hyb3dvar_estkf(collect_state_pdaf, &
                    init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prodRinvA_pdaf, &
                    cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
                    obs_op_pdaf, obs_op_adj_pdaf, init_obsvar_pdaf, &
                    prepoststep_ens_pdaf, status_pdaf)
         end select
   end select
   if (status_pdaf /= 0) stop 'assimilation_pdaf(): status_pdaf /= 0'
   if (verbosity >= debug) write(stderr,*) 'PDAF PUT STATUS= ',status_pdaf

   CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
        distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)
   if (verbosity >= debug) write(stderr,*) 'PDAF GET STATUS= ',status_pdaf
END SUBROUTINE assimilation_pdaf

!-----------------------------------------------------------------------

subroutine finish_pdaf() bind(c)

   !! Cleanup and finalize the EAT/PDAF component

!-----------------------------------------------------------------------
#ifdef _USE_PDAF__
   CALL finalize_pdaf(0) ! Basically CALL PDAF_deallocate()
#endif
end subroutine finish_pdaf

!-----------------------------------------------------------------------

! Below a number of callback routines implemented by the user.
! At present it is unclear how many - and which - are needed for the EAT implementation.

#define _CLASSICAL_OFFLINE_SERIAL_
#ifdef _CLASSICAL_OFFLINE_SERIAL_
SUBROUTINE init_ens_pdaf(filtertype,dim_p,dim_ens,state_p,Uinv,ens_p,flag)
   INTEGER, INTENT(in) :: filtertype
     !! Type of filter to initialize
   INTEGER, INTENT(in) :: dim_p
     !! PE-local state dimension
   INTEGER, INTENT(in) :: dim_ens
     !! Size of ensemble
   REAL(REAL64), INTENT(inout) :: state_p(dim_p)
     !! PE-local model state
     ! It is not necessary to initialize the array 'state_p' for SEIK.
     ! It is available here only for convenience and can be used freely.
   REAL(REAL64), INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1)
     !! Array not referenced for SEIK
   REAL(REAL64), INTENT(out)   :: ens_p(dim_p, dim_ens)
     !! PE-local state ensemble
   INTEGER, INTENT(inout) :: flag
     !! PDAF status flag

   if (verbosity >= debug) write(stderr,*) 'init_ens_pdaf() ',filtertype,dim_p,dim_ens
END SUBROUTINE init_ens_pdaf

!-----------------------------------------------------------------------

SUBROUTINE init_3dvar_pdaf()
   ! initialized via Python plugin
   if (verbosity >= debug) write(stderr,*) 'init_3dvar_pdaf()'
END SUBROUTINE init_3dvar_pdaf

!-----------------------------------------------------------------------

! Routine to collect a state vector from model fields
SUBROUTINE collect_state_pdaf(dim_p, state_p)
   INTEGER, INTENT(in) :: dim_p
   INTEGER, INTENT(inout) :: state_p(dim_p)

   ! can be empty as all states have already been collected
   if (verbosity >= debug) write(stderr,*) 'collect_state_pdaf()'
END SUBROUTINE collect_state_pdaf

!-----------------------------------------------------------------------

! Routine to distribute a state vector to model fields
SUBROUTINE distribute_state_pdaf(dim_p, state_p)
   INTEGER, INTENT(in) :: dim_p
   INTEGER, INTENT(inout) :: state_p(dim_p)

   ! can be empty as all states have already been distrubuted
   if (verbosity >= debug) write(stderr,*) 'distribute_state_pdaf()'
END SUBROUTINE distribute_state_pdaf

!-----------------------------------------------------------------------

! Initialize dimension of observation vector
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)
   INTEGER, INTENT(in)  :: step
     !! Current time step
   INTEGER, INTENT(out) :: dim_obs_p
     !! Dimension of observation vector

   if (verbosity >= debug) write(stderr,*) 'init_dim_obs_pdaf() ',step,dim_obs_p
   dim_obs_p = size(obs)
END SUBROUTINE init_dim_obs_pdaf

!-----------------------------------------------------------------------

! Implementation of the Observation operator
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)
   INTEGER, INTENT(in) :: step
     !! Currrent time step
   INTEGER, INTENT(in) :: dim_p
     !! PE-local dimension of state
   INTEGER, INTENT(in) :: dim_obs_p
     !! Dimension of observed state
   REAL(REAL64), INTENT(in) :: state_p(dim_p)
     !! PE-local model state
   REAL(REAL64), INTENT(out) :: m_state_p(dim_obs_p)
     !! PE-local observed state

   integer :: i

   if (verbosity >= debug) write(stderr,*) 'obs_op_pdaf() ',dim_p, dim_obs_p
   DO i = 1, dim_obs_p
      m_state_p(i) = state_p(iobs(i))
   END DO
END SUBROUTINE obs_op_pdaf

!-----------------------------------------------------------------------

! Routine to provide vector of measurements
SUBROUTINE init_obs_pdaf(step, dim_obs_f, observation_f)
   INTEGER, INTENT(in) :: step
     !! Current time step
   INTEGER, INTENT(in) :: dim_obs_f
     !! Dimension of full observation vector
   REAL(REAL64), INTENT(out) :: observation_f(dim_obs_f)
     !! Full observation vector

   if (verbosity >= debug) write(stderr,*) 'init_obs_pdaf() ',dim_obs_f,size(observation_f)
   observation_f = obs(1:dim_obs_f)
END SUBROUTINE init_obs_pdaf

!-----------------------------------------------------------------------

! ! Subroutine used in SEIK and ETKF
! User supplied pre/poststep routine
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p,state_p, Uinv, ens_p, flag)
   INTEGER, INTENT(in) :: step
     !! Current time step (not relevant for offline mode)
   INTEGER, INTENT(in) :: dim_p
     !! PE-local state dimension
   INTEGER, INTENT(in) :: dim_ens
     !! Size of state ensemble
   INTEGER, INTENT(in) :: dim_ens_p
     !! PE-local size of ensemble
   INTEGER, INTENT(in) :: dim_obs_p
     !! PE-local dimension of observation vector
   REAL(REAL64), INTENT(inout) :: state_p(dim_p)
     !! PE-local forecast/analysis state
     ! The array 'state_p' is not generally not initialized in the case of SEIK.
     ! It can be used freely here.
   REAL(REAL64), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1)
     !! Inverse of matrix U
   REAL(REAL64), INTENT(inout) :: ens_p(dim_p, dim_ens)
     !! PE-local state ensemble
   INTEGER, INTENT(in) :: flag
     !! PDAF status flag

   !! local variables
   INTEGER :: i, j, member
     !! Counters
   LOGICAL, SAVE :: firsttime = .TRUE.
     !! Routine is called for first time?
   REAL :: invdim_ens
     !! Inverse ensemble size
   REAL :: invdim_ensm1
     !! Inverse of ensemble size minus 1
   REAL :: rmserror_est
     !! estimated RMS error
   REAL, ALLOCATABLE :: variance(:)
     !! model state variances
   REAL, ALLOCATABLE :: field(:,:)
     !! global model field
   CHARACTER(len=2) :: ensstr
     !! String for ensemble member
   CHARACTER(len=2) :: stepstr
     !! String for time step
   CHARACTER(len=3) :: anastr
     !! String for call type (initial, forecast, analysis)

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
   if (verbosity >= debug) then
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

!-----------------------------------------------------------------------

SUBROUTINE prepoststep_3dvar_pdaf()
   if (verbosity >= debug) write(stderr,*) 'prepoststep_3dvar_pdaf()'
END SUBROUTINE prepoststep_3dvar_pdaf

!-----------------------------------------------------------------------

! Initialize mean observation error variance
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)
   INTEGER, INTENT(in) :: step
     !! Current time step
   INTEGER, INTENT(in) :: dim_obs_p
     !! PE-local dimension of observation vector
   REAL(REAL64), INTENT(in) :: obs_p(dim_obs_p)
     !! PE-local observation vector
   REAL(REAL64), INTENT(out) :: meanvar
     !! Mean observation error variance

   if (verbosity >= debug) write(stderr,*) 'init_obsvar_pdaf() ',rms_obs
   meanvar = rms_obs ** 2
END SUBROUTINE init_obsvar_pdaf

!-----------------------------------------------------------------------

SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)
   INTEGER, INTENT(in)  :: stepnow
     !! Number of the current time step
   INTEGER, INTENT(out) :: nsteps
     !! Number of time steps until next obs
   INTEGER, INTENT(out) :: doexit
     !! Whether to exit forecasting (1 for exit)
   REAL(REAL64), INTENT(out) :: time
     !! Current model (physical) time

   if (verbosity >= debug) write(stderr,*) 'next_observation_pdaf() '
   nsteps=1
   doexit=1
END SUBROUTINE next_observation_pdaf

!-----------------------------------------------------------------------

! ! Subroutines used in EnKF
! Add obs. error covariance R to HPH in EnKF
SUBROUTINE add_obs_error_pdaf()
   call abort('add_obs_error_pdaf')
END SUBROUTINE add_obs_error_pdaf

!-----------------------------------------------------------------------

! Initialize obs error covar R in EnKF
SUBROUTINE init_obscovar_pdaf()
   call abort('init_obscovar_pdaf')
END SUBROUTINE init_obscovar_pdaf

!-----------------------------------------------------------------------

! ! Subroutine used in SEIK and ETKF
! Provide product R^-1 A for some matrix A
SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p)
   INTEGER, INTENT(in) :: step
     !! Current time step
   INTEGER, INTENT(in) :: dim_obs_p
     !! PE-local dimension of obs. vector
   INTEGER, INTENT(in) :: rank
     !! Rank of initial covariance matrix
   REAL(REAL64), INTENT(in) :: obs_p(dim_obs_p)
     !! PE-local vector of observations
   REAL(REAL64), INTENT(in) :: A_p(dim_obs_p,rank)
     !! Input matrix from SEEK_ANALYSIS
   REAL(REAL64), INTENT(out) :: C_p(dim_obs_p,rank)
     !! Output matrix

   !! local variables
   integer :: i,j
   REAL(REAL64) :: ivariance_obs

   if (verbosity >= debug) write(stderr,*) 'prodRinvA_pdaf() ',dim_obs_p,rank,rms_obs
   ivariance_obs = 1.0 / rms_obs ** 2
   DO j = 1, rank
      DO i = 1, dim_obs_p
         C_p(i, j) = ivariance_obs * A_p(i, j)
      END DO
   END DO
END SUBROUTINE prodRinvA_pdaf

!-----------------------------------------------------------------------
! 3Dvar specific routines
!-----------------------------------------------------------------------

! ~PDAF_V2.0/tutorial/3dvar/online_2D_serialmodel/cvt_pdaf.F90
SUBROUTINE cvt_pdaf(iter, dim_p, dim_cvec, v_p, Vv_p)

   ! *** Arguments ***
   INTEGER, INTENT(in) :: iter
     !!< Iteration of optimization
   INTEGER, INTENT(in) :: dim_p
     !!< PE-local observation dimension
   INTEGER, INTENT(in) :: dim_cvec
     !!< Dimension of control vector
   REAL, TARGET, INTENT(in)    :: v_p(dim_cvec)
     !!< PE-local control vector
   REAL, TARGET, INTENT(inout) :: Vv_p(dim_p)
     !!< PE-local result vector

   IF (ASSOCIATED(pcvt_callback)) CALL pcvt_callback(1, iter, dim_p, -1, dim_cvec, &
       c_null_ptr, c_loc(v_p), c_loc(Vv_p))
END SUBROUTINE cvt_pdaf

!-----------------------------------------------------------------------

! ~PDAF_V2.0/tutorial/3dvar/online_2D_serialmodel/cvt_adj_pdaf.F90
SUBROUTINE cvt_adj_pdaf(iter, dim_p, dim_cvec, Vv_p, v_p)

!  *** Arguments ***
   INTEGER, INTENT(in) :: iter
     !!< Iteration of optimization
   INTEGER, INTENT(in) :: dim_p
     !!< PE-local observation dimension
   INTEGER, INTENT(in) :: dim_cvec
     !!< Dimension of control vector
   REAL, TARGET, INTENT(in)    :: Vv_p(dim_p)
     !!< PE-local input vector
   REAL, TARGET, INTENT(inout) :: v_p(dim_cvec)
     !!< PE-local result vector

   IF (ASSOCIATED(pcvt_callback)) CALL pcvt_callback(2, iter, dim_p, -1, dim_cvec, &
       c_null_ptr, c_loc(v_p), c_loc(Vv_p))
END SUBROUTINE cvt_adj_pdaf

!-----------------------------------------------------------------------

! ~PDAF_V2.0/tutorial/3dvar/online_2D_serialmodel/cvt_ens_pdaf.F90
SUBROUTINE cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, v_p, Vv_p)

   INTEGER, INTENT(in) :: iter
     !!< Iteration of optimization
   INTEGER, INTENT(in) :: dim_p
     !!< PE-local dimension of state
   INTEGER, INTENT(in) :: dim_ens
     !!< Ensemble size
   INTEGER, INTENT(in) :: dim_cvec_ens
     !!< Dimension of control vector
   REAL, TARGET, INTENT(in) :: ens_p(dim_p, dim_ens)
     !!< PE-local ensemble
   REAL, TARGET, INTENT(in) :: v_p(dim_cvec_ens)
     !!< PE-local control vector
   REAL, TARGET, INTENT(inout) :: Vv_p(dim_p)
     !!< PE-local state increment

   IF (ASSOCIATED(pcvt_callback)) CALL pcvt_callback(3, iter, dim_p, dim_ens, dim_cvec_ens, &
       c_loc(ens_p), c_loc(v_p), c_loc(Vv_p))
END SUBROUTINE cvt_ens_pdaf

!-----------------------------------------------------------------------

! ~PDAF_V2.0/tutorial/3dvar/online_2D_serialmodel/cvt_adj_ens_pdaf.F90
SUBROUTINE cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, Vv_p, v_p)

! *** Arguments ***
   INTEGER, INTENT(in) :: iter
     !!< Iteration of optimization
   INTEGER, INTENT(in) :: dim_p
     !!< PE-local dimension of state
   INTEGER, INTENT(in) :: dim_ens
     !!< Ensemble size
   INTEGER, INTENT(in) :: dim_cvec_ens
     !!< Number of columns in HV_p
   REAL, TARGET, INTENT(in) :: ens_p(dim_p, dim_ens)
     !!< PE-local ensemble
   REAL, TARGET, INTENT(in)    :: Vv_p(dim_p)
     !!< PE-local input vector
   REAL, TARGET, INTENT(inout) :: v_p(dim_cvec_ens)
     !!< PE-local result vector

   IF (ASSOCIATED(pcvt_callback)) CALL pcvt_callback(4, iter, dim_p, dim_ens, dim_cvec_ens, &
       c_loc(ens_p), c_loc(v_p), c_loc(Vv_p))
END SUBROUTINE cvt_adj_ens_pdaf

!-----------------------------------------------------------------------

! Implementation of the adjoint Observation operator
SUBROUTINE obs_op_adj_pdaf(step, dim_p, dim_obs_p, m_state_p, state_p)
   INTEGER, INTENT(in) :: step
     !! Currrent time step
   INTEGER, INTENT(in) :: dim_p
     !! PE-local dimension of state
   INTEGER, INTENT(in) :: dim_obs_p
     !! Dimension of observed state
   REAL(REAL64), INTENT(in) :: m_state_p(dim_obs_p)
     !! PE-local observed state
   REAL(REAL64), INTENT(out) :: state_p(dim_p)
     !! PE-local model state

   !! local variables
   integer :: i

   if (verbosity >= debug) write(stderr,*) 'obs_op_adj_pdaf() ',dim_p, dim_obs_p
   state_p = 0.0
   DO i = 1, dim_obs_p
      state_p(iobs(i)) = m_state_p(i)
   END DO
END SUBROUTINE obs_op_adj_pdaf

!-----------------------------------------------------------------------

SUBROUTINE abort(msg)
   character(len=*), intent(in) :: msg
   write (stderr, '(a)') msg
   stop 1
   ! call MPI_Abort(MPI_COMM_WORLD,-3,ierr)
END SUBROUTINE

#if 1
!KB - this routine needs content
SUBROUTINE obs_op_lin_pdaf()
   call abort('obs_op_lin_pdaf')
END SUBROUTINE obs_op_lin_pdaf

!KB - this routine needs content
SUBROUTINE prepoststep_pdaf()
   call abort('prepoststep_pdaf')
END SUBROUTINE prepoststep_pdaf

! ! Subroutines used in LSEIK and LETKF
! Provide number of local analysis domains
SUBROUTINE init_n_domains_pdaf()
   call abort('init_n_domains_pdaf')
END SUBROUTINE init_n_domains_pdaf

! Initialize state dimension for local ana. domain
SUBROUTINE init_dim_l_pdaf()
   call abort('init_dim_l_pdaf')
END SUBROUTINE init_dim_l_pdaf

! Initialize dim. of obs. vector for local ana. domain
SUBROUTINE init_dim_obs_l_pdaf()
   call abort('init_dim_obs_l_pdaf')
END SUBROUTINE init_dim_obs_l_pdaf

! Get state on local ana. domain from global state
SUBROUTINE g2l_state_pdaf()
   call abort('g2l_state_pdaf')
END SUBROUTINE g2l_state_pdaf

! Init global state from state on local analysis domain
SUBROUTINE g2l_obs_pdaf()
   call abort('g2l_obs_pdaf')
END SUBROUTINE g2l_obs_pdaf

! Restrict a global obs. vector to local analysis domain
SUBROUTINE l2g_state_pdaf()
   call abort('l2g_state_pdaf')
END SUBROUTINE l2g_state_pdaf

! Provide vector of measurements for local ana. domain
SUBROUTINE init_obs_l_pdaf()
   call abort('init_obs_l_pdaf')
END SUBROUTINE init_obs_l_pdaf

! Provide product R^-1 A for some matrix A (for LSEIK)
SUBROUTINE prodRinvA_l_pdaf()
   call abort('prodRinvA_l_pdaf')
END SUBROUTINE prodRinvA_l_pdaf

! Initialize local mean observation error variance
SUBROUTINE init_obsvar_l_pdaf()
   call abort('init_obsvar_l_pdaf')
END SUBROUTINE init_obsvar_l_pdaf

! Provide full vector of measurements for PE-local domain
SUBROUTINE init_obs_f_pdaf()
   call abort('init_obs_f_pdaf')
END SUBROUTINE init_obs_f_pdaf

! Obs. operator for full obs. vector for PE-local domain
SUBROUTINE obs_op_f_pdaf()
   call abort('obs_op_f_pdaf')
END SUBROUTINE obs_op_f_pdaf

! Get dimension of full obs. vector for PE-local domain
SUBROUTINE init_dim_obs_f_pdaf()
   call abort('init_dim_obs_f_pdaf')
END SUBROUTINE init_dim_obs_f_pdaf
#endif
#endif

end module pdaf_wrapper
