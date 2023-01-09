! Copyright (C) 2021 Bolding & Bruggeman

module ceat_filter_pdaf

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use iso_c_binding
   use pdaf_wrapper

   implicit none

   private

contains

   subroutine ceat_pdaf_set_observations(nobs_, iobs_, obs_, rms_obs_) bind(c)
      integer, value, intent(in) :: nobs_
      integer,      target, intent(inout) :: iobs_(nobs_)
      real(REAL64), target, intent(inout) :: obs_(nobs_)
      real(REAL64), target, intent(inout) :: rms_obs_(nobs_)
      iobs => iobs_
      obs => obs_
      rms_obs => rms_obs_
   end subroutine

   subroutine ceat_set_3dvar_callback(cb_3dvar) bind(c)
      type(c_funptr), intent(in), value :: cb_3dvar
      call c_f_procpointer(cb_3dvar, pcvt_callback)
   end subroutine

   subroutine ceat_pdaf_get_ensemble_state_pointer(p) bind(c)
      type(c_ptr), intent(out) :: p

      real(real64), pointer, contiguous :: model_states(:,:)

      call get_ensemble_state_pointer_pdaf(model_states)
      p = c_loc(model_states)
   end subroutine

end module
