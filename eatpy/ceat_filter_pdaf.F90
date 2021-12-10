! Copyright (C) 2021 Bolding & Bruggeman

module ceat_filter_pdaf

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use iso_c_binding
   use pdaf_wrapper

   implicit none

   private

contains

   subroutine ceat_pdaf_set_observations(nobs_, iobs_, obs_) bind(c)
      integer, value, intent(in) :: nobs_
      integer,      target, intent(inout) :: iobs_(nobs_)
      real(REAL64), target, intent(inout) :: obs_(nobs_)
      iobs => iobs_
      obs => obs_
   end subroutine

   subroutine ceat_pdaf_init(EAT_COMM_filter, state_size, ensemble_size, p, stat) bind(c)
      integer,     intent(in), value :: EAT_COMM_filter, state_size, ensemble_size
      type(c_ptr), intent(out)       :: p
      integer,     intent(out)       :: stat

      real(real64), pointer, contiguous :: model_states(:,:)
      call init_pdaf(EAT_COMM_filter, state_size, ensemble_size, model_states, stat)
      p = c_loc(model_states)
   end subroutine

end module
