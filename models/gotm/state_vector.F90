!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save the model results in file
!
! !INTERFACE:
   subroutine state_to_fields()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                  :: n,len
   type (type_field_set)                    :: field_set
   class (type_field_set_member),pointer    :: member
!EOP
!-------------------------------------------------------------------------
!BOC
   field_set = fm%get_state()
   member => field_set%first
   n = 1
   len = 0
   do while (associated(member))
      ! This field is part of the model state. Its name is member%field%name.
      if (associated(member%field%data_0d)) then
         ! Depth-independent variable with data pointed to by child%field%data_0d
         member%field%data_0d = state_vector(n)
         n = n + 1
      elseif (associated(member%field%data_1d)) then
         ! Depth-dependent variable with data pointed to by member%field%data_1d
         len = size(member%field%data_1d)
         member%field%data_1d = state_vector(n:n+len)
         n = n + len
      else
         stop 'no data assigned to state field'
      end if
      member => member%next
   end do

   return
   end subroutine state_to_fields
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Data from fields to state vector
!
! !INTERFACE:
   subroutine fields_to_state()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                  :: n,len
   type (type_field_set)                    :: field_set
   class (type_field_set_member),pointer    :: member
!EOP
!-------------------------------------------------------------------------
!BOC
   field_set = fm%get_state()
   member => field_set%first
   n = 1
   len = 0
   do while (associated(member))
      ! This field is part of the model state. Its name is member%field%name.
      if (associated(member%field%data_0d)) then
         ! Depth-independent variable with data pointed to by member%field%data_0d
         state_vector(n) = member%field%data_0d
         n = n+1
      elseif (associated(member%field%data_1d)) then
         ! Depth-dependent variable with data pointed to by member%field%data_1d
         len = size(member%field%data_1d)
         state_vector(n:n+len) = member%field%data_1d
         n = n+len
      else
         stop 'no data assigned to state field'
      end if
      member => member%next
   end do

   return
   end subroutine fields_to_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine state_vector_size(state_size)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: state_size
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   type (type_field_set)                    :: field_set
   class (type_field_set_member),pointer    :: member
!EOP
!-------------------------------------------------------------------------
!BOC
   field_set = fm%get_state()
   member => field_set%first
   state_size = 0
   do while (associated(member))
      if (associated(member%field%data_0d)) then
         state_size = state_size+1
      elseif (associated(member%field%data_1d)) then
         ! Depth-dependent variable with data pointed to by member%field%data_1d
         state_size = state_size+size(member%field%data_1d)
      else
         stop 'no data assigned to state field'
      end if
      member => member%next
   end do
   return
   end subroutine state_vector_size
!EOC


!-----------------------------------------------------------------------
