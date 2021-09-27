! Copyright (C) 2021 Bolding & Bruggeman

module fields_2d

   !! An re-implementation of the 'model' from here:
   !! http://pdaf.awi.de/files/pdaf_tutorial_onlineserial.pdf
   !! Used for testing the analysis step

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   IMPLICIT NONE

   character(len=256), parameter :: time_format='%Y-%m-%d %H:%M:%S'
   integer, parameter :: nx=36,ny=18
   real(real64) :: field(nx,ny),obsfield(nx,ny)
   integer :: nobsmin=10,nobsmax=40
   real(real64) :: obsstd=0.05
   real(real64), parameter :: pi=acos(-1._real64)

contains

!-----------------------------------------------------------------------

subroutine true_field(member,nmember)
   integer, intent(in) :: member
   integer, intent(in), optional :: nmember
   integer :: i,j
   real(real64) :: phi=0._real64

   if (present(nmember)) then
      phi=0.5_real64*(1._real64*member+5._real64)/nmember
   end if
   do j=1,ny
      do i=1,nx
         field(i,j)=sin(2*pi*(1._real64*i/nx+1._real64*j/ny+phi))
         !KB sin(2*pi*(i/18+j/36)+2*0.5*pi*(k+5)/dim_ens)
      end do
   end do
end subroutine true_field

subroutine update_field(nsteps)
   integer, intent(in), optional :: nsteps
   integer :: i,j,k,n=2
   real(real64) :: tmp(nx)

   if (present(nsteps)) n=nsteps

   ! *** Time step: Shift field vertically ***
   do k=1,n
      tmp=field(:,ny)
      do j = ny,2,-1
         field(:,j)=field(:,j-1)
      end do
      field(:,1)=tmp
   end do
end subroutine update_field

subroutine get_obs(n,obsstd,iobs,obs)
   integer, intent(in) :: n
   real(real64), intent(in) :: obsstd
   integer, intent(inout) :: iobs(:)
   real(real64), intent(inout) :: obs(:)

   real(real64), allocatable :: obserr(:)
   real(real64) :: x(2)
   integer :: i,j,k,nobs
   integer :: iseed(4)=(/1000,2034,0,3/)

   nobs=size(obs)
   if (n == 1) then
      call true_field(0)
   end if
   call update_field(nsteps=2)
   obsfield=-999._real64
   !KB PDAF-D_V1.13.2/tutorial/offline_2D_serial/init_dim_obs_pdaf.F90
   do k=1,nobs
      call random_number(x)
      i=1+FLOOR(nx*x(1))
      j=1+FLOOR(ny*x(2))
      obsfield(i,j)=field(i,j)
      iobs(k)=i+(j-1)*nx
      obs(k)=field(i,j)
   end do
   allocate(obserr(nobs))
   call dlarnv(2,iseed,nobs,obserr)
#if 0
   do k=1,nobs
      write(0,*) 'BBB ',obs(k),obsstd*obserr(k)
   end do
#endif
   obs=obs+obsstd*obserr
end subroutine get_obs

subroutine write_field(fn,f)
   character(len=*), intent(in) :: fn
   real(real64), intent(in) :: f(:,:)
   integer :: u
   integer :: j
   open(newunit=u,file=trim(fn))
   do j=ny,1,-1
      write(u,'(*(f16.6))') f(:,j)
   end do
   close(u)
end subroutine write_field

!-----------------------------------------------------------------------

end module fields_2d
