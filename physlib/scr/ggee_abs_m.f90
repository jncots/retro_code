module ggee_abs_m
!============================================================
! Integration the total probability of gamma-gamma
! absorbtion with density of soft photon field
! var%wabs(eg,res), eg in eV, res in [1/s]
!============================================================
 use cont_fun1d_m, only : cont_fun1d
 use gauss_kronrod, only : gk_adaptive
 
 implicit none
 private
 save
 
 
 public :: ggee_abs
 
 
 abstract interface
  subroutine fun1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type ggee_abs
  type(cont_fun1d) :: sph ! soft photon field
 contains 
  procedure :: wabs, set, del
 end type ggee_abs
 
 real(8), parameter :: me=0.510998928d6 ! electron rest mass in eV
 real(8), parameter :: me2=me**2! 2.611199044d11
 real(8), parameter :: ime2=1/me2 
 real(8), parameter :: re2c=2.380587754d-15 ! re2c=re^2*c, re is classical radius of electron
 real(8), parameter :: anorm=6.28d0*re2c
 real(8) :: egam, egami
 procedure(fun1d), pointer :: nden=>null()
 

 contains
 
subroutine set(this,fx,x1,x2)
 class(ggee_abs) :: this
 procedure(fun1d) :: fx
 real(8), intent(in) :: x1, x2
 call this%sph%set(fx,x1,x2)
end subroutine set

subroutine del(this)
 class(ggee_abs) :: this
 call this%sph%del
end subroutine del 
 
 
subroutine wabs(this,eg,res)
 class(ggee_abs) :: this
 real(8), intent(in) :: eg
 real(8), intent(out) :: res
 real(8) :: tmin, tmax, emin, emax, es
 real(8) :: eph1, eph2 ! min and max energy of soft photons in eV 
 
 eph1=this%sph%x1
 eph2=this%sph%x2
 nden=>this%sph%fx
 
 egam=eg
 es=me2/egam
 emin=eph1
 if (es>eph1) emin=es
 emax=eph2
 
 if (emin>=emax) then
  res=0d0
  return
 end if 
 
 egami=egam*ime2
 tmin=log(emin)
 tmax=log(emax)
 call gk_adaptive(iwgg,tmin,tmax,res)
 res=res*anorm
end subroutine wabs


subroutine iwgg(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eph, ww, fph 
! Calculations
 eph=exp(x)
 call w_abs(egam,eph,ww)
 call nden(eph,fph)
 res=fph*ww*eph
end subroutine iwgg

 
 
subroutine w_abs(egam,eph,res)
! Variables
 real(8), intent(in) :: egam, eph ! in eV
 real(8), intent(out) :: res ! in cm^3/s
 real(8) :: z, lz
! Calculations
 z=eph*egami
 if (z<=1d0) then
  res=0d0
  return
 end if
 lz=log(z)
 res=(z+0.41d0*lz)*log(0.541d0*z+1.406d0)
 res=res/(z-0.47d0*lz)
 res=res*((z-1)/z)**(3d0/2d0)
 res=res/z          ! * anorm
end subroutine w_abs 
 

end module ggee_abs_m
