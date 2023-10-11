module gam_den
!============================================================
! Calculation of gamma-ray density for the isotropic 
! gamma-ray production rate, which depends on radius
! (spherically symmetric case). Mathematically:
!
! f(r0)=int d^3 r1 Q(r1)/(4\pi)*exp(-lam*s)/s^2, where
! s^2=(\vec r1-\vec r0)^2, r0 is observation point,
! Q(r1) is a total production rate at point r1,
! lam is absorbtion coefficient [1/cm].
! Integration is performed over the volume d^3 r1.
!
! After manipulation it can be written:
!
! f(r0)=int_{0}^{R} dr1 Q(r1)*Phi, where
! Phi=r1/(2*r)*(E1(lam*(r1-r0))-E1(lam*(r1+r0))) if r0/=0, and
! Phi=exp(-lam*r1), if r0=0.
! This form is used for calculations in current module.
!
! Here f is a flux, to get density: g=f/c, where c is
! a speed of light
!
! Usage
! =====
! set_gamden(x1,x2,grate), x1, x2 are min and max radius
! for considered volume, grate(x) is production rate
!
! gamden(labs,x,res), labs is absorbtion coefficient,
! x is the radius
!
! Example of usage is below
!
! There is a possibility to calculate density outside the
! production region, i.e. r0 could be outside
! [rmin,rmax] region. But one should check it.
!============================================================
 use gauss_kronrod, only : gk_adaptive
 use peak_integration, only : peak_int
 use phys_const, only : c_light
 implicit none
 private
 save
 
 
 public :: gamden, set_gamden
 
 real(8) :: r1, r0, lam, rmin, rmax
 real(8), parameter :: tiny_log=700d0
 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 procedure (fun_prot), pointer :: prod_rate=>null()
 
 
 
 contains



subroutine set_gamden(x1,x2,grate)
! Variables
 real(8), intent(in) :: x1, x2
 procedure (fun_prot) :: grate
! Calculations  
 prod_rate=>grate
 rmin=x1
 rmax=x2
end subroutine set_gamden



subroutine gamden(labs,x,res)
! Variables
 real(8), intent(in) :: x, labs
 real(8), intent(out) :: res
 real(8) :: tmin, tmax, tmid
 real(8), parameter :: d1=3.5d0, d2=3.5d0
 real(8), parameter :: zero=1d-300
 save
! Calculations 
 
 r0=x ! distance from the centre
 lam=labs
 if (abs(lam)<zero) lam=zero
 
 tmin=(rmin) ! limits of intergration
 tmid=(r0)   ! to bypass a logarithmic divergence
 tmax=(rmax)
 
 
 if ((tmin<=tmid).and.(tmid<=tmax)) then
  call peak_int(igamden,tmin,tmax,tmid,d1,d2,res)
 else
  call gk_adaptive(igamden,tmin,tmax,res)
 end if 
 res=res/c_light
end subroutine gamden



 
subroutine igamden(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: qrate, ang
! Calculations
 r1=x
 call prod_rate(r1,qrate)
 call ang_part(ang)
 res=ang*qrate
end subroutine igamden
 
 
subroutine ang_part(res)
!=====================================================
! Analytical representation of the integration over
! all angles with absorbtion, which leads to
! exponential integral E1(x) (Slatec library is used)
!
! r0, r1, lam is input parameters (global variables)
!=====================================================
! Variables
 real(8), intent(out) :: res
 real(8), parameter :: rzero=1d7
 real(8) :: x1, x2, rel
 interface
! DE1(x) from Slatec library
  real(8) function de1(x)
  real(8) :: x
  end function de1
 end interface
! Calculations

! If arguments are large, e.g. lam is large 
 x2=lam*(r1+r0) 
 if (x2>tiny_log) x2=tiny_log
 x1=lam*abs(r1-r0)
 
 if (x1>x2) then
  res=0d0
  return
 end if 
 
! Take into account that r0 can be zero  
 if (abs(r0)>0d0) then
  rel=r1/r0
! If r0 becomes small, it is better to use approximation (which is more precise in this case)  
  if (rel>rzero) then
   res=exp(-lam*r1)
  else
   res=rel*(de1(x1)-de1(x2))/2
  end if
 
 else
  res=exp(-lam*r1)
 end if 
 
 
end subroutine ang_part


end module gam_den




!program main
! use gam_den, only : gamden, set_gamden
! implicit none
!! Variables
! integer :: i, nr, il
! integer, parameter :: nlam=5
! real(8) :: rmin, rmax, lambda
! real(8) :: r1, r2, r, res, res1
! real(8) :: lam(nlam)
! character(100) :: fname(nlam)
! 
! interface
!  subroutine g_rate(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine
! end interface
!
!! Calculations 
! rmin=0d0
! rmax=1d3
! call set_gamden(rmin,rmax,g_rate)
! 
! lam=[1d-8,1d0,1d1,1d2,1d3]
! fname(1)='den1d-5.dat'
! fname(2)='den1d0.dat'
! fname(3)='den1d1.dat'
! fname(4)='den1d2.dat'
! fname(5)='den1d3.dat'
! 
! do il=1,nlam
! lambda=lam(il)
! 
! r1=1d-2
! r2=1d3
! nr=100
! 
! open(1,file=fname(il))
! do i=0,nr
!  r=r1*(r2/r1)**(i*1d0/nr)
!  call gamden(lambda,r,res)
!  call g_rate(r,res1)  
!  write(1,*) r, res, res1
! end do
! close(1)
! 
! end do
!end program main
!
!
!subroutine g_rate(x,res)
!! Variables
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! real(8) :: a
!! Calculations  
! a=x
! res=1/(1d-2+x**2)
!
!end subroutine g_rate