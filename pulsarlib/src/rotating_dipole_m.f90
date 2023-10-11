module rotating_dipole_m
!==============================================================================
! Example of usage is in the end after the module
!
! Calculation of the magnetic field of the rotating dipole
!
! If magnetic moment at angle xi to rotation axis z, which rotates with
! angular frequency w, then
! Bx=Bx'*cos(wt)*cos(xi)-By'*sin(wt)+Bz'*cos(wt)*sin(xi)
! By=Bx'*sin(wt)*cos(xi)+By'*cos(wt)+Bz'*sin(wt)*sin(xi)
! Bz=-Bx'*sin(xi)+Bz'*cos(xi),
! where Bx', By', Bz' are components in the rotating frame
! Coordinats in the rotating frame are
! x'=x*cos(wt)*cos(xi)+y*sin(wt)*cos(xi)-z*sin(xi)
! y'=-x*sin(wt)+y*cos(wt)
! z'=x*cos(wt)*sin(xi)+y*sin(wt)*sin(xi)+z*cos(xi) 
!
! Calculation of dipole magnetic field:
! B'=B_{*}/2*(R_{*}/R0)^3*(3*x'*z' ,3*y'*z', 3*z'^2-rho^2)/rho^5,
! where 
! B_{*} is magnetic field at the pole, in Gauss
! R_{*} is radius of the star, in cm
! R0 is normalization scale, so that x=x0/R0, where x0 and R0 in cm
! rho=sqrt(x^2+y^2+z^2)
! x, y, z are position coordinates in R0 units
! x, y, z and Bx', By', Bz' are in the system with z along dipole direction
!==============================================================================
 implicit none
 private
 save
 
 public :: rotating_dipole
 
 real(8), parameter :: pi=3.14159265358979324d0 ! pi number
 real(8), parameter :: c_light=2.99792458d10
 real(8), parameter :: degrad=pi/180
 
 
 type rotating_dipole
!============================================================================== 
! bm_star - magnetic field at the pole in Gauss (B_*=1d12 G by default)
! r_star - radius of the star in cm (R_*=1d6 cm by default)
! r_unit is normalization scale, so that x'=x/r_unit, where x and r_unit in cm
! r_unit=1 cm by default
! dipole_angle - angle between magnetic moment and rotation axis in degrees
! (dipole_angle=90 degrees by default)
! omega - angular frequency in rad/s (omega=2*pi/33d-3 for P=33 ms by default)
!==============================================================================
  real(8) :: bm_star=1d12, r_star=1d6, r_unit=1d0
  real(8) :: dipole_angle=90, omega=2*pi/33d-3
  real(8), private :: mus, sxi, cxi, rm(3,3) 
 contains
  procedure :: set_dipole, set_time, calc, set_lc_units
 end type rotating_dipole
 
 
 
 contains

subroutine set_lc_units(this)
! Setting to calculate in units
! of light cylinder
 class(rotating_dipole) :: this

 this%r_unit=c_light/this%omega

end subroutine set_lc_units
 
subroutine set_dipole(this)
! Precalculation of the normalization and sxi, cxi
 class(rotating_dipole) :: this
 real(8) :: dang
 
 this%mus=(this%bm_star/2)*(this%r_star/this%r_unit)**3
 dang=degrad*this%dipole_angle
 this%sxi=sin(dang)
 this%cxi=cos(dang)
 
end subroutine set_dipole


subroutine set_time(this,time)
! Calculation of rotation matrix for
! particular 'time'
 class(rotating_dipole) :: this
 real(8), intent(in) :: time
 real(8) :: phase, cph, sph, cxi, sxi
! Setting of rotation matrix rm:
 phase=this%omega*time
 cph=cos(phase)
 sph=sin(phase)
 cxi=this%cxi
 sxi=this%sxi
 
 this%rm(:,1)=[cph*cxi,sph*cxi,-sxi]
 this%rm(:,2)=[-sph,cph,0d0]
 this%rm(:,3)=[cph*sxi,sph*sxi,cxi]

end subroutine set_time
 

subroutine calc(this,x,bm)
 class(rotating_dipole) :: this
 real(8), intent(in) ::  x(3)
 real(8), intent(out) :: bm(3)
 real(8) :: x0(3), bm0(3)
 real(8) :: rho2, rho, mu, z3
  
! Change coordinates to rotating reference frame:  
 x0(1)=dot_product(this%rm(:,1),x)
 x0(2)=dot_product(this%rm(:,2),x)
 x0(3)=dot_product(this%rm(:,3),x)
 
! Calculation of the dipole magnetic field in rotating reference frame: 
 rho2=dot_product(x0,x0)
 rho=sqrt(rho2)
 mu=this%mus/rho**5
 z3=3*x0(3)*mu ! 3*z*mu
 bm0=x0*z3
 bm0(3)=bm0(3)-rho2*mu
! Transformation of the magnetic field to the observer's reference frame:
 
 bm(1)=dot_product(this%rm(1,:),bm0)
 bm(2)=dot_product(this%rm(2,:),bm0)
 bm(3)=dot_product(this%rm(3,:),bm0)
 
end subroutine calc

end module rotating_dipole_m


! program main
!  use rotating_dipole_m, only : rotating_dipole
!  implicit none
!  type(rotating_dipole) :: rd
!  integer :: i
!  real(8) :: x(3), bm(3), t 
! 
!  t=0d0
!  x=[3d0,2d0,4d0]
!  call rd%set_lc_units
!  write(*,'(10Es14.6)') rd%r_unit
!!  rd%r_unit=rd%r_star
!  call rd%set_dipole
! 
!  call rd%set_time(t)
!  call rd%calc(x,bm)
! 
!  write(*,'(10Es14.6)') (bm(i),i=1,3)
! 

! end program main

