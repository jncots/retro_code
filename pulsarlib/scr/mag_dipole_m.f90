module mag_dipole_m

 implicit none
 private
 save
 
 public :: mag_dipole
 
 real(8), parameter :: pi=3.14159265358979324d0 ! pi number
 real(8), parameter :: c_light=2.99792458d10
 
 
 type mag_dipole
!============================================================================== 
! bm_star - magnetic field at the pole in Gauss (B_*=1d12 G by default)
! r_star - radius of the star in cm (R_*=1d6 cm by default)
! r_unit is normalization scale, so that x'=x/r_unit, where x and r_unit in cm
! r_unit=1 cm by default
! period=33 ms by default
! omega - angular frequency in rad/s
! rsu - star radius in r_unit
!==============================================================================
  real(8) :: bm_star=1d12, r_star=1d6, r_unit=1d0, period=33d-3
  real(8) :: omega, rsu
  real(8), private :: mus 
 contains
  procedure :: set, calc
 end type mag_dipole
 
 
 
 contains

subroutine set(this,unt)
!================================
! Setting to calculate in units
! unt='lc' of light cylinder
! unt='cm' of cantimeters
! unt='rs' of star radius
!================================
 class(mag_dipole) :: this
 character(2), optional, intent(in)  :: unt

 this%omega=2*pi/this%period

 if (present(unt)) then
  select case (unt)
   case('lc')
    this%r_unit=c_light/this%omega 
   case('cm')
    this%r_unit=1d0
   case('rs')
    this%r_unit=this%r_star
  case default
   write(*,*) 'mag_dipole: units set to unt=lc'
   this%r_unit=c_light/this%omega
  end select
 end if

 this%rsu=this%r_star/this%r_unit
 this%mus=(this%bm_star/2)*(this%rsu)**3

end subroutine set
 


subroutine calc(this,x,bm)
! Calculation of the dipole magnetic field in rotating reference frame 
 class(mag_dipole) :: this
 real(8), intent(in) ::  x(3)
 real(8), intent(out) :: bm(3)
 real(8) :: rho2, rho, mu, z3

 rho2=dot_product(x,x)
 rho=sqrt(rho2)
 mu=this%mus/rho**5
 z3=3*x(3)*mu ! 3*z*mu
 
 bm=x*z3
 bm(3)=bm(3)-rho2*mu
 
end subroutine calc

end module mag_dipole_m


!program main
! use mag_dipole_m, only : mag_dipole
! implicit none
! type(mag_dipole) :: md
! real(8) :: x(3), bm(3)
! 

! x=[0.1d0,0.5d0,1d0]

! call md%set('lc')
! call md%calc(x,bm)
! write(*,'(3Es14.6)') bm
! write(*,*) md%rsu, md%omega


!end program main
