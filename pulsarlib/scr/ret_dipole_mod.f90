module ret_dipole_mod

 implicit none
 private
 save


 public :: ret_dipole

 real(8), parameter :: pi=3.14159265358979324d0 ! pi number
 real(8), parameter :: c_light=2.99792458d10
 real(8), parameter :: pi_180=pi/180d0


 type ret_dipole
!============================================================================== 
! bm_star - magnetic field at the pole in Gauss (B_*=1d12 G by default)
! r_star - radius of the star in cm (R_*=1d6 cm by default)
! r_unit is normalization scale, so that x'=x/r_unit, where x and r_unit in cm
! r_unit=1 cm by default
! period=33 ms by default
! omega - angular frequency in rad/s
! rsu - star radius in r_unit
!==============================================================================
  real(8), private :: sp, cp, sa, ca
  logical, private :: seta=.true.
  real(8), private :: a0(3), b0(3), c0(3)
  real(8), private :: mu
  real(8) :: bm_star=1d12, r_star=1d6, r_unit=1d0, period=33d-3
  real(8) :: omega, rsu
 contains
  procedure :: phase, inclin, calc, set

 end type ret_dipole


 contains


subroutine set(this,unt)
 class(ret_dipole) :: this
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
 else
   this%r_unit=c_light/this%omega
 end if


 this%rsu=this%r_star/this%r_unit
 this%mu=(this%bm_star/2)*(this%rsu)**3
end subroutine set



subroutine phase(this,ph,deg)
 class(ret_dipole) :: this
 real(8), intent(in) :: ph
 character(3), intent(in), optional :: deg
 character(3) :: deg0

 if (present(deg)) then
  deg0=deg
 else
  deg0='rad'
 end if

 call sin_cos(ph,deg0,this%sp,this%cp)
 this%seta=.true.

end subroutine phase

subroutine inclin(this,alpha,deg)
 class(ret_dipole) :: this
 real(8), intent(in) :: alpha
 character(3), intent(in), optional :: deg
 character(3) :: deg0

 if (present(deg)) then
  deg0=deg
 else
  deg0='rad'
 end if

 call sin_cos(alpha,deg0,this%sa,this%ca)
 this%seta=.true.

end subroutine inclin


subroutine sin_cos(alpha,deg,sr,cr)
 real(8), intent(in) :: alpha
 character(3), intent(in) :: deg
 real(8), intent(out) :: sr, cr
 real(8) :: alpha0

 if (deg=='deg') then
  alpha0=alpha*pi_180
 else
  alpha0=alpha
 end if
 
 sr=sin(alpha0)
 cr=cos(alpha0)

end subroutine sin_cos


subroutine calc_abc(this)
 class(ret_dipole) :: this
 real(8) :: s0, c0

 s0=this%sa*this%sp
 c0=this%sa*this%cp

 this%a0=[c0,s0,this%ca]
 this%b0=[-s0,c0,0d0]
 this%c0=[-c0,-s0,0d0]
 this%seta=.false.
end subroutine calc_abc


subroutine calc(this,rv,bm)
 class(ret_dipole) :: this
 real(8), intent(in) :: rv(3)
 real(8), intent(out) :: bm(3)
 real(8) :: r2, r0, nr(3), a(3), b(3), c(3)

 if (this%seta) call calc_abc(this)
 
 r2=dot_product(rv,rv)
 r0=sqrt(r2)
 nr=rv/r0
 

 a=this%a0
 b=r0*this%b0
 c=r2*this%c0


 bm=3*(a+b)+c
 bm=dot_product(bm,nr)*nr

 bm=bm-(a+b+c)
 bm=bm*(this%mu/r0**3)
 

end subroutine calc


end module ret_dipole_mod


!program main
! use vector_ops_mod, only : vector_ops
! use ret_dipole_mod, only : ret_dipole
! implicit none
! type(ret_dipole) :: rdip
! type(vector_ops) :: vp
! real(8) :: rv(3), bm(3)

! call rdip%set
! call rdip%phase(30d0,'deg')
! call rdip%inclin(90d0,'deg')


! rv=vp%unit_vec(20d0,30d0)*0.1d0

! call rdip%calc(rv,bm)

! write(*,*) bm
! 

!end program main
