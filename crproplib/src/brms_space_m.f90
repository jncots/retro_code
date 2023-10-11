module brms_space_m
 use phys_const, only : c_light, pi, psec

 implicit none
 private
 save

 public :: brms_space !, main_calc

 type brms_space
  real(8) :: b0, z0, iz0, bcent
 contains
  procedure :: set, bm_exp, bm_box
 end type brms_space



 real(8), parameter :: kpsec=1d3*psec
 real(8), parameter :: rsun=8.5*kpsec, rbulge=3*kpsec
 real(8), parameter :: irsun=1/rsun, rgal=20*kpsec

 contains



subroutine set(this,b0,z0) 
 class(brms_space) :: this
 real(8), intent(in) :: b0, z0

 this%b0=b0
 this%z0=z0
 this%iz0=1/z0
 this%bcent=this%b0*exp(5.5d0/8.5d0)
 

end subroutine set



function bm_exp(this,x)
 class(brms_space) :: this
 real(8), intent(in) :: x(3)
 real(8) :: bm_exp
 real(8) :: r, z
 real(8) :: br


 r=sqrt(x(1)**2+x(2)**2)
 z=abs(x(3))

 br=this%bcent
 if (r>rbulge) br=br*exp(-(r-rbulge)*irsun)

 bm_exp=br*exp(-z*this%iz0)


end function bm_exp

function bm_box(this,x)
 class(brms_space) :: this
 real(8), intent(in) :: x(3)
 real(8) :: bm_box
 real(8) :: r, z

 r=sqrt(x(1)**2+x(2)**2)
 z=abs(x(3))

 if (z>this%z0) then
  bm_box=0d0
 else
  if (r>rgal) then
   bm_box=0d0
  else
   bm_box=this%b0
  end if
 end if

end function bm_box




!subroutine main_calc
! type(brms_space) :: bt
! real(8) :: xc(3)

! xc=[1d0,1d0,3d0]*kpsec

! call bt%set(1d-6,2*kpsec)
! write(*,*) bt%bm_box(xc)
! 
!end subroutine main_calc



end module brms_space_m


!program main
! use brms_space_m, only : main_calc

! call main_calc

!end program main
