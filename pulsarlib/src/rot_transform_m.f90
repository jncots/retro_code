module rot_transform_m
 implicit none
 private
 save
 
 public :: rot_transform
 
 real(8), parameter :: pi=3.14159265358979324d0 ! pi number
 real(8), parameter :: degrad=pi/180


 type rot_transform
  real(8) :: theta, phi ! in degrees
  real(8), private :: sxi, cxi, cph, sph, rm(3,3) 
 contains
  procedure :: set_theta, set_phi, calc_rm
  procedure :: to_rot, from_rot
 end type rot_transform

 contains

subroutine set_theta(this,theta)
 class(rot_transform) :: this
 real(8), intent(in) :: theta
 real(8) :: dtheta 

 this%theta=theta
 dtheta=degrad*theta
 this%sxi=sin(dtheta)
 this%cxi=cos(dtheta)
end subroutine set_theta

subroutine set_phi(this,phi)
 class(rot_transform) :: this
 real(8), intent(in) :: phi
 real(8) :: dphi  

 this%phi=phi
 dphi=degrad*phi
 this%cph=cos(dphi)
 this%sph=sin(dphi)
end subroutine set_phi

subroutine calc_rm(this)
 class(rot_transform) :: this
 this%rm(:,1)=[this%cph*this%cxi,this%sph*this%cxi,-this%sxi]
 this%rm(:,2)=[-this%sph,this%cph,0d0]
 this%rm(:,3)=[this%cph*this%sxi,this%sph*this%sxi,this%cxi]
end subroutine calc_rm


subroutine to_rot(this,x,x0)
 class(rot_transform) :: this
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: x0(3)
! Change coordinates to rotating reference frame:  
 x0(1)=dot_product(this%rm(:,1),x)
 x0(2)=dot_product(this%rm(:,2),x)
 x0(3)=dot_product(this%rm(:,3),x)
end subroutine to_rot


subroutine from_rot(this,bm0,bm)
 class(rot_transform) :: this
 real(8), intent(in) :: bm0(3)
 real(8), intent(out) :: bm(3)
! Transformation of the magnetic field to the observer's reference frame:
 
 bm(1)=dot_product(this%rm(1,:),bm0)
 bm(2)=dot_product(this%rm(2,:),bm0)
 bm(3)=dot_product(this%rm(3,:),bm0)
end subroutine from_rot

end module rot_transform_m



!program main
! use rot_transform_m, only : rot_transform
! type(rot_transform) :: rt
! real(8) :: x(3), x0(3), x1(3)

! call rt%set_theta(30d0)
! call rt%set_phi(120d0)
! call rt%calc_rm


! x=[1d0,2d0,3d0]

! call rt%to_rot(x,x0)
! call rt%from_rot(x0,x1)
! write(*,*) 'x=', x
! write(*,*) 'x0=', x0
! write(*,*) 'x1=', x1
! 


!end program main
