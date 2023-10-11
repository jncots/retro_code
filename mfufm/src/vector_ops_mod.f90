module vector_ops_mod
!========================================
! Contains shortcut functions for
! 3d vector operations
!========================================
 use phys_const, only : pi
 implicit none
 private
 save

 public :: vector_ops

 
 type vector_ops
 contains
  procedure, nopass :: cross_prod, unit_vec
  procedure, nopass :: vec2, vec_norm, vec_dir
  procedure, nopass :: ab_cos, ab_sin, ab_deg
 end type vector_ops

 real(8), parameter :: pi_deg=pi/180, deg_pi=180/pi


 contains


function cross_prod(a,b)
!=================================
! Cross product (axb) of vectors 
! a and b
!=================================
 real(8), intent(in) :: a(3), b(3)
 real(8) :: cross_prod(3)

 cross_prod(1)=a(2)*b(3)-a(3)*b(2)
 cross_prod(2)=a(3)*b(1)-a(1)*b(3)
 cross_prod(3)=a(1)*b(2)-a(2)*b(1)

end function cross_prod


function unit_vec(phi,theta)
!=====================================
! Return unit vector in the direction
! with spherical angles phi and theta
! given in degrees
!=====================================
 real(8), intent(in) :: phi, theta
 real(8) :: unit_vec(3)
 real(8) :: phi0, theta0, sth

 phi0=phi*pi_deg
 theta0=theta*pi_deg
 sth=sin(theta0)
  
 unit_vec(1)=cos(phi0)*sth
 unit_vec(2)=sin(phi0)*sth
 unit_vec(3)=cos(theta0)

end function unit_vec


function vec2(a)
 real(8), intent(in) :: a(3)
 real(8) :: vec2
 vec2=dot_product(a,a)
end function vec2

function vec_norm(a)
 real(8), intent(in) :: a(3)
 real(8) :: vec_norm
 vec_norm=sqrt(dot_product(a,a))
end function vec_norm


function vec_dir(a)
 real(8), intent(in) :: a(3)
 real(8) :: vec_dir(3)
 real(8) :: an 

 an=sqrt(dot_product(a,a))
 vec_dir=a/an

end function vec_dir

function ab_cos(a,b)
 real(8), intent(in) :: a(3), b(3)
 real(8) :: ab_cos
 real(8) :: ab, a2, b2
 
 ab=dot_product(a,b)
 a2=dot_product(a,a)
 b2=dot_product(b,b)

 ab_cos=ab/sqrt(a2*b2)

end function ab_cos


function ab_sin(a,b)
!==================================
! Alternative to sin=sqrt(1-cos^2)
!==================================
 real(8), intent(in) :: a(3), b(3)
 real(8) :: ab_sin
 real(8) :: axb(3), axb2, a2, b2

 axb=cross_prod(a,b)
 axb2=dot_product(axb,axb)
 a2=dot_product(a,a)
 b2=dot_product(b,b)
 ab_sin=sqrt(axb2/(a2*b2))
 
end function ab_sin


function ab_deg(a,b)
!==================================
! Angle between a and b in degrees
!==================================
 real(8), intent(in) :: a(3), b(3)
 real(8) :: ab_deg

 ab_deg=acos(ab_cos(a,b))*deg_pi

end function ab_deg

end module vector_ops_mod



! Test program
!program main
! use vector_ops_mod, only : vector_ops
! type(vector_ops) :: vp
! real(8), dimension(3) :: a, b, c, d

! a=vp%unit_vec(0d0,90d0)
! b=vp%unit_vec(0d0,30d0)
! 
! write(*,'(A,10Es30.17)') 'a=', a, 1-vp%vec_norm(vp%vec_dir(a*123d0))
! write(*,'(A,10Es30.17)') 'a2=', 1-vp%vec_norm(a)
! write(*,'(A,10Es30.17)') 'b=', b, vp%vec_norm(b)
! write(*,'(A,10Es14.6)') 'axb, bxa=', vp%cross_prod(a,b), vp%cross_prod(b,a)
! write(*,'(A,10Es14.6)') 'norm=', vp%vec_norm(a*2723), vp%vec_dir(a*123d0)
! write(*,'(A,10Es30.17)') 'sin=', vp%ab_sin(a,b), sqrt(1-(vp%ab_cos(a,b))**2)
! write(*,'(A,10Es30.17)') 'cos=', vp%ab_cos(a,b), vp%ab_cos(b,a), vp%ab_deg(a,b)


!end program main
