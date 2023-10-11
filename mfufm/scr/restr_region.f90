module restr_region
!======================================================
! Contains class 'r0region' with methods 'in', and
! 'shoot', which checks the conditions of getting of 
! a point in the cubic volume determined 
! by xb(2,3) array and crossing of that volume by the
! line drawn from a point in particular direction, 
! respectively. The details are in the description to 
! each subroutine.
!
! The example is after the module
!======================================================
 implicit none
 private
 save
 
 
 public :: r0region
 
 type r0region
! xb(1,1)=x_left, xb(2,1)=x_right, xb(1,2)=y_left, xb(2,2)=y_right
! xb(1,3)=z_left, xb(2,3)=z_right
  real(8) :: xb(2,3) 
 contains
  procedure :: inside=>check, shoot, mindl
 end type r0region
 
 
 contains


function mindl(this)
 class(r0region), intent(in) :: this
 real(8) :: mindl, mindl1

 mindl=this%xb(2,1)-this%xb(1,1)
 mindl1=this%xb(2,2)-this%xb(1,2)

 if (mindl1<mindl) mindl=mindl1
 mindl1=this%xb(2,3)-this%xb(1,3)
 
 if (mindl1<mindl) mindl=mindl1

end function mindl


function check(this,x)
!===========================================================
! Check whether point x(3) lies inside volume with borders
! xb. If it is , then check=.true., otherwise check=.false.
!===========================================================
 class(r0region), intent(in) :: this
 real(8), intent(in) :: x(3) 
 logical :: check
 
 check=.false.
 
 if ((this%xb(1,1)<=x(1)).and.(x(1)<=this%xb(2,1))) then
  if ((this%xb(1,2)<=x(2)).and.(x(2)<=this%xb(2,2))) then
   if ((this%xb(1,3)<=x(3)).and.(x(3)<=this%xb(2,3))) then
    check=.true.
   end if
  end if
 end if  
 
end function check


subroutine shoot(this,x0,nd,st,ti)
!===========================================================
! Check, wherther line drawn from point x0(3) in the 
! direction of vector nd(3) crosses the volume xb
! If the line x=x0+nd/|nd|*t crosses the volume then
! ti(2)=[t1,t2] is interval of t, when it occurs.
! Here nd/|nd| means normalization of vector nd, so
! vector nd/|nd| is a unit vector.
!===========================================================
 class(r0region) :: this
 real(8), intent(in) :: x0(3), nd(3) 
 real(8), intent(out) :: ti(2)
 real(8), parameter :: vsmall=1d-307, vlarge=1d307
 real(8) :: nn(3),  ti1(2), ti2(2), tx(2,3), tt
 integer :: i
 logical :: st, st1

 nn=nd/sqrt(dot_product(nd,nd)) 
 
 ti=0d0
 st=.false.
 
 do i=1,3
  if (abs(nn(i))<vsmall) then
   tx(1,i)=-vlarge
   tx(2,i)=vlarge
  else 
   tx(1,i)=(this%xb(1,i)-x0(i))/nn(i)
   tx(2,i)=(this%xb(2,i)-x0(i))/nn(i)
  end if
  
  if (tx(1,i)>tx(2,i)) then
   tt=tx(2,i)
   tx(2,i)=tx(1,i)
   tx(1,i)=tt
  end if 
 end do

 
 call isec(tx(:,1),tx(:,2),st1,ti1)
 if (st1) then
  call isec(ti1,tx(:,3),st1,ti2)
  if (st1) then
   st=.true.
   ti=ti2
  end if
 end if
 
end subroutine shoot


subroutine isec(x1,x2,st,x3)
!===========================================================
! Check for intersection between intervals x1=[x1(1),x1(2)]
! and x2=[x2(1),x2(2)], where xi(1)<xi(2). 
! If there is intersection then
! st=.true. and x3=[x3(1),x3(2)] is intersetion. If there is
! not then st=.false. and x3=[x3(1),x3(2)] interval between
! closest end of intervals x1 and x2.
!===========================================================
 real(8), intent(in) :: x1(2), x2(2) ! intervals for checking
 real(8), intent(out) :: x3(2)
 logical :: st
 
 st=.true.
 
 if ((x1(1)<=x2(1)).and.(x2(1)<=x1(2))) then
  x3(1)=x2(1)
 else if (x2(1)<x1(1)) then
  x3(1)=x1(1)
 else  
  x3(1)=x1(2)
  x3(2)=x2(1)
  st=.false.
  return
 end if
 
 if ((x1(1)<=x2(2)).and.(x2(2)<=x1(2))) then
  x3(2)=x2(2)
 else if (x1(2)<x2(2)) then
  x3(2)=x1(2)
 else  
  x3(1)=x2(2)
  x3(2)=x1(1)
  st=.false.
  return
 end if
 
end subroutine isec


end module restr_region


! program main
!  use restr_region, only : r0region
!  type(r0region) :: rempt
!  real(8) :: x(3), x1(3), bm(3), ti(2)
!  logical :: st
! 
!  rempt%xb(:,1)=[0d0,4d0] 
!  rempt%xb(:,2)=[0d0,1.15d0] 
!  rempt%xb(:,3)=[0d0,3d0] 
! 
!  x=[0d0,0d0,2.1d0]
! 
!  write(*,*) rempt%inside(x)
! 
!  x1=[2d0,0d0,0d0]
!  bm=[1d0,0d0,0d0]
!  call rempt%shoot(x1,bm,st,ti)
! 
!  write(*,*) st
!  write(*,'(10f14.6)') ti
!  write(*,*) rempt%mindl()

! end program main

