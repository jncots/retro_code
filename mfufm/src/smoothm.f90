module smoothm
!========================================================
! Smoothing of the drop transition of a function:
! Example: the function is given f(x)=y1, if x<a and
! f(x)=y2 if x>=a. To smooth transition we can use the 
! following:
! use smoothm, only : th_smooth
! type(th_smooth) :: a
! call a%set(x1,x2,y1,y2)! , where x1<a and x2>2
! call a%fx(x,res)!, where x1<x<x2 gives smooth
! transition with approximatly fx(x1)=y1 and fx(x2)=y2
! The example program is below
!========================================================
 implicit none
 private
 save
 
 public :: th_smooth
 
 type th_smooth
  logical, private :: nact=.true.
  real(8), private :: a0, d0, x0, y1
 contains 
  procedure :: set, fx
 end type th_smooth
 
 contains


subroutine set(this,x1,x2,y1,y2)
 class(th_smooth) :: this
 real(8), intent(in) :: x1, x2, y1, y2
 
 this%nact=.false.
 this%a0=(y2-y1)/2
 this%d0=6/(x2-x1)
 this%x0=(x2+x1)/2
 this%y1=y1
 
end subroutine set


subroutine fx(this,x,res)
 class(th_smooth) :: this
 real(8), intent(in) ::  x
 real(8), intent(out) :: res
 real(8) :: xa
 
 if (this%nact) then
  write(*,*) 'type th_smooth: The function has not been initiated'
  res=0d0
  return
 end if 
 
 xa=this%d0*(x-this%x0)
 res=this%y1+this%a0*(1+tanh(xa))

end subroutine fx

end module smoothm


! program main
 ! use smoothm, only : th_smooth
 ! use utools_mod, only : utools
 ! type(th_smooth) :: a
 ! type(utools) :: ut
 ! real(8), allocatable :: x(:)
 ! real(8) :: res
 
 ! call a%set(10d0,11d0,4.5d0,-8d0)
 ! call ut%grid(x,0d0,20d0,1000,'lin')
 
 ! open(1,file='sm_test.dat')
 ! do i=1,size(x)
  ! call a%fx(x(i),res)
  ! write(1,*) x(i), res
 ! end do
 ! close(1)
! end program main