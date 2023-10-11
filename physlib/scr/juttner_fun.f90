module juttner_fun
!===============================================
! Contains a propagator of the cosmic rays
! in spherically symmetric case, which includes
! both diffusion and ballistic regime.
! The formula is given in Aloisio et at (2009)
! (2009ApJ...693.1275A), Eq.19
!
! The module uses a function BesselK1 (dbesk1)
! from Slatec library
!===============================================
 implicit none
 private
 save
 
 real(8), parameter :: pi=3.14159265358979324d0
 real(8), parameter :: pi4=1/(4*pi)
 real(8), parameter :: pi2=sqrt(2/pi)*pi4
 real(8), parameter :: cs=2.99792458d10 ! cm/s
 
 
 public :: pjun
 
 contains


subroutine pjun(t,r,diff,res)
!===========================================================
! t is time in sec, r is a distance from the source in cm
! diff is the diffusion coefficient
!===========================================================
 real(8), intent(in) :: t, r, diff
 real(8), intent(out) :: res
 real(8) :: ct, y, x

 ct=cs*t
 y=r/ct
 
 if (y>=1d0) then
  res=0d0
  return
 end if
 
 x=ct*cs/(2*diff)
 call pjun0(x,y,res)
 res=res/ct**3

end subroutine pjun



subroutine pjun0(x,y,res)
!=======================================================
! Calculation of the function:
! Pj=1/(Z(x)*(1-y**2)**2)*exp(-x/sqrt(1-y**2))
!=======================================================
 real(8), intent(in) :: x, y
 real(8), intent(out) :: res
 real(8), parameter :: xmin=1d-10, xmax=700d0
 real(8) :: y2, t, t2
 
 interface
! dbesk1(x) is a BesselK1 from Slatec library
  real(8) function dbesk1(x)
  real(8) :: x
  end function dbesk1
 end interface
 
 y2=y**2
 t=1-y2
 t2=sqrt(t)
 
 if (x<xmin) then
  res=pi4*x*x*exp(-x/t2)
  res=res/t**2
  return
 end if
 
 if (x>xmax) then
  res=pi2*x*sqrt(x)
  res=res*exp(-x*y2/(t2*(1+t2)))
  res=res/t**2
  return
 end if 

 res=dbesk1(x)
 res=pi4*x/res
 res=res*exp(-x/t2)
 res=res/t**2
end subroutine pjun0



!subroutine test
! integer, parameter :: nt=7
! integer :: i, j, nx
! real(8) :: diff, x1, x2, tt, rr
! real(8) :: tarr(nt), res(nt)
! 
! tarr=[1d0,3d0,1d1,3d1,1d2,3d2,1d3]
! 
! diff=1d30
! 
! x1=1d-3*psec
! x2=1d3*psec
! nx=1000
! 
! open(1,file='junt_dist.dat')
! 
! do i=0,nx
!  rr=x1*(x2/x1)**(1d0*i/nx)
!  
!  do j=1,nt
!   tt=tarr(j)*year
!   call pjun(tt,rr,diff,res(j))
!  end do 
!  
!  write(1,*) rr/psec, (res(j),j=1,nt)
!  
! end do
!
! close(1)
!
!
!end subroutine test


end module juttner_fun

!program main
! use juttner_fun, only : test
! implicit none
! 
! 
!
! call test
!
!
!
!end program main