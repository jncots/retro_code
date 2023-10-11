module peak_integration
!============================================================
!
! Acceleration of integration of function with a sharp peak
!
! Usage:
! ======
! call peak_int(fun0,a,b,c,d1,d2,res,tol,nfun)
!
! tol and nfun are accuracy and number of calculated functions
! (optional)
! Calculation of the integral on a segment [a,b] with sharp
! peak at 'c', a<c<b (including logarithmic divergence).
! Acceleration of integration is reached by choosing
! appropriate change of variables with parameters of these
! changes d1, d2 (should be choosen by hand)
! The subroutine can be applied recursively up to 'nin' times
!
! Test program and example of usage are below
! Caution!
! ========
! There is risky 'if' condition (see below)
!============================================================
 use gauss_kronrod, only : gk_adaptive
 implicit none
 private
 save
 
 
 public :: peak_int
 
 abstract interface
  subroutine change_prot(t,x,dx)
   real(8), intent(in) :: t
   real(8), intent(out) :: x, dx
  end subroutine change_prot
 end interface
 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun_prot
 end interface
 
 type fun_set
  procedure (change_prot), pointer, nopass :: change=>null()
  procedure (fun_prot), pointer, nopass :: fun=>null()
 end type fun_set
 integer, parameter :: nin=100
 
 type(fun_set) :: ff(nin)
 
 
 integer :: ncur=0
 
 real(8) :: x1(nin), x2(nin), x0(nin), delta1(nin), delta2(nin)
 
 
 contains

 
recursive subroutine  peak_int(fun0,a,b,c,d1,d2,res,tol,nfun)
! Variables
 procedure (fun_prot) :: fun0
 real(8), intent(in) :: a, b, c, d1, d2
 real(8), intent(out) :: res
 real(8), optional, intent(in) :: tol
 integer, optional, intent(out) :: nfun
 real(8), parameter :: tol_def=1d-4 ! default value
 real(8) :: toler
 real(8) :: res1, res2, tmin, tmax
 integer :: nff1, nff2, nff
! Calculations
 ncur=ncur+1
 
 
 if (present(tol)) then
  toler=tol
 else
  toler=tol_def
 end if

 
 res1=0
 res2=0
 nff1=0
 nff2=0
 
 x1(ncur)=a
 x2(ncur)=b
 x0(ncur)=c
 delta1(ncur)=d1
 delta2(ncur)=d2
 
 tmin=0
 tmax=1
 if (c-a>0d0) then
  call set_change(change1,fun0)
  call gk_adaptive(fun_change,tmin,tmax,res1,tol=toler,nfun=nff1)
  res1=res1*d1*(c-a)
 end if
 if (b-c>0d0) then
  call set_change(change2,fun0)
  call gk_adaptive(fun_change,tmin,tmax,res2,tol=toler,nfun=nff2)
  res2=res2*d2*(b-c)
 end if 
 
 res=res1+res2
 nff=nff1+nff2
 if (present(nfun)) nfun=nff
 
 ff(ncur)%change=>null()
 ff(ncur)%fun=>null()
 
 ncur=ncur-1
end subroutine peak_int


subroutine change1(s,x,dx)
! Variables
 real(8), intent(in) :: s
 real(8), intent(out) :: x, dx
 real(8) :: sd
! Calculations 
 sd=s**delta1(ncur)
 x=x1(ncur)*sd+x0(ncur)*(1-sd)
 dx=sd/s
end subroutine change1

subroutine change2(s,x,dx)
! Variables
 real(8), intent(in) :: s
 real(8), intent(out) :: x, dx
 real(8) :: sd
! Calculations
 sd=s**delta2(ncur)
 x=x2(ncur)*sd+x0(ncur)*(1-sd)
 dx=sd/s
end subroutine change2


subroutine set_change(vchange,vfun)
! Variables 
 procedure(change_prot) :: vchange
 procedure(fun_prot) :: vfun
! Calculations 
 ff(ncur)%change=>vchange
 ff(ncur)%fun=>vfun
end subroutine set_change


subroutine fun_change(s,res)
! Variables
 real(8), intent(in) :: s
 real(8), intent(out) :: res
 real(8) :: x, dx
 real(8), parameter :: huge_num=1d307
! Calculations
 call ff(ncur)%change(s,x,dx)
 call ff(ncur)%fun(x,res)
 res=res*dx
! The following 'if' condition is quite risky
! and works only if the problem in small
! number of points, so it should be checked
! uncommenting the output
 if ((abs(res)>huge_num).or.(isnan(res))) then
  !write(*,*) 'nan', res
  !write(*,*) s, x, x0, x-x0
  !read(*,*)
  res=0d0
 end if
 
 
end subroutine fun_change


end module peak_integration


!=======================================================
!
! Test program and example of usage:
!
!=======================================================
!program main
! use peak_integration, only : peak_int
! use gauss_kronrod, only : gk_adaptive
! implicit none
! real(8) :: x1, x2, xpeak, d1, d2
! real(8) :: res1, res2, rel_err
! integer :: nff1, nff2
! 
! interface
!   subroutine test(x,res)
!    real(8), intent(in) :: x
!    real(8), intent(out) :: res
!   end subroutine
! end interface  
! 
! 
! 
! x1=0d0
! x2=1d0
! xpeak=0d0
! d1=0.5d0
! d2=10d0
! 
! call gk_adaptive(test,x1,x2,res1,tol=1d-4,eps=rel_err,nfun=nff1)
! call peak_int(test,x1,x2,xpeak,d1,d2,res2,tol=1d-8,nfun=nff2)
! 
! write(*,*) res1, res2, rel_err
! write(*,*) nff1, nff2
!
!end program main
!
!
!subroutine test(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! res=log(x)/x**0.9
!end subroutine test