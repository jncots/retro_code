module rkutta_m
 implicit none
 private
 save


 public :: solve_ode, sys_ode_res

 abstract interface
  subroutine fun_deriv(x,y,dydx)
   real(8), intent(in) :: x, y(:)
   real(8), intent(out) :: dydx(:)
  end subroutine
 end interface


 type sys_ode_res
  integer :: np=0
  integer :: nsys=0, ntot=0
  real(8), allocatable :: y(:,:)
 contains
  procedure :: del=>del_res, set=>set_res, add=>add_res
 end type sys_ode_res



 contains




subroutine del_res(this)
 class(sys_ode_res) :: this
 
 if (allocated(this%y)) deallocate(this%y)
 this%np=0
 this%nsys=0
 this%ntot=0
end subroutine del_res


subroutine set_res(this,nsys,ntot)
 class(sys_ode_res) :: this
 integer :: nsys, ntot

 call this%del
 allocate(this%y(0:nsys,ntot))
 this%nsys=nsys
 this%ntot=ntot

end subroutine set_res



subroutine add_res(this,x,y)
 class(sys_ode_res) :: this
 real(8), intent(in) :: x, y(:)

 this%np=this%np+1
 if (this%np<=this%ntot) then
  this%y(0,this%np)=x
  this%y(1:,this%np)=y
 else
  this%np=this%ntot
 end if


end subroutine add_res


subroutine solve_ode(ns,y,dy,x1,x2,eps,h0,xout,res)
 integer, intent(in) :: ns
 real(8) :: y(ns)
 procedure(fun_deriv) :: dy

 real(8), intent(in) :: x1, x2, eps
 real(8) :: h0 
 real(8), intent(out) :: xout
 type(sys_ode_res), optional :: res
 real(8) :: yeps(ns), dydx(ns), x

 integer :: i
 logical :: set_yeps
 
 x=x1
 if (present(res)) call res%add(x,y) 
 
 do
  call dy(x,y,dydx)
  yeps=y+h0*dydx
   
  set_yeps=.false.
  do i=1,ns
   if (abs(yeps(i))<1d-2) then
    set_yeps=.true.
    exit
   end if
  end do 
  if (set_yeps) yeps=1d-2

  call rk_adap_step(ns,y,dy,x,h0,yeps,eps)
!  write(*,'(A,10Es14.6)') 'h0,x,y,yeps=', h0, x, y, yeps
!  read(*,*)
  if (present(res)) call res%add(x,y)
  
  if (h0>0d0) then
   if (x+h0>x2) h0=x2-x
   if (x>=x2) exit
  else
   if (x+h0<x2) h0=x2-x
   if (x<=x2) exit
  end if
 end do

 xout=x

end subroutine solve_ode

      
subroutine rk_adap_step(ns,y,dy,t0,h0,yeps,eps)
!=======================================================================
! Stepper subroutine for the Embedded Runga-Kutta Method with 
! with Cash-Karp Parameters "rkutta_ck5"
!=======================================================================      
 integer, intent(in) :: ns
 real(8) :: y(ns)
 procedure(fun_deriv) :: dy
 real(8) :: t0, h0
 real(8), intent(in) :: yeps(ns), eps
 real(8) :: tnew, y_old(ns), yerr(ns), ypr(ns), rel_err
 logical :: stop_calc, non_zero, zero_true

 do
  
  call small_step(t0,h0,tnew,stop_calc)
  if (stop_calc) exit
   
  y_old=y
  call rkutta_ck5(ns,y,dy,t0,h0,yerr)
  
  ypr=abs(yeps/yerr)
  call check_zero(ns,yerr,zero_true)
  
!  write(*,'(A,10Es14.6)') 't0=', t0
!  write(*,*) 'yeps, yerr, ypr, zero_true:=', yeps, yerr, ypr, zero_true

  if (zero_true) then
    h0=2*h0
   exit
  end if
  
  call min_grzero(ns,ypr,rel_err,non_zero)

  if (non_zero) then
   rel_err=rel_err*eps
   call adjust_step(rel_err,h0,stop_calc)
   if (stop_calc) exit
  else
   exit
  end if

  y=y_old
 end do
 
 t0=tnew

end subroutine rk_adap_step


subroutine check_zero(ns,y,zero_true)
 integer, intent(in) :: ns
 real(8) :: y(ns)
 real(8), parameter :: zero=tiny(y(1))
 integer :: i 
 logical, intent(out) :: zero_true 


 zero_true=.true.
 do i=1,ns
  if (abs(y(i))>zero) then
   zero_true=.false.
   exit
  end if
 end do

end subroutine check_zero



subroutine min_grzero(ns,y,res,ires)
 integer, intent(in) :: ns 
 real(8), intent(in) :: y(ns)
 real(8), intent(out) :: res
 logical, intent(out) :: ires
 real(8), parameter ::  zero=tiny(res)
 integer :: i
 
 if (y(1)>zero) then
  res=y(1)
  ires=.true.
 else
  ires=.false.
 end if

 do i=2,ns
  
  if (ires) then
   if ((y(i)>zero).and.(y(i)<res)) res=y(i)
  else
   if (y(i)>zero) then
    res=y(i)
    ires=.true.
   end if
  end if
 
 end do

end subroutine min_grzero
   

subroutine small_step(t0,step,tnew,stop_calc)
 real(8), intent(in) :: t0, step
 real(8), intent(out) :: tnew
 logical, intent(out) :: stop_calc

 tnew=t0+step
 if (tnew==t0) then
  write(*,*) 'rkutta: too small step'
  stop_calc=.true.
 else
  stop_calc=.false.
 end if

end subroutine small_step


subroutine adjust_step(rel_err,step,stop_calc)
 real(8), intent(in) :: rel_err
 real(8) :: step
 logical, intent(out) :: stop_calc
 real(8), parameter :: fmin=1d-1, fmax=5d0 ! max decrease and increase of step
 real(8) :: fact


 if (rel_err<1d0) then
  
  stop_calc=.false.
  fact=0.9d0*rel_err**0.25d0
  if (fact>fmin) then
   step=step*fact
  else
   step=step*fmin
  end if

 else
  
  stop_calc=.true.
  fact=0.9d0*rel_err**0.2d0
  if (fact<fmax) then
   step=step*fact
  else
   step=step*fmax
  end if
 
 end if

end subroutine adjust_step


subroutine rkutta_ck5(ns,y,dy,t0,h0,yerr)
!=======================================================================
! Embedded Runga-Kutta Method with Cash-Karp Parameters 
! (see Numerical Recipes in FORTRAN 77, p.711 )
! Solution of ODE y'=f(x), y(i=1..n)
! On every step h 'rkutta_ck5' gives absolute error derr
! and solution y(i=1..n) at the point t with precision O(h^5)
! derives is a subroutine which calculates f(x)=dydx(i=1..n)
!=======================================================================
 integer, intent(in) :: ns
 real(8) :: y(ns)
 procedure(fun_deriv) :: dy
 real(8), intent(in) :: t0, h0
 real(8), intent(out) :: yerr(ns)
 real(8), dimension(ns) :: a1, a2, a3, a4, a5, a6, pc
 real(8) :: tc

 real(8), parameter :: d2=0.2d0, d3=0.3d0, d4=0.6d0, d6=0.875d0
 real(8), parameter :: b21=0.2d0, b31=0.075d0, b32=0.225d0
 real(8), parameter :: b41=0.3d0, b42=-0.9d0, b43=1.2d0
 real(8), parameter :: b51=-0.2037037037037037d0, b52=2.5d0
 real(8), parameter :: b53=-2.592592592592593d0, b54=1.296296296296296d0
 real(8), parameter :: b61=0.2949580439814815d-1, b62=0.341796875d0
 real(8), parameter :: b63=0.4159432870370370d-1, b64=0.4003454137731481d0
 real(8), parameter :: b65=0.6176757812500000d-1
 real(8), parameter :: c1=0.9788359788359788d-1, c3=0.4025764895330113d0
 real(8), parameter :: c4=0.2104377104377104d0,  c6=0.2891022021456804d0
 real(8), parameter :: e1=-0.429377480158732d-2, e3=0.186685860938579d-1
 real(8), parameter :: e4=-0.341550268308081d-1, e5=-0.1932198660714286d-1
 real(8), parameter :: e6=0.391022021456804d-1
 real(8), parameter :: c1z=0.1021773726851852d0, c3z=0.3839079034391534d0
 real(8), parameter :: c4z=0.2445927372685185d0, c5z=0.1932198660714286d-1
 real(8), parameter :: c6z=0.25d0


 call dy(t0,y,a1)
 a1=a1*h0

 tc=t0+d2*h0
 pc=y+b21*a1
 call dy(tc,pc,a2)
 a2=a2*h0
 
 tc=t0+d3*h0
 pc=y+b31*a1+b32*a2
 call dy(tc,pc,a3)
 a3=a3*h0
 
 tc=t0+d4*h0
 pc=y+b41*a1+b42*a2+b43*a3
 call dy(tc,pc,a4)
 a4=a4*h0
 
 tc=t0+h0
 pc=y+b51*a1+b52*a2+b53*a3+b54*a4
 call dy(tc,pc,a5)
 a5=a5*h0
 
 tc=t0+d6*h0       
 pc=y+b61*a1+b62*a2+b63*a3+b64*a4+b65*a5
 call dy(tc,pc,a6)
 a6=a6*h0
          
 y=y+c1*a1+c3*a3+c4*a4+c6*a6
 yerr=e1*a1+e3*a3+e4*a4+e5*a5+e6*a6

end subroutine rkutta_ck5



end module rkutta_m


!module test1
! implicit none
! private
! save

! public :: main_calc

! contains


!subroutine main_calc
! use rkutta_m, only : solve_ode, sys_ode_res
! integer, parameter :: ns=3
! integer :: i
! real(8) :: y(ns), yeps(ns), dydx(ns)
! real(8) :: x1, x2, eps, h, xout
! type(sys_ode_res) :: res
! 
! x1=0d0
! x2=3d1
! y=[1d0,0d0,0d0]
! eps=1d-6
! h=1d-1
! 
! call res%set(ns,6000)
! call solve_ode(ns,y,lorentz_sys,x1,x2,eps,h,xout,res)
! 
! write(*,*) res%np, xout
! open(1,file='ftest.dat')
! do i=1,res%np
!  write(1,*) res%y(:,i)
! end do
! close(1)

!end subroutine main_calc


!subroutine ode_sys(x,y,dydx)
! real(8), intent(in) :: x, y(:) 
! real(8), intent(out) :: dydx(:)
! 
! dydx(1)=y(2)
! dydx(2)=-y(1)

!end subroutine ode_sys


!subroutine lorentz_sys(x,y,dydx)
! real(8), intent(in) :: x, y(:) 
! real(8), intent(out) :: dydx(:)
! 
! dydx(1)=10*(y(2)-y(1))
! dydx(2)=y(1)*(27-y(3))-y(2)
! dydx(3)=y(1)*y(2)-8*y(3)/3

!end subroutine lorentz_sys



!end module test1




!program main
! use test1, only : main_calc

! call main_calc

!end program main
