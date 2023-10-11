module tab2fun
!================================================================
! Creates a continuous function of one argument from data arrays
! x(n), f(n), where x is argument, and f is value of the function.
! Linking
! =======
! use tab2fun, only : tabfun
! Initiation:
! ===========
! call ff%on(x,f) ! ff is a name of variable of type(tabfun)
! logarithmic interpolation is by default
! or
! call ff%on(x,f,'lin'), or call ff%on(x,f,'log'), where
! x(n), f(n) is data arrays
! Usage:
! ======
! ff%xmin, ff%xmax is minimum and maximum of the function:
! call ff%fun(x1,res) is a continuous function, where
! x1 is value(real(8)) and res is a result(real(8))
! ff%zero is variable determined by zero values of f array 
! ff%zero=0 - the function will given only 0
! ff%zero=1 - the ordinary operation
! Auxiliary function:
! ====================
! Assosiation of the function ff%fx(x,res)
! with bound-type function ff%fun(x,res) (max 25 assosiations)
! (useful for usage as argument of a subroutine with argument
! of 1d function) is at ff%on call
! call ff%off ! deassosiation of the ff%fx(x,res) 
! Deassosication should be done as the number of assosiations
! is restricted to 25. It can be increased by changing 'nfint'
! and adding the corresponding number of functions funX(x,res)
! where X is a number.
!
! Version history:
! ================
! 08.04.2016 Add the type procedure ff%fx, change name of
! type-bound functions to on and off
!================================================================
 use intpol_mod, only : zeros_cut, arr_ind_short, lin_int, log_int
 implicit none
 private
 save
 
 
 public :: tabfun 
 
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 
 type tabfun
  integer :: n1, n2, typ, zero
  real(8) :: xmin, xmax
  real(8), pointer :: x(:), f(:)
  integer, private :: nfx=0
  procedure(fun_1d), pointer, nopass :: fx=>null()
 contains
  procedure :: on=>tab, fun, off=>fx_off
 end type tabfun
 
 
 type tabfun_point
  class(tabfun), pointer :: p=>null()
  procedure(fun_1d), pointer, nopass :: f=>null()
 end type tabfun_point
 
  
 integer, parameter :: nfint=25
 type(tabfun_point) :: fint(nfint)  ! internal variable
 logical :: nfree(nfint)=.false. ! all functions are free
 logical :: connect=.true.
 
 contains


subroutine tab(fdat,xa,fa,typ)
! Variables 
 class(tabfun), target :: fdat
 real(8), intent(in), target :: xa(:), fa(:)
 character(3), optional, intent(in) :: typ ! 'lin' or 'log'
 character(3) :: typ1
 integer :: nfx
! Calculations
 if (connect) then    ! At first usage of the procedure
  call set_point_fun
  connect=.false. 
 end if
 
 if (fdat%nfx==0) then ! Connection of fx
  call free_fun(nfx)
  nfree(nfx)=.true.
  fdat%nfx=nfx
  fint(nfx)%p=>fdat
  fdat%fx=>fint(nfx)%f
 end if
 
 if (size(xa)/=size(fa)) then
  write(*,*) 'subroutine tab (type tabfun): disparity of argument and function array sizes'
  return
 end if
 
 fdat%x=>xa
 fdat%f=>fa
 if (present(typ)) then
  typ1=typ
 else
  typ1='log'
 end if 
 if (typ1=='lin') fdat%typ=0
 if (typ1=='log') fdat%typ=1
 
 call zeros_cut(fdat%f,fdat%n1,fdat%n2)
 if (fdat%n1==fdat%n2) then
  fdat%xmin=0d0
  fdat%xmax=0d0
  fdat%zero=0
 else 
  fdat%xmin=fdat%x(fdat%n1)
  fdat%xmax=fdat%x(fdat%n2)
  fdat%zero=1
 end if
 
 
end subroutine tab


subroutine fun(fdat,x,res)
! Variables 
 class(tabfun), intent(in) :: fdat
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 integer :: n1, n2
! Calculations
  
 if (fdat%zero==0) then
  res=0d0
  return
 end if
 
 if ((x<fdat%xmin).or.(x>fdat%xmax)) then
  res=0d0
  return
 end if
 
 call arr_ind_short(fdat%n1,fdat%n2,fdat%x,x,n1,n2) 
 if (fdat%typ==0) then
  call lin_int(fdat%x(n1),fdat%x(n2),fdat%f(n1),fdat%f(n2),x,res)
 else
  call log_int(fdat%x(n1),fdat%x(n2),fdat%f(n1),fdat%f(n2),x,res)
 end if 
 
 if (isnan(res)) call lin_int(fdat%x(n1),fdat%x(n2),fdat%f(n1),fdat%f(n2),x,res)
 
end subroutine fun



subroutine fx_off(fdat)
! Variables 
 class(tabfun) :: fdat
! Calculations
 if (fdat%nfx==0) return
 fint(fdat%nfx)%p=>null()
 fdat%fx=>null()
 nfree(fdat%nfx)=.false.
 fdat%nfx=0
end subroutine fx_off


subroutine free_fun(i)
 integer, intent(out) :: i 
 do i=1,nfint
  if (.not.nfree(i)) return
 end do 
 i=nfint+1
 write(*,*) 'free_fun: All functions are occupied'
end subroutine free_fun

subroutine set_point_fun
 fint(1)%f=>fun1
 fint(2)%f=>fun2
 fint(3)%f=>fun3
 fint(4)%f=>fun4
 fint(5)%f=>fun5
 fint(6)%f=>fun6
 fint(7)%f=>fun7
 fint(8)%f=>fun8
 fint(9)%f=>fun9
 fint(10)%f=>fun10
 fint(11)%f=>fun11
 fint(12)%f=>fun12
 fint(13)%f=>fun13
 fint(14)%f=>fun14
 fint(15)%f=>fun15
 fint(16)%f=>fun16
 fint(17)%f=>fun17
 fint(18)%f=>fun18
 fint(19)%f=>fun19
 fint(20)%f=>fun20
 fint(21)%f=>fun21
 fint(22)%f=>fun22
 fint(23)%f=>fun23
 fint(24)%f=>fun24
 fint(25)%f=>fun25
end subroutine set_point_fun



subroutine fun1(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(1)%p%fun(x,res)
end subroutine fun1

subroutine fun2(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(2)%p%fun(x,res)
end subroutine fun2

subroutine fun3(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(3)%p%fun(x,res)
end subroutine fun3

subroutine fun4(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(4)%p%fun(x,res)
end subroutine fun4

subroutine fun5(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(5)%p%fun(x,res)
end subroutine fun5

subroutine fun6(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(6)%p%fun(x,res)
end subroutine fun6

subroutine fun7(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(7)%p%fun(x,res)
end subroutine fun7

subroutine fun8(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(8)%p%fun(x,res)
end subroutine fun8

subroutine fun9(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(9)%p%fun(x,res)
end subroutine fun9

subroutine fun10(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(10)%p%fun(x,res)
end subroutine fun10

subroutine fun11(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(11)%p%fun(x,res)
end subroutine fun11

subroutine fun12(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(12)%p%fun(x,res)
end subroutine fun12

subroutine fun13(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(13)%p%fun(x,res)
end subroutine fun13

subroutine fun14(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(14)%p%fun(x,res)
end subroutine fun14

subroutine fun15(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(15)%p%fun(x,res)
end subroutine fun15

subroutine fun16(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(16)%p%fun(x,res)
end subroutine fun16

subroutine fun17(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(17)%p%fun(x,res)
end subroutine fun17

subroutine fun18(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(18)%p%fun(x,res)
end subroutine fun18

subroutine fun19(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(19)%p%fun(x,res)
end subroutine fun19

subroutine fun20(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(20)%p%fun(x,res)
end subroutine fun20


subroutine fun21(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(21)%p%fun(x,res)
end subroutine fun21

subroutine fun22(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(22)%p%fun(x,res)
end subroutine fun22

subroutine fun23(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(23)%p%fun(x,res)
end subroutine fun23

subroutine fun24(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(24)%p%fun(x,res)
end subroutine fun24

subroutine fun25(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 call fint(25)%p%fun(x,res)
end subroutine fun25


end module tab2fun

!! Test program and example of usage
!program main
!! Variables
! use tab2fun, only : tabfun
! use gauss_kronrod, only : gk_adaptive
! implicit none
! type(tabfun) :: ff, ff1, ff2
! integer :: n, i
! real(8), allocatable :: x(:), f(:), f1(:), f2(:)
! real(8) :: x1, x2, xx, res, res1, tmin, tmax
!! Calculations 
! 
! n=100
! allocate(x(n))
! allocate(f(n))
! allocate(f1(n))
! allocate(f2(n))
! x1=1d1
! x2=1d3
! do i=1,n
!  x(i)=x1*(x2/x1)**(i*1d0/n)
!  call test(x(i),f(i))
!  if ((i<8).or.(i>88)) f(i)=0d0
! end do
! 
! call ff%on(x,f)
! f1=10*f
! call ff1%on(x,f1)
! f2=100*f
! call ff2%on(x,f2)
! write(*,*) ff%zero, ff%xmin, ff%xmax
! 
! xx=35
! call test(xx,res)
! call ff%fx(xx,res1)
! write(*,*) res, res1, abs(res-res1)/res
! call gk_adaptive(ff%fx,ff%xmin,ff%xmax,res1)
! call gk_adaptive(test,ff%xmin,ff%xmax,res)
! write(*,*) res, res1, abs(res-res1)/res
! call gk_adaptive(ff1%fx,ff1%xmin,ff1%xmax,res1)
! call ff1%off
! call ff%off
! write(*,*) res1
! call gk_adaptive(ff2%fx,ff2%xmin,ff2%xmax,res1)
! write(*,*) res1
!
!end program main
!
!
!subroutine test(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
!
! res=x**2/(x**2+1+x*log(x)+x)
!
!end subroutine test