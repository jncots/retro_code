module int_3d_mod
!============================================================
! Triple integral is calculated using object of type(int_3dt)
! Example program is in the end of this file
! Usage:
! =====
! 
! type(int_3dt) :: fg
! 
! Settings:
! =========
! for example
! call fg%set(fun)
! call fg%set_l1(-2d0,2d0,2,tol=1d-8)
! call fg%set_l2(fx1,fx2,3,tol=1d-8)
! call fg%set_l3(fy1,fy2,1,2,3,tol=1d-8)
!
! The designations in a given case are following:
! int_l1 dx(2) int_l2 dx(3) int_l3 dx(1) fun(1,2,3)
! int_l1 gives the outer integration
! int_l2 is middle integration
! int_l3 is inner integration
! 3rd arguments of each functions is number of argument
! of fun(1,2,3) for corresponding integration.
!
! Thus, call fg%set_l1(-2d0,2d0,2,tol=1d-8) means that outer
! integration is performed along 2nd argument of function fun
! in the region from -2 to 2 with accuracy 1e-8
!
! call fg%set_l2(fx1,fx2,3,tol=1d-8) means that middle integration
! is performed along 3rd argument of function fun in the limits
! fx1(2) .. fx2(2), where 2 is 2nd argument function fun, as
! the outer integration is performed along 2nd argument
!
! call fg%set_l3(fy1,fy2,1,2,3,tol=1d-8) means that inner integration
! is performed along 1st argument of fun, in the limits
! fy1(2,3) ... fy2(2,3). Thus 4th and 5th arguments of set_l3 related
! to the functions fy1 and fy2. If we would have
! call fg%set_l3(fy1,fy2,1,3,2,tol=1d-8) then the integrations would be
! performed is the limits fy1(3,2)...fy2(3,2)
!
! Results:
! =========
! call fg%calc(res)
! write(*,*) res
!============================================================
 use gauss_kronrod, only : gk_adaptive
 implicit none
 private
 save
 
 public :: int_3dt
 
 
 abstract interface
  subroutine fun3d(x,y,z,res)
   real(8), intent(in) :: x, y, z
   real(8), intent(out) :: res
  end subroutine fun3d
  
  subroutine fun2d(x,y,res)
   real(8), intent(in) :: x, y
   real(8), intent(out) :: res
  end subroutine fun2d
  
  subroutine fun1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun1d
 end interface
 
 
 
 type int_3dt
  real(8) :: b01, b02
  procedure (fun1d), pointer, nopass :: b11=>null(), b12=>null()
  procedure (fun2d), pointer, nopass :: b121=>null(), b122=>null() 
  procedure (fun3d), pointer, nopass :: f3=>null()
  real(8) :: tol1, tol2, tol3
  real(8), private :: x0(3)
  integer :: i1, i2, i3, j1, j2
 contains
  procedure :: set, set_l1, set_l2, set_l3, calc
 end type int_3dt
 
 
 
 type(int_3dt), pointer :: lc
 real(8), parameter :: tol_def=1d-4 ! default value of error tolerance
 
 
 contains

subroutine set(this,fun)
 class(int_3dt) :: this
 procedure(fun3d) :: fun
! Calculations
 this%f3=>fun
end subroutine set
 
 
subroutine set_l1(this,x1,x2,i1,tol)
 class(int_3dt) :: this
 real(8), intent(in) :: x1, x2
 integer, intent(in) :: i1
 real(8), optional, intent(in) :: tol ! error tolerance
 
 this%b01=x1 
 this%b02=x2
 this%i1=i1
 
 if (present(tol)) then
  this%tol1=tol
 else
  this%tol1=tol_def
 end if
 
end subroutine set_l1
 
subroutine set_l2(this,f1,f2,i2,tol)
 class(int_3dt) :: this
 procedure(fun1d) :: f1, f2
 integer, intent(in) :: i2
 real(8), optional, intent(in) :: tol ! error tolerance
 
 this%b11=>f1
 this%b12=>f2
 this%i2=i2
 
 if (present(tol)) then
  this%tol2=tol
 else
  this%tol2=tol_def
 end if
 
end subroutine set_l2 

subroutine set_l3(this,f1,f2,i3,j1,j2,tol)
 class(int_3dt) :: this
 procedure(fun2d) :: f1, f2
 integer, intent(in) :: i3, j1, j2
 real(8), optional, intent(in) :: tol ! error tolerance
 
 this%b121=>f1
 this%b122=>f2
 this%i3=i3
 this%j1=j1
 this%j2=j2
 
 if (present(tol)) then
  this%tol3=tol
 else
  this%tol3=tol_def
 end if
 
end subroutine set_l3 
 
 


subroutine calc(this,res)
 class(int_3dt), target :: this
 real(8), intent(out) :: res
 integer :: sq(3), sr(3), kk

! Check of indices
 sq=[this%i1,this%i2,this%i3]
 sr=[this%i3,this%j1,this%j2]
 call check(sq,sr,kk)
 if (kk==0) then
  res=0d0
  return
 end if 
 
 lc=>this
 call gk_adaptive(int_1,lc%b01,lc%b02,res,tol=lc%tol1)
end subroutine calc


subroutine int_1(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x1, x2
 
 lc%x0(lc%i1)=x
 call lc%b11(lc%x0(lc%i1),x1)
 call lc%b12(lc%x0(lc%i1),x2)
 
 call gk_adaptive(int_2,x1,x2,res,tol=lc%tol2)
end subroutine int_1


subroutine int_2(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x1, x2
 
 lc%x0(lc%i2)=x 
 call lc%b121(lc%x0(lc%j1),lc%x0(lc%j2),x1)
 call lc%b122(lc%x0(lc%j1),lc%x0(lc%j2),x2)
 call gk_adaptive(int_3,x1,x2,res,tol=lc%tol3)
end subroutine int_2


subroutine int_3(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 lc%x0(lc%i3)=x
 call lc%f3(lc%x0(1),lc%x0(2),lc%x0(3),res)
end subroutine int_3



subroutine check(sq,sr,res)
 integer, intent(in) :: sq(3), sr(3)
 integer, intent(out) :: res
 integer :: kk
 
 res=1
 call check_arr(sq,kk)
 if (kk==0) then
  write(*,*) 'Setting of the order of integration is wrong:'
  write(*,*) 'set_l1(...,i1),set_l2(...,i2),set_l3(...,i3,j1,j2)'
  write(*,*) 'i1, i2, i3 should be equal to 1, 2, 3 and not equal to each other'
  res=0
  return
 end if
 
 call check_arr(sr,kk)
 if (kk==0) then
  write(*,*) 'Setting of the order of integration is wrong:'
  write(*,*) 'set_l3(...,i3,j1,j2) '
  write(*,*) 'i3, j1, i2 should be equal to 1, 2, 3 and not equal to each other'
  res=0
  return
 end if 
 
 
end subroutine check


subroutine check_arr(sq,res)
 integer, intent(in) :: sq(3)
 integer, intent(out) :: res
 integer :: i, i1 
 res=1
 do i=1,3
  if ((sq(i)<1).or.(sq(i)>3)) then
   res=0
   return
  end if
 end do
 
 
 do i=1,3
  i1=i+1
  if (i1>3) i1=1
  if (sq(i)==sq(i1)) then
   res=0
   return
  end if
 end do
 
end subroutine check_arr



end module int_3d_mod



! module test_int3d
 ! use int_3d_mod, only : int_3dt
 ! implicit none
 ! private
 ! save
 
 ! public :: main_calc
 
 ! contains
 
! subroutine main_calc
 ! type(int_3dt) :: fg
 ! real(8) :: res
 
 ! call fg%set(fun)
 ! call fg%set_l1(-2d0,2d0,2,tol=1d-8)
 ! call fg%set_l2(fx1,fx2,3,tol=1d-8)
 ! call fg%set_l3(fy1,fy2,1,2,3,tol=1d-8)

 ! call fg%calc(res)
 ! write(*,*) res

! end subroutine main_calc


! subroutine fx1(x,res)
 ! real(8), intent(in) :: x
 ! real(8), intent(out) :: res
 ! res=-sqrt(4-x**2)
! end subroutine fx1

! subroutine fx2(x,res)
 ! real(8), intent(in) :: x
 ! real(8), intent(out) :: res
 ! res=sqrt(4-x**2)
! end subroutine fx2

! subroutine fy1(x,y,res)
 ! real(8), intent(in) :: x, y
 ! real(8), intent(out) :: res
 ! res=-sqrt(12-x**2-y**2)
! end subroutine fy1


! subroutine fy2(x,y,res)
 ! real(8), intent(in) :: x, y
 ! real(8), intent(out) :: res
 ! res=sqrt(10-x**2-y**2)
! end subroutine fy2
 
! subroutine fun(x,y,z,res) 
 ! real(8), intent(in) :: x, y, z
 ! real(8), intent(out) :: res
 ! res=(x**2-(z+y)**2)**2
! end subroutine fun


! end module test_int3d



! program main
 ! use test_int3d, only : main_calc
 ! call main_calc
! end program main