module find_root_m
 implicit none
 private
 save


 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun_1d
 end interface


 public :: find_root

 type find_root
 contains
  procedure, nopass :: find=>find_sroot, set=>set_find_root
  procedure, nopass :: check=>check_roots, findm=>find_mroot
 end type find_root



 real(8) :: tol_c=1d-15
 integer :: max_eval_c=1000
 logical :: debug_mod=.false., silent_mod=.false.

 contains


subroutine find_mroot(x,fun,res,err_code)
 real(8), intent(in) :: x(:)
 procedure(fun_1d) :: fun
 real(8), allocatable, intent(out) :: res(:)
 integer, allocatable, intent(out), optional :: err_code(:)
 integer, allocatable :: err_code1(:)
 real(8), allocatable :: res1(:)
 integer :: i,  nroots, ntot
 real(8), allocatable :: xnew(:)


 call bracketing(x,fun,xnew)

 call check_roots(xnew,fun,res,err_code)
 ntot=size(res)
 allocate(res1(ntot))
 allocate(err_code1(ntot))

 nroots=0
 do i=1,size(res)
  if (err_code(i)/=3) then
   nroots=nroots+1
   res1(nroots)=res(i)
   err_code1(nroots)=err_code(i)
  end if
 end do
 
 deallocate(res)
 deallocate(err_code)


 if (nroots>0) then
  allocate(res(nroots))
  allocate(err_code(nroots))
  do i=1,nroots
   res(i)=res1(i)
   err_code(i)=err_code1(i)
  end do
 else
  allocate(res(1))
  allocate(err_code(1))
  res(1)=x(1)
  err_code(1)=3
 end if

 deallocate(res1)
 deallocate(err_code1)

end subroutine find_mroot


subroutine check_roots(x,fun,res,err_code)
 real(8), intent(in) :: x(:)
 procedure(fun_1d) :: fun
 real(8), allocatable, intent(out) :: res(:)
 integer, allocatable, intent(out), optional :: err_code(:)
 integer :: nseg, i

 nseg=size(x)-1
 
 if (allocated(res)) deallocate(res)
 if (allocated(err_code)) deallocate(err_code)

 allocate(res(nseg))
 allocate(err_code(nseg))

 do i=1,nseg
  call find_sroot(x(i),x(i+1),fun,res(i),err_code(i))
 end do
 
end subroutine check_roots



subroutine bracketing(x,fun,xnew)
 real(8), intent(in) :: x(:)
 procedure(fun_1d) :: fun
 real(8), allocatable, intent(out) :: xnew(:)
 integer :: nx, i, j
 real(8) :: fa, fb, res 
 real(8), allocatable :: xnew1(:)

 nx=size(x)
 
 allocate(xnew1(nx))

 j=1
 xnew1(1)=x(1)
 call fun(x(1),fa)
 
 do i=2,nx
  call fun(x(i),fb)
  if ((fb*fa)<=0d0) then
   j=j+1
   xnew1(j)=x(i)
   fa=fb
  end if
 end do

 if (allocated(xnew)) deallocate(xnew)
 allocate(xnew(j))

 do i=1,j
  xnew(i)=xnew1(i)
 end do
 
 deallocate(xnew1)

end subroutine bracketing



subroutine find_sroot(x1,x2,fun,res,err_code)
 real(8), intent(in) :: x1, x2
 procedure(fun_1d) :: fun
 real(8), intent(out) :: res
 integer, intent(out), optional :: err_code

 real(8) :: a, b, c
 real(8) :: fa, fb, fc
 integer :: stop_flag

 a=x1
 b=x2
 call fun(a,fa)
 call fun(b,fb)
 fc=fa

 do
  c=(a+b)/2
  call fun(c,fc)

  if (fc*fa<=0) then
   b=c
   fb=fc
  else
   a=c
   fa=fc
  end if


  call exit_cond(a,b,fa,fb,stop_flag)
  if (stop_flag>0) exit

 end do
 res=c
 if (present(err_code)) err_code=stop_flag
end subroutine find_sroot


subroutine exit_cond(a,b,fa,fb,stop_flag)
 real(8), intent(in) :: a, b, fa, fb
 integer, intent(out) :: stop_flag
 integer, save :: i=0
 real(8), parameter :: eps=epsilon(a)

 i=i+1

 if (debug_mod) write(*,*) i, ' a=', a, '  b=', b, '  fa=', fa, '  fb=', fb

 if ((b-a)<=eps) then
  if (.not.silent_mod) then
   write(*,*) 'find_root: b-a<=eps:'
   write(*,*) 'i=', i, ' b-a=', b-a, '  eps=', eps
   write(*,*) 'a=', a, '  b=', b, '  fa=', fa, '  fb=', fb
  end if
  stop_flag=2
  i=0
  return
 end if

 if (fa*fb>0d0) then
  if (.not.silent_mod) then
   write(*,*) 'find_root: fa*fb>0'
   write(*,*) 'a=', a, '  b=', b, '  fa=', fa, '  fb=', fb
  end if
  stop_flag=3
  i=0
  return
 end if

 if (i>=max_eval_c) then
  if (.not.silent_mod) then
   write(*,*) 'find_root: too many evaluations, limit of evaluations=', &
    max_eval_c
   write(*,*) 'a=', a, '  b=', b, '  fa=', fa, '  fb=', fb
  end if
  stop_flag=4
  i=0
  return
 end if

 if ((abs(fa)<tol_c).or.(abs(fb)<tol_c)) then
  stop_flag=1
  i=0
  return
 end if

 stop_flag=0

end subroutine exit_cond

subroutine set_find_root(cmode,tol,max_eval,debug,silent)
 integer, intent(in) :: cmode
 real(8), intent(in), optional :: tol
 integer, intent(in), optional :: max_eval
 logical, intent(in), optional :: debug, silent

 if (cmode==0) then ! default mode
  tol_c=1d-15
  max_eval_c=1000
  debug_mod=.false.
  silent_mod=.false.
  return
 end if

 if (present(tol)) tol_c=tol
 if (present(max_eval)) max_eval_c=max_eval
 if (present(debug)) debug_mod=debug
 if (present(silent)) silent_mod=silent

end subroutine set_find_root

end module find_root_m
