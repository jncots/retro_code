module lookup2d
!=====================================================================
! lookup2d is similar to lookup, but remembers the order of calling
! subroutine fun(x,res,id), where id is calling number
!
! Contains subroutine for function tabulation 'looktab'
! Usage:
! =====
! Minimal:
!
! use lookup, only : looktab
! real(8) :: x1, x2
! real(8), allocatable :: x(:), fx(:)
! interface
! subroutine fun(x,res)
!  real(8), intent(in) :: x
!  real(8), intent(out) :: res
! end subroutine fun
! end interface
! x1=1d-2 ! [x1,x2] is tabulation region of fun
! x2=50d0
! call looktab(fun,x1,x2,x,fx) ! x(:), fx(:) is result
!
! Additional options:
! 'iscale' determines initial scale of tabulation of
! 'inx' initial points:
! iscale='lin' or 'log', by defaut iscale=df_iscale='log'
! by default inx=df_inx=100
! nmax is maximum number of points, irrespective to accuracy
! nx is exact number of point, irrespective to accuracy
! tol is accuracy of tabulation = error in integral/integral,
! where integral is integral calculated using tabulated points
! tol=df_tol=1d-4 by default
! int_res is integral calculated using tabulated points
! int_err is attained accuracy = error in integral/integral
!=====================================================================
 implicit none
 private
 save
 
 public :: looktab2d
 
 abstract interface
  subroutine ifun1d(x,res,id)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
   integer, intent(out) :: id
  end subroutine ifun1d
 end interface
 
 
 
 type point
  integer :: n ! id of the point
  real(8) :: x, f
 end type point
 
 type segment
  type(point) :: p1, p2, pm
  real(8) :: a1, a2, dx, eps
 end type segment
 
 
 contains


subroutine looktab2d(fun,x1,x2,x,fx,idx,iscale,inx,nx,nmax,tol,&
 int_res,int_err)
 procedure(ifun1d) :: fun
 real(8), intent(in) :: x1, x2
 real(8), allocatable :: x(:), fx(:)
 integer, allocatable :: idx(:)
 type(point), allocatable :: pfx(:)
 character(3), optional :: iscale
 integer, intent(in), optional :: inx, nmax, nx
 real(8), intent(in), optional :: tol
 real(8), intent(out), optional :: int_res, int_err ! integral value and error
 character(3), parameter :: df_iscale='log'
 real(8), parameter :: df_tol=1d-4 ! default tol
 integer, parameter :: df_inx=100, df_nmax=10000
 
 character(3) :: im
 integer :: inx0, nmax0, nx0, inxp, isa, i, isw
 real(8) :: tol0, atot, eps
 
 type(segment) :: s0, s1, s2
 type(segment), allocatable :: sa(:)
 
 
! Options
 if (present(iscale)) then
  im=iscale
 else
  im=df_iscale
 end if 
 
 if (present(inx)) then
  inx0=inx
 else
  inx0=df_inx
 end if
 inx0=inx0-1 !  number of intervals=number of points-1
 
 if (present(nmax)) then
  nmax0=nmax
 else
  nmax0=df_nmax
 end if
 nmax0=nmax0-1
 
 
 if (present(nx)) then
  nx0=nx
  inxp=1
 else
  nx0=100
  inxp=0
 end if
 nx0=nx0-1
 
 if (present(tol)) then
  tol0=tol
 else
  tol0=df_tol
 end if 
 
 if (inxp==1) then
  if (inx0>nx0) inx0=nx0/2
 end if 
 
 if (inx0>nmax0) inx0=nmax0/2
 
 
! Calculation
 if (im=='lin') then
  call init_lin(fun,inx0,x1,x2,pfx)
 else
  call init_log(fun,inx0,x1,x2,pfx)
 end if 
 
 atot=0d0
 eps=0d0
 isa=0
 allocate(sa(nmax0+2))
 
 do i=1,inx0
  call initseg(fun,pfx(i),pfx(i+1),s0)
  atot=atot+s0%a1+s0%a2
  eps=eps+s0%eps
  call addseg(s0,sa,isa)
 end do
 
 deallocate(pfx)
 
 if (inxp==1) then
 if (isa<nx0) then
  do
   call divseg(fun,sa(isa),s1,s2)
   atot=atot+s1%eps+s2%eps
   eps=eps+s1%eps+s2%eps-sa(isa)%eps
   isa=isa-1
   call addseg(s1,sa,isa)
   call addseg(s2,sa,isa)
   if (isa>=nx0) exit
  end do
 end if
 else
 if ((abs(eps/atot)>tol0).and.(isa<nmax0)) then
  do
   call divseg(fun,sa(isa),s1,s2)
   atot=atot+s1%eps+s2%eps
   eps=eps+s1%eps+s2%eps-sa(isa)%eps
   isa=isa-1
   call addseg(s1,sa,isa)
   call addseg(s2,sa,isa)
   if ((abs(eps/atot)<tol0).or.(isa>nmax0)) exit
  end do
 end if 
 end if
 
 
 
 
 do
  isw=0
  do i=1,isa-1
   if (sa(i)%p1%x>sa(i+1)%p1%x) then
     s0=sa(i)
     sa(i)=sa(i+1)
     sa(i+1)=s0
     isw=1
   end if
  end do
  if (isw==0) exit
 end do
 
 if (allocated(x)) deallocate(x)
 if (allocated(fx)) deallocate(fx)
 if (allocated(idx)) deallocate(idx)
 allocate(x(isa+1))
 allocate(fx(isa+1))
 allocate(idx(isa+1))
 
 do i=1,isa
  x(i)=sa(i)%p1%x
  fx(i)=sa(i)%p1%f
  idx(i)=sa(i)%p1%n
 end do
 x(isa+1)=sa(isa)%p2%x
 fx(isa+1)=sa(isa)%p2%f
 idx(isa+1)=sa(isa)%p2%n
 deallocate(sa)
 
 if (present(int_res)) then
  int_res=atot
 end if
 
 if (present(int_err)) then
  int_err=abs(eps/atot)
 end if
 
end subroutine looktab2d



subroutine addseg(s0,sa,isa)
 type(segment), intent(in) :: s0
 type(segment), allocatable :: sa(:)
 integer :: isa
 real(8) :: a, b, c, x
 integer :: i, j, ii
 integer :: n1, n2, nmid
 
 
 a=abs(sa(1)%eps)
 b=abs(sa(isa)%eps)
 x=abs(s0%eps)
 
 if (x>=b) then
  isa=isa+1
  sa(isa)=s0
  return
 end if
 
 if (x<=a) then
  ii=1
 else
  n1=1
  n2=isa
  do
   nmid=(n1+n2)/2
   if (nmid==n1) exit
   c=abs(sa(nmid)%eps) 
   if (x<c) then
    n2=nmid
   else
    n1=nmid
   end if  
  end do 
  ii=n2
 end if

 
 do i=1,isa-ii+1
  j=isa-i+1
  sa(j+1)=sa(j)
 end do
 sa(ii)=s0
 isa=isa+1

end subroutine addseg



subroutine init_lin(fun,nx,x1,x2,pfx)
 procedure(ifun1d) :: fun
 integer, intent(in) :: nx
 real(8), intent(in) :: x1, x2
 type(point), allocatable :: pfx(:)
 integer :: i
 
 if (allocated(pfx)) deallocate(pfx)
 allocate(pfx(nx+1))
 
 do i=0,nx  
  pfx(i+1)%x=x1+(x2-x1)*i/nx
  call fun(pfx(i+1)%x,pfx(i+1)%f,pfx(i+1)%n)
 end do
 
end subroutine init_lin


subroutine init_log(fun,nx,x1,x2,pfx)
 procedure(ifun1d) :: fun
 integer, intent(in) :: nx
 real(8), intent(in) :: x1, x2
 type(point), allocatable :: pfx(:)
 integer :: i
 
 if (allocated(pfx)) deallocate(pfx)
 allocate(pfx(nx+1))
 
 do i=0,nx  
  pfx(i+1)%x=x1*(x2/x1)**(i*1d0/nx)
  call fun(pfx(i+1)%x,pfx(i+1)%f,pfx(i+1)%n)
 end do
 
end subroutine init_log





subroutine initseg(fun,p1,p2,s0)
 procedure(ifun1d) :: fun
 type(point), intent(in) :: p1, p2
 type(segment), intent(out) :: s0
 real(8) :: atot, dx
 
 s0%p1=p1
 s0%p2=p2
 dx=(s0%p2%x-s0%p1%x)/2
 
 atot=(s0%p1%f+s0%p2%f)*dx
 
 s0%pm%x=(s0%p1%x+s0%p2%x)/2
 call fun(s0%pm%x,s0%pm%f,s0%pm%n)
 
 s0%dx=dx/2
 s0%a1=(s0%p1%f+s0%pm%f)*s0%dx
 s0%a2=(s0%p2%f+s0%pm%f)*s0%dx
 s0%eps=s0%a1+s0%a2-atot
end subroutine initseg



subroutine divseg(fun,s0,s1,s2)
 procedure(ifun1d) :: fun
 type(segment), intent(in) :: s0
 type(segment), intent(out) :: s1, s2

 call segval(fun,s1,s0%p1,s0%pm,s0%a1,s0%dx)
 call segval(fun,s2,s0%pm,s0%p2,s0%a2,s0%dx)
 
end subroutine divseg



subroutine segval(fun,s0,p1,p2,a0,dx)
 procedure(ifun1d) :: fun
 type(segment) :: s0
 type(point), intent(in) :: p1, p2
 real(8), intent(in) :: a0, dx
 s0%p1=p1
 s0%p2=p2
 s0%pm%x=p1%x+dx
 call fun(s0%pm%x,s0%pm%f,s0%pm%n)
 s0%dx=dx/2
 s0%a1=(s0%p1%f+s0%pm%f)*s0%dx
 s0%a2=(s0%p2%f+s0%pm%f)*s0%dx
 s0%eps=s0%a1+s0%a2-a0
end subroutine segval


end module lookup2d



!program main
! use lookup2d, only : looktab2d
! implicit none
! real(8) :: x1, x2, eps
! real(8), allocatable :: x(:), fx(:)
! integer, allocatable :: idx(:)
! integer :: i
! 
! interface
! subroutine fun(x,res,id)
!  real(8), intent(in) :: x
!  real(8), intent(out) :: res
!  integer, intent(out) :: id
! end subroutine fun
! end interface
! 
! x1=1d-2
! x2=50d0
! call looktab2d(fun,x1,x2,x,fx,idx,inx=10,int_err=eps)
! 
! write(*,*) eps, size(x)
! open(1,file='fx.dat')
! do i=1,size(x)
!  write(1,*) x(i), fx(i), idx(i)
! end do
! close(1)
! 
!
!end program main
!
!
!subroutine fun(x,res,id)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! integer, intent(out) :: id
! integer, save :: idc=0
!  
!!  res=x**8*exp(-(x))
!! res=x**2
!  res=x**8*exp(-(x))/(1d-3+(x-20)**2)
!  res=res/(1d-2+(x-10)**2)
!  idc=idc+1
!  id=idc
!! res=x**2*exp(-(x))/(1d-1+(x-20)**2)*exp(-1/x)
!! res=res/(0.5d0+(x-10)**2)
!
!end subroutine fun