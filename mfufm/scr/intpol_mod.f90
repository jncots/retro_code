module intpol_mod
!================================================================
! History of versions:
! v2: change in arr_ind_short: addition of tolerance for out of
! boundary range (see the part with 'tolerance' variables in use)
!
! 01.04.2016 Add 'arr_ind_breal' - real version of arr_ind_b,
! which is real(8), and 'lin_int_real' - real version of lin_int
! 21.04.2016 Add simple bisectional root finding procedure
! bs_root
!================================================================
 implicit none
 private
 save
 
 
 public :: arr_ind_b, arr_ind_c, log_int, lin_int
 public :: zeros_cut, arr_ind_short, remap
 public :: plaw_index, get_grid, arr_ind_breal, lin_int_real
 public :: bs_root, logl_int
 
 contains




subroutine bs_root(x1,x2,fun,res,err_code)
! Variables
 real(8), intent(in) :: x1, x2
 real(8), intent(out) :: res
 interface
   subroutine fun(x,res)
    real(8), intent(in) :: x
    real(8), intent(out) :: res
   end subroutine
 end interface
 integer, intent(out) :: err_code
 
 integer :: i
 real(8) :: a, b, c
 real(8) :: fa, fb, fc, fc0
! Calculations
 
 err_code=0
 a=x1
 b=x2
 call fun(a,fa)
 call fun(b,fb)
 
 fc0=fa
 if (fa*fb<0d0) then
  i=0
  do 
   c=(a+b)/2
   fc0=fc
   call fun(c,fc)
   if ((abs((fc-fc0)/fc)<1d-15).or.(fc<1d-16)) exit
   if (fc*fa<0) then
    b=c
    fb=fc
   else
    a=c
    fa=fc
   end if
   i=i+1
   if (i>1000) then
    write(*,*) 'bs_root: limit is exceeded, err=', abs((fc-fc0)/fc)
    exit
   end if
  end do
 else
  write(*,*) 'bs_root: f(x1) and f(x2) are of the same sign'
  err_code=1
 end if
 res=c
end subroutine bs_root


 
subroutine arr_ind_breal(n,array,value,n1,n2)
!
! Binary search of indices n1 and n2 of two elements
! between which the 'value' is located
!
! Variables
 integer, intent(in) :: n
 real, intent(in) :: array(n), value
 integer, intent(out) :: n1, n2
 integer :: nmid
 real :: a, b, c
! Calculations

 n1=1
 n2=n
  
 a=array(n1)
 b=array(n2)

 if (value<a) then
  n2=n1+1
  write(*,*) 'arr_ind_breal: value<array(1)'
  return
 end if
 
 if (value>b) then
  n1=n2-1
  write(*,*) 'arr_ind_breal: value>array(nmax)'
  return
 end if
 
 do
  nmid=(n1+n2)/2
  if (nmid==n1) return
  c=array(nmid)
  
  if (value<c) then
   n2=nmid
  else
   n1=nmid
  end if  
 end do  

end subroutine arr_ind_breal 
 

subroutine lin_int_real(x1,x2,y1,y2,x,res)
!
! Linear interpolation in log scale
!
 real, intent(in) :: x1, x2, y1, y2, x
 real, intent(out) :: res
 
 res=y1+(y2-y1)*(x-x1)/(x2-x1)
 
end subroutine lin_int_real

 
 

 
subroutine arr_ind_b(n,array,value,n1,n2)
!
! Binary search of indices n1 and n2 of two elements
! between which the 'value' is located
!
! Variables
 integer, intent(in) :: n
 real(8), intent(in) :: array(n), value
 integer, intent(out) :: n1, n2
 integer :: nmid
 real(8) :: a, b, c
! Calculations

 n1=1
 n2=n
  
 a=array(n1)
 b=array(n2)

 if (value<a) then
  n2=n1+1
  write(*,*) 'arr_ind_b: value<array(1)'
  return
 end if
 
 if (value>b) then
  n1=n2-1
  write(*,*) 'arr_ind_b: value>array(nmax)'
  return
 end if
 
 do
  nmid=(n1+n2)/2
  if (nmid==n1) return
  c=array(nmid)
  
  if (value<c) then
   n2=nmid
  else
   n1=nmid
  end if  
 end do  

end subroutine arr_ind_b



subroutine arr_ind_c(n,array,value,number)
!
! Binary search of index 'number' of the element closest to 'value'
!
! Variables
 integer, intent(in) :: n
 real(8), intent(in) :: array(n), value
 integer, intent(out) :: number
 integer :: nmin, nmax, nmid
 real(8) :: a, b, c
! Calculations

 nmin=1
 nmax=n
  
 a=array(nmin)
 b=array(nmax)

 if (value<a) then
  number=1
  write(*,*) 'arr_ind_c: value<array(1)'
  return
 end if
 
 if (value>b) then
  number=n
  write(*,*) 'arr_ind_c: value>array(nmax)'
  return
 end if
 
 
 do
  nmid=(nmax+nmin)/2
  
  if (nmid==nmin) then  
   a=array(nmin)
   b=array(nmax)
   
   if ((value-a)>(b-value)) then
    number=nmax
   else
    number=nmin
   end if
   
   return
  end if

  c=array(nmid)
  
  if (value<c) then
   nmax=nmid
  else
   nmin=nmid
  end if
  
 end do  

end subroutine arr_ind_c



subroutine log_int(x1,x2,y1,y2,x,res)
!
! Linear interpolation in log scale
!
 real(8), intent(in) :: x1, x2, y1, y2, x
 real(8), intent(out) :: res
 real(8) :: xp 
 
 xp=log(x/x1)/log(x2/x1)
 res=y1*(y2/y1)**xp
 
end subroutine log_int




subroutine logl_int(x1,x2,y1,y2,x,res)
!
! Linear interpolation in log scale
!
 real(8), intent(in) :: x1, x2, y1, y2, x
 real(8), intent(out) :: res
 real(8) :: xp 
 
 xp=log(x/x1)/log(x2/x1)
 res=y1*(y2/y1)**xp
 if (isnan(res)) res=y1+(y2-y1)*(x-x1)/(x2-x1)
 
end subroutine logl_int




subroutine lin_int(x1,x2,y1,y2,x,res)
!
! Linear interpolation in log scale
!
 real(8), intent(in) :: x1, x2, y1, y2, x
 real(8), intent(out) :: res
 
 res=y1+(y2-y1)*(x-x1)/(x2-x1)
 
end subroutine lin_int



subroutine zeros_cut(fdat,nmin,nmax)
!
! Find zeros (fdat<zero=1d-307) in the begining and the
! end of the array fdat and gives the range [nmin,nmax],
! where fdat>zero, although zeros can be still presented
! inside of this range 
!
 real(8), intent(in) :: fdat(:)
 integer, intent(out) :: nmin, nmax
 real(8), parameter :: zero=1d-307
 integer :: n, i, j
 
 n=size(fdat)
 nmin=0
 nmax=0
 
 do i=1,n
  if (abs(fdat(i))>zero) then
   nmin=i
   exit
  end if
 end do
 
 if (nmin==0) return
 
 do i=0,n-1
  j=n-i 
   if (abs(fdat(j))>zero) then
   nmax=j
   exit
  end if
 end do
end subroutine zeros_cut





subroutine arr_ind_short(nmin,nmax,xa,x,n1,n2,nst)
!
! Binary search of indices n1 and n2 of two elements
! between which the 'x' is located.
! Search is performed between nmin and nmax.
!
 integer, intent(in) :: nmin, nmax
 real(8), intent(in), target :: xa(:)
 real(8), intent(in) :: x
 integer, intent(out) :: n1, n2
 integer, optional, intent(in) :: nst
 real(8), pointer :: a, b, c
 integer :: nmid
 real(8) :: ob
 real(8), parameter :: tol=1d-14 ! usually 1d-15 for real(8), 1d-14 is taken with a margin
 ! tol is relative error of how far 'x' can go beyond defined range [a,b]
 real(8), parameter :: zero=1d-307
 integer :: i, nstart

 if (present(nst)) then
   nstart=nst-1
 else
   nstart=0
 end if
 

 n1=nmin
 n2=nmax

 a=>xa(n1-nstart)
 b=>xa(n2-nstart)


 if (x<a) then
! Output of left border information
  n2=n1+1
   
  if (abs(a)>zero) then
   ob=(a-x)/a
  else 
   ob=a-x
  end if 
! Warning of exceeding of tolerance   
  if (ob<tol) write(*,'(A,10Es14.6)') 'arr_ind_short: value<array(1)', x, a, ob
  return
 end if
 
 if (x>b) then
! Output of right border information
  n1=n2-1

  if (abs(b)>zero) then
   ob=(x-b)/b
  else 
   ob=x-b
  end if 
! Warning of exceeding of tolerance 
  if (ob>tol) write(*,'(A,10Es14.6)') 'arr_ind_short: value>array(nmax)', x, b, ob 
  return
 end if 

 do
  nmid=(n1+n2)/2
  if (nmid==n1) return
  c=>xa(nmid-nstart)
  
  if (x<c) then
   n2=nmid
  else
   n1=nmid
  end if
 end do 

end subroutine arr_ind_short




subroutine remap(x1,f1,x2,f2)
!
! Remaping of function defined at point x1 to
! interpolated values at points x2
!
! Variables
 real(8), intent(in) :: x1(:), f1(:)
 real(8), intent(out) :: x2(:), f2(:)
 integer :: nx1, nx2, nlow, nhigh, n1, n2, i
 real(8) :: xmin1, xmax1, xmin2, xmax2, vx
! Calculations

 nx1=size(x1)
 xmin1=x1(1)
 xmax1=x1(nx1)
 
 nx2=size(x2)
 xmin2=x2(1)
 xmax2=x2(nx2)
 
 f2=0d0
 nlow=1
 nhigh=nx2
! if [x2(1):x2(nx2)] greater than [x1(1):x1(nx1)]
 if (xmin2<xmin1) call arr_ind_b(nx2,x2,xmin1,n1,nlow)
 if (xmax2>xmax1) call arr_ind_b(nx2,x2,xmax1,nhigh,n2)
 !write(*,'(I5,10Es14.6)') nx1, xmin1, xmax1
 !write(*,'(I5,10Es14.6)') nx2, xmin2, xmax2
 !write(*,*) nlow, nhigh
 !read(*,*)
 
 
 
 do i=nlow,nhigh
  
  vx=x2(i)
  call arr_ind_b(nx1,x1,vx,n1,n2)
  call log_int(x1(n1),x1(n2),f1(n1),f1(n2),vx,f2(i))
  if (isnan(f2(i))) call lin_int(x1(n1),x1(n2),f1(n1),f1(n2),vx,f2(i))
  
 end do

end subroutine remap


subroutine plaw_index(xarg,fun,xarg_a,fun_a)
!
! Calculation of power-law index d ln(f(x))/d ln(x) using central finite difference
!
! Variables
 real(8), intent(in) :: xarg(:), fun(:)
 real(8), intent(out) :: xarg_a(:), fun_a(:)
 integer :: i, ns
 real(8) :: res, x0, f0, xx
 real(8), parameter :: zero=1d-307
! Calculations

! Check of input data
 
 if (size(xarg)/=size(fun)) then
  write(*,*) 'plaw_index: size(xarg)/=size(fun)'
  write(*,*) 'size(xarg)=', size(xarg)
  write(*,*) 'size(fun)=', size(fun)
  return
 end if 
 
 ns=size(xarg)-1
 
 do i=1,ns
  x0=(xarg(i)+xarg(i+1))/2
  xarg_a(i)=x0
  if ((abs(fun(i))<zero).or.(abs(fun(i+1))<zero)) then
   xx=0d0
   f0=(xx-xx)/(xx-xx)
  else
   res=(fun(i+1)-fun(i))/(xarg(i+1)-xarg(i))
   f0=fun(i)+res*(x0-xarg(i))
  end if
  fun_a(i)=res*x0/f0
 end do
 

end subroutine plaw_index


subroutine get_grid(typ,nx,x1,x2,gridx)
!
! Produce of array gridx(0:nx) with values distibuted
! lineary (typ='lin') or logarithmically (typ='log')
! in the interval [x1,x2]
!
! Variables
 character(3) :: typ
 integer, intent(in) :: nx
 real(8) :: x1, x2
 real(8), allocatable, intent(out) :: gridx(:)
 integer :: i
 real(8) :: dx
! Calculations
 if (allocated(gridx)) deallocate(gridx)
 allocate(gridx(0:nx))
 if (typ=='lin') then
  dx=(x2-x1)/nx
  do i=0,nx
   gridx(i)=x1+dx*i
  end do
  return
 end if

 if (typ=='log') then
  dx=(x2/x1)**(1d0/nx)
  gridx(0)=x1
  gridx(nx)=x2
  do i=1,nx-1
   gridx(i)=x1*dx**i
  end do
  return
 end if
end subroutine get_grid



end module intpol_mod
