module utools_mod
 implicit none
 private
 
 
 public :: utools
 
 type utools
 contains
  procedure, nopass :: grid, get_unit, bslash
  procedure, nopass :: unit_vec, cross_prod
  procedure, nopass :: uniform_reduction, frac_rep
 end type utools
 
 
 contains


subroutine grid(gx,x1,x2,nx,typ,nst)
!====================================
! create gx(nst:nx) from x1 to x2
! in log scale (default) or lin scale
! typ='lin'
!====================================
 real(8), allocatable :: gx(:)
 real(8), intent(in) :: x1, x2
 integer, optional, intent(in) :: nx
 character(3), optional, intent(in) :: typ
 integer, optional, intent(in) :: nst
 integer, parameter :: nst_def=1
 integer :: nn, i, n1
 real(8) :: dx
 character(3), parameter :: typ_def='log'
 character(3) :: typ0

 if (present(nst)) then
  n1=nst
 else
  n1=nst_def
 end if


 if (present(nx)) then
  if (allocated(gx)) then
   if (size(gx)/=(nx-n1+1)) then
    deallocate(gx)
    allocate(gx(n1:nx))
   end if
  else
   allocate(gx(n1:nx))
  end if  
  nn=nx
 else
  if (allocated(gx)) then
   nn=size(gx)
  else
   write(*,*) 'utools%grid: gx is not allocated and nx is not given'
   return
  end if
 end if
 
 if (present(typ)) then
  typ0=typ
 else 
  typ0=typ_def
 end if
 
 if (typ0=='log') then
  dx=(x2/x1)**(1d0/(nn-n1))
  gx(n1)=x1
  gx(nn)=x2
  do i=n1,nn-2
   gx(i+1)=x1*dx**(i-n1+1)
  end do
  return
 end if
 
 if (typ0=='lin') then
  dx=(x2-x1)/(nn-n1)
  gx(n1)=x1
  gx(nn)=x2
  do i=n1,nn-2
   gx(i+1)=x1+dx*(i-n1+1)
  end do
  return
 end if 
 

end subroutine grid



subroutine get_unit(unit)
!====================================
! Search for free 'unit' from (1..99)
!====================================
! Variables
 integer, intent(out) :: unit
 integer :: i
! Calculations
 unit=0
 do i=1,99
  if ((i/=5).and.(i/=6).and.(i/=9)) then
   open(unit=i,err=10,status='scratch')
   close(unit=i)
   unit=i
   return
  end if
10  continue
 end do
end subroutine get_unit


subroutine bslash(ns,str)
!==========================================================
! Conversion of the string 'str' of length ns to the 
! string without escape characters
!==========================================================
 integer, intent(in) :: ns
 character(ns) :: str, str1, str2
 integer, parameter :: nsym=8
 character(3), parameter :: symb(nsym)=['\a','\b','\f','\n','\r','\t','\v','\0']
 character(3), parameter :: symb1(nsym)=['\\a','\\b','\\f','\\n','\\r','\\t','\\v','\\0']
 integer :: nn, i, j
 
 nn=len_trim(str)
 i=0
 do
  i=i+1
  do j=1,nsym
   if (str(i:i)==symb(j)) then
    str1=trim(str(:i-1))
	str2=trim(str(i+1:))
	str=trim(str1)//trim(symb1(j))//trim(str2)
	i=i+1
	nn=nn+1
   end if
  end do
  if ((i==nn).or.(i==ns)) exit
 end do
 
end subroutine bslash


subroutine unit_vec(phi,theta,nv)
 use phys_const, only :  pi
 real(8), parameter :: deg=pi/180
 real(8), intent(in) :: phi, theta
 real(8), intent(out) :: nv(3)
 real(8) :: phi1, theta1
 
 phi1=phi*deg
 theta1=theta*deg
  
 nv(1)=cos(phi1)*sin(theta1)
 nv(2)=sin(phi1)*sin(theta1)
 nv(3)=cos(theta1)
end subroutine unit_vec

subroutine cross_prod(a,b,res)
! Cross product
 real(8), intent(in) :: a(3), b(3)
 real(8), intent(out) :: res(3)
 
 res(1)=a(2)*b(3)-a(3)*b(2) 
 res(2)=a(3)*b(1)-a(1)*b(3)
 res(3)=a(1)*b(2)-a(2)*b(1)
end subroutine cross_prod



subroutine uniform_reduction(nold,nnew,anum)
!================================================================
! It is auxiliary function for reduction of large arrays.
! Instead of all "nold" points we need choose only "nnew" points
! in uniform way, i.e. each part of array should be represented
! equally. It is done by using expansion of the fraction
! nold/nsec=n(0)+1/n(1)+1/n(2)+1/n(3)+1/n(4)+...
! Then we should take only every n(0) elements adding or 
! subtracting unity every n(1), n(2) ... similar like in the
! gregorian year = 365+1/4-1/100+1/400 we add one day per 4 years,
! not add 1700, 1800, 1900 but still add in 1600 and 2000 year.
! Result: anum(i), i.e.
! if we want to reduce array a(i=1,nold) to b(i=1,nnew), then
! b(i)=a(anum(i))
!================================================================
 integer, intent(in) :: nold, nnew
 integer, allocatable, intent(out)  :: anum(:)
 integer, allocatable :: nfrac(:)
 integer :: nsec, jump0, jump, i, j, k, nfracs
 integer :: fc1, fc2 ! fraction counters
 integer, parameter :: uunit=1 ! just unit
 integer :: nnew1

 if (allocated(anum)) deallocate(anum)
 allocate(anum(nnew))

 anum(1)=1         ! include first
 anum(nnew)=nold   ! and last points of array
 
 nsec=nnew-1       ! nsec is number of sections between nnew points, e.g. for 3 points there are 2 sections

 call frac_rep(nold,nsec,nfrac) ! nold/nsec=n(0)+1/n(1)+1/n(2)+1/n(3)+1/n(4)+...


! Setting basic jump

 nfracs=size(nfrac)
 if (nfracs>1) then
  if (nfrac(1)==1) then
   jump0=nfrac(0)+1
   fc1=2
  else
   jump0=nfrac(0)
   fc1=1
  end if

  fc2=nfracs-1
  do i=fc1,fc2
   if (abs(nfrac(i))>nsec) exit ! if nfrac(i)>nsec, then accounting 
  end do                        ! of the addition 1/nfrac(i) will never be realized
  fc2=i-1
  if (fc2<fc1) fc2=fc1
 else
  jump0=nfrac(0)
 end if


 i=1 ! i=1 and j=1 are already accounted in anum(1) 
 j=1

 nnew1=nnew-1

 do
  if (j>=nnew1) exit
  jump=jump0
  if (nfracs>1) then
   do k=fc1,fc2
    if (mod(j,nfrac(k))==0) jump=jump+sign(uunit,nfrac(k)) ! change of step
   end do
  end if

  i=i+jump
  if (i>=nold) then
   do 
    j=j+1
    anum(j)=anum(j-1)+1
    if (j>=nnew1) exit
   end do
   exit
  end if

  j=j+1
  anum(j)=i
   
 end do

! write(*,*) 'anum(1)=', anum(1)
! write(*,*) 'anum(end)=', anum(j)

! do i=1,nnew
!  write(*,*) 'i, anum(i)=', i, anum(i)
!  read(*,*)
! end do

! read(*,*)

 deallocate(nfrac)

! fixing, if there is equal elements

! do i=3,nnew
!  if (anum(i-1)>=anum(i)) then
!   anum(i-1)=floor((anum(i)+anum(i-2))/2d0)
!  end if
! end do


end subroutine uniform_reduction



subroutine frac_rep(num,den,nres)
!==============================================
! Representation of the fraction in 
! the form: 
! num/den=n(0)+1/n(1)+1/n(2)+1/n(3)+1/n(4)+...
! In this notation n(i)=nres(i)
!==============================================
 integer, intent(in) :: num, den
 integer, allocatable, intent(out) :: nres(:)
 integer, parameter :: ntot=100 ! actually, it is much smaller (typically smaller than 10)
 integer :: delta(0:ntot), rem(0:ntot), nw(0:ntot)
 integer :: i, j

 i=0
 nw(i)=num/den
 rem(i)=mod(num,den)
 delta(i)=den
 

 if (rem(0)>0) then
  do
   i=i+1
   nw(i)=delta(i-1)/rem(i-1)  
   rem(i)=mod(delta(i-1),rem(i-1))
   delta(i)=nw(i)*(rem(i)+nw(i)*rem(i-1))

   if (delta(i)<0) exit ! because of overflow 
   if (rem(i)==0) exit  ! because of end of process
   if (i>ntot) exit     ! should never happen 
  end do

  if (allocated(nres)) deallocate(nres)
  allocate(nres(0:i))
 
  nres(0)=nw(0)
  do j=1,i
   nres(j)=-((-1)**j)*nw(j)
  end do
 else
  if (allocated(nres)) deallocate(nres)
  allocate(nres(0:0))
 
  nres(0)=nw(0)
 
 end if

end subroutine frac_rep




end module utools_mod

! program main
!  use utools_mod, only : utools
!  type(utools) :: ut
!  real(8), allocatable :: x(:)
!  integer :: nn
! 
!  call ut%grid(x,1d0,100d0,10,'log',-5)
!  
!  write(*,*) size(x)
!  do i=-5,10
!   write(*,*) i, x(i)
!  end do

! 
! 
!  call ut%get_unit(nn)
!  write(*,*) nn

! end program main

