module fsamplem
!====================================================================
! Class 'fsample' samples a function of one argument with controlled
! accuracy. Output are x(:), f(:) arrays, where x(:) are points
! and f(:) are values at these points.
!
! Description:
! ===========
! use fsamplem, only : fsample
! type(fsample) :: a
!# 1
! call a%fun(fun,x1,x2) - setting function 'fun' in the range 
! from  x1 to x2. 'fun' in the form:
! subroutine fun(x,res)
!  real(8), intent(in) :: x
!  real(8), intent(out) :: res
!   ...
! end subroutine fun
!# 2
! call a%igrid(npoint,iscale) - setting initial grid
! iscale='log' - logarithmic scale
! iscale='lin' - linear scale
! for precalculation of npoints
! both npoint and isclale are optional parameters
! if they are not given then:
! npoint=50, iscale='log' are default values
!# 3
! call a%want(nsample,relerr,maxfun) - setting of calculations
! nsample is the total number of points which will be output,
! irrespective of relerr - relative error.
! if nsample (which is optinal, by default =100) is not given, then
! calculations will continue up to reaching whether relerr 
! (which is optional,by default =1d-4) or maxfun (by default =10000).
! maxfun is the maximum number of evaluations of fun.
!# 4
! call a%sample(x,f,call_order) - calculation and output of sampling
! integer :: call_order(:) is array with serial number of calling 
! function fun in the process of sampling
!# 5
! call a%print_status - auxiliary function, which prints the status
! of sampling (not sampled, reached maxfun, reached relerr or 
! reached required nsample)
!# 6
! call a%clean - clean variable and calculation returning to default
! parameters
!
! Usage:
! =====
! use fsamplem, only : fsample
! type(fsample) :: a
!
! call a%fun(fun,1d-2,50d0)
! call a%igrid(100,'log')
! call a%want(300)      ! or call a%want(relerr=1d-5,maxfun=20000)
! call a%sample(x,f,nc)
! write(*,*) a%error, a%integral, a%relerr, a%nsample, a%nfun
! call a%print_status
! call a%clean
!
!! a%error is total error of sampling (error of integral)
!! a%integral is value of integral calculated using sampled function
!! a%relerr is reached relative error
!! a%nsample is number of sampled points
!! a%nfun is total number of function evaluations
!
! Test program is in the end of the file
!====================================================================
 implicit none
 private
 save
 
 public :: fsample
 
 abstract interface
  subroutine fun1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun1d
 end interface


 type segdat
  integer :: i
  real(8) :: a1, a2, dx, eps
 end type segdat
 

 type fsample
  real(8), allocatable :: x(:), f(:)
  integer, allocatable :: next(:), nseg(:)
  integer(1), allocatable :: nip(:)
  type(segdat), allocatable :: sg(:)
  integer :: nc, isg
  integer(1) :: nxp=0, excode=0
  real(8) :: x1, x2, eps, atot
  procedure(fun1d), pointer, nopass :: fx=>null()
  integer :: nmax=10000, inx=50, nx=100    ! default parameters
  real(8) :: tol=1d-4                      ! default parameters
  character(3) :: iscale='log'             ! default parameter
  real(8) :: error, integral, relerr       ! interface values
  integer :: nsample, nfun  
 contains
  procedure :: fun=>setfun, igrid=>set_igrid, want
  procedure :: print_status=>cstatus, sample=>samplefun, clean=>cleandat
  procedure :: calc, outfx, outnc
  procedure :: del, set 
  procedure :: fcall, fcall0
  procedure :: stseg, initsamp, refine
 end type fsample
 

 contains
 
! Interface functions

subroutine setfun(this,fun,x1,x2)
 class(fsample) :: this
 procedure(fun1d) :: fun
 real(8), intent(in) :: x1, x2
 
 this%fx=>fun
 this%x1=x1
 this%x2=x2
 
end subroutine setfun


subroutine set_igrid(this,npoint,iscale)
 class(fsample) :: this
 integer, intent(in), optional :: npoint
 character(3), intent(in), optional :: iscale
 
 if (present(npoint)) this%inx=npoint
 if (present(iscale)) this%iscale=iscale

end subroutine set_igrid



subroutine want(this,nsample,relerr,maxfun)
 class(fsample) :: this
 integer, intent(in), optional :: nsample
 real(8), intent(in), optional :: relerr
 integer, intent(in), optional :: maxfun
 
 if (present(nsample)) then
  this%nxp=1
  this%nx=nsample
 end if
 if (present(relerr)) this%tol=abs(relerr)
 if (present(maxfun)) this%nmax=maxfun

end subroutine want
 
subroutine cstatus(this)
!
! print status of calculation
!
 class(fsample) :: this
 
 if (this%excode==0) then
  write(*,*) 'fsample: Calculation has not been executed'
 end if
 
 if (this%excode==1) then
  write(*,*) 'fsample: Required accuracy is attained'
 end if
 
 if (this%excode==2) then
  write(*,*) 'fsample: Required number of values is calculated'
 end if

 if (this%excode==3) then
  write(*,*) 'fsample: Calculation has been interrupted because &
          the maximum number of function evaluations is reached'
 end if
end subroutine cstatus


subroutine cleandat(this)
!
! Set to default
!
 class(fsample) :: this  
 
 call this%del
 this%nxp=0
 this%excode=0
 this%nmax=10000
 this%inx=50 
 this%nx=100
 this%tol=1d-4
 this%iscale='log'
 this%fx=>null()  
 
end subroutine cleandat  
 
 

subroutine samplefun(this,x,f,call_order)
 class(fsample) :: this
 real(8), allocatable :: x(:), f(:)
 integer, allocatable, optional :: call_order(:)
! Consistency check
 
 if (this%nxp==1) then ! nsample is given
  if (.not.((this%inx<this%nx).and.(this%nx<this%nmax))) then
   write(*,*) 'fsample: npoint<nsample<maxfun is not fulfilled'
   return
  end if 
 else
  if (.not.(this%inx<this%nmax)) then
   write(*,*) 'fsample: npoint<maxfun is not fulfilled'
   return
  end if 
 end if 

! Calculations 
 this%nx=this%nx-1 ! number of intervals instead of number of points
 call this%calc
 this%nmax=this%nmax+1 ! because it was changed in this%calc
 call this%outfx(x,f)
 
 this%error=this%eps
 this%integral=this%atot
 this%relerr=abs(this%eps/this%atot)
 this%nsample=this%isg+1
 this%nfun=this%nc
  
 if (present(call_order)) call this%outnc(call_order)
end subroutine samplefun
 
 
! End of interface functions 
 
 
subroutine calc(this)
!
! Main calculation
!
 use utools_mod, only : utools
 class(fsample) :: this
 type(utools) :: ut
 real(8), allocatable :: x(:)
 real(8) :: rel
 integer :: ocode
 
 
 call this%set(this%nmax)
 call ut%grid(x,this%x1,this%x2,this%inx,this%iscale)
 this%nmax=this%nmax-1 ! for not overflowing
 call this%initsamp(x)
 deallocate(x)
 
 if (this%nxp==1) then
  do
   call this%refine
   
   if (this%isg>=this%nx) then
    this%excode=2
    exit
   end if
   
   if (this%nc>=this%nmax) then
    this%excode=3
    exit
   end if
  
  end do
 else
  do
   call this%refine
   rel=abs(this%eps/this%atot)
   
   if (rel<=this%tol) then
    this%excode=1
    exit
   end if	
   
   if (this%nc>=this%nmax) then
    this%excode=3
    exit
   end if
   
 end do
 
 end if
 
end subroutine calc


subroutine outfx(this,x,f)
!
! Order and output the calculated points
!
 class(fsample) :: this
 real(8), allocatable :: x(:), f(:)
 integer :: ic, i, nt, ns, i1
 
 if (allocated(x)) deallocate(x)
 if (allocated(f)) deallocate(f)
 
 nt=this%nc
 ns=this%isg+1
 
 allocate(x(ns))
 allocate(f(ns))
 
 ic=1
 i1=1
 do i=1,nt
  if (this%nip(ic)==1) then
   x(i1)=this%x(ic)
   f(i1)=this%f(ic)
   i1=i1+1
  end if 
  ic=this%next(ic)
 end do

end subroutine outfx


subroutine outnc(this,nc)
!
! Output array with calling order of calculated points
!
 class(fsample) :: this
 integer, allocatable :: nc(:)
 integer :: i, ic, nt, ns, i1
 
 if (allocated(nc)) deallocate(nc)
 nt=this%nc
 ns=this%isg+1
 allocate(nc(ns))
 ic=1
 i1=1
 do i=1,nt
  if (this%nip(ic)==1) then
   nc(i1)=ic
   i1=i1+1
  end if 
  ic=this%next(ic)
 end do
end subroutine outnc
  

subroutine del(this)
 class(fsample) :: this
 if (allocated(this%x)) deallocate(this%x)
 if (allocated(this%f)) deallocate(this%f)
 if (allocated(this%next)) deallocate(this%next)
 if (allocated(this%nseg)) deallocate(this%nseg)
 if (allocated(this%nip)) deallocate(this%nip)
 if (allocated(this%sg)) deallocate(this%sg)
end subroutine del

subroutine set(this,n)
 class(fsample) :: this
 integer, intent(in) :: n
 
 call this%del
 allocate(this%x(n))
 allocate(this%f(n))
 allocate(this%next(n))
 allocate(this%nseg(n))
 allocate(this%nip(n))
 allocate(this%sg(n))
 
end subroutine set

 

subroutine fcall0(this,x)
 class(fsample) :: this
 real(8), intent(in) :: x
 real(8) :: res
 integer :: i
 
 call this%fx(x,res)
 this%nc=this%nc+1
 i=this%nc
 this%x(i)=x
 this%f(i)=res
end subroutine fcall0


 
subroutine fcall(this,x,res,i)
!
! Call sampled function and count number of calls (in this%nc)
!
 class(fsample) :: this
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 integer, intent(out) :: i
 
 call this%fx(x,res)
 this%nc=this%nc+1
 i=this%nc
 this%x(i)=x
 this%f(i)=res
end subroutine fcall


 
subroutine initsamp(this,x)
 class(fsample) :: this
 real(8) :: x(:)
 integer :: i, nst
 
 this%atot=0d0
 this%eps=0d0
 this%nc=0
 this%isg=0
 this%nip=0
 nst=size(x)
 do i=1,nst
  call this%fcall0(x(i))
  this%next(i)=i+1
 end do
 this%next(nst)=1
 this%nip(nst)=1

 call this%stseg(nst-1)
end subroutine initsamp


subroutine stseg(this,n)
 class(fsample) :: this
 integer, intent(in) :: n
 integer :: i, nput
 type(segdat) :: s0
 
 i=1
 call initseg(this,i,s0)
 this%atot=this%atot+s0%a1+s0%a2
 this%eps=this%eps+s0%eps
 this%isg=1
 this%sg(1)=s0
 this%nseg(1)=1
 
 
 do i=2,n
  call initseg(this,i,s0)
  this%atot=this%atot+s0%a1+s0%a2
  this%eps=this%eps+s0%eps
  nput=this%isg+1
  call addseg(s0,nput,this%isg,this%nseg,this%sg)
 end do

end subroutine stseg



subroutine initseg(this,i,s0)
 class(fsample) :: this
 integer, intent(in) :: i
 type(segdat), intent(out) :: s0
 integer :: i1, i2
 real(8) :: dx, a0
 i1=i
 i2=i+1
 dx=(this%x(i2)-this%x(i1))/2
 a0=(this%f(i1)+this%f(i2))*dx
 call segval(this,s0,i1,i2,a0,dx)
end subroutine initseg


subroutine refine(this)
 class(fsample) :: this
 type(segdat) :: s1, s2
 integer :: nsg, nput
  
 nsg=this%nseg(this%isg)

 call divseg(this,this%sg(nsg),s1,s2)
 this%atot=this%atot+s1%eps+s2%eps
 this%eps=this%eps+s1%eps+s2%eps-this%sg(nsg)%eps
 this%isg=this%isg-1
 
 nput=nsg
 call addseg(s1,nput,this%isg,this%nseg,this%sg)
 nput=this%isg+1
 call addseg(s2,nput,this%isg,this%nseg,this%sg)
  
end subroutine refine
 
 
 
subroutine addseg(s0,nput,nmax,nsa,sa)
 type(segdat), intent(in) :: s0
 integer :: nmax, nput
 integer :: nsa(:)
 type(segdat) :: sa(:)
 real(8) :: a, b, c, x
 integer :: i, j, ii
 integer :: n1, n2, nmid
 
 
 sa(nput)=s0
 a=abs(sa(nsa(1))%eps)
 b=abs(sa(nsa(nmax))%eps)
 x=abs(s0%eps)
 
 if (x>=b) then
  nmax=nmax+1
  nsa(nmax)=nput
  return
 end if
 
 if (x<=a) then
  ii=1
 else
  n1=1
  n2=nmax
  do
   nmid=(n1+n2)/2
   if (nmid==n1) exit
   c=abs(sa(nsa(nmid))%eps) 
   if (x<c) then
    n2=nmid
   else
    n1=nmid
   end if  
  end do 
  ii=n2
 end if


 do i=1,nmax-ii+1
  j=nmax-i+1
  nsa(j+1)=nsa(j)
 end do
 nsa(ii)=nput
 nmax=nmax+1
 
end subroutine addseg

 

subroutine divseg(this,s0,s1,s2)
 class(fsample) :: this
 type(segdat), intent(in) :: s0
 type(segdat), intent(out) :: s1, s2
 integer :: i1, i2, i3
 
 i1=s0%i
 i2=this%next(i1)
 i3=this%next(i2)
 
 call segval(this,s1,i1,i2,s0%a1,s0%dx)
 call segval(this,s2,i2,i3,s0%a2,s0%dx)
 
end subroutine divseg

 
 
subroutine segval(this,s0,i1,i2,a0,dx)
 class(fsample) :: this
 type(segdat), intent(out) :: s0
 integer, intent(in) :: i1, i2
 real(8), intent(in) :: a0, dx
 integer :: i12
 real(8) :: x12, f12, dx2
 
  
 x12=this%x(i1)+dx
 call this%fcall(x12,f12,i12)
 this%nip(i1)=1
 this%next(i1)=i12 ! instead of i2
 this%next(i12)=i2
 dx2=dx/2
 
 s0%i=i1 
 s0%dx=dx2
 s0%a1=(this%f(i1)+f12)*dx2
 s0%a2=(this%f(i2)+f12)*dx2
 s0%eps=s0%a1+s0%a2-a0
end subroutine segval
 
 
 
end module fsamplem



! program main
 ! call test
! end program main



! subroutine test
 ! use fsamplem, only : fsample
 ! use utools_mod, only : utools
 ! type(fsample) :: a
 ! type(utools) :: ut
 ! real(8), allocatable :: x(:), f(:)
 ! integer, allocatable :: nc(:)
 ! integer :: i, nst, ic, nn
 ! real(8) :: res
 ! real :: p1, p2
 
 ! interface
  ! subroutine fund(x,res)
   ! real(8), intent(in) :: x
   ! real(8), intent(out) :: res
  ! end subroutine fund
 ! end interface
 

 ! call a%fun(fun,1d-2,50d0)
 ! call a%igrid(100,'log')
! ! call a%want(nsample=1000,maxfun=10000,relerr=1d-3)
 ! call a%want(relerr=1d-4)
 
 ! call cpu_time(p1) 
 ! call a%sample(x,f,nc)
 ! call cpu_time(p2)
 ! write(*,*) p2-p1
 ! call a%print_status
 ! write(*,*) a%error, a%integral, a%relerr, a%nsample, a%nfun
 ! call a%clean
 
 ! open(1,file='fun_samp.dat')
 ! do i=1,size(x)
  ! write(1,*) x(i), f(i), nc(i)
 ! end do
 ! close(1)
 
 ! call ut%grid(x,1d-2,50d0,300)
 ! open(1,file='fun_log.dat')
 ! do i=1,size(x)
  ! call fun(x(i),res)
  ! write(1,*) x(i), res
 ! end do
 ! close(1)

 ! end subroutine test


! subroutine fun(x,res)
 ! real(8), intent(in) :: x
 ! real(8), intent(out) :: res
 
 
 ! res=x**8*exp(-(x))/(1d-3+(x-20)**2)
 ! res=res/(1d-2+(x-10)**2)
 
! end subroutine fun
