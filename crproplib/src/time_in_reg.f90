module time_in_reg_m
!===============================================================
! Module contains the class time_in_reg, which records 
! the time spent in the specified region
!===============================================================
 use ode_solver_m, only : ode_solution, ode_solver
 implicit none
 private
 save


 public :: time_in_reg

 abstract interface
  function inreg_fun(x)
   real(8), intent(in) :: x(3)
   logical :: inreg_fun
  end function inreg_fun
 end interface


 type seg_inside
  real(8) :: dt, p1(7), p2(7)
 end type seg_inside


! Class for the record of the entry and exiting the region
 type time_in_reg
  integer :: narr=0, nseg=0
  real(8) :: ttime=0
  type(seg_inside), allocatable :: seg(:)
 contains
  procedure :: del=>del_td, set=>set_td, double=>double_td
  procedure :: enter_pt, exit_pt
  procedure :: wander_time, wander_time_border
 end type time_in_reg



 contains


subroutine del_td(this)
 class(time_in_reg) :: this

 this%narr=0
 this%nseg=0
 this%ttime=0
 if (allocated(this%seg)) deallocate(this%seg)

end subroutine del_td

subroutine set_td(this,narr)
 class(time_in_reg) :: this
 integer, intent(in), optional :: narr
 integer :: narr0

 if (present(narr)) then
  narr0=narr
 else
  narr0=5
 end if
 
 call this%del
 allocate(this%seg(narr0))
 this%narr=narr0

end subroutine set_td

subroutine double_td(this)
 class(time_in_reg) :: this
 type(time_in_reg) :: that
 integer :: i, narr2

 narr2=2*this%narr

 call that%set(this%narr)
 that%nseg=this%nseg
 that%ttime=this%ttime

 do i=1,this%nseg
  that%seg(i)=this%seg(i)
 end do

 call this%set(narr2)
 this%nseg=that%nseg
 this%ttime=that%ttime

 do i=1,that%nseg
  this%seg(i)=that%seg(i)
 end do
 
 call that%del

end subroutine double_td


subroutine enter_pt(this,pt)
! Writes down entry point
 class(time_in_reg) :: this
 real(8) :: pt(7)
 integer :: nss

 if (this%narr==0) call this%set
 nss=this%nseg+1
 if (nss>this%narr) call this%double

 this%nseg=nss
 this%seg(nss)%p1=pt
 this%seg(nss)%dt=0d0

end subroutine enter_pt


subroutine exit_pt(this,pt)
! Writes down exit point
 class(time_in_reg) :: this
 real(8) :: pt(7)
 integer :: nss

 nss=this%nseg
 this%seg(nss)%p2=pt
 this%seg(nss)%dt=this%seg(nss)%p2(1)-this%seg(nss)%p1(1)
 this%ttime=this%ttime+this%seg(nss)%dt

end subroutine exit_pt




subroutine wander_time(this,path,inreg)
!==============================================================
! Calculates the time spent in the region "inreg",
! recording entrance and exit points and time between 
! these events. Do not take into account crossing of 
! the region. For more accurate treatment of the border 
! crossing see "wander_time_border" function
!==============================================================
 class(time_in_reg) :: this
 type(ode_solution), intent(in) :: path
 procedure(inreg_fun) :: inreg
 logical :: prev_in, now_in
 integer :: i
 real(8) :: xc(3)

 xc=path%y(1:3,1)

 call this%del
 
 if (inreg(xc)) then
  prev_in=.true.
  call this%enter_pt(path%y(0:6,1))
 else
  prev_in=.false.
 end if

 do i=2,path%np-1
  xc=path%y(1:3,i)
  if (inreg(xc)) then
   now_in=.true.
   if (.not.prev_in) call this%enter_pt(path%y(0:6,i))
  else
   now_in=.false.
   if (prev_in) call this%exit_pt(path%y(0:6,i-1))
  end if
  prev_in=now_in
 end do


 i=path%np
 xc=path%y(1:3,i)
 
 if (inreg(xc)) then
  if (prev_in) call this%exit_pt(path%y(0:6,i))
 else
  if (prev_in) call this%exit_pt(path%y(0:6,i-1))
 end if  

end subroutine wander_time



subroutine wander_time_border(this,path,inreg,lmin)
!==============================================================
! More accurace version of function "wander_time"
!
! Calculates the time spent in the region "inreg",
! recording entrance and exit points and time between 
! these events, taking into account accurate calculation
! of border crossing with accuracy lmin.
!==============================================================
 class(time_in_reg) :: this
 type(ode_solution), intent(in) :: path
 procedure(inreg_fun) :: inreg
 real(8), intent(in) :: lmin
 logical :: in_before, in_now
 integer :: i, np
 real(8) :: yprev(0:6), ynow(0:6), yc1(0:6), yc2(0:6)
 real(8) :: lsout
 integer :: nmax=30, stat, nevl


 call this%del

 ynow=path%y(:,1)
 in_now=inreg(ynow(1:3))
 if (in_now) call this%enter_pt(ynow) 
 yprev=ynow
 in_before=in_now

 np=path%np

 do i=2,np-1
  ynow=path%y(:,i)
  in_now=inreg(ynow(1:3))

  if (in_before) then
   if (.not.in_now) then
     call find_border_points(inreg,yprev,ynow,lmin,nmax,yc1,yc2,stat,nevl,lsout)
     call this%exit_pt(yc1)

!     write(*,*) '1:', i, stat, nevl
!     write(*,*) 'lsout=', lsout/lmin
!     read(*,*) 

   end if    
  else
   if (in_now) then
     call find_border_points(inreg,yprev,ynow,lmin,nmax,yc1,yc2,stat,nevl,lsout)
     call this%enter_pt(yc2)

!    write(*,*) '2:', i, stat, nevl
!    write(*,*) 'lsout=', lsout/lmin
!    read(*,*)
 
   end if 
  end if
  
  yprev=ynow
  in_before=in_now
 end do

 ynow=path%y(:,np)
 in_now=inreg(ynow(1:3))

 if (in_before) then
   if (.not.in_now) then
     call find_border_points(inreg,yprev,ynow,lmin,nmax,yc1,yc2,stat,nevl,lsout) 
     call this%exit_pt(yc1)
   else
     call this%exit_pt(ynow)
   end if    
  else
   if (in_now) then
     call find_border_points(inreg,yprev,ynow,lmin,nmax,yc1,yc2,stat,nevl,lsout)
     call this%enter_pt(yc2)
     call this%exit_pt(ynow)
   end if 
  end if

end subroutine wander_time_border


subroutine find_time(y1,y2,y)
 real(8), intent(in) :: y1(4), y2(4)
 real(8) :: y(4)
 real(8) :: t1, t2, x1(3), x2(3), x(3)
 real(8) :: dx(3), lx, dx1(3), lx1

 t1=y1(1)
 t2=y2(1)
 
 x1=y1(2:4)
 x2=y2(2:4)

 x=y(2:4)

 dx=x2-x1
 lx=sqrt(dot_product(dx,dx))

 if (lx<1d-307) then
  y(1)=t2
  return
 end if

 dx1=x-x1
 lx1=sqrt(dot_product(dx1,dx1))
 
 y(1)=t1+(t2-t1)*lx1/lx
end subroutine find_time


subroutine find_border_points(inreg,y1,y2,acc,nmax,yc1,yc2,stat,nevl,ls)
!==============================================================
! Envelope function for "find_border" and "find_time" to 
! deal with result vectors y(0:6) as whole
!==============================================================
 procedure(inreg_fun) :: inreg
 real(8), intent(in) :: y1(0:6), y2(0:6), acc
 integer, intent(in) :: nmax
 real(8), intent(out) :: yc1(0:6), yc2(0:6), ls
 integer, intent(out) :: stat, nevl

 call find_border(inreg,y1(1:3),y2(1:3),acc,nmax,yc1(1:3),yc2(1:3),stat,nevl,ls)
 call find_time(y1(0:3),y2(0:3),yc1(0:3))
 call find_time(y1(0:3),y2(0:3),yc2(0:3))

 yc1(4:6)=y1(4:6)
 yc2(4:6)=y1(4:6)  ! should be y1 !!!


end subroutine find_border_points




subroutine find_border(inreg,x1,x2,acc,nmax,t1,t2,stat,nevl,ls)
!==============================================================
! Calculation of the points t1 and t2 located inside and 
! outside of the specified region given by function 'inreg'.
! The border is considered simple, i.e. the line between
! x1 and x2 crosses the region border only once.
! 
! Input:
! ======
! inreg - function, defining the region
! x1 and x2 - the initial points, between which the border
! is located
! acc - the maximal distance between t1 and t2, which 
! is tolerated
! nmax - maximal number of evaluations of function 'inreg'
!
! Output:
! ======
! t1 and t2 - resulting points
! stat - number of exit point from the subroutine. There are
! 4 reasons to finish calculations. stat=3 means success 
! 
! nevl - actuall number of evaluations
! ls - actuall distance between t1 and t2
!==============================================================
 procedure(inreg_fun) :: inreg
 real(8), intent(in) :: x1(3), x2(3), acc
 integer, intent(in) :: nmax
 real(8), intent(out) :: t1(3), t2(3), ls
 integer, intent(out) :: stat, nevl
 real(8), parameter :: mindx=1d-307
 real(8) :: dx(3), lx, nx(3), tc(3)
 logical :: f1, f2, fc
 

 ls=0d0
 nevl=0
 dx=x2-x1
 lx=sqrt(dot_product(dx,dx))

 t1=x1
 t2=x2
 
 if (lx<mindx) then
  stat=1
  return
 end if 

 nx=dx/lx


 f1=inreg(t1)
 f2=inreg(t2)
 nevl=2

 if (f1.eqv.f2) then
  stat=2
  return
 end if

 ls=lx
 do 
  ls=ls/2
  tc=t1+nx*ls 
  fc=inreg(tc)
  nevl=nevl+1

  if (fc.eqv.f1) then
   t1=tc
  else
   t2=tc
  end if
  
!  write(*,*) 'f1=', f1, t1
!  write(*,*) 'f2=', f2, t2
!  write(*,*) 'fc=', fc, tc
!  read(*,*)

  if (ls<acc) then
   stat=3
   exit
  end if
 
  if (nevl==nmax) then
   stat=4
   exit
  end if
 
 end do

end subroutine find_border



end module time_in_reg_m


