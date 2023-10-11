module ode_solver_m
 use utools_mod, only : utools
 use fifo_rec_m, only : fifo_rec
 implicit none
 private
 save


 public :: ode_solution, ode_solver

 real(8), parameter :: def_eps=1d-6

 type ode_solution
  integer :: np=0
  integer :: nsys=0, ntot=1000
  integer :: nsh_old=1000, nsh_new=100
  integer :: new_piece=0
  logical :: shrink_tr=.false.
  integer(8) :: max_size=0
  real(8), allocatable :: y(:,:)
 contains
  procedure :: del=>del_res, set=>set_res, add_res
  procedure :: double=>double_res, reduce_to_every, reduce_to_npoints
  procedure :: cut_np=>cut_np_res, write=>write_res, read=>read_res
  procedure :: shrink_array, add=>add_shr, shrink_end
  procedure :: shrink_on, shrink_off
  procedure :: data_max_size
  procedure :: copy_from
 end type ode_solution


 type ode_solver
  integer :: ns
  real(8) :: xc, hc
  real(8), allocatable :: yc(:), yeps(:)
  real(8) :: xstart, xfinish, eps=def_eps
  logical :: stop_flag=.false., out_flag=.false., yeps_flag=.false.
  logical :: hc_pos=.true., hc_init=.false.
  logical :: cond1_flag=.false., cond2_flag=.false.
  logical :: rec_last=.false.
  integer :: out_reason=0
  integer :: adop_step_count
  procedure(fderiv), pointer, nopass :: fdydx=>null()
  procedure(fcond), pointer, nopass :: cond1=>null(), cond2=>null()
  type(ode_solution), pointer :: ysol=>null()
  type(fifo_rec), pointer :: ylast=>null()
 contains
  procedure :: ode_sys, init_cond, toler, init_step, solve_to, add_cond
  procedure :: sol_out_to, ode_sys_del, solve_ode
  procedure :: keep_last_steps
 end type ode_solver


 abstract interface
  subroutine fderiv(x,y,dydx)
   real(8), intent(in) :: x, y(:)
   real(8), intent(out) :: dydx(:)
  end subroutine
 end interface

 abstract interface
  subroutine fcond(this)
   import ode_solver
   type(ode_solver) :: this
  end subroutine
 end interface




 contains

!====================================================================
! Type ode_solution subroutines
!====================================================================

subroutine del_res(this)
 class(ode_solution) :: this
 
 if (allocated(this%y)) deallocate(this%y)
 this%np=0
 this%nsys=0
 this%ntot=1000
end subroutine del_res


subroutine set_res(this,nsys,ntot,np)
 class(ode_solution) :: this
 integer, intent(in) :: nsys
 integer, intent(in), optional :: ntot, np

 call this%del
 this%nsys=nsys
 if (this%max_size==0) call this%data_max_size

 if (present(ntot)) this%ntot=ntot
 if (this%ntot>this%max_size) this%ntot=this%max_size
 if (present(np)) this%np=np

 allocate(this%y(0:nsys,this%ntot))

end subroutine set_res


subroutine copy_from(this,that)
 class(ode_solution) :: this
 type(ode_solution), intent(in) :: that

 call this%set(that%nsys,that%ntot,that%np)
 this%nsh_old=that%nsh_old
 this%nsh_new=that%nsh_new
 this%new_piece=that%new_piece
 this%shrink_tr=that%shrink_tr
 this%max_size=that%max_size
 this%y=that%y

end subroutine copy_from



subroutine data_max_size(this,sizebyte)
 class(ode_solution) :: this
 real(8), optional :: sizebyte
 real(8) :: one_rec, sb
 real(8), parameter :: real8=8d0 ! 8 bytes for real(8)
 real(8), parameter :: gbyte=1024**3
 real(8), parameter :: sizebyte_def=2*gbyte

 if (present(sizebyte)) then
  sb=sizebyte
 else
  sb=sizebyte_def
 end if
  
 one_rec=(this%nsys+1)*real8
 this%max_size=int8(sb/one_rec)

end subroutine data_max_size


recursive subroutine double_res(this)
 class(ode_solution) :: this
 real(8), allocatable :: y1(:,:)
 integer :: i, ntot2, np, np2

 ntot2=2*this%ntot
 np=this%np
 if (ntot2>this%max_size) ntot2=this%max_size
 
 if (np>ntot2) then
  np2=np/2.5d0
  call shrink_array(this,1,np,np2)
  return
 end if
 
 allocate(y1(0:this%nsys,ntot2))
 do i=1,this%ntot
   y1(:,i)=this%y(:,i)
 end do
 
 deallocate(this%y)
 allocate(this%y(0:this%nsys,ntot2))
 do i=1,this%ntot
   this%y(:,i)=y1(:,i)
 end do
 deallocate(y1) 
 
 this%ntot=ntot2
 
end subroutine double_res



subroutine add_res(this,x,y)
 class(ode_solution) :: this
 real(8), intent(in) :: x, y(:)
 integer :: np

 np=this%np+1
 
 if (np>this%ntot) call this%double
 if (this%np+1<np) np=this%np+1
 
 this%y(0,np)=x
 this%y(1:,np)=y
 this%np=np
 
end subroutine add_res


subroutine cut_np_res(this)
 class(ode_solution) :: this
 type(ode_solution) :: that
 integer :: i


 call that%set(this%nsys,this%np,this%np)

 do i=1,this%np
  that%y(:,i)=this%y(:,i)
 end do
 
 call this%set(that%nsys,that%np,that%np)
 do i=1,that%np
  this%y(:,i)=that%y(:,i)
 end do
 call that%del

end subroutine cut_np_res


subroutine write_res(this,ndev,tp)
 class(ode_solution) :: this
 integer, intent(in) :: ndev
 character(2), optional :: tp ! r4, r8
 integer :: rtype, i
 real, allocatable :: y(:,:)
 
 call this%cut_np
 rtype=4
 if (present(tp)) then
  if (tp=='r8') rtype=8
 end if

 write(ndev) rtype, this%nsys, this%ntot, this%np


 if (rtype==4) then
  allocate(y(0:this%nsys,this%ntot))
  do i=1,this%ntot
   y(:,i)=real(this%y(:,i))
  end do
  write(ndev) y
  deallocate(y)
 else 
  write(ndev) this%y 
 end if
 
 
end subroutine write_res


subroutine read_res(this,ndev)
 class(ode_solution) :: this
 integer, intent(in) :: ndev
 integer :: rtype, i
 real, allocatable :: y(:,:)
 integer :: nsys, ntot, np, rstat
 integer :: nlast
 

 read(ndev,iostat=rstat) rtype, nsys, ntot, np
 if (rstat<0) return

 if (rtype==4) then
  allocate(y(0:nsys,ntot))
  nlast=ntot
  do i=1,ntot
   read(ndev,iostat=rstat) y(:,i)
   if (rstat<0) then
    nlast=i-1
    exit
   end if
  end do
  rstat=0
  
  call this%set(nsys,nlast,nlast)

  
  if (rstat<0) then
   deallocate(y)
   call this%del
   return 
  end if

  do i=1,this%ntot
   this%y(:,i)=real(y(:,i),8)
  end do
  deallocate(y)
 else
  read(ndev,iostat=rstat) this%y
  
  if (rstat<0) then
   call this%del
   return 
  end if
 
 end if

end subroutine read_res


!====================================================================
! Type ode_solution subroutines related to reduction
! of data
!====================================================================
subroutine reduce_to_every(this,dx,that)
!==========================================
! Reduce this to that, so that indepedent
! variable changes with the step not 
! smaller than dx
!==========================================
 class(ode_solution) :: this
 real(8), intent(in) :: dx
 type(ode_solution) :: that
 real(8) :: x
 integer :: i
 
 call that%set(this%nsys)
 call that%add_res(this%y(0,1),this%y(1:,1)) ! start point
 
 x=this%y(0,1)
 do i=2,this%np-1
  if (abs(this%y(0,i)-x)>dx) then
   x=this%y(0,i)
   call that%add_res(this%y(0,i),this%y(1:,i))
  end if
 end do
 
 call that%add_res(this%y(0,this%np),this%y(1:,this%np)) ! end point
end subroutine reduce_to_every


subroutine reduce_to_npoints(this,nnew,that)
!==========================================
! Reduce this to that with that containing
! only nnew points using uniform reduction.
!==========================================
 class(ode_solution) :: this
 integer, intent(in) :: nnew
 type(ode_solution) :: that
 type(utools) :: ut
 integer, allocatable :: anum(:)
 integer :: i, j, nold

 nold=this%np
 call that%set(this%nsys)
 
 if (nold<=1.6*nnew) then
  do i=1,nold
   call that%add_res(this%y(0,i),this%y(1:,i))
  end do
  
 else
  call ut%uniform_reduction(this%np,nnew,anum)
  
  do i=1,nnew
   j=anum(i)
   call that%add_res(this%y(0,j),this%y(1:,j))
  end do
  deallocate(anum)
 end if

end subroutine reduce_to_npoints


subroutine shrink_array(this,i1,i2,nnew,in1,in2)
!==========================================
! Reduce number of points in the segment 
! from i1 to i2
! (i1 and i2 are included in the segment)
! to nnew, if i2-i1+1>1.6*nnew
! It gives in1 and in2, which are new 
! indecies of array elements, which before 
! reduction were i1 and i2
!==========================================
 class(ode_solution) :: this
 integer, intent(in) :: i1, i2, nnew
 integer, intent(out), optional :: in1, in2
 type(ode_solution) :: that
 type(utools) :: ut
 integer, allocatable :: anum(:)
 integer :: i, j, nold, i3

 nold=i2-i1+1

 if (present(in1)) in1=0
 if (present(in2)) in2=0
! if nold/nnew<1.6 ignore reduction
! because reduction work not always correct
 if (nold<=1.6*nnew) return
  

 call ut%uniform_reduction(nold,nnew,anum)
 call that%set(this%nsys)

 do i=1,nnew
  j=anum(i)+i1-1
  call that%add_res(this%y(0,j),this%y(1:,j))
 end do


 deallocate(anum)

 i3=i2+1
 do i=i3,this%np
  call that%add_res(this%y(0,i),this%y(1:,i))
 end do


 do i=1,nnew
  j=i+i1-1
  this%y(:,j)=that%y(:,i)
 end do

 if (present(in1)) in1=i1
 if (present(in2)) in2=j

 i3=nnew+1
 do i=i3,that%np
  j=j+1
  this%y(:,j)=that%y(:,i)
 end do
 
 this%np=j
 call that%del
 call this%cut_np
 call this%double
 
end subroutine shrink_array


subroutine add_shr(this,x,y)
!==========================================
! Add x and y to this and reduce
! number of points every this%nsh_old points
! to this%nsh_new, if required
!==========================================
 class(ode_solution) :: this
 real(8), intent(in) :: x, y(:)
 integer :: i1, i2, in1, in2

 this%np=this%np+1

 if (this%np>this%ntot) call this%double
 this%y(0,this%np)=x
 this%y(1:,this%np)=y

 
 if (this%shrink_tr) then
  if ((this%np-this%new_piece)==this%nsh_old) then

   i1=this%new_piece+1
   i2=this%np

   call this%shrink_array(i1,i2,this%nsh_new,in1,in2)
   this%new_piece=in2

  end if
 end if
  

end subroutine add_shr


subroutine shrink_end(this)
!==========================================
! Reduce cut the end of the array
! according to the policy used in
! add_shr
!==========================================
 class(ode_solution) :: this
 integer :: i1, i2, in1, in2, nnew

 if (this%shrink_tr) then
   if (this%new_piece>=this%nsh_new) then
   nnew=(this%nsh_new*1d0/this%nsh_old)*(this%np-this%new_piece)     
  
   i1=this%new_piece+1
   i2=this%np

!   write(*,*) 'nsh_new=', this%nsh_new
!   write(*,*) 'nsh_old=', this%nsh_old
!   write(*,*) 'np=',  this%np
!   write(*,*) 'new_piece=',  this%new_piece
!   write(*,*) 'nnew, i1, i2=',  nnew, i1, i2

   call this%shrink_array(i1,i2,nnew,in1,in2)
   this%new_piece=in2
  end if
 end if
 
end subroutine shrink_end


subroutine shrink_on(this,nold,nnew,npiece)
!==========================================
! Set parameters of reduction for 
! add_shr and shrink_end functions
!==========================================
 class(ode_solution) :: this
 integer, intent(in) :: nold, nnew
 integer, intent(in), optional :: npiece
 
 this%shrink_tr=.true.
 this%nsh_old=nold
 this%nsh_new=nnew
 if (present(npiece)) this%new_piece=npiece

end subroutine shrink_on


subroutine shrink_off(this)
!==========================================
! Swicht off the reduction of data,
! which is by default
!==========================================
 class(ode_solution) :: this
 
 this%shrink_tr=.false.
 this%nsh_old=1000
 this%nsh_new=100
 this%new_piece=0

end subroutine shrink_off


!====================================================================
! Type ode_solver interface subroutines
!====================================================================


subroutine ode_sys_del(this)
 class(ode_solver) :: this

 this%ns=0
 this%xc=0d0
 this%hc=0d0
 if (allocated(this%yc)) deallocate(this%yc)
 if (allocated(this%yeps)) deallocate(this%yeps)
 this%xstart=0d0
 this%xfinish=0d0
 this%eps=def_eps
 this%stop_flag=.false. 
 this%out_flag=.false. 
 this%yeps_flag=.false.
 this%yeps_flag=.false.
 this%hc_pos=.true.
 this%hc_init=.false. 
 this%cond1_flag=.false. 
 this%cond2_flag=.false.
 this%rec_last=.false. 
 this%out_reason=0
 this%fdydx=>null()
 this%cond1=>null()
 this%cond2=>null()
 this%ysol=>null()
 this%ylast=>null()

end subroutine ode_sys_del



subroutine ode_sys(this,ns,fdydx)
 class(ode_solver) :: this
 integer, intent(in) :: ns
 procedure(fderiv) :: fdydx
 
 call this%ode_sys_del
 
 this%ns=ns
 this%fdydx=>fdydx
 allocate(this%yc(ns))
 allocate(this%yeps(ns))


end subroutine ode_sys

subroutine init_cond(this,xc,yc)
 class(ode_solver) :: this
 real(8), intent(in) :: xc, yc(:)
 this%xc=xc
 this%xstart=xc
 if (.not.(allocated(this%yc))) then
  write(*,*) 'ode_solver: Please initialise ode with ode_sys function'
  return
 end if

 this%yc=yc
end subroutine init_cond


subroutine toler(this,eps)
 class(ode_solver) :: this
 real(8), intent(in) :: eps
 this%eps=eps
end subroutine toler

subroutine init_step(this,hc)
 class(ode_solver) :: this
 real(8), intent(in) :: hc
 this%hc=hc
 this%hc_init=.true.
end subroutine init_step

subroutine solve_to(this,xfinish)
 class(ode_solver) :: this
 real(8), intent(in) :: xfinish

 this%xfinish=xfinish
end subroutine solve_to
 

subroutine add_cond(this,ncond,cond)
 class(ode_solver) :: this
 integer :: ncond
 procedure(fcond) :: cond
  
 select case(ncond)
  case (1)
   this%cond1=>cond
   this%cond1_flag=.true.
  case (2)
   this%cond2=>cond
   this%cond2_flag=.true.
  case default
   write(*,*) 'Condition number ', ncond, ' >', 2, ' -possible number of conditions'
 end select  

end subroutine add_cond


subroutine sol_out_to(this,ysol)
 class(ode_solver) :: this
 type(ode_solution), target :: ysol
 
 call ysol%set(this%ns)  !  record the integration
 this%ysol=>ysol
 this%out_flag=.true.
end subroutine sol_out_to


subroutine keep_last_steps(this,ylast,nstep)
 class(ode_solver) :: this
 type(fifo_rec), target :: ylast
 integer, intent(in), optional :: nstep

 if (present(nstep)) call ylast%set_ntot(nstep)
 this%ylast=>ylast
 this%rec_last=.true.

end subroutine keep_last_steps



subroutine auto_step(this)
 type(ode_solver) :: this 
 if (this%hc_init) then
  if (this%xfinish>=this%xstart) then
   this%hc=abs(this%hc)
   this%hc_pos=.true.
  else
   this%hc=-abs(this%hc)
   this%hc_pos=.false.
  end if
 else
  this%hc=(this%xfinish-this%xstart)/100
  if (this%hc>=0d0) then
   this%hc_pos=.true.
  else
   this%hc_pos=.false.
  end if
 end if
end subroutine auto_step


!====================================================================
! Solution subroutines
!====================================================================


subroutine solve_ode(this)
 class(ode_solver) :: this

 call auto_step(this)
 if (this%out_flag) call this%ysol%add(this%xc,this%yc)

 this%stop_flag=.false.
 do

  if (this%cond1_flag) then    
   call this%cond1(this)
   if (this%stop_flag) then
    this%out_reason=1
    exit
   end if
  end if
  
  call rk_adap_step(this%ns,this)
  if (this%stop_flag) then
   this%out_reason=2
   exit
  end if
  
  if (this%out_flag) call this%ysol%add(this%xc,this%yc)
  if (this%rec_last) call this%ylast%add(this%xc,this%yc)

  if (this%cond2_flag) then    
   call this%cond2(this)
  
   if (this%stop_flag) then
    this%out_reason=3
    exit
   end if
  end if

  call fin_cond(this)

  if (this%stop_flag) then
   this%out_reason=4
   exit
  end if
 end do
 

end subroutine solve_ode



subroutine fin_cond(this)
 type(ode_solver) :: this

 if (this%hc_pos) then
  if (this%xc+this%hc>this%xfinish) this%hc=this%xfinish-this%xc
  if (this%xc>=this%xfinish) this%stop_flag=.true.
 else
  if (this%xc+this%hc<this%xfinish) this%hc=this%xfinish-this%xc
  if (this%xc<=this%xfinish) this%stop_flag=.true.
 end if

end subroutine fin_cond




subroutine rk_adap_step(ns,this)
!=======================================================================
! Stepper subroutine for the Embedded Runga-Kutta Method with 
! with Cash-Karp Parameters "rkutta_ck5"
!======================================================================= 
 integer, intent(in) :: ns
 type(ode_solver) :: this
 real(8) :: xnew, yold(ns), yerr(ns), yre(ns), rel_err
 real(8) :: hc_min
 logical :: stop_calc, reach_mstep
 integer, parameter :: max_adopt=1000000

 this%adop_step_count=0
 do
  if (this%adop_step_count>max_adopt) then
    this%stop_flag=.true.
    exit
  end if
   
!  if (this%adop_step_count>100) then
!   write(*,'(A,10I6)') 'adop_step_count=', this%adop_step_count 
!   write(*,'(A,10Es14.6)') '1this%hc=', this%hc
!  end if
! Check for too small step
  call small_step(this%xc,this%hc,xnew,reach_mstep)
  
!  if (this%adop_step_count>100) then
!   write(*,'(A,10Es14.6)') '2this%hc=', this%hc
!   read(*,*)
!  end if  

 
  yold=this%yc
  call rkutta_ck5(ns,this%yc,this%fdydx,this%xc,this%hc,yerr)

! Check for yerr=0, which means that derivative is zero at (xc,yc)  
  call zero_array(ns,yerr,stop_calc) 
  if (stop_calc) then
   this%hc=2*this%hc
   exit
  end if

! The way of evaluation of inverse relative error  
  if (this%yeps_flag) then
   yre=abs(this%yeps/yerr)
  else
   yre=abs(this%yc/yerr)
  end if

!  if (this%adop_step_count>100) then
!   write(*,*)'======================================================'
!   write(*,*) 'hstep=', this%hc
!   write(*,*) 'yeps_flag=', this%yeps_flag
!   write(*,'(A,10Es14.6)') 'yerr=', yerr
!   write(*,'(A,10Es14.6)') 'yc=', this%yc
!   write(*,'(A,10Es14.6)') 'yeps=', this%yeps
!   write(*,*) 'vel=', sqrt(this%yc(4)**2+this%yc(5)**2+this%yc(6)**2)
!   write(*,*)
!   write(*,'(A,10Es14.6)') 'yre=', yre
!   write(*,'(A,10Es14.6)') 'minval(yre), eps, m*eps=', minval(yre), this%eps, this%eps*minval(yre)
!   write(*,*)'======================================================'
!  end if


! Check for yre=0, which means yc=0 or yeps
  call zero_array(ns,yre,stop_calc)
  if (stop_calc) exit

! Step change
  rel_err=this%eps*minval(yre)
  if (reach_mstep) then
   hc_min=this%hc
   call adjust_step(rel_err,this%hc,stop_calc)
   if (abs(this%hc)<abs(hc_min)) stop_calc=.true.
  else
   call adjust_step(rel_err,this%hc,stop_calc)
  end if

!  call adjust_step(rel_err,this%hc,stop_calc)
!  if (this%adop_step_count>100) then
!   write(*,*)'======================================================'
!   write(*,*) 'hstep1=', this%hc
!   write(*,*)'======================================================'
!  end if

  if (stop_calc) exit

  this%yc=yold ! in case of new iteration with new smaller step
  this%adop_step_count=this%adop_step_count+1
 end do
 
 this%xc=xnew
  
end subroutine rk_adap_step


subroutine zero_array(ns,y,zero_true)
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

end subroutine zero_array


subroutine small_step(t0,step,tnew,reach_mstep)
 real(8), intent(in) :: t0
 real(8), intent(inout) :: step
 real(8), intent(out) :: tnew
 logical, intent(out) :: reach_mstep
 real(8), parameter :: mstep=1d-15
 real(8) :: min_step

 reach_mstep=.false.
 min_step=abs(t0*mstep)
 if (abs(step)<min_step) then
  reach_mstep=.true.
  step=sign(min_step,step)
 end if
 tnew=t0+step

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
 procedure(fderiv) :: dy
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



end module ode_solver_m



!module test2
! use ode_solver_m, only : ode_solution
! implicit none
! private
! save

! public :: main_calc

! contains 


!subroutine main_calc
! type(ode_solution) :: sol
! integer, parameter :: ns=3
! integer :: i
! 
! 
! call sol%set(ns)
! call sol%double
! call sol%double
! call sol%double

! do i=1,sol%ntot
!  write(*,'(I5,10Es14.6)') i, sol%y(:,i)
! end do

!end subroutine main_calc



!end module test2



!module test1
! use ode_solver_m, only : ode_solution, ode_solver
! implicit none
! private
! save

! public :: main_calc

! contains


!subroutine main_calc
! type(ode_solver) :: ode
! type(ode_solution) :: sol
! integer, parameter :: ns=3
! integer :: i
! real(8) :: y(ns), yeps(ns), dydx(ns)
! real(8) :: x1, x2, eps, h, xout
! 
! x1=0d0
! x2=1d2
! y=[1d0,0d0,0d0]
! eps=1d-10
! 
! call sol%set(ns,6000)

! call ode%ode_sys(ns,lorentz_sys)
! call ode%init_cond(x1,y)
! call ode%solve_to(x2)
! call ode%sol_out_to(sol)
! call ode%toler(eps)
!! call ode%add_cond(1,cond1)
!! call ode%add_cond(2,cond2)
! call ode%solve_ode

! 
! write(*,*) sol%np, ode%xc, ode%yc
! open(1,file='ftest.dat')
! do i=1,sol%np
!!  write(*,'(I5,10Es14.6)') i, sol%y(:,i)
!  write(1,*) sol%y(:,i)
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
! dydx(2)=y(1)*(1000-y(3))-y(2)
! dydx(3)=y(1)*y(2)-8*y(3)/3

!end subroutine lorentz_sys

!subroutine cond1(this)
! type(ode_solver) :: this
! real(8) :: dydx(3)
! 
! call this%fdydx(this%xc,this%yc,dydx)
! this%yeps=this%yc+dydx*this%hc

!end subroutine cond1

!subroutine cond2(this)
! type(ode_solver) :: this

! if (this%yc(3)>30d0) this%stop_flag=.true.

!end subroutine cond2

!end module test1




!program main
! use test2, only : main_calc

! call main_calc

!end program main
