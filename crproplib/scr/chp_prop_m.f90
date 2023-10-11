module chp_prop_m
 use ode_solver_m, only : ode_solution, ode_solver
 use fifo_rec_m, only : fifo_rec
 use vector_ops_mod, only : vector_ops
 use phys_const, only : e_charge, erg_eV, c_light, pi, psec, year

 implicit none
 private
 save

 public :: chp_prop


 abstract interface
  function fun_1d(x)
   real(8), intent(in) :: x(3)
   real(8) :: fun_1d(3)
  end function fun_1d
  
  subroutine fcond(this)
   import ode_solver
   type(ode_solver) :: this
  end subroutine
 
 end interface


 type chp_prop
  integer :: zcharge=1 ! charge of the particle
  real(8) :: energy, ipos(3), idir(3) ! energy, initial position and direction
  real(8) :: accuracy=1d-6, prop_time=1d6*year
  integer :: np_tot=1000
  real(8) :: u1_par, u2_par, rg
  integer :: nsteps=0
  integer :: nsh_old, nsh_new
  logical :: shrink=.false.
  type(ode_solution) :: sol, sol_last
  integer :: out_reason=0
  procedure(fun_1d), pointer, nopass :: mfield=>null()
  procedure(fcond), pointer, nopass :: cond1=>null()
 contains
  procedure :: set=>set_chp, calc=>calc_chp, larm_rad
  procedure :: reduc_out, reduc_out_off
 end type chp_prop

 
 real(8), parameter :: emcg=c_light*e_charge/erg_eV


 real(8) :: u1_par, u2_par
 procedure(fun_1d), pointer :: mfield=>null()
 procedure(fcond), pointer :: cond1=>null()
 
 contains


subroutine set_chp(this,mfield,energy,ipos,idir,zch)
 class(chp_prop) :: this
 procedure(fun_1d) :: mfield
 real(8), intent(in) :: energy, ipos(3), idir(3)
 integer, optional, intent(in) :: zch

 this%energy=energy
 this%ipos=ipos
 this%idir=idir
 this%mfield=>mfield

 call this%reduc_out_off
 
 if (present(zch)) this%zcharge=zch
 
 this%u1_par=c_light
 this%u2_par=this%zcharge*emcg/this%energy
 this%rg=c_light/(this%u2_par*psec)

 
end subroutine set_chp


subroutine reduc_out(this,nsh_old,nsh_new)
 class(chp_prop) :: this
 integer, intent(in) :: nsh_old, nsh_new

 this%shrink=.true.
 this%nsh_old=nsh_old
 this%nsh_new=nsh_new

end subroutine reduc_out


subroutine reduc_out_off(this)
 class(chp_prop) :: this

 this%shrink=.false.
 this%nsh_old=0
 this%nsh_new=0

end subroutine reduc_out_off


subroutine calc_chp(this,cond,prop_time,accuracy,np_tot,nlast_steps)
 class(chp_prop) :: this
 procedure(fcond) :: cond
 real(8), optional, intent(in) :: prop_time, accuracy
 integer, optional, intent(in) :: np_tot
 integer, optional, intent(in) :: nlast_steps
 type(ode_solver) :: ode
 type(ode_solution) :: sol
 type(fifo_rec) :: ylast
 integer, parameter :: ns=6
 real(8) :: y(ns)
 real(8) :: t1, t2, eps
 integer :: i, nstep


 if (present(nlast_steps)) then
  nstep=nlast_steps
 else
  nstep=1000
 end if

 this%cond1=>cond
 if (present(prop_time)) this%prop_time=prop_time
 if (present(accuracy)) this%accuracy=accuracy
 if (present(np_tot)) this%np_tot=np_tot

 
 y(1:3)=this%ipos ! coordinates
 y(4:6)=this%idir ! velocity
 
 t1=0d0
 t2=this%prop_time
 eps=this%accuracy

! external variables
 u1_par=this%u1_par
 u2_par=this%u2_par
 mfield=>this%mfield
 cond1=>this%cond1
! end of external variables


 if (this%shrink) call sol%shrink_on(this%nsh_old,this%nsh_new)

! call sol%data_max_size(300d0)
! write(*,'(A,10Es14.6)') 'sol%max_size=', sol%max_size*1d0

 call ode%ode_sys(ns,motion_ode)
 call ode%init_cond(t1,y)
 call ode%solve_to(t2)
 call ode%sol_out_to(sol)
 call ode%toler(eps)
 call ode%init_step(1d0)
 call ode%add_cond(1,cond1)
 call ode%keep_last_steps(ylast,nstep)
 call ode%solve_ode

 this%out_reason=ode%out_reason


 if (this%shrink) call sol%shrink_end

 call sol%reduce_to_npoints(this%np_tot,this%sol)


 this%nsteps=sol%np
 call sol%del

 call this%sol_last%set(ns)
 call ylast%arrange

 do i=1,ylast%ntot
  call this%sol_last%add_res(ylast%y(0,i),ylast%y(1:,i))
 end do

 call this%sol_last%cut_np
 call ylast%del

end subroutine calc_chp


function larm_rad(this,bm)
 class(chp_prop) :: this
 real(8), intent(in) :: bm ! in Gauss
 real(8) :: larm_rad ! in pc

 larm_rad=this%rg/bm
 
end function larm_rad


subroutine motion_ode(x,y,dydx)
 real(8), intent(in) :: x, y(:)
 real(8), intent(out) :: dydx(:)
 type(vector_ops) :: vop
  
 dydx(1:3)=u1_par*y(4:6)
 dydx(4:6)=u2_par*vop%cross_prod(y(4:6),mfield(y(1:3)))

end subroutine motion_ode


end module chp_prop_m



!module test_chp_prop
! use ode_solver_m, only : ode_solution, ode_solver
! use vector_ops_mod, only : vector_ops
! use phys_const, only : me_ev, e_charge, erg_eV, c_light, pi, psec
! use phys_const, only : mprot, year
! use utools_mod, only : utools
! use turb_realiz_m, only : turb_realiz
! use turbmf_tab_m, only : turbmf_tab
! use chp_prop_m, only : chp_prop

! implicit none
! private
! save

! public :: main_calc

! type(turb_realiz) :: bt1
! type(turbmf_tab) :: btab
! integer :: nstep
! 

! contains



!subroutine cond1(this)
!!================================================
!! Condition for fast calculation
!!================================================
! type(ode_solver) :: this
! real(8) :: rad, vel
! real(8) :: dydx(6), min_vel=1d-2, min_dist=1d0*psec
! integer :: i
! 
! nstep=nstep+1

! call this%fdydx(this%xc,this%yc,dydx)
!! this%yeps=this%yc+dydx*this%hc
!! rad=sqrt(this%yeps(1)**2+this%yeps(2)**2+this%yeps(3)**2)
!! vel=sqrt(this%yeps(4)**2+this%yeps(5)**2+this%yeps(6)**2)
! 
! this%yeps_flag=.false.
! 
! do i=1,3
!  if (abs(this%yeps(i))<min_dist) then
!   this%yeps(i)=min_dist
!   this%yeps_flag=.true.
!  end if
! end do


! do i=4,6
!  if (abs(this%yeps(i))<min_vel) then
!   this%yeps(i)=1d0
!   this%yeps_flag=.true.
!  end if
! end do

! if (mod(nstep,10000)==0) then
!  write(*,*) nstep, this%hc/year, this%xc/year, this%yc(1:3)/psec
! end if
!! read(*,*)
!! cur=this%yc(1)
!! eps=abs(cur-prev)/cur


!end subroutine cond1


!subroutine sch_prop
! integer :: ndev, i
! type(chp_prop) :: pr
! type(vector_ops) :: vop
! real(8) :: ipos(3), idir(3)
! character(500) :: fdir, fname

! ipos=[5d0,1d0,9d0]
! idir=vop%unit_vec(0d0,92d0)

! call pr%set(mfield,1d17,ipos,idir)
! write(*,*) pr%larm_rad(1d-6)
! nstep=0
! call pr%calc(cond1,1d6*year,1d-6,10000)


! 
! fdir='/home/anton/work/project/halo/chp_prop/traj_test/'
! fname=trim(fdir)//'traj1.dat'

! open(newunit=ndev,file=fname)
! 
! do i=1,pr%sol%np
!!  write(*,*) pr%sol%y(1:3,i)/psec
!  write(ndev,*) pr%sol%y(1:3,i)/psec
! end do

! close(ndev)
! 
! 

!end subroutine sch_prop



!subroutine set_mfield_tab
! character(500) :: fdir, fname

! fdir='/home/anton/work/project/halo/chp_prop/nest_tabs1/'
! fname=trim(fdir)//'tbt_info.dat'

! btab%fname=fname
! call btab%read(1)
! 
!end subroutine set_mfield_tab


!subroutine set_mfield
! integer :: nmodes
! real(8) :: lmin, lmax, bms 


! nmodes=1000
! lmin=20*psec
! lmax=200*psec
! bms=1d-6

! call bt1%set(nmodes,lmin,lmax,bms)

!end subroutine set_mfield


!function mfield(xc)
! real(8), intent(in) :: xc(3)
! real(8):: mfield(3)

!! mfield=[0d0,0d0,1d-6]
! call bt1%calc(xc,mfield)
!! call btab%cfield(xc,mfield)

! 


!! write(*,*) 'rg, psec=', c_light/(psec*u2_par*sqrt(dot_product(bm,bm)))
!! write(*,*) 'tg, year=', 1/(year*u2_par*sqrt(dot_product(bm,bm)))
!! read(*,*)

!! bm=bm*(1+xc(3)/(1*psec))

!end function mfield


!subroutine main_calc

! call set_mfield
! call set_mfield_tab
! call sch_prop
! 
!end subroutine main_calc



!end module test_chp_prop


!program main
! use test_chp_prop, only : main_calc

! call main_calc

!end program main




