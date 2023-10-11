module traj_analysis_m
 use cent_grida_m, only : cent_grida
 use restr_region, only : r0region
 implicit none
 private
 save


 public :: traj_analysis

 type traj_analysis
 contains
  procedure, nopass :: rec_time2grid
 end type traj_analysis



 contains


subroutine rec_time2grid(t1,t2,x1,x2,cgrid)
!============================================
! Record of the time spent by a particle
! in the cells of the grid "cgrid" moving
! along straight line between x1 and x2
! t1 is the time, when particle in x1
! and t2 is the time for x2 point
!============================================ 
 real(8), intent(in) :: t1, t2, x1(3), x2(3)
 type(cent_grida), intent(in) :: cgrid
 type(r0region) :: rcube
 real(8) :: dtm, dx(3), lx, nx(3), dtime, xs(3)
 real(8) :: ltot, ti(2), ls, dt
 logical :: fr, st
 integer :: nc(2,3)
 real(8), parameter :: shift_out=1d-6

 dtm=t2-t1
 dx=x2-x1
 lx=sqrt(dot_product(dx,dx))
 xs=x1 ! start point


! In case of motionless particle
 if (lx<=0d0) then
  call cgrid%add_val(xs,dtm)
  return
 end if

 nx=dx/lx
 dtime=dtm/lx
 ltot=0d0

 do 
  call cgrid%find_cell(xs,fr,nc,rcube%xb)
  if (.not.fr) return ! outside of grid
  if (rcube%inside(x2)) exit
  
  call rcube%shoot(xs,nx,st,ti)
  if (ti(2)<=0d0) ti(2)=shift_out*rcube%mindl()

  ls=ti(2)*(1d0+shift_out)
  ltot=ltot+ls
  
  dt=dtime*ls    ! determine time interval spent in the cell
  call cgrid%add_val(xs,dt)
  xs=xs+nx*ls
 end do
 
 dt=dtime*(lx-ltot)
 call cgrid%add_val(xs,dt)

end subroutine rec_time2grid




end module traj_analysis_m
