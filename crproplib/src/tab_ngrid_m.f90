module tab_ngrid_m
 use utools_mod, only : utools
 use phys_const, only : c_light, pi, psec
 use turb_realiz_m, only : turb_realiz
 use field_3dlm, only : field_3dl, reind
 use timer_module, only : timer_class
 use nest_grid_m, only : nest_grid

 implicit none
 private
 save

 public :: main_calc, start_old_nest, start_new_nest


 contains


subroutine print_nest(tg)
 type(nest_grid), intent(in) :: tg
 integer :: i

 write(*,'(A,10Es14.6)') 'lmin, lmax= ', tg%lmin/psec, tg%lmax/psec
 do i=1,tg%nlev
  write(*,*) 
  write(*,*) 'i, ng=', i, tg%ngrid(i)
  write(*,'(A,10Es14.6)') 'dx, n*dx bmi = ', tg%dx(i)/psec, tg%ngrid(i)*tg%dx(i)/psec,  tg%bmi(i)
  write(*,'(A,10Es14.6)') 'l1, l2= ', tg%lk1(i)/psec, tg%lk2(i)/psec
  write(*,*) trim(tg%fnd(i))
  write(*,*) 'seed=', tg%seed(:,i)
 end do

end subroutine print_nest



subroutine start_old_nest(fname,tg,bt)
 character(500), intent(in) :: fname
 type(nest_grid):: tg
 type(turb_realiz), allocatable :: bt(:)
 integer :: i

 tg%fname=trim(fname)
 call tg%read
 allocate(bt(tg%nlev))
 do i=1,tg%nlev
!  bt(i)%manual_seed=.true.
!  allocate(bt(i)%seed(size(tg%seed(:,i),1)))
!  bt(i)%seed=tg%seed(:,i)
  call bt(i)%set(tg%nmodes,tg%lk1(i),tg%lk2(i),tg%bmi(i))
!  tg%seed(:,i)=bt(i)%seed
 end do

end subroutine start_old_nest



subroutine start_new_nest(nlev,ngrid,alpha,lmax,bm,fnd,fname,tg,bt)
 integer, intent(in) :: nlev, ngrid
 real(8), intent(in) :: alpha, lmax, bm
 character(500) :: fnd(nlev)
 character(500), intent(in) :: fname
 type(nest_grid):: tg
 type(turb_realiz), allocatable :: bt(:)
 integer :: i


 tg%ngrid_def=ngrid
 tg%fname=trim(fname)
 call tg%set(alpha,lmax,bm,nlev)
 tg%fnd=fnd 


 allocate(bt(tg%nlev))
 do i=1,tg%nlev
  call bt(i)%set(tg%nmodes,tg%lk1(i),tg%lk2(i),tg%bmi(i))
  tg%seed(:,i)=bt(i)%seed
 end do

 call tg%write


end subroutine start_new_nest




subroutine calc_turb_lev(ilev,tg,bt)
 integer, intent(in) :: ilev
 type(nest_grid):: tg
 type(turb_realiz) :: bt(:)
 type(utools) :: ut
 type(field_3dl) :: tbm
 type(timer_class) :: tmr
 integer :: i, j, k
 integer :: nx, num
 real(8) :: dx, xmin, xmax
 real(8) :: xc(3), bm(3)
 real(8), allocatable :: xg(:)



 dx=tg%dx(ilev)
 nx=tg%ngrid(ilev)
 xmin=0d0
 xmax=nx*dx
 nx=nx+1
 
 call ut%grid(xg,xmin,xmax,nx,'lin')

 call tbm%set(nx,nx,nx)
 tbm%x=xg
 tbm%y=xg
 tbm%z=xg
 tbm%fname=trim(tg%fnd(ilev))


 call tmr%start(nx*nx*nx,1.0)
 do k=1,nx
  do j=1,nx
   do i=1,nx
    xc=[xg(i),xg(j),xg(k)]
    call bt(ilev)%calc(xc,bm)
    call reind(i,j,k,num)
    tbm%f(num)%bm=bm
    call tmr%loop
   end do
  end do
 end do

 call tbm%write
 

end subroutine calc_turb_lev



subroutine calc_tab
 integer, parameter :: nlev=5
 character(500) :: fdir, fname
 character(500) :: fnd(nlev)
 type(nest_grid):: tg
 type(turb_realiz), allocatable :: bt(:)
 integer :: ilev


 fdir='/home/anton/work/project/halo/chp_prop/nest_tabs_test/'
 fname=trim(fdir)//'tbt_info.dat'
 fnd(1)=trim(fdir)//'tb_tab1.dat'
 fnd(2)=trim(fdir)//'tb_tab2.dat'
 fnd(3)=trim(fdir)//'tb_tab3.dat'
 fnd(4)=trim(fdir)//'tb_tab4.dat'
 fnd(5)=trim(fdir)//'tb_tab5.dat'


! start_new_nest(nlev,255,5d0/3d0,200*psec,1d-6,fnd,fname,tg,bt)
 call start_old_nest(fname,tg,bt)
 call print_nest(tg)

 ilev=3
 call calc_turb_lev(ilev,tg,bt)
 

end subroutine calc_tab



subroutine main_calc

 call calc_tab

end subroutine main_calc



end module tab_ngrid_m


!program main
! use tab_ngrid_m, only : main_calc

! call main_calc

!end program main
