module turbmf_tab_m
 use phys_const, only : c_light, pi, psec
 use turb_realiz_m, only : turb_realiz
 use field_3dlm, only : field_3dl, reind
 use nest_grid_m, only : nest_grid
 use tab_ngrid_m, only : start_old_nest 

 implicit none
 private
 save

 public :: turbmf_tab, main_calc


 type turbmf_tab
  character(500) :: fname
  integer :: nlev
  type(field_3dl), allocatable :: tbm(:)
 contains
  procedure :: read=>read_ttab, cfield=>cfield_ttab
 end type turbmf_tab
 


 contains


subroutine read_ttab(this,nlev)
 class(turbmf_tab) :: this
 integer, intent(in) :: nlev
 integer :: i
 type(nest_grid):: tg
 type(turb_realiz), allocatable :: bt(:)


 this%nlev=nlev
 call start_old_nest(this%fname,tg,bt)

 if (allocated(this%tbm)) deallocate(this%tbm)
 allocate(this%tbm(nlev))

 do i=1,nlev
  this%tbm(i)%fname=trim(tg%fnd(i))
  call this%tbm(i)%read
 end do
 

end subroutine read_ttab


subroutine cfield_ttab(this,x,bm)
 class(turbmf_tab) :: this
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: bm(3)
 real(8) :: b0(3), xc(3), xmax
 integer :: i, ix(3), nx

 bm=0d0
 xc=x

! Periodic boundary conditions
 nx=this%tbm(1)%nx
 xmax=this%tbm(1)%x(nx)
 xc(1)=xc(1)-floor(x(1)/xmax)*xmax
 
 nx=this%tbm(1)%ny
 xmax=this%tbm(1)%y(nx)
 xc(2)=xc(2)-floor(x(2)/xmax)*xmax
 
 nx=this%tbm(1)%nz
 xmax=this%tbm(1)%z(nx)
 xc(3)=xc(3)-floor(x(3)/xmax)*xmax


 do i=1,this%nlev
  call this%tbm(i)%cfield(xc(1),xc(2),xc(3),b0)
  bm=bm+b0
  xc=this%tbm(i)%xc1
 end do

end subroutine cfield_ttab




subroutine read_nest_tabs
 type(turbmf_tab) :: bs
 character(500) :: fdir, fname
 type(nest_grid):: tg
 type(turb_realiz), allocatable :: bt(:)
 type(field_3dl), allocatable :: tbm(:)
 integer :: i, nlev
 integer :: nx(2), ny(2), nz(2)
 real(8) :: xc(3), bm(3), btot(3)

 fdir='/home/anton/work/project/halo/chp_prop/nest_tabs_test/'
 fname=trim(fdir)//'tbt_info.dat'

 bs%fname=fname
 call bs%read(2)
 
 xc=[-1d0,-10d0,-436d0]*psec

 call bs%cfield(xc,bm)
 
 write(*,*) 'btot=', sqrt(dot_product(bm,bm))


end subroutine read_nest_tabs




subroutine main_calc

 call read_nest_tabs

end subroutine main_calc



end module turbmf_tab_m


!program main
! use turbmf_tab_m, only : main_calc

! call main_calc

!end program main
