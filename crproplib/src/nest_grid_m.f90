module nest_grid_m
 use phys_const, only : c_light, pi, psec
 use turb_realiz_m, only : turb_realiz

 implicit none
 private
 save

 public :: nest_grid, main_calc


 type nest_grid
  integer :: ngrid_def=255, nmodes=1000
  integer, private :: n_oscill=4
  integer :: nlev
  integer, allocatable :: ngrid(:)
  real(8) :: lmax, lmin, alpha, brms
  real(8), allocatable :: dx(:)
  real(8), allocatable :: lk1(:), lk2(:), bmi(:)
  character(500), allocatable :: fnd(:)
  character(500) :: fname='dat.dat'
  integer, allocatable :: seed(:,:)
 contains
  procedure :: del=>del_ngrid, set=>set_ngrid
  procedure :: write=>write_ngrid, read=>read_ngrid
 end type nest_grid



 contains

subroutine del_ngrid(this)
 class(nest_grid) :: this

 if (allocated(this%ngrid)) deallocate(this%ngrid)
 if (allocated(this%dx)) deallocate(this%dx)
 if (allocated(this%lk1)) deallocate(this%lk1)
 if (allocated(this%lk2)) deallocate(this%lk2)
 if (allocated(this%bmi)) deallocate(this%bmi)
 if (allocated(this%fnd)) deallocate(this%fnd)
 if (allocated(this%seed)) deallocate(this%seed)

end subroutine del_ngrid



subroutine set_ngrid(this,alpha,lmax,brms,nlev)
 class(nest_grid) :: this
 real(8), intent(in) :: alpha, lmax, brms
 integer, intent(in) :: nlev
 real(8) :: al1, bnorm
 integer :: i, nss
 
 call this%del

 this%alpha=alpha
 this%lmax=lmax
 this%brms=brms
 this%nlev=nlev

 allocate(this%ngrid(nlev))
 allocate(this%dx(nlev))
 allocate(this%lk1(nlev))
 allocate(this%lk2(nlev))
 allocate(this%bmi(nlev))
 allocate(this%fnd(nlev))

! call random_seed(size = nss)
 nss=12 ! should be changed later
 allocate(this%seed(nss,nlev))

 
 this%ngrid(1)=this%ngrid_def
 this%dx(1)=this%lmax*this%n_oscill/this%ngrid(1)
 this%lk1(1)=2*this%dx(1)
 this%lk2(1)=lmax
 
 

 do i=2,nlev
  this%ngrid(i)=this%ngrid_def
  this%dx(i)=this%dx(i-1)/this%ngrid(i)
  this%lk1(i)=2*this%dx(i)
  this%lk2(i)=this%dx(i-1)
 end do

 this%lmin=this%lk1(nlev)

 al1=alpha-1
 bnorm=brms**2/(this%lmax**al1-this%lmin**al1)

 do i=1,nlev
  this%bmi(i)=sqrt((this%lk2(i)**al1-this%dx(i)**al1)*bnorm)
 end do
 

end subroutine set_ngrid



subroutine write_ngrid(this)
 class(nest_grid) :: this
 integer :: ndev
 
 if (this%fname=='dat.dat') then
  write(*,*) 'nest_grid: The name of the output file is not set'
  return
 end if


 open(newunit=ndev,file=this%fname,form='unformatted')
 
 write(ndev) this%ngrid_def
 write(ndev) this%n_oscill
 write(ndev) this%nmodes
 write(ndev) this%alpha
 write(ndev) this%lmax
 write(ndev) this%brms
 write(ndev) this%nlev
 write(ndev) this%seed
 write(ndev) this%fnd
 
 close(ndev)



end subroutine write_ngrid

subroutine read_ngrid(this)
 class(nest_grid) :: this
 real(8) :: alpha, lmax, brms
 integer :: ndev, nlev
 
 if (this%fname=='dat.dat') then
  write(*,*) 'nest_grid: The name of the input file is not set'
  return
 end if


 open(newunit=ndev,file=this%fname,form='unformatted')
 
 read(ndev) this%ngrid_def
 read(ndev) this%n_oscill
 read(ndev) this%nmodes
 read(ndev) alpha
 read(ndev) lmax
 read(ndev) brms
 read(ndev) nlev
 call this%set(alpha,lmax,brms,nlev)
 read(ndev) this%seed
 read(ndev) this%fnd
 
 close(ndev)


end subroutine read_ngrid



subroutine main_calc
 type(nest_grid):: tg
 type(turb_realiz), allocatable :: bt(:)
 integer :: i, nlev
 character(500) :: fdir, fname


 fdir='/home/anton/work/project/halo/chp_prop/nest_tabs_test11/'

 nlev=5
 fname=trim(fdir)//'tab_info.dat'
 tg%fname=trim(fname)
 call tg%set(5/3d0,200*psec,1d-6,nlev)
 tg%fnd(1)=trim(fdir)//'tb_tab1.dat'
 tg%fnd(2)=trim(fdir)//'tb_tab2.dat'
 tg%fnd(3)=trim(fdir)//'tb_tab3.dat'
 tg%fnd(4)=trim(fdir)//'tb_tab4.dat'
 tg%fnd(5)=trim(fdir)//'tb_tab5.dat'


 allocate(bt(tg%nlev))
 do i=1,tg%nlev
  call bt(i)%set(tg%nmodes,tg%lk1(i),tg%lk2(i),tg%bmi(i))
  tg%seed(:,i)=bt(i)%seed
 end do

 call tg%write


end subroutine main_calc


end module nest_grid_m


!program main
! use nest_grid_m, only : main_calc

! call main_calc

!end program main
