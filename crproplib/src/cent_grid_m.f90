module cent_grid_m
 use intpol_mod, only : arr_ind_short, lin_int
 implicit none
 private
 save
 
 public :: cent_grid, reind


 type cent_val
  real :: val
 end type cent_val


 
 type cent_grid
  character(500) :: fname='dat.dat'
  integer :: nx, ny, nz, nf
  real :: xc1(3)
  real, allocatable :: x(:), y(:), z(:)
  type(cent_val), allocatable :: f(:)
 contains
  procedure :: del, set, write=>write_f, read=>read_f
  procedure :: find_cell, add_val, get_val
 end type cent_grid
 
 contains

subroutine del(this) 
 class(cent_grid) :: this
 if (allocated(this%x)) deallocate(this%x)
 if (allocated(this%y)) deallocate(this%y)
 if (allocated(this%z)) deallocate(this%z)
 if (allocated(this%f)) deallocate(this%f)
 this%nx=0
 this%ny=0
 this%nz=0 
end subroutine del 
 
subroutine set(this,nx,ny,nz)
 class(cent_grid) :: this
 integer, intent(in) :: nx, ny, nz
 integer, parameter :: zero=0
 integer :: nmax, bit0, i
 call this%del
 this%nx=nx
 this%ny=ny
 this%nz=nz
 allocate(this%x(0:nx))
 allocate(this%y(0:ny))
 allocate(this%z(0:nz))
 
 nmax=nx
 if (nmax<ny) nmax=ny
 if (nmax<nz) nmax=nz
 bit0=0
 do
  if (ibset(zero,bit0)>=nmax) exit
  bit0=bit0+1  
 end do 
 this%nf=(ibset(zero,bit0))**3
 
 allocate(this%f(this%nf))
 do i=1,this%nf
  this%f(i)%val=0d0
 end do
end subroutine set


subroutine write_f(this)
 integer, parameter :: ndev=25
 class(cent_grid) :: this
 integer :: i, j, k, ib
 
 if (this%fname=='dat.dat') then
  write(*,*) 'The name of the output file is not set'
  return
 end if

 open(ndev,file=this%fname,form='unformatted')
 write(ndev) this%nx
 write(ndev) this%ny
 write(ndev) this%nz
 write(ndev) this%x
 write(ndev) this%y
 write(ndev) this%z
 write(ndev) this%f
 close(ndev)
 
end subroutine write_f



subroutine read_f(this)
 integer, parameter :: ndev=21
 class(cent_grid) :: this
 integer :: i, j, k, ib
 integer :: nx, ny, nz
 
 if (this%fname=='dat.dat') then
  write(*,*) 'The name of the output file is not set'
  return
 end if

 open(ndev,file=this%fname,form='unformatted')
 read(ndev) nx
 read(ndev) ny
 read(ndev) nz
 
 call this%set(nx,ny,nz)
 
 read(ndev) this%x
 read(ndev) this%y
 read(ndev) this%z
 read(ndev) this%f
 close(ndev)
 
end subroutine read_f


subroutine reind(nx0,ny0,nz0,num)
 integer, intent(in) :: nx0, ny0, nz0
 integer, intent(out) :: num
 integer, parameter :: zero=0
 integer :: nmax, bitn, bit0, i
 integer :: nx, ny, nz
 
 nx=nx0-1
 ny=ny0-1
 nz=nz0-1
 
 nmax=nx
 if (nmax<ny) nmax=ny
 if (nmax<nz) nmax=nz
 
 bitn=0
 bit0=0
 num=0
 do
  if (ibset(zero,bit0)>=nmax) exit
  call mvbits(nx,bit0,1,num,bitn)
  bitn=bitn+1
  call mvbits(ny,bit0,1,num,bitn)
  bitn=bitn+1
  call mvbits(nz,bit0,1,num,bitn)
  bitn=bitn+1
  bit0=bit0+1  
 end do 
 
 num=num+1
end subroutine reind


subroutine find_cell(this,xc,fr,nc,xcube)
 class(cent_grid) :: this
 real(8), intent(in) :: xc(3)
 logical, intent(out) :: fr
 integer, intent(out) :: nc(2,3)
 real(8), optional, intent(out) :: xcube(2,3)
 real(8) :: x, y, z
 integer :: i

 x=xc(1)
 y=xc(2)
 z=xc(3)

 if ((x<this%x(0)).or.(x>this%x(this%nx))) then
  fr=.false.
  return
 end if
 
 if ((y<this%y(0)).or.(y>this%y(this%ny))) then
  fr=.false.
  return
 end if
 
 if ((z<this%z(0)).or.(z>this%z(this%nz))) then
  fr=.false.
  return
 end if

 fr=.true.

 call arr_ind_short(0,this%nx,dble(this%x),x,nc(1,1),nc(2,1),0)
 call arr_ind_short(0,this%ny,dble(this%y),y,nc(1,2),nc(2,2),0)
 call arr_ind_short(0,this%nz,dble(this%z),z,nc(1,3),nc(2,3),0)

 if (present(xcube)) then
  xcube(1,1)=this%x(nc(1,1))
  xcube(2,1)=this%x(nc(2,1))

  xcube(1,2)=this%y(nc(1,2))
  xcube(2,2)=this%y(nc(2,2))

  xcube(1,3)=this%z(nc(1,3))
  xcube(2,3)=this%z(nc(2,3))
 end if 


end subroutine find_cell




subroutine add_val(this,xc,val)
 class(cent_grid) :: this
 real(8), intent(in) :: xc(3), val
 integer :: nc(2,3), num
 logical :: fr

 call this%find_cell(xc,fr,nc)
 if (fr) then
  call reind(nc(2,1),nc(2,2),nc(2,3),num)
  this%f(num)%val=this%f(num)%val+val
 end if

end subroutine add_val


subroutine get_val(this,xc,val)
 class(cent_grid) :: this
 real(8), intent(in) :: xc(3)
 real(8), intent(out) :: val
 integer :: nc(2,3), num
 logical :: fr
 
 call this%find_cell(xc,fr,nc)
 if (fr) then
  call reind(nc(2,1),nc(2,2),nc(2,3),num)
  val=this%f(num)%val
 else
  val=0d0
 end if

end subroutine get_val



!subroutine test_cent_grid
! use utools_mod, only : utools
! type(utools) :: ut
! type(cent_grid) :: gr
! real(8), allocatable :: x(:)
! real(8) :: xc(3), val, val1
! integer :: i, nx, ny, nz


! nx=20
! ny=15
! nz=10

! call gr%set(nx,ny,nz)

! call ut%grid(x,0d0,20d0,nx,'lin',0)
! gr%x=x
! call ut%grid(x,0d0,20d0,ny,'lin',0)
! gr%y=x
! call ut%grid(x,0d0,20d0,nz,'lin',0)
! gr%z=x
! 
! xc=[3d0,12d0,10d0]

! val=17d0
! call gr%add_val(xc,1d0)
! call gr%add_val(xc,2d0)
! call gr%add_val(xc,3d0) 

!! xc=[1d0,1d0,1d0]
! call gr%get_val(xc,val1)
! write(*,*) val1
! 

!end subroutine test_cent_grid



!subroutine main_calc

! call test_cent_grid
! 
!end subroutine main_calc



end module cent_grid_m


!program main
! use cent_grid_m, only : main_calc

! call main_calc

!end program main
