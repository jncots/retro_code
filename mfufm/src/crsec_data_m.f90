module crsec_data_m
 implicit none
 private
 save

 public :: crsec_data

 type crs1d
  integer :: a, n
  real(8), allocatable :: v(:)
 end type crs1d

 type crsec_data
  type(crs1d) :: x(0:3)
  integer :: imax
  character(500) :: fname
 contains
  procedure :: set, del, pos, write=>write_real, read=>read_real
 end type crsec_data


 contains


subroutine read_real(this)
 class(crsec_data) :: this
 integer, parameter :: ndev=35
 integer :: nx, ny, nz, ax(3)
 real, allocatable :: x(:), y(:), z(:), f(:)
  
 open(ndev,file=this%fname,form='unformatted',access='stream',status='old')

 read(ndev) nx, ny, nz
 read(ndev) ax

 allocate(x(nx))
 allocate(y(ny))
 allocate(z(nz))
 allocate(f(nx*ny*nz))

 read(ndev) x, y, z
 read(ndev) f
 close(ndev)

 call this%set(real(x,8),real(y,8),real(z,8),ax)
 this%x(0)%v=real(f,8)
 
  
 deallocate(x)
 deallocate(y)
 deallocate(z)
 deallocate(f)

end subroutine read_real


subroutine write_real(this)
 class(crsec_data) :: this
 integer, parameter :: ndev=35
  
 open(ndev,file=this%fname,form='unformatted',access='stream',status='replace')

 write(ndev) this%x(1)%n, this%x(2)%n, this%x(3)%n
 write(ndev) this%x(1)%a, this%x(2)%a, this%x(3)%a
 write(ndev) real(this%x(1)%v), real(this%x(2)%v), real(this%x(3)%v)
 write(ndev) real(this%x(0)%v)
  
 close(ndev)

end subroutine write_real



subroutine pos(this,iter,x)
 class(crsec_data) :: this
 integer, intent(in) :: iter
 real(8), intent(out) :: x(3)
 integer :: i, nn(3)
 
 call ind_calc3d(iter,this%x(1)%n,this%x(2)%n,this%x(3)%n,nn(1),nn(2),nn(3))
 do i=1,3
  x(this%x(i)%a)=this%x(i)%v(nn(i))
 end do


end subroutine pos

subroutine ind_calc3d(n0,nx,ny,nz,i,j,k)
 integer, intent(in) :: n0, nx, ny, nz
 integer, intent(out) :: i, j, k
 integer :: nc0, nc, nn
 
 nc0=n0
 nn=nx*ny
 call ind_calc(nc0,nn,k,nc)
 nc0=nc
 nn=nx
 call ind_calc(nc0,nn,j,nc)
 nc0=nc
 nn=1
 call ind_calc(nc0,nn,i,nc)

end subroutine ind_calc3d



subroutine ind_calc(n0,nx,i,nr)
 integer, intent(in) :: n0, nx
 integer, intent(out) :: i, nr

 i=n0/nx
 nr=n0-i*nx
 if (nr==0) then
  nr=nx
 else
  i=i+1 
 end if
 
end subroutine ind_calc




subroutine set(this,x1,x2,x3,ax)
 class(crsec_data) :: this
 integer, intent(in) :: ax(3)
 real(8) :: x1(:), x2(:), x3(:)
 integer :: i, ax1(0:3)
 
 call this%del

 this%x(1)%n=size(x1)
 this%x(2)%n=size(x2)
 this%x(3)%n=size(x3)
 this%x(0)%n=this%x(1)%n*this%x(2)%n*this%x(3)%n
 ax1(1:3)=ax
 ax1(0)=0

 do i=0,3
  this%x(i)%a=ax1(i)
  allocate(this%x(i)%v(this%x(i)%n))
 end do
 
 this%x(1)%v=x1
 this%x(2)%v=x2
 this%x(3)%v=x3
 this%imax=this%x(0)%n

end subroutine set


subroutine del(this)
 class(crsec_data) :: this
 integer :: i

 do i=0,3
  if (allocated(this%x(i)%v)) deallocate(this%x(i)%v)
 end do

end subroutine del




end module crsec_data_m

!program main
! use crsec_data_m, only : crsec_data
! use utools_mod, only : utools
! implicit none
! type(crsec_data) :: cs, cs1
! type(utools) :: ut
! real(8), allocatable :: x1(:), x2(:), x3(:)
! real(8) :: xc(3)
! integer :: ax(3), i

!! call ut%grid(x1,0d0,10d0,11,'lin')
!! call ut%grid(x2,0d0,10d0,11,'lin')
!! call ut%grid(x3,0d0,10d0,11,'lin')
!! ax=[3,1,2]

!! call cs%set(x1,x2,x3,ax)
!! 
!! do i=1,cs%imax
!!  call cs%pos(i,xc)
!!  cs%x(0)%v(i)=i
!! end do

!! cs%fname='test.dat'
!! call cs%write
!! call cs%del 

! cs1%fname='test.dat'
! call cs1%read

! do i=1,cs1%imax
!  call cs1%pos(i,xc)
!  write(*,'(10Es14.6)') xc, cs1%x(0)%v(i)
! end do


!end program main
