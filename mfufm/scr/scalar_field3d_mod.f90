module scalar_field3d_mod
 implicit none
 private
 save

 public :: scalar_field3d

 type scalar_field3d
  integer :: nx, ny, nz
  real, allocatable :: x(:), y(:), z(:), f(:,:,:)
  character(500) :: fname
 contains
  procedure :: del, set, write=>write_data, read=>read_data
  procedure :: sfield
 end type scalar_field3d




 contains



subroutine del(this)
 class(scalar_field3d) :: this

 this%nx=0
 this%ny=0
 this%nz=0
 if (allocated(this%x)) deallocate(this%x)
 if (allocated(this%y)) deallocate(this%y)
 if (allocated(this%z)) deallocate(this%z)
 if (allocated(this%f)) deallocate(this%f)
 
end subroutine del

subroutine set(this,nx,ny,nz)
 class(scalar_field3d) :: this
 integer, intent(in) :: nx, ny, nz

 call this%del
 this%nx=nx
 this%ny=ny
 this%nz=nz

 allocate(this%x(nx))
 allocate(this%y(ny))
 allocate(this%z(nz))
 allocate(this%f(nx,ny,nz))
 
end subroutine set


subroutine write_data(this)
 class(scalar_field3d) :: this
 integer, parameter :: ndev=25

 open(ndev,file=this%fname,form='unformatted',access='stream',status='replace')
 write(ndev) this%nx, this%ny, this%nz
 write(ndev) this%x, this%y, this%z
 write(ndev) this%f
 close(ndev)

end subroutine write_data


subroutine read_data(this)
 class(scalar_field3d) :: this
 integer, parameter :: ndev=25
 integer :: nx, ny, nz

 open(ndev,file=this%fname,form='unformatted',access='stream',status='old')

 read(ndev) nx, ny, nz
 call this%set(nx,ny,nz)
 read(ndev) this%x, this%y, this%z
 read(ndev) this%f
 close(ndev)

end subroutine read_data



subroutine sfield(this,x3,res)
 use intpol_mod, only : arr_ind_short, lin_int
 class(scalar_field3d) :: this
 real(8), intent(in) :: x3(3)
 real(8), intent(out) :: res
 real(8) :: x, y, z
 real(8) :: fpp1(2), resx(2,2), resy(2)
 real(8) :: x1, x2, y1, y2, z1, z2
 integer :: i, j, k
 integer :: nx(2), ny(2), nz(2), numx(2)
 

 x=x3(1)
 y=x3(2)
 z=x3(3)
 if ((x<this%x(1)).or.(x>this%x(this%nx))) then
  res=0d0
  return
 end if
 
 if ((y<this%y(1)).or.(y>this%y(this%ny))) then
  res=0d0
  return
 end if
 
 if ((z<this%z(1)).or.(z>this%z(this%nz))) then
  res=0d0
  return
 end if
 
 
 call arr_ind_short(1,this%nx,dble(this%x),x,nx(1),nx(2))
 call arr_ind_short(1,this%ny,dble(this%y),y,ny(1),ny(2))
 call arr_ind_short(1,this%nz,dble(this%z),z,nz(1),nz(2))
 
 x1=this%x(nx(1))
 x2=this%x(nx(2))
 y1=this%y(ny(1))
 y2=this%y(ny(2))
 z1=this%z(nz(1))
 z2=this%z(nz(2))
 

 do k=1,2
  do j=1,2
   do i=1,2   
    fpp1(i)=this%f(nx(i),ny(j),nz(k))
   end do
   call lin_int(x1,x2,fpp1(1),fpp1(2),x,resx(j,k))
  end do
  
  call lin_int(y1,y2,resx(1,k),resx(2,k),y,resy(k))
 end do 
  
 call lin_int(z1,z2,resy(1),resy(2),z,res)

end subroutine sfield


end module scalar_field3d_mod



!program main
! use scalar_field3d_mod, only : scalar_field3d
! type(scalar_field3d) :: cfl
! character(500) :: fdir, fname
! real(8) :: x(3), res

! fdir='/home/anton/work/project/pulsar_escape&
! /mag_field/field_lines1/test_check/flines3d/'

! fname=trim(fdir)//'last_fline_300.dat'

! cfl%fname=fname
! call cfl%read

! x=[0.3d0,0.5d0,1.4d0]

! call cfl%sfield(x,res)


! write(*,*) res

! 
!end program main
