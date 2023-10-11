module cent_grida_m
 use intpol_mod, only : arr_ind_short, lin_int
 implicit none
 private
 save
 
 public :: cent_grida


 type cent_val
  real :: val
 end type cent_val


 
 type cent_grida
  character(500) :: fname='dat.dat'
  integer :: nx, ny, nz
  real :: xc1(3)
  real, allocatable :: x(:), y(:), z(:)
  type(cent_val), allocatable :: f(:,:,:)
 contains
  procedure :: del, set, write=>write_f, read=>read_f
  procedure :: find_cell, add_val, get_val
  procedure :: get_den, norm_den
 end type cent_grida
 
 contains

subroutine del(this) 
 class(cent_grida) :: this
 if (allocated(this%x)) deallocate(this%x)
 if (allocated(this%y)) deallocate(this%y)
 if (allocated(this%z)) deallocate(this%z)
 if (allocated(this%f)) deallocate(this%f)
 this%nx=0
 this%ny=0
 this%nz=0 
end subroutine del 
 
subroutine set(this,nx,ny,nz)
 class(cent_grida) :: this
 integer, intent(in) :: nx, ny, nz
 integer, parameter :: zero=0
 integer :: i, j, k
 call this%del
 this%nx=nx
 this%ny=ny
 this%nz=nz
 allocate(this%x(0:nx))
 allocate(this%y(0:ny))
 allocate(this%z(0:nz))

 allocate(this%f(this%nx,this%ny,this%nz))

 do k=1,this%nz
  do j=1,this%ny
   do i=1,this%nx
    this%f(i,j,k)%val=0d0
   end do
  end do
 end do
end subroutine set


subroutine write_f(this)
 integer, parameter :: ndev=25
 class(cent_grida) :: this
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
 class(cent_grida) :: this
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


subroutine find_cell(this,xc,fr,nc,xcube)
 class(cent_grida) :: this
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
 class(cent_grida) :: this
 real(8), intent(in) :: xc(3), val
 integer :: nc(2,3), i, j, k
 logical :: fr

 call this%find_cell(xc,fr,nc)
 if (fr) then
  i=nc(2,1)
  j=nc(2,2)
  k=nc(2,3)
  this%f(i,j,k)%val=this%f(i,j,k)%val+val
 end if

end subroutine add_val


subroutine get_val(this,xc,val)
 class(cent_grida) :: this
 real(8), intent(in) :: xc(3)
 real(8), intent(out) :: val
 integer :: nc(2,3), i, j, k
 logical :: fr
 
 call this%find_cell(xc,fr,nc)
 if (fr) then
  i=nc(2,1)
  j=nc(2,2)
  k=nc(2,3)
  val=this%f(i,j,k)%val
 else
  val=0d0
 end if

end subroutine get_val


subroutine get_den(this,xc,val,units,norm)
!================================================
! Calculation of the density (data normalized by
! volume of the cell). The data could be 
! additionally normalized using 'norm' parameter
!================================================
 class(cent_grida) :: this
 real(8), intent(in) :: xc(3) 
 real(8), intent(out) :: val
 real(8), intent(in), optional :: units, norm
 integer :: nc(2,3), i, j, k
 real(8) :: xcube(2,3), dx, dy, dz, vol, anorm
 logical :: fr

 if (present(units)) then
  anorm=units**3
 else
  anorm=1d0
 end if

 if (present(norm)) anorm=anorm*norm


 call this%find_cell(xc,fr,nc,xcube)
 if (fr) then
  i=nc(2,1)
  j=nc(2,2)
  k=nc(2,3)
  
  dx=xcube(2,1)-xcube(1,1)
  dy=xcube(2,2)-xcube(1,2)
  dz=xcube(2,3)-xcube(1,3)
  vol=dx*dy*dz
  val=this%f(i,j,k)%val*(anorm/vol)

  if (val<0d0) then
    val=0d0
!    write(*,*) 'this%nx, this%ny, this%nz', this%nx, this%ny, this%nz
!    write(*,*) 'i, j, k', i, j, k
!    write(*,*) 'val=', this%f(i,j,k)%val
  end if
 else
  val=0d0
 end if

end subroutine get_den

subroutine norm_den(this,units,norm)
!================================================
! Transformation of data to the density 
! (data normalized by volume of the cell). 
! The data could be additionally normalized 
! using 'norm' parameter
!================================================
 class(cent_grida) :: this
 real(8), intent(in) :: units, norm
 real(8) :: anorm, dx, dy, dz, vol
 integer :: i, j, k

 anorm=norm*units**3
 do k=1,this%nz
  dz=this%z(k)-this%z(k-1)
  do j=1,this%ny
   dy=this%y(j)-this%y(j-1)
   do i=1,this%nx
    dx=this%x(i)-this%x(i-1)
    vol=dx*dy*dz
    this%f(i,j,k)%val=this%f(i,j,k)%val*(anorm/vol)
   end do
  end do
 end do
end subroutine norm_den



!subroutine test_cent_grid
! use utools_mod, only : utools
! type(utools) :: ut
! type(cent_grida) :: gr
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



end module cent_grida_m


!program main
! use cent_grida_m, only : main_calc

! call main_calc

!end program main
