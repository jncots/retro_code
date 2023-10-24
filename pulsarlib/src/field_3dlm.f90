module em_fieldm
 type em_field
  real :: bm(3)
 end type em_field

end module em_fieldm


module field_3dlm
 use em_fieldm, only : em_field
 use intpol_mod, only : arr_ind_short, lin_int
 implicit none
 private
 save
 
 public :: field_3dl, reind
 
 type field_3dl
  character(500) :: fname='dat.dat'
  integer :: nx, ny, nz, nf
  real :: xc1(3)
  real, allocatable :: x(:), y(:), z(:)
  type(em_field), allocatable :: f(:)
 contains
  procedure :: del, set, write=>write_f, read=>read_f
  procedure :: cfield
 end type field_3dl
 
 contains

subroutine del(this) 
 class(field_3dl) :: this
 if (allocated(this%x)) deallocate(this%x)
 if (allocated(this%y)) deallocate(this%y)
 if (allocated(this%z)) deallocate(this%z)
 if (allocated(this%f)) deallocate(this%f)
 this%nx=0
 this%ny=0
 this%nz=0 
end subroutine del 
 
subroutine set(this,nx,ny,nz)
 class(field_3dl) :: this
 integer, intent(in) :: nx, ny, nz
 integer, parameter :: zero=0
 integer :: nmax, bit0
 call this%del
 this%nx=nx
 this%ny=ny
 this%nz=nz
 allocate(this%x(nx))
 allocate(this%y(ny))
 allocate(this%z(nz))
 
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
end subroutine set


subroutine write_f(this)
 integer, parameter :: ndev=25
 class(field_3dl) :: this
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
 class(field_3dl) :: this
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


subroutine cfield(this,x,y,z,fl)
 class(field_3dl) :: this
 real(8), intent(in) :: x, y, z
 real(8), intent(out) :: fl(3)
 real(8) :: fpp1(3), fpp2(3), resx(3,2,2), resy(3,2)
 real(8) :: x1, x2, y1, y2, z1, z2
 integer :: i, j, k, im
 integer :: nx(2), ny(2), nz(2), numx(2)
 
 if ((x<this%x(1)).or.(x>this%x(this%nx))) then
  fl=0d0
  return
 end if
 
 if ((y<this%y(1)).or.(y>this%y(this%ny))) then
  fl=0d0
  return
 end if
 
 if ((z<this%z(1)).or.(z>this%z(this%nz))) then
  fl=0d0
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

 this%xc1=[x-x1,y-y1,z-z1]
 

 do k=1,2
  do j=1,2
   do i=1,2
    call reind(nx(i),ny(j),nz(k),numx(i))
   end do
   
   fpp1=this%f(numx(1))%bm(:)
   fpp2=this%f(numx(2))%bm(:)
   
   do im=1,3  
    call lin_int(x1,x2,fpp1(im),fpp2(im),x,resx(im,j,k))
   end do
  end do  
  
  do im=1,3  
    call lin_int(y1,y2,resx(im,1,k),resx(im,2,k),y,resy(im,k))
  end do
 end do
  
 do im=1,3  
  call lin_int(z1,z2,resy(im,1),resy(im,2),z,fl(im))
 end do

end subroutine cfield


subroutine cfield_cubic(this,x,y,z,fl)
  class(field_3dl) :: this
  real(8), intent(in) :: x, y, z
  real(8), intent(out) :: fl(3)
  real(8) :: fpp1(3), fpp2(3), resx(3,2,2), resy(3,2)
  real(8) :: x1, x2, y1, y2, z1, z2
  integer :: i, j, k, im
  integer :: nx(2), ny(2), nz(2), numx(2)
  
  if ((x<this%x(2)).or.(x>this%x(this%nx-1))) then
   fl=0d0
   return
  end if
  
  if ((y<this%y(2)).or.(y>this%y(this%ny-1))) then
   fl=0d0
   return
  end if
  
  if ((z<this%z(2)).or.(z>this%z(this%nz-1))) then
   fl=0d0
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
 
  this%xc1=[x-x1,y-y1,z-z1]
  
 
  do k=1,2
   do j=1,2
    do i=1,2
     call reind(nx(i),ny(j),nz(k),numx(i))
    end do
    
    fpp1=this%f(numx(1))%bm(:)
    fpp2=this%f(numx(2))%bm(:)
    
    do im=1,3  
     call lin_int(x1,x2,fpp1(im),fpp2(im),x,resx(im,j,k))
    end do
   end do  
   
   do im=1,3  
     call lin_int(y1,y2,resx(im,1,k),resx(im,2,k),y,resy(im,k))
   end do
  end do
   
  do im=1,3  
   call lin_int(z1,z2,resy(im,1),resy(im,2),z,fl(im))
  end do
 
 end subroutine cfield_cubic


end module field_3dlm
