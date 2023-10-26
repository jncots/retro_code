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
  real(8) :: t3c_coef(4, 4)
  real, allocatable :: x(:), y(:), z(:)
  type(em_field), allocatable :: f(:)
 contains
  procedure :: del, set, write=>write_f, read=>read_f
  procedure :: cfield, cfield_cubic
  procedure :: init_t3c_coef, t3cint, tcint_3d
 end type field_3dl
 
 contains


subroutine init_t3c_coef(this)
  class(field_3dl) :: this

  this%t3c_coef = 0.5d0*transpose(reshape((/0, 2, 0, 0,&
                                     -1, 0, 1, 0,&
                                     2, -5, 4, -1,&
                                     -1, 3, -3, 1/),&
                                     shape(this%t3c_coef)))

end subroutine init_t3c_coef

function t3cint(this, u, p)
  class(field_3dl) :: this
  real(8), intent(in) :: u, p(4)
  real(8) :: t3cint
  real(8) :: u2, u3, uvec(4)

  u2 = u*u
  u3 = u2*u
  uvec = [1d0, u, u2, u3]

  uvec = matmul(uvec, this%t3c_coef)
  t3cint = dot_product(uvec, p)
end function t3cint


function tcint_3d(this, upoint, fun_vals)
  class(field_3dl) :: this
  real(8), intent(in) :: upoint(3), fun_vals(4, 4, 4)
  real(8) :: tcint_3d
  real(8) :: tz(4, 4), ty(4)
  integer :: i, j

  do i=1,4
       do j=1,4
       tz(i, j) = this%t3cint(upoint(3),fun_vals(i, j, :))
       end do
  end do

  do i=1,4
       ty(i) = this%t3cint(upoint(2), tz(i, :))
  end do

  tcint_3d = this%t3cint(upoint(1), ty)
end function tcint_3d

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

 call this%init_t3c_coef
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
  real(8) :: x1, x2, y1, y2, z1, z2
  real(8) :: fun_vals(3, 4, 4, 4), upoint(3)
  integer :: i, j, k, im
  integer :: ix1, ix2, iy1, iy2, iz1, iz2, num
  
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
  
  
  call arr_ind_short(1,this%nx,dble(this%x),x,ix1,ix2)
  call arr_ind_short(1,this%ny,dble(this%y),y,iy1,iy2)
  call arr_ind_short(1,this%nz,dble(this%z),z,iz1,iz2)
  
  x1=this%x(ix1)
  x2=this%x(ix2)
  y1=this%y(iy1)
  y2=this%y(iy2)
  z1=this%z(iz1)
  z2=this%z(iz2)

  upoint = [(x-x1)/(x2-x1), (y-y1)/(y2-y1), (z-z1)/(z2-z1)]

  do i = 1, 4
    do j = 1, 4
         do k = 1, 4
           do im=1,3
            call reind(ix1-2+i,iy1-2+j,iz1-2+k, num)
            fun_vals(im, i, j, k) = &
            this%f(num)%bm(im)
           end do  
         end do
    end do
  end do
  
  do im=1,3
   fl(im) = this%tcint_3d(upoint, fun_vals(im,:,:,:))
  end do 
  
 end subroutine cfield_cubic


end module field_3dlm
