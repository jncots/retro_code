module tabfun_2d_mod
 use intpol_mod, only : arr_ind_short, logl_int
 implicit none
 private
 save
 
 public :: tabfun_2d
 
 
 abstract interface
  subroutine fun_2d(x,y,res)
   real(8), intent(in) :: x, y
   real(8), intent(out) :: res
  end subroutine fun_2d
 end interface
 
 
 
 type tabfun_2d
  integer :: nx, ny
  integer, private :: i1=1, i2=2
! interpolation (wraping up) goes in order: first i1, then i2 
  real(8) :: x1, x2, y1, y2
  real(8), pointer :: xa(:), ya(:), fa(:,:)
 contains
  procedure :: set, fxy
 end type tabfun_2d
 
 
 contains


subroutine set(this,x,y,f)
 class(tabfun_2d) :: this
 real(8), intent(in), target :: x(:), y(:), f(:,:)
 integer :: nxf, nyf
 
 this%nx=size(x)
 this%ny=size(y)
 nxf=size(f,1)
 nyf=size(f,2)
 
 if (this%nx/=nxf) then
  write(*,*) 'tabfun_2d: 1st dimention of data array f is not consistent with 1st grid'
  return
 end if
 
 if (this%ny/=nyf) then
  write(*,*) 'tabfun_2d: 2nd dimention of data array f is not consistent with 2nd grid'
  return
 end if
 
 this%xa=>x
 this%ya=>y
 this%fa=>f
 
 this%x1=this%xa(1)
 this%x2=this%xa(this%nx)
 this%y1=this%ya(1)
 this%y2=this%ya(this%ny)
end subroutine set


! subroutine order(this,i1,i2)
 ! class(tabfun_2d) :: this
 ! integer, intent(in) :: i1, i2
 ! integer :: ss(2), i
 
 ! ss=[i1,i2]
 
 ! if (i1==i2) then
  ! write(*,*) 'tabfun_2d: in order(i1,i2), i1, i2 should be different'
  ! return
 ! end if
 
 ! do i=1,2
  ! if ((ss(i)<1).or.(ss(i)>2)) then
   ! write(*,*) 'tabfun_2d: in order(i1,i2), i1, i2 should be 1 or 2'
   ! return
  ! end if 
 ! end do
 
 ! this%i1=i1
 ! this%i2=i2
! end subroutine order 
 

 
subroutine fxy(this,x,y,res)
 class(tabfun_2d) :: this
 real(8), intent(in) :: x, y
 real(8), intent(out) :: res
 real(8) :: va(2,2), v(2), res1(2)
 integer :: na(2,2), r1(2), r2(2)
 integer :: i1, i2, k
 
 if ((x<this%x1).or.(x>this%x2)) then
  res=0d0
  return
 end if 
 
 if ((y<this%y1).or.(y>this%y2)) then
  res=0d0
  return
 end if 
 
 
 v=[x,y]
 
 call arr_ind_short(1,this%nx,this%xa,x,na(1,1),na(2,1))
 va(1,1)=this%xa(na(1,1))
 va(2,1)=this%xa(na(2,1))
 
 call arr_ind_short(1,this%ny,this%ya,y,na(1,2),na(2,2))
 va(1,2)=this%ya(na(1,2))
 va(2,2)=this%ya(na(2,2))
 
 i1=this%i1
 i2=this%i2
 
 do k=1,2
  r1(i1)=na(1,i1)
  r2(i1)=na(2,i1)
  r1(i2)=na(k,i2)
  r2(i2)=na(k,i2)
  
  call logl_int(va(1,i1),va(2,i1),&
	   this%fa(r1(1),r1(2)),this%fa(r2(1),r2(2)),v(i1),res1(k))	   
 end do
 
 call logl_int(va(1,i2),va(2,i2),&
	res1(1),res1(2),v(i2),res)

 if (isnan(res)) then
  write(*,*) 'tabfun_2d%fxy, res=NaN'
  write(*,'(3(A,Es14.6))') 'xmin=', this%xa(1), '  xmax=', this%xa(this%nx),'  x=', x
  write(*,'(3(A,Es14.6))') 'ymin=', this%ya(1), '  xmax=', this%ya(this%ny),'  y=', y
 end if 

end subroutine fxy


end module tabfun_2d_mod


! program main
 ! use taber, only : tab_er
 ! use tabfun_2d_mod, only : tabfun_2d
 ! implicit none
 ! type(tab_er) :: tt
 ! type(tabfun_2d) :: ft
 ! integer :: ie, ir
 ! real(8) :: res
 
 ! call tt%set(105,234)
 ! call tt%set_ea(1d0,100d0)
 ! call tt%set_ra(1d-2,10d0)
 
 ! call srand(309)
 ! do ir=1,tt%nr
  ! do ie=1,tt%ne
   ! tt%fa(ie,ir)=rand()
  ! end do
 ! end do 

 
 ! call ft%set(tt%ea,tt%ra,tt%fa)
 ! call ft%fxy(4.3d0,3.34d0,res) 
 ! write(*,*) res
 
 

! end program main