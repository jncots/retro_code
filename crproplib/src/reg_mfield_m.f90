module reg_mfield_m
 use phys_const, only : c_light, pi, psec
 use utools_mod, only : utools

 implicit none
 private
 save

 public :: reg_mfield !,  main_calc


 type reg_mfield
  real(8) :: itanp, phi_par, rcent, br_val, iz0, ppar
  real(8) :: bt0(4), rt0(4), ht(4)
  character(3) :: ftype ! ass, bss
 contains
  procedure :: set, bm_str, bm=>bm_field
 end type reg_mfield
 
 
 real(8), parameter :: kpsec=1d3*psec
 real(8), parameter :: rsun=8.5*kpsec
 real(8), parameter :: degree=pi/180, pi2=2*pi
 real(8), parameter :: p_ass=-5*degree, p_bss=-6*degree


 contains


subroutine set(this,ftype)
 class(reg_mfield) :: this
 character(3) :: ftype
 
 this%ftype=ftype
 
 if (ftype=='ass') then
  call calc_par(this,p_ass)
 else if (ftype=='bss') then
  call calc_par(this,p_bss)
 else
  write(*,*) 'There is no such model as '//ftype
 end if
 
end subroutine set

subroutine calc_par(this,ppar)
 class(reg_mfield) :: this
 real(8), intent(in) :: ppar
 real(8) :: itanp, phi_par, rcent, br_val, iz0
 real(8) :: dpar, br0

 dpar=-0.6*kpsec
 itanp=1/tan(ppar)
 phi_par=log(1+dpar/rsun)*itanp-pi/2

 br0=2d-6 ! Gauss
 rcent=5*kpsec
 br_val=br0*rsun/(rcent*cos(phi_par))
 iz0=1/kpsec


 this%itanp=itanp
 this%phi_par=phi_par
 this%rcent=rcent
 this%br_val=br_val
 this%iz0=iz0
 this%ppar=ppar


 this%bt0=[4d0,2d0,4d0,4d0]*1d-6 
 this%rt0=[6d0,6d0,6d0,5d0]*kpsec
 this%ht=[1.3d0,1.3d0,1.5d0,1.5d0]*kpsec

end subroutine calc_par



function bm_str(this,rc,theta,zc)
 class(reg_mfield), intent(in) :: this
 real(8), intent(in) :: rc, theta, zc
 real(8) :: bm_str, cph, bm_rad

 if (rc<this%rcent) then
   bm_rad=this%br_val
 else
   bm_rad=this%br_val*this%rcent/rc
 end if

 cph=cos(theta-this%itanp*log(rc/rsun)+this%phi_par)
 if (this%ftype=='ass') cph=abs(cph)
 bm_str=bm_rad*cph*exp(-abs(zc)*this%iz0)

end function bm_str


function bm_field(this,x)
 class(reg_mfield) :: this
 real(8), intent(in) :: x(3)
 real(8) :: bm_field(3)
 real(8) :: r, theta, z
 real(8) :: ir, ct, st, b0, br, bt

 r=sqrt(x(1)**2+x(2)**2)
 ir=1/r
 ct=x(2)*ir
 st=x(1)*ir

 theta=acos(ct)
 if (x(1)<0d0) theta=pi2-theta
 z=x(3)

 b0=this%bm_str(r,theta,z)
 br=b0*sin(this%ppar)
 bt=b0*cos(this%ppar)

 bm_field(1)=br*st+bt*ct
 bm_field(2)=br*ct-bt*st
 bm_field(3)=0d0

 bm_field=bm_field+halo_field(this,r,ct,st,z)

end function bm_field



function halo_field(this,r,ct,st,z)
 class(reg_mfield) :: this
 real(8), intent(in) :: r, ct, st, z
 real(8) :: halo_field(3)
 integer :: im
 real(8) :: az, rrt, bt, ht, wt, ww


 if (this%ftype=='ass') then
  if (z>0d0) then
   im=1
  else
   im=2
  end if
 else if (this%ftype=='bss') then
  if (z>0d0) then
   im=3
  else
   im=4
  end if
 else
  write(*,*) 'There is no such model as '//this%ftype
 end if
 
 az=abs(z)

 rrt=r/this%rt0(im)
 bt=this%bt0(im)*rrt*exp(1-rrt)
 ht=this%ht(im) 

 if (az<=ht) then
  wt=0.25*kpsec
 else
  wt=0.4*kpsec
 end if

 ww=(az-ht)/wt
 bt=bt/(1+ww**2)

 bt=sign(bt,z)

 halo_field(1)=-bt*ct
 halo_field(2)=bt*st
 halo_field(3)=0d0

end function halo_field


!subroutine main_calc
! 
! call plot_field
!end subroutine main_calc


!subroutine plot_field
! type(reg_mfield) :: rb
! type(utools) :: ut
! integer :: nx, ny, nreg, i, j
! real(8), allocatable :: x(:), y(:)
! real(8) :: reg, xc(3)
! character(500) :: fdir, fname

! nreg=100
! nx=nreg
! ny=nreg
! reg=10*kpsec
! call ut%grid(x,-reg,reg,nx,'lin')
! call ut%grid(y,-reg,reg,ny,'lin')

! call rb%set('ass')

! fdir='/home/anton/work/project/&
! halo/chp_prop/reg_field/'

! fname=trim(fdir)//'ass_field.dat'

! open(1,file=fname)

! do i=1,nx
!  do j=1,ny
!   xc=[x(i),y(j),0d0]
!   write(1,*) rb%bm(xc)
!  end do
! end do

! close(1)

! fname=trim(fdir)//'xarr.dat'
! open(1,file=fname)
! do i=1,size(x)
!  write(1,*) x(i)/kpsec
! end do
! close(1)

! fname=trim(fdir)//'yarr.dat'
! open(1,file=fname)
! do i=1,size(y)
!  write(1,*) y(i)/kpsec
! end do
! close(1)


!end subroutine plot_field




end module reg_mfield_m


!program main
! use reg_mfield_m, only : main_calc

! call main_calc

!end program main
