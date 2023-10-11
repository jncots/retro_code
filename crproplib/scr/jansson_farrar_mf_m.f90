module jansson_farrar_mf_m
!==========================================================
! Galactic magnetic field model developed by
! Jansson&Farrar(2012):
!
! Jansson, Ronnie, and Glennys R. Farrar. 
! “A New Model of the Galactic Magnetic Field.” 
! The Astrophysical Journal 757 (2012): 14. 
! https://doi.org/10.1088/0004-637X/757/1/14.
!
! Jansson, Ronnie, and Glennys R. Farrar. 
! “The Galactic Magnetic Field.” 
! The Astrophysical Journal Letters 761 (2012): L11. 
! https://doi.org/10.1088/2041-8205/761/1/L11.
!
! For the instructions of usage see
! "subroutine usage_example" in the end of the file
!==========================================================
 use phys_const, only : pi, psec
 use intpol_mod, only : arr_ind_short

 implicit none
 private
 save


 public :: jansson_farrar_mf
 public :: main_calc        ! for test purposes


 real(8), parameter :: kpc=1d3*psec
 real(8), parameter :: pi_180=pi/180, pi2=2*pi
 real(8), parameter :: ezcyl(3)=[0d0,0d0,1d0]


 type jansson_farrar_mf
  logical :: not_init_model=.true.
  real(8) :: rc, zc, azc, irc, phi, ercyl(3), ephi(3) ! see set_coord subroutine
  real(8) :: rpi(11), sin_ia, cos_ia, itan_ia, rpim(40) ! see init_disk subroutine
  real(8) :: reg_b5(11) ! see init_regmf_disk subroutine
  real(8) :: b_north, b_south, r_north, r_south ! see init_regmf_halo subroutine
  real(8) :: iwrad, izhalo, hdisk, iwdisk              ! see init_regmf_halo subroutine
  real(8) :: bxf, cos_tx0, sin_tx0, itan_tx0, rxc, irx     ! see init_regmf_xhalo subroutine
  real(8) :: turb_b5(11), iz0turb_disk                     ! see init_turbmf subroutine
  real(8) :: turbb0_halo, ir0turb_halo,  iz0turb_halo
  real(8) :: beta_str, beta_iso_str
 contains
  procedure :: init_model, set_coord, disk_region
  procedure :: regl_bdisk, regl_bhalo, regl_bxhalo ! components of regular field
  procedure :: regl_bvec, turb_bval, stri_bvec, total_bvec
 end type jansson_farrar_mf


 contains


subroutine set_coord(this,x)
!=====================================================
! Convertion of the cartesian coordinate x(3)=[x,y,z]
! to cylindrical one [rc,phi,zc] and calculation of
! radial and azimuthal unit direction vectors
! ercyl and ephi
!=====================================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in) :: x(3)

 this%rc=sqrt(x(1)**2+x(2)**2)        ! cylindrical radius
 this%zc=x(3)                         ! z coordinate
 this%azc=abs(this%zc)                ! absolute value of z
 this%irc=1/this%rc                   ! inverse of cylindrical radius
 this%phi=acos(x(1)*this%irc)         ! angle 0<phi<pi
 if (x(2)<0d0) this%phi=pi2-this%phi  ! angle 0<phi<2*pi

 this%ercyl=[x(1),x(2),0d0]*this%irc  ! unit vector in radial (cylindrical) direction
 this%ephi=[-x(2),x(1),0d0]*this%irc  ! unit vector in azimuthal direction

 if (this%not_init_model) call this%init_model

end subroutine set_coord


subroutine init_model(this)
!=====================================================
! Initiation of all parameters of the model
!=====================================================
 class(jansson_farrar_mf) :: this

 call init_disk(this)
 call init_regmf_disk(this)
 call init_regmf_halo(this)
 call init_regmf_xhalo(this)
 call init_turbmf(this)


 this%not_init_model=.false.

end subroutine init_model



subroutine init_disk(this)
!=====================================================
! Initiation of all 11 regions in the disk 
! (see comments)
!=====================================================
 class(jansson_farrar_mf) :: this
 real(8) :: ia, ia1, ia2
 integer :: nrot1, nrot2, nspiral, nrot, ig, i1, i2
 real(8) :: rsp

! Disk region borders
! Azimuthally symmetric
 this%rpi(1)=3d0     ! central region: circle with rc<3 kpc
 this%rpi(2)=5d0     ! molecular ring: ring 3 kpc < rc < 5 kpc
! 8 spiral regions 3..10
 this%rpi(3)=5.1d0   ! spiral region border: at phi=pi 5 kpc<rc<5.1 kpc 
 this%rpi(4)=6.3d0   ! spiral region border: at phi=pi 5.1 kpc<rc<6.3 kpc 
 this%rpi(5)=7.1d0   ! and so on ...
 this%rpi(6)=8.3d0
 this%rpi(7)=9.8d0 
 this%rpi(8)=11.4d0
 this%rpi(9)=12.7d0
 this%rpi(10)=15.5d0 ! last spiral border: at phi=pi 12.7 kpc<rc< 15.5 kpc
! Border of the disk (end of the spiral region)
 this%rpi(11)=20d0   ! outside Galaxy 
 this%rpi=this%rpi*kpc ! in kpc


! Spiral parameters
 ia=11.5 ! in degree, spiral opening angle
 ia1=ia*pi_180
 ia2=(90-ia)*pi_180

 this%sin_ia=sin(ia1)
 this%cos_ia=cos(ia1)
 this%itan_ia=1/tan(ia2) ! (i)nverse tan(ia)

! Multiplication of the spirals on phi=pi axis
! for several periods of rotation.
! It is required for the determination 
! of the spiral region in disk_region subroutine
 nrot1=-2 ! from 2*pi*nrot1
 nrot2=2  ! to 2*pi*nrot2
 nspiral=8 ! number of spirals

 do nrot=nrot1,nrot2
  rsp=exp(pi2*nrot*this%itan_ia)
  ig=(nrot-nrot1)*nspiral
  i1=ig+1
  i2=ig+nspiral
  this%rpim(i1:i2)=this%rpi(3:10)*rsp ! multiplication on phi=pi axis
 end do

end subroutine init_disk

function disk_region(this)
!==========================================================
! Determines from the spherical coordinates rc [in cm],
! and phi [counted counterclockwise from x axis]
! number of the spiral (3..10), of gives identifing number
! for outside disc region (11), if r>20 kpc,
! for molecular ring (2), if 3 kpc<r<5 kpc, or center(1)
! if r<3 kpc
! 
! Input:
! ======
! rc [in cm] is the (cylindrical) radius in the Galactic disk
! phi [in radians] is the azimuthal angle counted counterclockwise 
! from x axis 
!
! Output:
! ======
! disk_region=1..11 is the number of the region
! 1 is central region (r<3 kpc)
! 2 is molecular ring (3 kpc<r<5 kpc)
! 3..10 is spiral number (5 kpc<r<20 kpc)
! 11 is the region outside of disk (r>20 kpc)
!
!==========================================================
 class(jansson_farrar_mf), intent(in) :: this
 real(8) :: rc, phi
 integer :: disk_region
 integer, parameter :: nspiral=8 ! number of spirals
 integer :: n1, n2
 real(8) :: r0

 rc=this%rc
 phi=this%phi

 if (rc<this%rpi(1)) then
  disk_region=1
  return
 end if
 
 if (rc<this%rpi(2)) then
  disk_region=2
  return
 end if

 if (rc>this%rpi(11)) then
  disk_region=11
  return
 end if

 r0=rc*exp((pi-phi)*this%itan_ia)
 call arr_ind_short(1,size(this%rpim),this%rpim,r0,n1,n2)
 disk_region=mod(n2,nspiral)
 if (disk_region==0) disk_region=nspiral
 disk_region=disk_region+2    ! take into account two inside regions

end function disk_region


subroutine init_regmf_disk(this)
!=====================================================
! Initiation of the strength of the regular 
! magnetic field in the 11 disk regions normalized 
! at distance 5 kpc (see comments)
!=====================================================
 class(jansson_farrar_mf) :: this

 this%reg_b5(1)=0.1d0      ! b_centre ?
 this%reg_b5(2)=0.1d0      ! b_ring
 this%reg_b5(3)=0.1d0      ! in spirals 3 .. 10
 this%reg_b5(4)=3.0d0
 this%reg_b5(5)=-0.9d0
 this%reg_b5(6)=-0.8d0
 this%reg_b5(7)=-2.0d0
 this%reg_b5(8)=-4.2d0
 this%reg_b5(9)=0d0
 this%reg_b5(10)=2.7d0
 this%reg_b5(11)=0d0       ! outside disk
 this%reg_b5=this%reg_b5*1d-6 ! in microGauss

end subroutine init_regmf_disk


subroutine init_regmf_halo(this)
!=====================================================
! Initiation of the parameters describing the 
! structure of the strength of regular magnetic
! field in the halo (toroidal part)
!=====================================================
 class(jansson_farrar_mf) :: this

 this%b_north=1.4d-6   ! Gauss, northern halo
 this%b_south=-1.1d-6  ! Gauss, southern halo
 this%r_north=9.22*kpc ! transition radius, north
 this%r_south=16.7*kpc ! transition radius, south
 this%iwrad=1/(0.2*kpc)     ! inverse of transition width of the radius
 this%izhalo=1/(5.3*kpc)    ! inverse of vertical scale height of halo

 this%hdisk=0.4*kpc         ! disk/halo transition
 this%iwdisk=1/(0.27*kpc)   ! inverse of transition width

end subroutine init_regmf_halo


subroutine init_regmf_xhalo(this)
!=====================================================
! Initiation of the parameters describing the 
! structure of the strength of regular magnetic
! field in the halo (x-field part)
!=====================================================
 class(jansson_farrar_mf) :: this
 real(8) :: theta_x0

 this%bxf=4.6d-6               ! field strength at origin
 theta_x0=49*pi_180            ! elevated angle at z=0, r>rxc
 this%cos_tx0=cos(theta_x0)
 this%sin_tx0=sin(theta_x0)
  
 this%itan_tx0=1/tan(theta_x0) ! 1/tan of elevated angle at z=0, r>rxc
 this%rxc=4.8*kpc              ! radius where theta_x=tx0
 this%irx=1/(2.9*kpc)          ! inverse of exponential scale length

 
end subroutine init_regmf_xhalo


subroutine init_turbmf(this)
!=====================================================
! Initiation of the parameters describing the 
! structure of the strength of turbulent magnetic
! field in the disk and halo
!=====================================================
 class(jansson_farrar_mf) :: this

! Magnetic field of the disk regions
 this%turb_b5(1)=7.63d0    ! b_centre
 this%turb_b5(2)=7.63d0    ! b_ring
 this%turb_b5(3)=10.81d0
 this%turb_b5(4)=6.96d0
 this%turb_b5(5)=9.59d0
 this%turb_b5(6)=6.96d0
 this%turb_b5(7)=1.96d0
 this%turb_b5(8)=16.34d0
 this%turb_b5(9)=37.29d0
 this%turb_b5(10)=10.35d0
 this%turb_b5(11)=0d0
 this%turb_b5=this%turb_b5*1d-6 ! in microGauss

 this%iz0turb_disk=1/(0.61*kpc)
 
 this%turbb0_halo=4.68d-6 ! random field, halo component
 this%ir0turb_halo=1/(10.97*kpc)
 this%iz0turb_halo=1/(2.84*kpc)

 this%beta_str=1.36d0  ! (B_stri)**2=beta_str*(B_reg)**2
 this%beta_iso_str=sqrt(this%beta_str*3) ! 3 to take into account isotropy of isotropic random field

end subroutine init_turbmf




function flog(zd,hd,iwd)
!==========================================================
! Logistic function flog
! flog approx=0 for zd<hd
! flog approx=1 for zd>hd 
! the smaller wd the faster transition from 0 to 1
!==========================================================
 real(8), intent(in) :: zd, hd, iwd
 real(8)  :: flog

 flog=2*(hd-abs(zd))*iwd ! iwd=1/wd
 flog=1/(1+exp(flog))

end function flog


function regl_bdisk(this,x)
!===========================================
! Calculation of disk component of 
! the regular magnetic field (vector)
!===========================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in), optional :: x(3)
 real(8) :: regl_bdisk(3), res(3), zext
 real(8), parameter :: r5=5*kpc
 integer :: nreg

 if  (present(x)) then
   call this%set_coord(x)
   zext=1-flog(this%zc,this%hdisk,this%iwdisk)      ! extent along z axis
 end if
 
 nreg=disk_region(this)
 if (nreg<3) then
  res=this%reg_b5(nreg)*this%ephi ! constant purely azimuthal field
 else
  res=this%sin_ia*this%ercyl+this%cos_ia*this%ephi ! unit vector along spirals
  res=(this%reg_b5(nreg)*r5*this%irc)*res             ! b5*(r5/rc) changes as 1/r
 end if

 if  (present(x)) then
   regl_bdisk=zext*res
 else
   regl_bdisk=res
 end if
 
end function regl_bdisk


function regl_bhalo(this,x)
!===========================================
! Calculation of halo toroidal component of 
! the regular magnetic field (vector)
!===========================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in), optional :: x(3)
 real(8) :: regl_bhalo(3), zext, btor

 if  (present(x)) then
   call this%set_coord(x)
   zext=flog(this%zc,this%hdisk,this%iwdisk)      ! extent along z axis
 end if

 btor=exp(-this%azc*this%izhalo)
 
 if (this%zc>=0d0) then
  btor=btor*this%b_north*(1-flog(this%rc,this%r_north,this%iwrad))
 else
  btor=btor*this%b_south*(1-flog(this%rc,this%r_south,this%iwrad))
 end if

 regl_bhalo=btor*this%ephi ! purely azimuthal field

 if  (present(x)) regl_bhalo=zext*regl_bhalo
 
end function regl_bhalo


function regl_bxhalo(this,x)
!===========================================
! Calculation of halo poloidal (X) component 
! of the regular magnetic field (vector)
!===========================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in), optional :: x(3)
 real(8) :: regl_bxhalo(3), res(3)
 real(8) :: rad_plus, rxc_at_z
 real(8) :: rpp, bx_val, rpa

 if  (present(x)) call this%set_coord(x)

 rad_plus=this%azc*this%itan_tx0
 rxc_at_z=this%rxc+rad_plus

 if (this%rc>=rxc_at_z) then
  rpp=this%rc-rad_plus
  bx_val=this%bxf*exp(-rpp*this%irx) ! in the disk
  bx_val=bx_val*rpp*this%irc         ! in the halo change as 1/r

  if (this%zc>0) then
   res=this%cos_tx0*this%ercyl+this%sin_tx0*ezcyl ! poloidal field
  else if (this%zc==0) then
   res=ezcyl
  else 
   res=-this%cos_tx0*this%ercyl+this%sin_tx0*ezcyl
  end if
 else
! rc/rxc_at_z find where rc would be on unit segment
! and then lineary project on the segment rxc
  rpp=this%rxc*this%rc/rxc_at_z
  bx_val=this%bxf*exp(-rpp*this%irx)
  bx_val=bx_val*(rpp*this%irc)**2

  rpa=this%rc-rpp
  if (this%zc>0) then
   bx_val=bx_val/sqrt(rpa**2+(this%azc)**2) ! normalization for vector (rpa, zca)
   res=rpa*this%ercyl+this%azc*ezcyl
  else if (this%zc==0) then
   res=ezcyl
  else 
   bx_val=bx_val/sqrt(rpa**2+(this%azc)**2) ! normalization for vector (rpa, zca)
   res=-rpa*this%ercyl+this%azc*ezcyl
  end if
 end if

 regl_bxhalo=bx_val*res

end function regl_bxhalo


function regl_bvec(this,x)
!===========================================
! Calculation of total regular field
!(vector)
!===========================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in), optional :: x(3)
 real(8) :: regl_bvec(3), zext, res(3)

 if  (present(x)) call this%set_coord(x)

 zext=flog(this%zc,this%hdisk,this%iwdisk)      ! extent along z axis
 res=(1-zext)*this%regl_bdisk()+zext*this%regl_bhalo()
 regl_bvec=res+this%regl_bxhalo()
 
end function regl_bvec


function turb_bval(this,x)
!===========================================
! Calculation of strength of turbulent 
! magnetic field
!===========================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in), optional :: x(3)
 real(8) :: turb_bval,  bdisk, bhalo
 real(8), parameter :: r5=5*kpc
 integer :: nreg

 if  (present(x)) call this%set_coord(x)

! Disk component 
 nreg=disk_region(this)
 if (nreg<3) then
  bdisk=this%turb_b5(nreg)
 else
  bdisk=this%turb_b5(nreg)*r5*this%irc
 end if
 bdisk=bdisk*exp(-(this%zc*this%iz0turb_disk)**2) ! z extension

! Halo component
 bhalo=this%turbb0_halo*exp(-this%rc*this%ir0turb_halo) ! along radius
 bhalo=bhalo*exp(-0.5d0*(this%zc*this%iz0turb_halo)**2) ! along z direction

! Total field
 turb_bval=sqrt(bdisk**2+bhalo**2)
 
end function turb_bval


function stri_bvec(this,iso_turb,reg_field)
!============================================================
! Calculation of striated field which is along regular 
! galactic magnetic field but change sign 
! (along or opposite to regular magnetic field). To model 
! anisotropic turbulence we use isotropic turbulence and 
! project is along direction of regular magnetic field
! Then iso_turb is bm/sqrt(<bm^2>). Because we 
! take projection and thus loosing strength, we multiply
! result by sqrt(3) (which is included in this%beta_iso_str)
! variable)
!============================================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in) :: iso_turb(3)   ! bm/sqrt(<bm^2>)
 real(8), intent(in) :: reg_field(3)
 real(8) :: stri_bvec(3)
 real(8) :: regbm, str_val

 regbm=sqrt(dot_product(reg_field,reg_field))

 str_val=dot_product(iso_turb,reg_field)/regbm
 str_val=this%beta_iso_str*str_val
 
 stri_bvec=str_val*reg_field

end function stri_bvec


function total_bvec(this,iso_turb,x)
!===========================================
! Calculation of total magnetic field for
! the jansson_farrar model
!===========================================
 class(jansson_farrar_mf) :: this
 real(8), intent(in) :: iso_turb(3) ! normalized turbulence
 real(8), intent(in), optional :: x(3)
 real(8) :: total_bvec(3)
 real(8) :: reg(3), str(3), turb(3)

 if  (present(x)) call this%set_coord(x)

 reg=this%regl_bvec()
 str=this%stri_bvec(iso_turb,reg)
 turb=iso_turb*this%turb_bval()
 total_bvec=reg+turb+str

end function total_bvec




subroutine usage_example
 type(jansson_farrar_mf) :: jff
 real(8) :: x(3), iso_turb(3), tot(3)
 real(8) :: reg(3), str(3), dir, reg_field(3)

! Coordinates for the calculation:
 x=[-8.5d0,6d0,0d0]*kpc

! First calculate normalized isotropic turbulence bm/sqrt(<bm^2>)
 iso_turb=[-3d0,-2d0,-1d0]*1d-1

! The total field of the model can be calculated as:
 tot=jff%total_bvec(iso_turb,x)
! where iso_turb is vector of isotropic turbulence (normalized) calculated
! at the same point x
 write(*,*) 'total_field=', tot
 write(*,*) 'with strength=', sqrt(dot_product(tot,tot))



! Calculation of components:
 call jff%set_coord(x) ! one can determine position for magnetic field calculation
! just for test
! Then one can find magnetic field in that point separatly for different components
 write(*,*) 'regl_bdisk=', jff%regl_bdisk(x)
 write(*,*) 'regl_bhalo=', jff%regl_bhalo(x)
 write(*,*) 'regl_bxhalo=', jff%regl_bxhalo(x)


! But in real calculations we need only the total regular field
! which can be used in two forms of call:

! The way number 1:
 call jff%set_coord(x) ! position for magnetic field calculation
 write(*,*) 'regular_field=', jff%regl_bvec() ! gives vector of magnetic field

! The way number 2:
  write(*,*) 'regular_field=', jff%regl_bvec(x) ! just put x as argument

! The same way for the strength of turbulent field:
 write(*,*) 'turbulent_field=', jff%turb_bval(x) ! gives root mean square of the strength of 
! turbulent magnetic field

! For striated field, which is along or in opposite direction to regular field we should
! supply iso_turb=bm(x)/sqrt(<bm^2>) - normalized (to rms) vector of isotropic turbulent magnetic field 

 iso_turb=[-3d0,-2d0,-1d0]*1d-1
 reg_field=jff%regl_bvec(x)
 str=jff%stri_bvec(iso_turb,reg_field) ! this is vector of striated field
 write(*,*) 'striated_field=', str

! Now check whether it is along or opposite to regular field:
 reg=jff%regl_bvec(x)
 dir=dot_product(reg,str)/sqrt(dot_product(reg,reg)*dot_product(str,str))
 write(*,*) 'dir=', dir ! should be 1 or -1.


end subroutine usage_example


subroutine main_calc
 call usage_example
end subroutine main_calc



end module jansson_farrar_mf_m


!program main
! use jansson_farrar_mf_m, only : main_calc

! call main_calc

!end program main
