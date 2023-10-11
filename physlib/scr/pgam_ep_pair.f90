module pgam_ep_pair
 use gauss_kronrod, only : gk_adaptive
 use intpol_mod, only : bs_root
 use timer_module, only : timer_class
 use description_line, only : descr_line, infoline
 use word_processor, only : line_break
 use house_keeping, only : get_unit
 implicit none
 private
 save
 
 
 public :: beheit_tab
 
 real(8), parameter :: mp=0.938272046d9 ! proton mass in eV
 real(8), parameter :: me=0.51099891d6 ! electron mass in eV
 real(8), parameter :: mep=me/mp
 real(8) :: omega, xi, rho, rho1, rho2, gp
 real(8) :: elen, eps
! norm=(cs*hbar/(me*cs^2))^3*(me*cs^2/hbar), [cm^3/s]
 real(8) :: norm=4.4704745616869275d-11
 real(8) :: eph1, eph2
 
 integer :: ndes
 character(len=:), allocatable :: description
 
 
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type beheit_tab
  integer :: ne, np
  real(8) :: eph1, eph2
  real(8) :: ee1, ee2
  real(8) :: ep1, ep2
  real(8), allocatable :: eea(:), epa(:,:), rate(:,:)
  procedure(fun_1d), pointer, nopass :: sph=>null()
  type(infoline) :: info
  type(infoline), allocatable :: inf(:)
  character(500) :: file_out
 contains
  procedure :: phot=>sph_set, elec=>elec_set, prot=>prot_set
  procedure :: calc=>rate_per_proton, add_info
 end type beheit_tab
 
 
 procedure(fun_1d), pointer :: sph_dist
 
 
 contains


!================================================
! Procedures for writing information about
! calculation in the output file
!================================================
subroutine add_info(this)
! Variables 
 class(beheit_tab) :: this
 type(infoline), allocatable :: inf1(:)
 integer :: nsize, i
! Calculations 
 if (allocated(this%inf)) then
  nsize=size(this%inf)
  allocate(inf1(nsize))
  inf1=this%inf
  deallocate(this%inf)
  allocate(this%inf(nsize+1))
  do i=1,nsize
   this%inf(i)=inf1(i)
  end do
  this%inf(nsize+1)=this%info
  deallocate(inf1)
 else
  allocate(this%inf(1))
  this%inf(1)=this%info
 end if
end subroutine add_info


subroutine write_info(inf)
!
! Write information of calculation into description(ndes) variable
!
! Variables
 type(infoline), intent(in) :: inf(:)
 type(infoline), allocatable :: note(:)
 integer :: nleng, ninf, nnote, i
 character(:), allocatable :: line
! Calculations
 nleng=70
 line='This file contains the rate of electron-positron pair production in interactions &
 of high-energy protons with a specific soft photon field (see description below).'
 call line_break(nleng,line)
 
 ninf=size(inf)
 nnote=ninf+1
 
 allocate(note(nnote))
 
 note(1)%typ='des'
 note(1)%intro='Introduction'
 note(1)%info=line
 
 do i=2,nnote
  note(i)=inf(i-1)
 end do
 
 call descr_line(nnote,note,ndes,description)

end subroutine write_info


!================================================
! Procedures for setting parameters of
! calculation
!================================================
subroutine sph_set(this,eph1,eph2,sph)
! Variables 
 class(beheit_tab) :: this
 real(8), intent(in) :: eph1, eph2
 procedure(fun_1d) :: sph
! Calculations
 this%eph1=eph1
 this%eph2=eph2
 this%sph=>sph
end subroutine sph_set

subroutine elec_set(this,ee1,ee2,ne)
! Variables 
 class(beheit_tab) :: this
 real(8), intent(in) :: ee1, ee2
 integer, intent(in) :: ne
! Calculations
 this%ee1=ee1
 this%ee2=ee2
 this%ne=ne
end subroutine elec_set

subroutine prot_set(this,ep1,ep2,np)
! Variables 
 class(beheit_tab) :: this
 real(8), intent(in) :: ep1, ep2
 integer, intent(in) :: np
! Calculations
 this%ep1=ep1
 this%ep2=ep2
 this%np=np
end subroutine prot_set

 
 
!================================================
! Procedures for calculation
!================================================ 
subroutine rate_per_proton(this)
! Variables
 class(beheit_tab) :: this
 integer :: ndev
 integer :: i, j, nee, nep, ne, np, npl
 real(8) :: ee01, ee02, ee1, ee2, ee_min,  ep1, ep2
 real(8), allocatable :: eea(:), epa(:,:), rate(:,:)
 type(timer_class) :: tmr
! Calculations
 
 
 call write_info(this%inf)
 write(*,*) description
 
 sph_dist=>this%sph
 eph1=this%eph1
 eph2=this%eph2
! Definition of proton energy range
 nep=this%np
 ep2=this%ep2
 
! Definition of electron energy range 
 nee=this%ne
 ee1=this%ee1
 ee2=this%ee2
 call el_min_ph(eph2,ee_min)
 if (ee1<ee_min) ee1=ee_min
 if (ee1>=ee2) then
  write(*,'(A,10Es14.6)') 'ee1>=ee2, ee1, ee2=',ee1, ee2
 end if
 
 ne=nee+1
 np=nep+1
 allocate(this%eea(ne))
 allocate(this%epa(np,ne))
 allocate(this%rate(np,ne))
 this%ne=ne
 this%np=np
 
 
 do i=0,nee
  this%eea(i+1)=ee1*(ee2/ee1)**(i*1d0/nee)
  call prot_min_en(eph2,this%eea(i+1),ep1)
  ep1=2*ep1 
  do j=0,nep
   this%epa(j+1,i+1)=ep1*(ep2/ep1)**(j*1d0/nep)
  end do 
 end do
 
 
 tmr%proc_name='rate_per_proton'
 call tmr%start(ne*np,1.0)
 
 do i=1,ne
  do j=1,np
   call ee_rate(this%eea(i),this%epa(j,i),this%rate(j,i))
!   write(*,'(A,2I5,10Es14.6)') 'i, j, eea(i),epa(j,i),rate(j,i)', i, j, eea(i),epa(j,i),rate(j,i)
   call tmr%loop
  end do
 end do
 
 npl=0
 call get_unit(ndev)
 open(ndev,file=this%file_out,form='unformatted')
 write(ndev) ndes
 write(ndev) description
 write(ndev) npl 
 write(ndev) ne
 write(ndev) np
 write(ndev) this%eea
 write(ndev) this%epa
 write(ndev) this%rate
 close(ndev) 
 
end subroutine rate_per_proton



subroutine ee_rate(ee,ep,res)
!
!==============================================
! Production rate of electrons+positrons
! (see explanation below, after !!!)
! ee, ep electron and proton energies in eV
! res is rate in 1/(eV*cm^3)
!==============================================
! Variables
 real(8), intent(in) :: ee, ep
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations

 call eph_min_en(ee,ep,gp,xi,rho)
 rho1=2*gp*rho ! lower limit for integration over omega ('isph' subroutine)
 rho2=rho1-1   ! lower limit for integration over emin  ('iemin' subroutine)
 
 tmin=eph1/me
 if (tmin<rho) tmin=rho
 tmax=eph2/me

 if (tmin>=tmax) then
  res=0d0
  write(*,*) 'ee_rate: res=0'
 else
  tmin=log(tmin)
  tmax=log(tmax)
  call gk_adaptive(isph,tmin,tmax,res)
  res=norm*res/gp**3 !!! intead of 1/(2*gp**3) since we multiply it by 2 to take
                     !!! both electrons and positrons into account (electrons+positrons)
 end if

end subroutine ee_rate

subroutine isph(x,res)
! me, rho1, gp, sph_dist - parameters
! Variables
 real(8), intent(in) :: x ! log(eph)
 real(8), intent(out) :: res
 real(8) :: eph, fsph, tmin, tmax
! Calculations
 eph=exp(x)
 call sph_dist(eph*me,fsph)
! if (fsph<1d-100) then
!  write(*,*) eph*me, fsph
! end if
 tmin=rho1
 tmax=2*gp*eph
 call gk_adaptive(iomega,tmin,tmax,res)
 res=res*fsph/eph
end subroutine isph

subroutine iomega(x,res)
! rho2 - parameter
! Variables
 real(8), intent(in) :: x ! omega
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations
 omega=x
 tmin=rho2
 tmax=omega-1
 call gk_adaptive(iemin,tmin,tmax,res)
 res=omega*res
end subroutine iomega
      
      
subroutine iemin(x,res)
! xi, omega - parameters
! Variables
 real(8), intent(in) :: x ! is the energy of electron in proton reference frame
 real(8), intent(out) :: res
 real(8) :: cos0
! Calculations
 cos0=(x-xi)/sqrt(x**2-1)
 call wcross(omega,x,cos0,res)
end subroutine iemin


subroutine wcross(omega,Em,co,res)
 implicit real(8) (a-h,o-y)
! Variables
 real(8), intent(in) :: omega, em, co
 real(8), intent(out) :: res
! Calculations
 E2 = omega-Em
 t1 = Em**2
 t2 = t1-1
 P1 = dsqrt(t2)
 t3 = omega**2
 BBB=E2**2-1
 if (BBB.lt.0d0) then
  P2 = 0d0
 else
  P2 = dsqrt(BBB)
 end if
 xx=1d0/t1
 if (xx.lt.1d-4) then
  ff=xx/2d0+xx**2/8d0+xx**3/16d0
  ff1=1d0-co
  if (ff1.lt.0d0) ff1=0d0 
  fff=ff1+ff*co
 else
  fff=1d0-dsqrt(1d0-xx)*co
 end if
 DD = Em*fff
 www = 2d0*omega*Em*fff
 t12 = P2**2+www
 TT = dsqrt(t12)
 vvv=(TT+P2)**2
 del = dlog(vvv/www)
 t17 = 1/t2
 t18 = E2*Em
 t19 = P2*P1
 t23 = dlog((t18+t19+1)/omega)
 YY = 2*t17*t23
 t25 = 1/P2
 t27 = dlog(E2+P2)
 t28 = t25*t27
 Y1 = 2*t28
 t32 = 1-co**2
 t36 = DD**2
 t37 = t36**2
 t38 = 1/t37
 t46 = 1/t36
 t48 = t1-1-t3
 t49 = 1/t12
 t53 = 1/DD
 t66 = E2**2
 SS = t19/t3/omega*(-4*t32*(2*t1+1)*t17*t38+(5*t1-2*t18+3)*t17*t46+&
 t48*t49*t46+2*E2*t17*t53+2/P1/t2*t23*t25*(2*Em*t32*(3*omega+t2*E2)&
 *t38+(2*t1*(t1+t66)-7*t1-3*t18-t66+1)*t46+omega*(t1-t18-1)*t53)-&
 del*t25/TT*(2*t46-3*omega*t53-omega*t48*t49*t53)-4*t28*t53) !  /2 two is taken into account in (\alpha^3)/2
 SS = SS/P1*1.943D-7 ! cross section is devided on P1=Pminus  and * (\alpha^3)/2
 res=SS
end subroutine wcross

!================================================
! Procedures for min/max energy calculations
!================================================
subroutine el_min_ph(eph,ee_min)
!================================================
! Calculation of minimal energy that can be
! produced by photon with energy eph
!================================================
! Variables
 real(8), intent(in) :: eph
 real(8), intent(out) :: ee_min
! Calculations
  ee_min=me**2/((1+mep)*eph)
end subroutine el_min_ph


subroutine eph_min_en(ee,ep,gp,xi,rho)
!================================================
! Calculations of minimal photon energy 'rho'
! (in units of me) for production of electron with
! energy 'ee' by proton with energy 'ep'
! All energies in eV
!
! xi=ge/gp, rho=eph_min/me
!================================================
! Variables
 real(8), intent(in) :: ee, ep
 real(8), intent(out) :: gp, xi, rho
 real(8) :: ge
! Calculations
 ge=ee/me
 gp=ep/mp
 xi=ge/gp
 rho=4*mep*(1+1/xi)+1
 rho=(sqrt(rho)+1)*xi
 rho=1+2/rho
 rho=(1+xi)/(1-xi*mep)*rho
 rho=rho/(4*gp)
end subroutine eph_min_en



subroutine el_min_max_en(eph,ep,ee_min,ee_max)
!=======================================================
! Calculation of min and max energy of electrons
! 'ee_min' and 'ee_max',
! which can be produced by proton with energy 'ep', and
! soft photon with energy 'eph'
! All energies in eV
!=======================================================
! Variables
 real(8), intent(in) :: eph, ep
 real(8), intent(out) :: ee_min, ee_max
 real(8) :: eps, gp, x, tp, res
! Calculations
 eps=eph/me
 gp=ep/mp
 x=gp*eps
 res=me*gp/(1+4*x*mep)
 tp=(sqrt(x)+sqrt(x-1))**2
 ee_min=res/tp
 ee_max=res*tp
end subroutine el_min_max_en


subroutine prot_min_en(eph_max,ee,ep_min)
!=========================================================
! Calculation of minimal energy of proton 'ep_min', which 
! can produce electrons with energy 'ee' in the soft
! photon field with maximum energy 'eph_max'
! All energies are in eV
!=========================================================
! Variables
 real(8), intent(in) :: eph_max, ee
 real(8), intent(out) :: ep_min
 real(8) :: xmin, xmax
 integer :: err_code
! Calculations
 
 eps=eph_max/me
 elen=ee/me
 xmin=ee/mp
 xmax=1d23/mp
 if (xmin*eps<1) xmin=(1+1d-14)/eps
 if (xmin>=xmax) xmax=xmin*1d5
 call bs_root(xmin,xmax,elen_max,ep_min,err_code)
 ep_min=ep_min*mp
 
end subroutine prot_min_en


subroutine elen_max(gp,res)
!=========================================================
! For subroutine 'prot_min_en'
!=========================================================
! Variables
 real(8), intent(in) :: gp
 real(8), intent(out) :: res
 real(8) :: x
! Calculations
 x=gp*eps
 res=gp/(1+4*x*mep)
 res=res*(sqrt(x)+sqrt(x-1))**2
 res=res-elen
end subroutine elen_max



end module pgam_ep_pair



!program main
! call tab_rate
!end program main
!
!
!
!subroutine tab_rate
! use pgam_ep_pair, only : beheit_tab
! type(beheit_tab) :: fc
! 
! interface
! subroutine sph_dist(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! end subroutine sph_dist
! end interface
! 
! 
! 
! fc%info%typ='sub'
! fc%info%intro='Parameters of tabulation'
! fc%info%sub_name='tab_rate'
! fc%info%file_name='pgam_ep_pair.f90'
! call fc%add_info
! 
! fc%info%typ='sub'
! fc%info%intro='Parameters of soft photon field'
! fc%info%sub_name='sph_dist'
! fc%info%file_name='pgam_ep_pair.f90'
! call fc%add_info
! 
! call fc%phot(1d-10,150d0,sph_dist)
! call fc%elec(1d12,1d20,100)
! call fc%prot(1d12,1d21,100)
! fc%file_out='bethe_heitler_sph.dat'
! call fc%calc
!
!
!end subroutine tab_rate
!
!
!subroutine sph_dist(x,res)
! use thermal_photon_field, only : therm_phf, set_thnden, thnden
!! Variables 
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! type(therm_phf) :: sph1
!! Calculations
! call sph1%init(4)
! sph1%temp=[3d0,3d-1,6d-3,2.35d-4]
! sph1%uden=[5d4,5d4,5d4,0.26d0]
! call set_thnden(sph1)
! call thnden(x,res)
!end subroutine sph_dist




module beheit_module
!====================================================================
! Calculation of bethe-heitler production rate for specific
! proton spectrum
!====================================================================
 use gauss_kronrod, only : gk_adaptive
 use house_keeping, only : get_unit
 use intpol_mod, only : arr_ind_short, logl_int
 use tab2fun, only : tabfun
 implicit none
 private
 
 
 public :: beheit
 
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 
 type beheit 
  integer :: npl=0, ne=0, np=0
  real(8) :: ep1, ep2
  real(8), allocatable :: ee(:), ep(:,:), rate(:,:)
  real(8), allocatable :: epl(:), loss(:), int_len(:)
  procedure(fun_1d), pointer, nopass :: protsp=>null()
  character(500) :: tab_file=''
 contains
  procedure :: read_tab, set_protsp, bh_rate
 end type beheit


 procedure (fun_1d), pointer :: pr_spec=>null(), phi_fun=>null()

 contains


subroutine read_tab(this)
!=============================================
! Reading of production rate for single
! proton
!=============================================
 class(beheit) :: this
 integer :: ndev, ndes
 character(len=:), allocatable :: description
 
 
 if (allocated(this%ee)) deallocate(this%ee)
 if (allocated(this%ep)) deallocate(this%ep)
 if (allocated(this%rate)) deallocate(this%rate)
 if (allocated(this%epl)) deallocate(this%epl)
 if (allocated(this%loss)) deallocate(this%loss)
 if (allocated(this%int_len)) deallocate(this%int_len)
 this%npl=0
 this%ne=0
 this%np=0
  
 call get_unit(ndev)
 open(ndev,file=this%tab_file,form='unformatted')
 read(ndev) ndes
 allocate(character(len=ndes) :: description)
 read(ndev) description
 read(ndev) this%npl
 read(ndev) this%ne
 read(ndev) this%np 
 allocate(this%ee(this%ne))
 allocate(this%ep(this%np,this%ne))
 allocate(this%rate(this%np,this%ne))
 read(ndev) this%ee
 read(ndev) this%ep
 read(ndev) this%rate
 
 if (this%npl>0) then
  allocate(this%epl(this%npl))
  allocate(this%loss(this%npl))
  allocate(this%int_len(this%npl))
  read(ndev) this%epl
  read(ndev) this%loss
  read(ndev) this%int_len
 end if 
 close(ndev)
end subroutine read_tab

subroutine set_protsp(this,ep1,ep2,protsp)
!===================================================
! Setting of proton spectrum
!===================================================
 class(beheit) :: this
 real(8), intent(in) :: ep1, ep2  ! in eV
 procedure(fun_1d) :: protsp      ! in 1/eV

 this%ep1=ep1
 this%ep2=ep2
 this%protsp=>protsp
end subroutine set_protsp



subroutine bh_rate(this,x,res)
!===================================================
! Calculation of Bethe-Heitler production rate
! for given proton spectrum
!===================================================
! Variables
 class(beheit), intent(in) :: this
 real(8), intent(in)  :: x   ! energy of electron in eV
 real(8), intent(out) :: res ! production rate of electrons in 1/(eV sec)
 integer :: n1, n2
 real(8) :: res1, res2
 type(tabfun) :: phi1, phi2
! Calculations
 if ((x<this%ee(1)).or.(x>this%ee(this%ne))) then
  res=0d0
  return
 end if
 call arr_ind_short(1,this%ne,this%ee,x,n1,n2) 
 call phi1%on(this%ep(:,n1),this%rate(:,n1))
 call conv_fp(this,phi1,res1)
 call phi2%on(this%ep(:,n2),this%rate(:,n2))
 call conv_fp(this,phi2,res2)
 call logl_int(this%ee(n1),this%ee(n2),res1,res2,x,res)
 
 
 call phi1%off
 call phi2%off
end subroutine bh_rate


subroutine conv_fp(this,phic,res)
! Variables
 class(beheit), intent(in) :: this
 type(tabfun), intent(in) :: phic
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations

 if (phic%zero==0) then
  res=0d0
  return
 end if
 
 tmin=phic%xmin
 tmax=phic%xmax
 if (tmin<this%ep1) tmin=this%ep1
 if (tmax>this%ep2) tmax=this%ep2
 
 if (tmin>tmax) then
  res=0d0
  return
 end if
 
 phi_fun=>phic%fx
 pr_spec=>this%protsp
 tmin=log(tmin)
 tmax=log(tmax)
 call gk_adaptive(iconv_fp,tmin,tmax,res)
 phi_fun=>null()
 pr_spec=>null()

end subroutine conv_fp


subroutine iconv_fp(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: ep, fp, phi
! Calculations
 ep=exp(x)
 call pr_spec(ep,fp)
 call phi_fun(ep,phi)
 res=fp*phi*ep
end subroutine iconv_fp


end module beheit_module


!program main
! use beheit_module, only : beheit
! type(beheit) :: bc
! real(8) :: e1, e2, ee, res
! integer :: i, ne
! 
! interface
!  subroutine psp(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine psp 
! end interface
! 
! 
! bc%tab_file='bethe_heitler_sph.dat'
! call bc%read_tab
! call bc%set_protsp(1d9,1d19,psp)
! 
! e1=1d11
! e2=1d18
! ne=100
! 
! 
! open(1,file='bh_test.dat')
! do i=0,ne
!  ee=e1*(e2/e1)**(1d0*i/ne)
!  call bc%bh_rate(ee,res)
!  write(1,*) ee, res*ee**2
! end do
! close(1)
!
!end program main
!
!
!subroutine psp(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! 
! res=1/x**2
!
!end subroutine psp


