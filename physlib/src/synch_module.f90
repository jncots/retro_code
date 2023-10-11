module synch_module
!=====================================================================
! Example of usage is in the end of this file.
!
! Contains object of synch_class for calculation of synchrotron
! radiation:
!
! use synch_module, only : synch_class
! type(synch_class) :: sc
! call sc%set_elsp(e1,e2,elsp) ! settings for electron spectrum,
! e1, e1 in eV, elsp in 1/eV (electron-volt)
! call sc%set_bm(bm,turb) ! setting of magnetic field, bm in Gauss
! turb(optional) is 0 or 1, determines sc%turb
! sc%turb=0 ! calculations for homogeneous field, sc%turb=1 by default 
!
! call sc%calc_en(eg) ! eg is optional
! sc%ph1 and sc%ph2 are minimum and maximum energy of synchrotron
! radiation for given spectrum and magnetic field (always)
! if eg is present: calculations of sc%el_min and sc%el_peak
! are minimum energy of electron which produce photon with energy eg
! and energy of electron which produce synchrotron with peak at eg
! 
! call sc%syn_rad(eg,res) ! eg is energy of photon in eV,
! res is production rate of synchrotron radiation in 1/(eV s)
! call sc%syn_loss(ee,res) ! ee is energy of electron in eV,
! res is energy losses dE/dt in eV/s
!=====================================================================
 use gauss_kronrod, only : gk_adaptive
 implicit none
 private
 
 
 public :: synch_class
 
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun_1d
 end interface
 

 type synch_class
  real(8) :: bm
  real(8), private :: norm, iwc, asyn, bloss
  real(8) :: e1, e2
  real(8) :: ph1, ph2
  real(8) :: el_peak, el_min
  integer :: turb=1 ! turb=1 for turbulent (chaotic) magnetic field, turb=0 for homogeneous
  procedure(fun_1d), pointer, nopass :: elsp=>null() ! electron spectrum
  procedure(fun_1d), pointer, nopass :: em_fun=>null() ! emissivity function
 contains
  procedure :: set_elsp, set_bm, calc_en, syn_rad, syn_loss
 end type synch_class
 
 type(synch_class), pointer :: loc_syn
 

 real(8), parameter :: syn_norm=3.538078273d4 ! 1/(s)
! syn_norm=sqrt(3)/(2*Pi)*(ez^2)/(hbar*cs)*ez*Bm/(me*cs), for Bm=1 G
! real(8), parameter :: e_cycl=1.736514534d-8 ! eV
! e_cycl=3/2*e/(me*c)*hbar/eV, cyclotron energy for B=1 G and gamma=1
 real(8), parameter :: e_cev=6.650257238d-20
! e_cev=3/2*e/(me*c)*hbar/eV*1/(me_eV)**2, cyclotron energy for B=1 G and E=1 eV 
 real(8), parameter :: xmax=2d-2 ! 1/x, where max x=50 for G(x), min energy of electrons
 real(8), parameter :: aloss=2.529014647511923d-15
! aloss=4/9*e^4/(me^2*c^3)/eV/(meV)^2 
 
 contains


subroutine set_elsp(this,e1,e2,elsp)
 class(synch_class) :: this
 real(8), intent(in) :: e1, e2 ! min and max energy of electron spectrum eV
 procedure(fun_1d) :: elsp     ! electron spectrum in 1/eV
 this%e1=e1
 this%e2=e2
 this%elsp=>elsp
end subroutine set_elsp


subroutine set_bm(this,bm,turb)
 class(synch_class) :: this
 integer, intent(in), optional :: turb
 real(8), intent(in) :: bm
 real(8) :: aa
 
 this%bm=bm    ! magnetic field in Gauss
 this%norm=syn_norm*this%bm
 this%iwc=1d0/(e_cev*this%bm)
 this%bloss=aloss*(this%bm)**2
 
 if (present(turb)) this%turb=turb 
 
 if (this%turb==0) this%bloss=1.5d0*this%bloss
 if (this%turb==0) then
  this%em_fun=>f_synch
 else
  this%em_fun=>g_synch
 end if
end subroutine set_bm


subroutine syn_loss(this,ee,res)
 class(synch_class) :: this
 real(8), intent(in) :: ee
 real(8), intent(out) :: res
 res=this%bloss*ee**2 ! res in eV/sec
end subroutine syn_loss



subroutine calc_en(this,egam)
!=======================================================
! Calculation of ph1 and ph2 which are minimum and
! maximum energy of synchrotron radiation for given
! spectrum and magnetic field (always)
!
! Calculates energies el_peak, el_min which gives egam
! photon energy for given magnetic field bm (optional)
!========================================================
 class(synch_class) :: this
 real(8), intent(in), optional :: egam
 real(8), parameter :: xpeak=0.2292d0, xmax=50d0
 real(8) :: aa
 
 aa=1/this%iwc
 this%ph1=1d-4*aa*this%e1**2 ! min and max energy of radiation
 this%ph2=50*aa*this%e2**2
 
 if (present(egam)) then
  aa=egam*this%iwc
  this%el_min=sqrt(aa/xmax)
  this%el_peak=sqrt(aa/xpeak)
 end if
end subroutine calc_en



subroutine syn_rad(this,egam,res)
 class(synch_class), target :: this
 real(8), intent(in) :: egam
 real(8), intent(out) :: res
 real(8) :: tmin, tmax

 this%asyn=egam*this%iwc
 tmin=sqrt(xmax*this%asyn)
 if (tmin<this%e1) tmin=this%e1
 tmin=log(tmin)
 tmax=log(this%e2)
 if (tmin>=tmax) then
  res=0d0
  return
 end if
 
 loc_syn=>this 
 call gk_adaptive(isynch,tmin,tmax,res,tol=1d-4)
 res=this%norm*res/egam ! 1/(s eV)
end subroutine syn_rad


subroutine isynch(x,res)
!
! Integrand for convolutions of electron spectrum with
! synchrotron generation function
!
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: ee, xsyn, gsyn, fe
 
 ee=exp(x)
 xsyn=loc_syn%asyn/ee**2
 call loc_syn%em_fun(xsyn,gsyn)
 call loc_syn%elsp(ee,fe)
 res=gsyn*fe*ee
end subroutine isynch



subroutine g_synch(x,res)
!
! Emmisivity function in turbulent magnetic field
! (the directions of magnetic field is averaged over
! all directions)
!
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x13, x23, x43, p1, p2, p3

 x13=x**(1d0/3d0)
 x23=x13**2
 x43=x23**2
 
 p1=1d0+2.21d0*x23+0.347d0*x43
 p2=1d0+1.353d0*x23+0.217d0*x43
 p3=1.808d0*x13/sqrt(1d0+3.4d0*x23)
 
 res=p3*(p1/p2)*exp(-x)
end subroutine g_synch


subroutine f_synch(x,res)
!
! Emmisivity function in homogeneous magnetic field
!
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x13, x23, x43, p1, p2, p3
 
 x13=x**(1d0/3d0)
 x23=x13**2
 x43=x23**2
 
 p1=1d0+0.884d0*x23+0.471d0*x43
 p2=1d0+1.64d0*x23+0.974d0*x43
 p3=2.15d0*x13*(1d0+3.06*x)**(1d0/6d0)
 
 res=p3*(p1/p2)*exp(-x)
end subroutine f_synch

end module synch_module


!program test
! use synch_module, only : synch_class
! use synch_rad, only : set_synch_bm, set_synch_ee, synch
! use phys_const, only : year
! implicit none
! type(synch_class) :: sc
! real(8) :: eg1, eg2, eg, res, spnorm, rel
! integer :: i, ng
! 
! 
! interface
!  subroutine elsp(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine elsp
! end interface
! 
!! sc%turb=0
! call sc%set_elsp(1d7,1d17,elsp)
! call sc%set_bm(1d-4)
! call sc%calc_en(1d7)
! write(*,'(10Es14.6)') sc%ph1,  sc%ph2
! write(*,'(10Es14.6)') sc%el_min,  sc%el_peak
! 
! call set_synch_bm(1d-4,elsp,spnorm)
! call set_synch_ee(1d7,1d17)
! 
! 
! eg1=1d5
! eg2=1d15
! ng=100
! 
! open(1,file='test_syn.dat')
! do i=0,ng
!  eg=eg1*(eg2/eg1)**(1d0*i/ng)
!  call sc%syn_rad(eg,res)
!  call synch(eg,rel)
!  rel=rel*spnorm
!  write(1,*) eg, res*eg**2, rel*eg**2, (res-rel)/rel
! end do
! close(1)
! 
!  
! eg1=1d6
! eg2=1d18
! ng=100
! 
! open(1,file='loss_syn.dat')
! do i=0,ng
!  eg=eg1*(eg2/eg1)**(1d0*i/ng)
!  call sc%syn_loss(eg,res)
!  write(1,*) eg, 1/(year*res/eg)
! end do
! close(1)
!
!
!
!end program test
!
!
!subroutine elsp(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! res=1/x**2
!end subroutine elsp

