module synch_cooling
!============================================================
! Example of usage see below
!============================================================ 
 use gauss_kronrod, only : gk_adaptive
 implicit none
 private
 save
 
 public :: set_synch_cool, synch_cool, syn_cool_data
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type syn_cool_data
  logical :: time_depended ! time (.true.) or stationary (.false.) type of problem
  logical :: impulsive ! impulsive (.true.) or continious (.false.) type of problem
  real(8) :: bm ! magnetic field in Gauss, [G]
  real(8) :: asyn ! precalculated coefficient for b(Ee)=asyn*Ee**2
  real(8), private :: itl ! precalculated coefficient for asyn*time
  real(8) :: time ! time of cooling in seconds, [sec]
  real(8) :: emin, emax ! minimum and maximum energy of electron spectrum of the source, [eV]
  procedure (fun_1d), pointer, nopass :: q_sp=>null()! function of electron spectrum in [1/eV]
 contains
  procedure :: init_spec, init_bm, init_time, clean=>clean_data
 end type syn_cool_data

 
 class(syn_cool_data), pointer :: fcalc
 
 contains

!============================================================
! Subroutines for type(syn_cool_data)
!============================================================ 

subroutine clean_data(fdat)
! Variables
 class(syn_cool_data) :: fdat
! Calculations
 fdat%time_depended=.false.
 fdat%impulsive=.false.
 fdat%bm=0d0
 fdat%asyn=0d0
 fdat%itl=0d0
 fdat%time=0d0
 fdat%emin=0d0
 fdat%emax=0d0
 fdat%q_sp=>null()
end subroutine clean_data

 
subroutine init_spec(fdat,qs)
! Variables
 class(syn_cool_data) :: fdat
 interface
  subroutine qs(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
! Calculations
 fdat%q_sp=>qs
end subroutine init_spec

subroutine init_bm(fdat,bm)
!
! Initiation of magnetic field
!
! Variables
 class(syn_cool_data) :: fdat
 real(8), intent(in) :: bm
 real(8), parameter :: asyn=2.5290146475119223d-15
! asyn=4/9*(e^4/m^2c^3)*1/eV*(eV/mc^2)^2
! for turbulent magnetic field !!! because multiplied by 2/3
! asyn*Bm^2(in Gauss)*E^2(in eV)=I[eV/s] 
! Calculations
 fdat%bm=bm
 fdat%asyn=asyn*(fdat%bm)**2
 fdat%itl=fdat%asyn*fdat%time
end subroutine init_bm

subroutine init_time(fdat,time)
!
! Initiation of cooling time
!
! Variables
 class(syn_cool_data) :: fdat
 real(8), intent(in) :: time
! Calculations
 fdat%time=time
 fdat%itl=fdat%asyn*fdat%time
end subroutine init_time

!============================================================
! Subroutines for synchrotron cooling
!============================================================ 

subroutine set_synch_cool(fdat)
! Variables 
 class(syn_cool_data), target, intent(in) :: fdat
! Calculations
 fcalc=>fdat
end subroutine set_synch_cool


subroutine synch_cool(Ee,res)
! Variables
 real(8), intent(in) :: Ee
 real(8), intent(out) :: res
! Calculations
 if (fcalc%impulsive) then
  call synch_cool_imp(Ee,res)
 else
  call synch_cool_cont(Ee,res)
 end if
end subroutine synch_cool

subroutine synch_cool_imp(Ee,res)
! Cooled spectrum of impulse source
! Variables
 real(8), intent(in) :: Ee
 real(8), intent(out) :: res
 real(8) :: itl, E_eff
! Calculations


! If energy is too large 
 if (Ee>=fcalc%emax) then
  res=0d0
  return
 end if 
 
 itl=fcalc%itl*Ee   ! inverse of time loss
 if (itl>=1d0) then
  res=0d0
  return
 else
  E_eff=Ee/(1d0-itl)
  
  if (E_eff>=fcalc%emax) then
   res=0d0
   return
  end if
  
  if (E_eff<=fcalc%emin) then
   res=0d0
   return
  end if 
 
 end if
 
 call fcalc%q_sp(E_eff,res)
 res=res*(E_eff/Ee)**2
end subroutine synch_cool_imp


subroutine synch_cool_cont(Ee,res)
! Cooled spectrum of electrons
 real(8), intent(in) :: Ee
 real(8), intent(out) :: res
 real(8) :: E0, E_eff, tmin, tmax, itl
! Calculations
 
! If energy is too large 
 if (Ee>=fcalc%emax) then
  res=0d0
  return
 end if

! Lower limit of integration
 if (Ee>fcalc%emin) then
  E0=Ee
 else
  E0=fcalc%emin
 end if 
! Upper limit
 if (fcalc%time_depended) then
  itl=fcalc%itl*Ee   ! inverse of time loss
  if (itl>=1d0) then
   E_eff=fcalc%emax
  else
   E_eff=Ee/(1d0-itl)
   if (E_eff>fcalc%emax) E_eff=fcalc%emax
   if (E_eff<=fcalc%emin) then
    res=0d0
    return
   end if 
  end if
 else
  E_eff=fcalc%emax ! For stationary problem t>>tcool
 end if
 
 
 tmin=log(E0)
 tmax=log(E_eff)
 !write(*,*) 'synch_cool_cont:'
 !write(*,'(10Es14.6)') fcalc%emin, fcalc%emax, Ee, E0, E_eff
 !write(*,'(10Es14.6)') exp(tmin), exp(tmax), exp(tmin)-exp(tmax),&
 !(exp(tmin)-E0)/E0, (exp(tmax)-E_eff)/E_eff
 !read(*,*)
 call gk_adaptive(icool,tmin,tmax,res,tol=1d-4)
 res=res/(fcalc%asyn*Ee**2)
end subroutine synch_cool_cont


subroutine icool(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: Ee
! Calculations
 Ee=exp(x)
 call fcalc%q_sp(Ee,res)
 res=res*Ee
end subroutine icool


end module synch_cooling


!============================================================
! Example of usage (test program):
!============================================================ 
!module test_synch 
! use synch_cooling, only : set_synch_cool, synch_cool, syn_cool_data
! 
! use phys_const, only : year
! 
! implicit none
! private
! save
! 
! public :: main_test
! 
! 
! contains
!
!subroutine main_test
! type(syn_cool_data) :: calc, calc1
! integer :: ne, i
! real(8) :: Ee1, Ee2, Bm, Ee, res
! real(8) :: emin, emax, tc
!
!
! Ee1=1d9
! Ee2=1d17
! ne=1000
! 
! tc=1d0*year
! Bm=1d-4
! emin=1d14
! emax=1d17
! 
! 
! calc%emin=emin
! calc%emax=emax
! calc%time_depended=.true.
! calc%impulsive=.true.
! call calc%init_spec(f1spec)
! call calc%init_bm(Bm)
! call calc%init_time(tc)
! 
! call set_synch_cool(calc)
! open(1,file='test1.dat')
! do i=0,ne
!  Ee=Ee1*(Ee2/Ee1)**(1d0*i/ne) 
!  call synch_cool(Ee,res)
!  write(1,*) Ee, res*Ee**2
! end do
! close(1)
! 
! 
! calc1%emin=emin
! calc1%emax=emax
! calc1%time_depended=.true.
! calc1%impulsive=.false.
! call calc1%init_spec(f2spec)
! call calc1%init_bm(Bm)
! call calc1%init_time(tc)
! 
! call set_synch_cool(calc1)
! open(1,file='test_cont.dat')
! do i=0,ne
!  Ee=Ee1*(Ee2/Ee1)**(1d0*i/ne) 
!  call synch_cool(Ee,res)
!  write(1,*) Ee, res*Ee**2 
! end do
! close(1) 
! 
! 
! 
!
!
!end subroutine main_test
!
!
!subroutine f1spec(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! res=(1d10/x)**1d0
!end subroutine f1spec
!
!
!subroutine f2spec(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! res=(1d10/x)**2d0
!end subroutine f2spec
!
!
!end module test_synch
!
!program main
! use test_synch, only : main_test
! 
! call main_test
!
!end program main


module synch_rad
 use gauss_kronrod, only : gk_adaptive
 implicit none
 private
 save
 
 public :: set_synch_bm, set_synch_ee, synch, synch_info, synch_el_en
 
! real(8), parameter :: e_cycl=1.736514534d-8 ! eV
! e_cycl=3*e*hbar/(2*me*c*eV), cyclotron energy for B=1 G and gamma=1
 real(8), parameter :: syn_norm=3.538078273d4 ! 1/(s)
! syn_norm=sqrt(3)/(2*Pi)*(ez^2)/(hbar*cs)*ez*Bm/(me*cs), for Bm=1 G
! real(8), parameter :: e_cev=e_cycl/(me_ev)**2
 real(8), parameter :: e_cev=6.650257238d-20
 real(8), parameter :: xmax=2d-2 ! 1/x, where max x=50 for G(x), min energy of electrons
 
 real(8) :: eemin, eemax, iwc, asyn
 

 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 procedure (fun_prot), pointer :: el_sp=>null()
 
 contains
 
subroutine synch_info(bb,x1,x2,emin,emax)
!
! Calculates energy range for given magnetic field bb
! bb - magnetic field in Gauss
! x1, x2 - min and max energy of electrons in eV
! emin, emax - min and max energy of photons in eV
!
! Variables
 real(8), intent(in) :: bb, x1, x2
 real(8), intent(out) :: emin, emax ! eV, energy range of radiation
 real(8) :: aa
! Calculations
 aa=e_cev*bb
 emin=1d-4*aa*x1**2
 emax=50*aa*x2**2
end subroutine synch_info

subroutine synch_el_en(bb,xg,emin,epeak)
!
! Calculates energies epeak, emin which gives xg photon energy for given bb
!
! Variables
 real(8), intent(in) :: bb, xg
 real(8), intent(out) :: epeak, emin
 real(8), parameter :: xpeak=0.2292d0, xmax=50d0
 real(8) :: aa
! Calculations
 aa=e_cev*bb
 emin=sqrt(xg/(xmax*aa))
 epeak=sqrt(xg/(xpeak*aa))
 
end subroutine synch_el_en


 
subroutine set_synch_bm(Bm,elsp,norm)
! Variables
 real(8), intent(in) :: Bm ! magnetic field strength in Gauss
 real(8), intent(out) :: norm
 
 interface
  subroutine elsp(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
! Calculations

 el_sp=>elsp
 norm=syn_norm*Bm 
 iwc=1d0/(e_cev*Bm)
end subroutine set_synch_bm



subroutine set_synch_ee(x1,x2)
! Variables
 real(8), intent(in) :: x1, x2 ! min and max energy of electrons in spectrum elsp
! Calculations
 eemin=x1
 eemax=x2
end subroutine set_synch_ee







subroutine synch(Egam,res)
! Variables
 real(8), intent(in) :: Egam
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations
 asyn=Egam*iwc
 tmin=sqrt(xmax*asyn)
 if (tmin<eemin) tmin=eemin
 tmin=log(tmin)
 tmax=log(eemax)
 if (tmin>=tmax) then
  res=0d0
  return
 end if
  
 call gk_adaptive(isynch,tmin,tmax,res,tol=1d-4)
 res=res/Egam
end subroutine synch
 
 
subroutine isynch(x,res)
!
! Integrand for convolutions of electron spectrum with
! synchrotron generation function
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: Ee, xsyn, gsyn, fe
! Calculations
 
 Ee=exp(x)
 xsyn=asyn/Ee**2
 call g_synch(xsyn,gsyn)
 call el_sp(Ee,fe)
 res=gsyn*fe*Ee
 
end subroutine isynch


subroutine g_synch(x,res)
!
! Emmisivity function in turbulent magnetic field
! (the directions of magnetic field is averaged over
! all directions)
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x13, x23, x43, p1, p2, p3
! Calculations
 x13=x**(1d0/3d0)
 x23=x13**2
 x43=x23**2
 
 p1=1d0+2.21d0*x23+0.347d0*x43
 p2=1d0+1.353d0*x23+0.217d0*x43
 p3=1.808d0*x13/sqrt(1d0+3.4d0*x23)
 
 res=p3*(p1/p2)*exp(-x)

end subroutine g_synch


end module synch_rad




!program test
! use synch_rad, only : synch_info, synch_el_en
! implicit none
! real(8) :: Bm, Ee1, Ee2, Eg1, Eg2
! 
! Bm=1d-4
! 
! Ee1=1d9
! Ee2=1d10
! call synch_info(Bm,Ee1,Ee2,Eg1,Eg2)
! 
! write(*,'(10Es14.6)') Eg1, Eg2
!
!
!
!end program test