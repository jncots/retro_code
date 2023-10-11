module tab_ggee
!===========================================================
! Module integrates the produciton rate of the process
! gamma+soft photon->electron+positron
! with soft photon field and output tabulated result
! as function of electron energy and gamma ray energy.
! To calculate the production rate, one need one further
! integration with gamma-ray spectrum.
!
! It also contains the integration the total probability of
! gamma-gamma absorbtion with density of soft photon field
! and output the result as tabulated function of
! absorbtion probabililty of gamma ray with energy egam
!
! Usage:
! =====
! ng=1000 ! size of grid for total interaction probability
! ne=600  ! size of grid of pair production rate
! (ne of electron energies over ne of gamma energies)
! egam1=1d8   ! minimum energy of gamma-ray spectrum in eV
! egam2=1d18  ! maximum energy of gamma-ray spectrum in eV
! eph1=1d-10  ! minimum energy of soft photon field spectrum in eV
! eph2=1d3    ! maximum energy of soft photon field spectrum in eV
! nden(x,res) ! soft photon field spectrum in [1/eV]
! fname='ggee.dat' ! path and name of the file for saving of results
! inf is the input information about calculation which should be
! written in the form of subroutines in the file
! call ree_gam(ng,ne,egam1,egam2,eph1,eph2,nden,fname)
!
! Example of usage (test program) is in the end of the file
!
!===========================================================
 use gauss_kronrod, only : gk_adaptive
 use description_line, only : descr_line, infoline
 use word_processor, only : line_break
 use phys_const, only : c_light
 implicit none
 private
 save
 

 public :: ree_gam
 
 
 real(8), parameter :: app=2.9915347994420107d-14
! app=2*2*Pi*re^2*c [cm^3/s], where re=e^2/(me*c^2) is the classical radius of electron
! first 2 takes into account both electron and positrion contribution
 real(8), parameter :: me=0.510998928d6 ! electron rest mass in eV
 real(8), parameter :: me2=me**2! 2.611199044d11
 real(8), parameter :: ime2=1/me2 
 real(8), parameter :: re2c=2.380587754d-15 ! re2c=re^2*c, re is classical radius of electron
 
 real(8) :: eel ! electron energy
 real(8) :: egam ! gamma ray energy
 real(8) :: u1, z1 ! part of u and z variables which do not depend on eph
 real(8) :: eph1, eph2 ! min and max energy of soft photons in eV 
 real(8) :: egam1, egam2 ! min and max energy of gamma rays in eV
 
 integer :: ndes
 character(len=:), allocatable :: description
 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 procedure (fun_prot), pointer :: nden=>null()
 
  contains
  
  
 
subroutine ree_gam(neg,ne,eg1,eg2,es1,es2,sph_den,fname,inf)
! Variables
 real(8), intent(in) :: eg1,eg2,es1,es2
 procedure (fun_prot) :: sph_den
 character(100) :: fname
 integer, intent(in) :: neg, ne ! array size of absorbtion and rate
 integer :: ng, ie, ig, i
 real(8) :: eg22, zeta, zeta1, eel1, eel2
 real(8) :: egam_min, egam_max, egam, eg, res
 real(8), allocatable :: fres(:,:), ee_arg(:), gam_arg(:,:)
 real(8), allocatable :: xgam(:), wgam(:)
 type(infoline) :: inf(:)
! Calculations 
 
 call info(inf)
 
 nden=>sph_den ! soft photon field energy distribution
 eph1=es1 ! soft photon field boundaries
 eph2=es2
 
 egam1=eg1 ! gamma-ray boundaries
 egam2=eg2


! Calculation min and max energy of produced electrons for a given
! gamma-ray distribution and photon field. min and max is determined
! by egam2 and eph2, which are maximum energies of both distributions
 eg22=egam2/2
 zeta=me**2/(egam2*eph2)
 zeta1=1+sqrt(1-zeta)
 eel2=eg22*zeta1
 eel1=eg22*zeta/zeta1
 
 ng=ne ! ne=600
 
 allocate(fres(ng+1,ne+1))
 allocate(ee_arg(ne+1))
 allocate(gam_arg(ng+1,ne+1))
 
 do ie=0,ne
  eel=eel1*(eel2/eel1)**(1d0*ie/ne)
  
  ee_arg(ie+1)=eel !
  
  egam_min=eel
  if (egam_min<egam1) egam_min=egam1
  egam_max=egam2
  
  do ig=0,ng
   egam=egam_min*(egam_max/egam_min)**(1d0*ig/ng)
   call soft_ph(egam,res)
   gam_arg(ig+1,ie+1)=egam !
   fres(ig+1,ie+1)=res   !
  end do
 end do
 
 
 allocate(xgam(neg+1))
 allocate(wgam(neg+1))
 do i=0,neg
  eg=egam1*(egam2/egam1)**(1d0*i/neg)
  call wgg(eg,res)
  xgam(i+1)=eg
  wgam(i+1)=res/c_light ! 1/cm
 end do
 
 
 open(1,file=fname,form='unformatted')
 write(1) ndes
 write(1) description
 write(1) size(ee_arg)
 write(1) ee_arg
 write(1) size(gam_arg,1), size(gam_arg,2)
 write(1) gam_arg
 write(1) fres
 write(1) size(xgam)
 write(1) xgam
 write(1) wgam
 close(1)
 
 deallocate(fres)
 deallocate(ee_arg)
 deallocate(gam_arg)
 deallocate(xgam)
 deallocate(wgam)
 
 
end subroutine ree_gam



subroutine info(inf)
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
 line='This file contains the rate of pair production in gamma-gamma absorbtion &
 for specific soft photon field (see description below), the range of gamma-ray, &
 and electron energies, as well as the tabulated function of the absorbtion probability &
 of gamma ray with energy egam for the same soft photon field'
 call line_break(nleng,line)
 
 ninf=size(inf)
 nnote=ninf+1
 
 allocate(note(nnote))
 
 note(1)%typ='des'
 note(1)%intro='Introduction'
 note(1)%info=line
 
 do i=2,nnote
  note(i)=inf(i-1)
 
 !note(2)%typ='sub'
 !note(2)%intro='Parameters of tabulation'
 !note(2)%sub_name='tab_param'
 !note(2)%file_name='tab_ggee.f90'
 !
 !
 !note(3)%typ='sub'
 !note(3)%intro='Setup of the calculations'
 !note(3)%sub_name='ree_gam'
 !note(3)%file_name='tab_ggee.f90'
 !
 !note(4)%typ='sub'
 !note(4)%intro='Soft photon field'
 !note(4)%sub_name='nden'
 !note(4)%file_name='tab_ggee.f90'
 end do
 
 call descr_line(nnote,note,ndes,description)


end subroutine info





subroutine soft_ph(egam,res)
! Variables
 real(8), intent(in) :: egam
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
 real(8) :: elg, eph_min, eph_max
! Calculations
 z1=egam*ime2
 elg=eel/egam
 if (elg>=1d0) then
  res=0d0
  return
 end if 
 u1=4*elg*(1-elg)
 
 eph_min=1/(z1*u1)
 if (eph_min<eph1) eph_min=eph1
 eph_max=eph2
 if (eph_min>=eph_max) then
  res=0d0
  return
 end if
 tmin=log(eph_min)
 tmax=log(eph_max)
 call gk_adaptive(isoft_ph,tmin,tmax,res)
 res=app*res/egam
 
end subroutine soft_ph



subroutine isoft_ph(x,res)
!==========================================================
! Integrand of the integration of pair production rate with
! soft photon field 'nden(eph,res)', where eph is energy of
! soft photons
!==========================================================
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eph, z, u, phi, fph
! Calculations
 eph=exp(x)
 z=eph*z1
 u=z*u1
 phi=2*u*log(u)+(1-u)*(2+u-2*z)
 phi=phi/(z*u**2)
 call nden(eph,fph)
 res=fph*phi*eph
end subroutine isoft_ph


!============================================================
! Integration the total probability of gamma-gamma
! absorbtion with density of soft photon field
! and output the result as tabulated function of
! absorbtion probabililty of gamma ray with energy egam
!============================================================

subroutine wgg(eg,res)
! Variables
 real(8), intent(in) :: eg
 real(8), intent(out) :: res
 real(8) :: tmin, tmax, emin, emax, es
! Calculations
 egam=eg
 
 es=me2/egam
 emin=eph1
 if (es>eph1) emin=es
 emax=eph2
 
 if (emin>=emax) then
  res=0d0
  return
 end if 
 
 tmin=log(emin)
 tmax=log(emax)
 call gk_adaptive(iwgg,tmin,tmax,res)
end subroutine wgg


subroutine iwgg(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eph, ww, fph 
! Calculations
 eph=exp(x)
 call w_abs(egam,eph,ww)
 call nden(eph,fph)
 res=fph*ww*eph
end subroutine iwgg

 
 
subroutine w_abs(egam,eph,res)
! Variables
 real(8), intent(in) :: egam, eph ! in eV
 real(8), intent(out) :: res ! in cm^3/s
 real(8) :: z, lz
! Calculations
 z=egam*eph*ime2
 if (z<=1d0) then
  res=0d0
  return
 end if
 lz=log(z)
 res=(z+0.41d0*lz)*log(0.541d0*z+1.406d0)
 res=res/(z-0.47d0*lz)
 res=res*((z-1)/z)**(3d0/2d0)
 res=res*6.28d0*re2c/z
end subroutine w_abs



end module tab_ggee


!============================================================
! Example of usage (test program)
!============================================================
!program main
!
! call tab_param
!
!end program main
!
!
!subroutine tab_param
! use tab_ggee, only : ree_gam
! implicit none
!! Variables
! integer :: ng, ne
! real(8) :: egam1, egam2, eph1, eph2
! character(100) :: fname
! 
! interface
!  subroutine nden(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine
! end interface
! 
! 
!! Calculations 
! egam1=1d8
! egam2=1d18
! eph1=1d-10
! eph2=1d3
! fname='ggee.dat'
! 
! ng=1000
! ne=600
! call ree_gam(ng,ne,egam1,egam2,eph1,eph2,nden,fname)
!
!end subroutine tab_param
!
!subroutine nden(x,res)
!! Variables
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! integer, parameter :: n=4
! real(8) :: t(n), u(n), nd
! integer :: i
!! Calculations
! 
! t=[3d0,3d-1,6d-3,2.35d-4]
! u=[5d3,4d4,4d4,0.26d0]
! 
! res=0d0
! do i=1,n
!  call therm_den(u(i),t(i),x,nd)
!  res=res+nd
! end do 
! 
!end subroutine nden
!
! 
! 
!subroutine therm_den(uden,temp,x,res)
!! Variables
!! energy density [eV/cm^3], temperature [eV], photon energy [eV]
! real(8), intent(in) :: uden, temp, x
! real(8), intent(out) :: res ! number density [1/cm^3 eV] in all directions
!! inorm=1/a, where a=int(x^3/(exp(x)-1),x=0..infinity)=Pi^4/15
! real(8), parameter :: inorm=0.153989733820265028d0 
!! Calculations
! res=x**2/(exp(x/temp)-1)
! res=res*uden*inorm/temp**4
!end subroutine therm_den

