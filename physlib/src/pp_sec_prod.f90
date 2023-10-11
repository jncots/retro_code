! written by Andrew M Taylor, Oct 2014
! modified for fortran 90 by Anton Prosekin, Feb 2015
!
! For usage put the line 'use ppgam_cross, only : MinProtonEn, ppGamMinMax, ppGamCrossSec'
! before any other commands (before implicit none)
! Only subroutines MinProtonEn, ppGamMinMax, ppGamCrossSec are accessible outside
! (see description for each one in the head of corresponding subroutine)
!
module ppgam_cross
!
! Only subroutines MinProtonEn, ppGamMinMax, ppGamCrossSec are accessible outside
! (see description for each one in the head of corresponding subroutine)
! For usage put the line 'use ppgam_cross, only : ppGamMinMax, ppGamCrossSec'
! before any other commands (before implicit none)
!
 implicit none
 private
 public :: MinProtonEn, ppGamMinMax, ppGamCrossSec
 real(8), parameter :: mpi=0.134976d0 ! [GeV] pi-zero mass
 real(8), parameter :: mp=0.938272d0 ! [GeV] proton mass
 real(8), parameter :: mres=1.1883d0 ! [GeV] resonance mass
 real(8), parameter :: gres=0.2264d0 ! [GeV] resonance width
 real(8), parameter :: sigma0=7.66d-3 ! [mb] (millibarn)  
 real(8), parameter :: Tpth=0.279660549640189625d0 ! Tpth=2.0*mpi+mpi**2/(2.0*mp)
! gamma=sqrt(Mres**2*(Mres**2+Gres**2))
! k_param=sqrt(8.0)*Mres*Gres*gamma/(pi*sqrt(Mres**2+gamma))
 real(8), parameter :: k=0.206256239756209630d0
 real(8), parameter :: mpip=mpi/mp
! Parameters a
 real(8), parameter :: a1(5)=[0.728d0,0.596d0,0.491d0,0.2503d0,0.117d0]
 real(8), parameter :: a2(5)=[0.652d0,0.0016d0,0.488d0,0.1928d0,0.483d0]
 real(8), parameter :: a3(5)=[5.436d0,0.254d0,0.072d0,0.075d0,0.166d0]
 real(8), parameter :: a4(5)=[0.908d0,9.0d-4,6.089d0,0.176d0,0.448d0]
 real(8), parameter :: a(5,4)=reshape([a1,a2,a3,a4],[5,4])
! Parameters b 
 real(8), parameter :: b0=5.9d0, b1low(3)=[9.53d0,0.52d0,0.054d0], b1high(3)=[9.13d0,0.35d0,9.7d-3]
 real(8), parameter :: b2(3)=[9.06d0,0.3795d0,0.01105d0], b3(3)=[10.77d0,0.412d0,0.01264d0]
 real(8), parameter :: b4(3)=[13.16d0,0.4419d0,0.01439d0]
 real(8), parameter :: b(3,0:4)=reshape([b1low,b1high,b2,b3,b4],[3,5])
 real(8), parameter :: GeV=1d-9, eV=1d9 ! [GeV/eV]
 real(8), parameter, public :: millibarn=1d-27 ! [cm^2/millibarn]
 real(8), parameter, public :: millibarnGeV=millibarn*GeV
 
 integer, save :: m
 real(8), save :: Tp, s, Tptilda
  
 contains

 
subroutine MinProtonEn(Eg0,Epmin)
!
! Calculates the minimal energy of proton Epmin[in eV] which can produce the photon
! with energy Eg0 [in eV].
!
 implicit none
! Arguments
 real(8), intent(in) :: Eg0
 real(8), intent(out) :: Epmin
! Variables
 real(8) :: Ega, g, x, p, gp, x2, g2
 real(8) :: p1, p2, p3, p4
! Calculations

 Ega=2*Eg0*GeV/mpi
 g=(1+Ega**2)/(2*Ega)
 x=mpip 
 
 p=(2*g-x)*x
 gp=g*x-2
 x2=x*x
 g2=g*g-1
 
 p1=p*gp+x2
 p2=sqrt(x2*g2*(p+8)*p)
 p3=2*(p-1)
 p4=2+(p1+p2)/p3 
 
 Epmin=(p4-1)*mp*eV
 
end subroutine MinProtonEn


subroutine ppGamMinMax(Ep,Eph1,Eph2,Epth)
!
! Calculation of minimal and maximal photon energy for a given energy of proton
! and output of threshold energy.
! All energies should be given in [eV]
!
 implicit none
! Arguments
 real(8), intent(in) :: Ep
 real(8), intent(out) :: Eph1, Eph2, Epth
! Variables
 real(8) :: Epim_LAB, Egmin, Egmax
! Calculation
 Tp=Ep*GeV-mp
 
 if (Tp<Tpth) then
  Eph1=0d0
  Eph2=0d0
  Epth=(Tpth+mp)*eV
  write(*,*) 'GamMinMax: The output gives zeros.'
  write(*,'(A,Es14.6)') 'The energy of proton should be higher that threshold E_th [in eV]:', Epth
  return
 end if
 
! if (Tp+mp>1d6) then
!  write(*,*) 'Too high energy. The calculation could be inaccurate.'
! end if
 
 call get_Eg_minmax(Epim_LAB,Egmin,Egmax)
 Eph1=Egmin*eV
 Eph2=Egmax*eV
 Epth=(Tpth+mp)*eV
end subroutine ppGamMinMax

subroutine ppGamCrossSec(model,Ep,Eg1,res)
!
! Calculation of d\sigma/dEg(Ep,Eg) for model=1..4
! Geant4(1), pythia8(2), SIBYLL(3), QGSJET(4)
! proton energy Ep and photon energy Eg1 in [eV]
! res is in [cm^2/eV]
!
 implicit none
! Arguments
 integer, intent(in) :: model
 real(8), intent(out) :: Ep, Eg1 ! in eV
 real(8), intent(out) :: res
! Variables
 real(8) :: Eg
 real(8) :: Epim_LAB, Egmin,Egmax
 real(8) :: sigmapi0, Amax, Xgtildam,C,C1,C2,C3,C4,alpha,beta, ds_dEg
 
 m=model
 Tp=Ep*GeV-mp
 Tptilda=Tp/mp
 Eg=Eg1*GeV
 
 if (Tp<Tpth) then
  res=0d0
  return
 end if 
  
! tracking=0
 
 call get_Eg_minmax(Epim_LAB,Egmin,Egmax)
 
 if ((Eg<Egmin).or.(Eg>Egmax)) then
  res=0d0
  return
 end if
 
 call get_sigmapi0(sigmapi0)
 call get_Amax(Epim_LAB,sigmapi0,Amax)
 call get_parameters3(Egmax,Xgtildam,C,C1,C2,C3,C4,alpha,beta)
 call get_fX(Eg,Xgtildam,Amax,beta,alpha,C,C1,C2,C3,C4,ds_dEg)
! if ((Tp>=Tpth).and.(Eg>=Egmin).and.(Eg<Egmax)) then
!  if (isnan(ds_dEg)) print*,'Eg:',Eg,' Tp:',Tp
! endif

  if ((Tp>=Tpth).and.(Eg>=Egmin).and.(Eg<Egmax)) then
   if (isnan(ds_dEg)) ds_dEg=0d0
  endif
 
! res=ds_dEg*millibarn*GeV
  res=ds_dEg ! in millibarn/GeV
 
! if (EgdNdE<=1.0e-40) EgdNdE=0d0

end subroutine ppGamCrossSec
 
subroutine get_Eg_minmax(Epim_LAB,Egmin,Egmax)
!
! Probably calculates min and max energy of photons Egmin, Egmax
! for kinetic energy of proton Tp
!
 implicit none
! Arguments 
 real(8), intent(out) :: Epim_LAB, Egmin, Egmax
! Variables 
 real(8) :: Epim_CM,ppim_CM,gamma_CM,beta_CM
 real(8) :: gamma_LAB,beta_LAB, gp
! Calculations
 
 s=2.0*mp*(Tp+2.0*mp)
 Epim_CM=(s-4.0*mp**2+mpi**2)/(2.0*sqrt(s))
 ppim_CM=sqrt(Epim_CM**2-mpi**2)

 gamma_CM=(Tp+2.0*mp)/(sqrt(s))
 beta_CM=sqrt(1.0-gamma_CM**(-2))

 Epim_LAB=gamma_CM*(Epim_CM+beta_CM*ppim_CM)
      
 gamma_LAB=Epim_LAB/mpi
 beta_LAB=sqrt(1.0-gamma_LAB**(-2))

!      Egmax=Tp-Tpth ! approximate expression

 gp=gamma_LAB*(1.0+beta_LAB)
 Egmin=(mpi/2.0)/gp
 Egmax=(mpi/2.0)*gp

end subroutine get_Eg_minmax


subroutine get_sigmapi0(sigmapi0)
!
!
!
 implicit none
! Arguments
 real(8), intent(out) :: sigmapi0
! Variables
 real(8) :: s_inel, eta, fBW, sigma1pi, sigma2pi, Q, npi0, Q2
! Calculations

 eta=sqrt((s-mpi**2-4.0*mp**2)**2-16.0*mpi**2*mp**2)/(2.0*mpi*sqrt(s))

! fBW=K/(((sqrt(s)-mp)**2-Mres**2)**2+(Mres**2)*(Gres**2))
 fBW=mp*K/(((sqrt(s)-mp)**2-Mres**2)**2+(Mres**2)*(Gres**2))

 if (Tp.le.2.0) then
  sigma1pi=sigma0*eta**(1.95)*(1.0+eta+eta**5)*(fBW**1.86)
 endif

 s_inel=(30.7-0.96*log(Tp/Tpth)+0.18*(log(Tp/Tpth)**2))*((1.0-(Tpth/Tp)**1.9)**3)

 if (Tp.lt.0.56) then
  sigma2pi=0.0
 elseif (Tp.ge.0.56.and.Tp.le.2.0) then
  sigma2pi=5.7/(1.0+exp(-9.3*(Tp-1.4)))
 endif

! if (tracking.eq.1) then
!  print*,'s_inel:',s_inel
! endif

 if (Tp.lt.2.0) then
  sigmapi0=sigma1pi+sigma2pi
 elseif (Tp.ge.2.0.and.Tp.lt.5.0) then
  Q=(Tp-Tpth)/mp
  npi0=-6.0*10**(-3.0)+0.237*Q-0.023*Q**2
  sigmapi0=npi0*s_inel
 elseif (Tp.ge.5.0.and.Tp.le.100.0) then
  Q2=(Tp-3.0)/mp
  npi0=a(1,1)*(Q2**a(4,1))*(1.0+exp(-a(2,1)*(Q2**a(5,1))))*(1.0-exp(-a(3,1)*(Q2**0.25)))
  sigmapi0=npi0*s_inel
 elseif (Tp.gt.100.0) then
  Q2=(Tp-3.0)/mp
  npi0=a(1,M)*(Q2**a(4,M))*(1.0+exp(-a(2,M)*(Q2**a(5,M))))*(1.0-exp(-a(3,M)*(Q2**0.25)))
  sigmapi0=npi0*s_inel
 endif

 if (Tp.gt.50.0.and.M.eq.2) then
  Q2=(Tp-3.0)/mp
  npi0=a(1,M)*(Q2**a(4,M))*(1.0+exp(-a(2,M)*(Q2**a(5,M))))*(1.0-exp(-a(3,M)*(Q2**0.25)))
  sigmapi0=npi0*s_inel
 endif

end subroutine get_sigmapi0

subroutine get_Amax(Epim_LAB,sigmapi0,Amax)
!
!
!
 implicit none
! Arguments 
 real(8), intent(in) :: Epim_LAB, sigmapi0
 real(8), intent(out) :: Amax
! Variables
 integer :: m1
! Calculations
      
 if (Tp.ge.1.0.and.Tp.lt.5.0) then
  m1=0
 elseif (Tp.ge.5.0) then
  m1=1
 endif 
      
! low energy maximum
 if (Tp.lt.1.0) then
  Amax=b0*sigmapi0/Epim_LAB
! Geant4 (default)
 elseif (Tp.ge.1.0.and.Tp.le.100.0) then
  Amax=b(1,m1)*Tptilda**(-b(2,m1))*exp(b(3,m1)*(log(Tptilda))**2)*sigmapi0/mp      
 elseif (Tp.gt.100.0) then
  Amax=b(1,M)*Tptilda**(-b(2,M))*exp(b(3,M)*(log(Tptilda))**2)*sigmapi0/mp
 end if

! Pythia is exception since it can take over at 50 GeV
 if (Tp.gt.50.0.and.M.eq.2) then
  Amax=b(1,M)*Tptilda**(-b(2,M))*exp(b(3,M)*(log(Tptilda))**2)*sigmapi0/mp
 end if
end subroutine get_Amax


subroutine get_parameters3(Egmax,Xgtildam,C,C1,C2,C3,C4,alpha,beta)
!
!
!   
 implicit none
 real(8) :: Egmax,Q,Q2,Xgtildam
 real(8) :: C,C1,C2,C3,C4,mu,alpha,beta

 Q=(Tp-0.3)/mp
 Q2=(Tp-1.0)/mp
            
 Xgtildam=Egmax+(mpi**2)/(4.0*Egmax)

! if (tracking.eq.1) then
!  print*,'Q:',Q
! endif

 C=5.0*mpi/Xgtildam
 C1=3.0*mpi/Xgtildam
 C2=3.5*mpi/Xgtildam
 C3=3.55*mpi/Xgtildam
 C4=3.55*mpi/Xgtildam

 mu=1.2*Q2**1.2*exp(-1.2*Q2)

 if (Tp.lt.1.0) then
  alpha=0.0
  beta=3.29-0.2*Tptilda**(-1.5)
! Geant4 (default)
 elseif (Tp.ge.1.0) then
  if (Tp.le.4.0) then
   alpha=mu+1.45
   beta=mu+2.45
  elseif (Tp.gt.4.0.and.Tp.le.20.0) then
   alpha=mu+1.5
   beta=1.5*mu+4.95
  elseif (Tp.gt.20.0.and.Tp.le.100.0) then
   alpha=1.0
!  beta=4.4
   beta=4.2
  elseif (Tp.gt.100.0) then
   alpha=1.0
   beta=4.9
  endif

 endif

! pythia8
 if (M.eq.2.and.Tp.gt.50.0) then
  alpha=1.0
  beta=4.00001
! SIBYLL
 elseif (M.eq.3.and.Tp.gt.100.0) then
  alpha=1.0
  beta=3.6
! QGSJET
 elseif (M.eq.4.and.Tp.gt.100.0) then
  alpha=1.0
  beta=4.5
 endif

! if (tracking.eq.1) then
!  print*,'C:',C,'Xgtildam:',Xgtildam
!  print*,'alpha:',alpha,' beta:',beta
! endif

end subroutine get_parameters3


subroutine get_fX(Eg,Xgtildam,Amax,beta,alpha,C,C1,C2,C3,C4,ds_dEg)
!
!
!
 implicit none
! Arguments
 real(8), intent(in) :: Eg, Xgtildam, Amax, beta ,alpha, C,C1,C2,C3,C4
 real(8), intent(out) :: ds_dEg
! Variables
 real(8) :: Xgtilda, Xg
  
 Xgtilda=Eg+(mpi**2.0)/(4.0*Eg)
 Xg=(Xgtilda-mpi)/(Xgtildam-mpi)
! low energies
 if (Tp.lt.1.0) then
  ds_dEg=Amax*(1.0-Xg)**beta/((1.0+Xg/C)**(alpha))
! Geant4 (default)
 elseif (Tp.ge.1.0) then
  if (Tp.le.20.0) then
   ds_dEg=Amax*(1.0-Xg)**beta/((1.0+Xg/C1)**(alpha))
  elseif (Tp.gt.20.0) then
   ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C1)**(alpha))
  endif
 endif

! pythia8
 if (M.eq.2.and.Tp.gt.50.0) then
  ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C2)**(alpha))
! SIBYLL
 elseif (M.eq.3.and.Tp.gt.100.0) then
  ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C3)**(alpha))
! QGSJET
 elseif (M.eq.4.and.Tp.gt.100.0) then
  ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C4)**(alpha))
 endif

end subroutine get_fX
 

end module ppgam_cross



module pp_sec_prod
 use gauss_kronrod, only : gk_adaptive
 use phys_const, only : millibarn, c_light
 use ppgam_cross, only : MinProtonEn, ppGamCrossSec
 implicit none
 private
 save
 
 public :: set_protsp, pp_el_dist, pp_gam_dist, set_protsp_min_max
 
 real(8), parameter :: kpi=1d0/0.17d0
 real(8), parameter :: mprot=0.938272046d9 ! eV
 real(8), parameter :: mpion=0.1349766d9 ! eV
 real(8), parameter :: mpi2=mpion**2
 real(8) :: emin_dfa
 real(8) :: lep, xgm, Ee_pp
 real(8) :: epmin, epmax
 real(8) :: Egam
 integer :: isw
 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 procedure (fun_prot), pointer :: prot_sp=>null()
  
 contains




subroutine set_protsp(nden,pdist,norm)
! Variables 
 real(8), intent(in) :: nden ! in 1/cm^3
 real(8), intent(out) :: norm
 
 interface
  subroutine pdist(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
! Calculations 
 
 prot_sp=>pdist
 norm=millibarn*c_light*nden
 isw=0
end subroutine set_protsp


subroutine set_protsp_min_max(x1,x2)
! Variables
 real(8), intent(in) :: x1, x2 ! min and max of the spectrum
 epmin=x1
 epmax=x2
end subroutine set_protsp_min_max





subroutine pp_gam_dist(x,res)
!
! Calculation of the gamma-ray spectrum where x is gamma ray energy
! and res is production rate using Kafexiu et al
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8), parameter :: igev=1d-9 ! take into account that cross section is millibarn/GeV -> millibarn/eV
 real(8) :: ep1, tmin, tmax
! Calculations
 
 if (abs(epmax-epmin)<1d-307) then
  res=0d0
  return
 end if 
 
 
 Egam=x
 call MinProtonEn(Egam,ep1)
 if (ep1<epmin) ep1=epmin
 
 tmin=log(ep1)
 tmax=log(epmax)
 
 call gk_adaptive(igam,tmin,tmax,res,tol=1d-4)
 res=res*igev
end subroutine pp_gam_dist


subroutine igam(x,res)
!
! Intergrand for intergration with proton spectrum
!
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 integer, parameter :: model=3 ! pythia8
 real(8) :: Ep, sig, fp
! Calculations

 Ep=exp(x)
 call ppGamCrossSec(model,Ep,Egam,sig)
 call prot_sp(Ep,fp)
 res=sig*fp*Ep
 
end subroutine igam


subroutine pp_el_dist(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8), parameter :: Ejoin=1d11
! integer :: isw=0
 real(8) :: res_k, res_d, cadj
 save
! Calculations 
 if (isw==0) then
  call kel_sp(Ejoin,res_k)
  call dfa_sp(Ejoin,res_d)
  cadj=res_k/res_d
  if (isnan(cadj)) cadj=1d0
  isw=1
 end if
  
 if (x>Ejoin) then
  call kel_sp(x,res)
 else
  call dfa_sp(x,res)
  res=cadj*res
 end if
end subroutine pp_el_dist



subroutine kel_sp(x,res)
!
! Energy spectrum of electrons in pp with E>100 GeV, according Kelner et al 2006
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: tmin, tmax, fe
! Calculations
 Ee_pp=x
 
 tmin=Ee_pp/epmax ! integration of proton spectrum over interval (epmin,epmax)
 tmax=Ee_pp/epmin 
 
 if (tmax>1d0) tmax=1d0
 if (tmin<1d-3) tmin=1d-3
 if (tmin>=tmax) then
  res=0d0
  return
 end if

 call gk_adaptive(ielec,tmin,tmax,fe,tol=1d-4)
 res=fe
 
end subroutine kel_sp


subroutine ielec(x,res)
!
! Integrand
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: Ep, fe, fp 
! Calculations

 Ep=Ee_pp/x
 lep=log(Ep/1d12) ! log(Ep/1TeV) 
 xgm=x
 
 call pp_el(fe)
 call prot_sp(Ep,fp) 
 res=fe*fp/x
  
end subroutine ielec




subroutine pp_el(res)
!
! Distribution function for electrons(positrons) produced in
! pp interactsions (Kelner et al 2006, eq. 62) * cross section
!
! Variables
 real(8), intent(out) :: res
 real(8) :: l2, b_e, beta_e, k_e, lnx, sig
! Calculations 
 l2=lep**2
 b_e=1d0/(69.5d0+2.65d0*lep+0.3d0*l2)
 beta_e=1d0/(0.201d0+6.2d-2*lep+4.2d-4*l2)**(0.25d0)
 k_e=(0.279d0+0.141d0*lep+1.72d-2*l2)
 k_e=k_e/(5.59d0+4.6d0*lep+l2)
 
 lnx=-dlog(xgm)
 
 res=(1d0+k_e*lnx**2)**3
 res=res*b_e*lnx**5
 res=res/(xgm*(1+0.3d0/(xgm**beta_e)))
! Inelastic cross section 
 sig=34.3d0+1.88*lep+0.25d0*l2
 res=res*sig 
end subroutine pp_el


subroutine dfa_sp(x,res)
!
! Spectrum of electons in delta-functional approximation
!
! Variables
 real(8), intent(in) :: x ! energy of secondary particle, in eV
 real(8), intent(out) :: res ! dN/dE, energy spectrum
 real(8) :: tmin, tmax, em_kp
! Calculations
! tmin=0d0
! tmax=1d0
 emin_dfa=x+mpi2/(4*x)
 
 em_kp=emin_dfa*kpi
 tmin=em_kp/(epmax-mprot)
 tmax=em_kp/(epmin-mprot)
  
 if ((tmax>1d0).or.(tmax<0d0)) tmax=1d0
 if (tmin>=tmax) then
  res=0d0
  return
 end if
  
 call gk_adaptive(pi_prod_rate,tmin,tmax,res,tol=1d-4)
 res=2*emin_dfa*kpi*res
end subroutine dfa_sp


subroutine pi_prod_rate(x,res)
!
! Production rate of pi mesons in delta-functional approximation
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8), parameter :: eth=1.22d9, ceth=1d11 ! eV, threshold energy
 real(8) :: epi, xp, lp, sig, fp
! Calculations
 epi=emin_dfa/x
 xp=mprot+epi*kpi

! Inelastic cross section 
 if (xp<eth) then
  res=0d0
  return
 end if
 
 lp=log(xp*1d-12)
 sig=34.3d0+(1.88d0+0.25d0*lp)*lp
 if (xp<ceth) sig=sig*(1-(eth/xp)**4)**2
! Production rate 
 call prot_sp(xp,fp)
 res=sig*fp
 res=res/sqrt((epi-mpion)*(epi+mpion))
 res=res/x**2 
end subroutine pi_prod_rate


end module pp_sec_prod