!==================================================================
! Calculation of absorption of the photon in the strong magnetic
! and electric fields which are perpendicular to each other.
!
! This file contains 2 modules:
! eb_par_abs
! semf_phabs_m
!
! The 'semf_phabs_m' uses 'eb_par_abs'
!==================================================================


module eb_par_abs
!==================================================================
! Calculation of absorption of photon in a strong parallel!!!
! electric and magnetic fields following prescription of Eq.65 in:
!
! Urrutia, L. F. “Vacuum Polarization in Parallel Homogeneous
! Electric and Magnetic Fields.” Physical Review D 17, no. 8
! (April 15, 1978): 1977–84. doi:10.1103/PhysRevD.17.1977.
!
! It could also be used for calculations in magnetic field alone
!==================================================================
 use gauss_kronrod, only : gk_adaptive
 use phys_const, only : e_charge, hbar, c_light, me_ev,&
 erg_eV, alpha_fs, pi

 implicit none
 private
 save

 public :: k_abs, landau_limit

 real(8), parameter :: me=me_ev*erg_eV/c_light**2
 real(8), parameter :: emc=e_charge/(me*c_light)
 real(8), parameter :: hb_mc2=hbar/(me*c_light**2) ! in sec
 real(8), parameter :: alpha_k=sqrt(3d0)*alpha_fs/(pi*hb_mc2) ! normalization factor for k_abs
 real(8), parameter :: alam=3*emc*hb_mc2/2

 real(8) :: lam

 contains


subroutine landau_limit(eb2,wsin,res,res1)
 real(8), intent(in) :: eb2, wsin
 real(8), intent(out) :: res, res1 ! lower and upper limit
 real(8) :: chi, norm, norm1


 chi=eb2*wsin*alam*2/3

 norm=alpha_fs*emc*eb2
 norm1=sqrt(3**3/(1d0*2**9))
 norm1=norm*norm1
 res=norm1*exp(-8/(3*chi))  ! lower limit

 res1=0.3796123083d0*norm/chi**(1d0/3d0) ! upper limit


end subroutine landau_limit



subroutine k_abs(eb2,we,wsin,res)
!==================================================================
! Input:
! ======
!
! eb2=sqrt(em**2+bm**2), where em and bm are strength of electric
! and magnetic field in cgs units, respectively.
! if only magnetic field is present in a given reference frame
! then eb2=bm
!
! wsin - is the energy of the photon in units of electron rest mass
! multiplied by sin(theta), where theta is the angle between
! direction of photon and a common direction of electric and
! magnetic fields
! we  - the same as wsin but without multiplication by sin(theta)
! Output:
! =======
! res is total probability of the photon to annihilate into
! electron-positron pair per second
!==================================================================
 real(8), intent(in) :: eb2, we, wsin
 real(8), intent(out) :: res
 real(8) :: tmin, tmax

 lam=alam*eb2*wsin

 tmin=0d0
 tmax=1d0
 call gk_adaptive(ik_abs,tmin,tmax,res,tol=1d-4)
 res=res*alpha_k/we

end subroutine k_abs

subroutine ik_abs(x,res)
!==================================================================
! integrand for k_abs
!==================================================================
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: v, v2, iv, rho

 v=x
 v2=v**2
 iv=1/(1-v2)

 rho=4*iv/lam
 call bessel_k23(rho,res)

 res=res*iv*(1-v2/9)

end subroutine ik_abs


subroutine bessel_k23(x,res)
!==================================================================
! BesselK(2/3,x): computation of the modified (hyperbolic) Bessel
! function of the 2/3 order
!==================================================================
 real(8), intent(in) :: x
 real(8), intent(out) :: res

 integer :: kode, n, nz
 real(8) :: fnu, y(1)

 interface
! dbesk(x,fnu,kode,n,y,nz) is a BesselK1 from Slatec library
  subroutine dbesk(x,fnu,kode,n,y,nz)
   integer :: kode, n, nz
   real(8) :: x, fnu, y(n)
  end subroutine dbesk
 end interface

 fnu=2d0/3
 kode=2
 n=1
 call dbesk(x,fnu,kode,n,y,nz)
 res=y(1)*exp(-x)

end subroutine bessel_k23

end module eb_par_abs




module semf_phabs_m
!========================================================================
! Module semf_phabs_m (strong electromagnetic field photon absorption)
! Calculation of the absorption in electric and magnetic field which are
! perpendicular each other.
!
! One can find a reference frame where electric field disappears.
! The velosity of this reference frame is perpendicular to electric and
! magnetic field and is expressed as
!
! v=(ExB)/B^2, where E and B are vectors of electric and magnetic field
! Effectively it lead to magnetic field in the new reference frame
! B'=B*sqrt(1-E^2/B^2), i.e. it decreases its value
! All of these done in 'trans_mf' subroutine.
!
! The calculation of the final result in 'ebgam_abs', see description
! after subroutine declaration
!========================================================================
 implicit none
 private
 save



 public ::  ebgam_abs, ebgam_abs_test, ebgam_abs_trans


 contains



subroutine ebgam_abs_test(em,bm,nk,ek,res,res_min,res_max)
 use phys_const, only : me_ev
 use eb_par_abs, only : k_abs, landau_limit
!=====================================================================
! Test version of ebgam_abs with res_min and res_max small energy
! and large energy limits of the general formula
! Input:
! ======
! em, bm - perpendicular to each other electric and magnetic fields
! nk - photon direction
! ek - photon energy
! Output:
! =======
! res is total probability of the photon to annihilate into
! electron-positron pair per second in current (observer's) ref frame
!=====================================================================
 real(8), intent(in) :: em(3), bm(3), nk(3), ek
 real(8), intent(out) :: res, res_min, res_max
 real(8) :: grf, gvrf, nrf(3), bm1(3)
 real(8) :: nk1(3), ek1
 real(8) :: bm1a, we, wsin


 call trans_mf(em,bm,grf,gvrf,nrf,bm1)

 if (grf>1d0) then
  call trans_photon(grf,gvrf,nrf,nk,ek,nk1,ek1)
 else
  ek1=ek
  bm1=bm
  nk1=nk
  grf=1d0
 end if

 bm1a=sqrt(dot_product(bm1,bm1))

 wsin=dot_product(nk1,bm1)/bm1a ! cos(theta)
 wsin=1-wsin**2
 if (wsin<0d0) wsin=0d0
 wsin=sqrt(wsin)
 we=ek1/me_ev
 wsin=wsin*we

 call k_abs(bm1a,we,wsin,res)
 res=res/grf                 ! time dilation

 call landau_limit(bm1a,wsin,res_min,res_max)

 res_min=res_min/grf
 res_max=res_max/grf

end subroutine ebgam_abs_test



subroutine ebgam_abs(em,bm,nk,ek,res)
! old determine R=R'/Gamma, although it should be R=R'/delta,
! where Gamma is Lorenz factor, delta=1/(Gamma(1-V*n)) is Doppler factor
 use phys_const, only : me_ev
 use eb_par_abs, only : k_abs
!=====================================================================
! Input:
! ======
! em, bm - perpendicular to each other electric and magnetic fields
! nk - photon direction
! ek - photon energy
! Output:
! =======
! res is total probability of the photon to annihilate into
! electron-positron pair per second in current (observer's) ref frame
!=====================================================================
 real(8), intent(in) :: em(3), bm(3), nk(3), ek
 real(8), intent(out) :: res
 real(8) :: grf, gvrf, nrf(3), bm1(3)
 real(8) :: nk1(3), ek1
 real(8) :: bm1a, we, wsin


 call trans_mf(em,bm,grf,gvrf,nrf,bm1)

 if (grf>1d0) then
  call trans_photon(grf,gvrf,nrf,nk,ek,nk1,ek1)
 else
  ek1=ek
  bm1=bm
  nk1=nk
  grf=1d0
 end if

 bm1a=sqrt(dot_product(bm1,bm1))

 wsin=dot_product(nk1,bm1)/bm1a ! cos(theta)

 wsin=1-wsin**2
 if (wsin<0d0) wsin=0d0
 wsin=sqrt(wsin)
 we=ek1/me_ev
 wsin=wsin*we


 call k_abs(bm1a,we,wsin,res)
 res=res*ek1/ek                 ! R=R'/delta, ek1/ek=1/delta


! if (abs(1-ek1/(ek*grf))>5d-1) then
! write(*,'(A,10Es20.12)') 'new, old, new/old=', ek1/ek, 1/grf, abs(1-ek1/(ek*grf))
! read(*,*)
! end if

end subroutine ebgam_abs


subroutine ebgam_abs_trans(em,bm,nk,ek,abs_prob,prim_sin,obs_sin,prim_bm,obs_bm,inv_doppler_fact,lor_fact)
! old determine R=R'/Gamma, although it should be R=R'/delta,
! where Gamma is Lorenz factor, delta=1/(Gamma(1-V*n)) is Doppler factor
 use phys_const, only : me_ev
 use eb_par_abs, only : k_abs
 use vector_ops_mod, only : vector_ops
!=====================================================================
! Input:
! ======
! em, bm - perpendicular to each other electric and magnetic fields
! nk - photon direction
! ek - photon energy
! Output:
! =======
! abs_prob is total probability of the photon to annihilate into
! electron-positron pair per second in current (observer's) ref frame
!=====================================================================
 real(8), intent(in) :: em(3), bm(3), nk(3), ek
 real(8), intent(out) :: abs_prob, prim_sin, obs_sin, prim_bm, obs_bm, inv_doppler_fact, lor_fact
 real(8) :: grf, gvrf, nrf(3), bm1(3)
 real(8) :: nk1(3), ek1
 real(8) :: bm1a, we, wsin
 type(vector_ops) :: vop


 call trans_mf(em,bm,grf,gvrf,nrf,bm1)

 if (grf>1d0) then
  call trans_photon(grf,gvrf,nrf,nk,ek,nk1,ek1)
 else
  ek1=ek
  bm1=bm
  nk1=nk
  grf=1d0
 end if

 bm1a=sqrt(dot_product(bm1,bm1))


 wsin=dot_product(nk1,bm1)/bm1a ! cos(theta)

 wsin=1-wsin**2
 if (wsin<0d0) wsin=0d0
 wsin=sqrt(wsin)

 we=ek1/me_ev
 wsin=wsin*we


 call k_abs(bm1a,we,wsin,abs_prob)
 inv_doppler_fact=ek1/ek 
 abs_prob=abs_prob*inv_doppler_fact                ! R=R'/delta, ek1/ek=1/delta

 inv_doppler_fact=ek1/ek
 lor_fact=grf


 prim_sin=vop%ab_sin(bm1,nk1)
 obs_sin=vop%ab_sin(bm,nk)

 prim_bm=vop%vec_norm(bm1)
 obs_bm=vop%vec_norm(bm)



! if (abs(1-ek1/(ek*grf))>5d-1) then
! write(*,'(A,10Es20.12)') 'new, old, new/old=', ek1/ek, 1/grf, abs(1-ek1/(ek*grf))
! read(*,*)
! end if

end subroutine ebgam_abs_trans





subroutine ebgam_abs_old(em,bm,nk,ek,res)
! old determine R=R'/Gamma, although it should be R=R'/delta,
! where Gamma is Lorenz factor, delta=1/(Gamma(1-V*n)) is Doppler factor
 use phys_const, only : me_ev
 use eb_par_abs, only : k_abs
!=====================================================================
! Input:
! ======
! em, bm - perpendicular to each other electric and magnetic fields
! nk - photon direction
! ek - photon energy
! Output:
! =======
! res is total probability of the photon to annihilate into
! electron-positron pair per second in current (observer's) ref frame
!=====================================================================
 real(8), intent(in) :: em(3), bm(3), nk(3), ek
 real(8), intent(out) :: res
 real(8) :: grf, gvrf, nrf(3), bm1(3)
 real(8) :: nk1(3), ek1
 real(8) :: bm1a, we, wsin


 call trans_mf(em,bm,grf,gvrf,nrf,bm1)

 if (grf>1d0) then
  call trans_photon(grf,gvrf,nrf,nk,ek,nk1,ek1)
 else
  ek1=ek
  bm1=bm
  nk1=nk
  grf=1d0
 end if

 bm1a=sqrt(dot_product(bm1,bm1))

 wsin=dot_product(nk1,bm1)/bm1a ! cos(theta)
 wsin=1-wsin**2
 if (wsin<0d0) wsin=0d0
 wsin=sqrt(wsin)
 we=ek1/me_ev
 wsin=wsin*we


 call k_abs(bm1a,we,wsin,res)
 res=res/grf                 ! time dilation

end subroutine ebgam_abs_old



subroutine trans_mf(em,bm,grf,gvrf,nrf,bm1)
!====================================================================
! Input:
! ======
! em, bm - electric and magnetic fields in CGS units
! Output:
! =======
! grf - Lorentz factor of moving (new) reference frame
! gvrf - Lorentz factor * velosity magnitude
! nrf - direction of the new reference frame
! bm1 - magnetic field in the new ref frame
!====================================================================
 use utools_mod, only : utools
 real(8), intent(in) :: em(3), bm(3)
 real(8), intent(out) :: grf, gvrf, nrf(3)
 real(8), intent(out) :: bm1(3) ! velosity of the reference frame
 real(8) :: vrf(3), vrf2, vrfa
 type(utools) :: ut

 call ut%cross_prod(em,bm,vrf)
 vrf=vrf/dot_product(bm,bm)
 vrf2=dot_product(vrf,vrf)
 vrfa=sqrt(vrf2)
 nrf=vrf/vrfa          ! unit direction vector of ref frame velosity

 if (vrf2>=1d0) then
  vrf2=vrf2-1d-8
  vrfa=sqrt(vrf2)
 end if

 grf=1/sqrt(1-vrf2)   ! Lorentz factor of ref frame
 gvrf=grf*vrfa        ! Lorentz factor * velosity magnitude

 bm1=bm/grf           ! magnetic field in the new ref frame

! write(*,*) 'vrf=', vrf

end subroutine trans_mf



subroutine trans_photon(glf,gv,nv,nk,ek,nk1,ek1)
!====================================================================
! Calculation of the energy and direction of the photon in
! the new reference frame moving in direction nv with velosity v
! Input:
! ======
! glf - Lorentz factor glf=1/sqrt(1-v^2), where v in speed of light units
! gv=glf*v - relativistic velosity
! nv(3)=[nvx,nvy,nvz] - direction of the velosity in initial
! reference frame (unit vector nv^2=1)
! nk(3) - direction of the photon in initial reference frame
! ek - energy of the photon, in any energy units
! Output:
! =======
! nk1 - direction of the photon in the new reference frame
! ek1 - energy in the new reference frame, in the same units as ek
!====================================================================
 real(8), intent(in) :: glf, gv, nv(3)
 real(8), intent(in) :: nk(3), ek
 real(8), intent(out) :: nk1(3), ek1
 real(8) :: xi, nka

 xi=dot_product(nv,nk) ! angle between velocity and photon momentum
 nk1=nk+gv*(gv*xi/(1+glf)-1)*nv
 nka=dot_product(nk1,nk1)
 nk1=nk1/sqrt(nka)     ! direction of in the new reference frame
 ek1=ek*(gv*(1-xi)+1/(glf+gv)) ! energy in the new reference frame

end subroutine trans_photon


end module semf_phabs_m

! Example of usage
!
! program main
!  use utools_mod, only : utools
!  use semf_phabs_m, only : ebgam_abs, ebgam_abs_test
!  type(utools) :: ut
!  real(8) :: bm(3), vr(3), em(3)
!  real(8) :: nk(3), ek, res, resmin, resmax
!  real(8), allocatable :: ea(:)
!
!
!  call ut%unit_vec(0d0,0d0,bm)
!  bm=bm*1d6
!  call ut%unit_vec(0d0,10d0,vr)
!  vr=1d-2*vr
!  call ut%cross_prod(vr,bm,em)
!  em=-em
!  call ut%unit_vec(0d0,90d0,nk)
!
!  call ut%grid(ea,1d9,1d15,100)
!
!
!  open(1,file='abs_test.dat')
!  do i=1,size(ea)
!   ek=ea(i)
!   call ebgam_abs_test(em,bm,nk,ek,res,resmin,resmax)
!   write(1,*) ea(i), res, resmin, resmax
!  end do
!  close(1)
!
! end program main
