module eb_par_abs
!==================================================================
! Calculation of absorption of photon in a strong parallel!!! 
! electric and magnetic fields following prescription of Eq.65 in:
!
! Urrutia, L. F. “Vacuum Polarization in Parallel Homogeneous 
! Electric and Magnetic Fields.” Physical Review D 17, no. 8 
! (April 15, 1978): 1977–84. doi:10.1103/PhysRevD.17.1977.
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
 res=norm1*exp(-8/(3*chi))
 
 res1=0.3796123083d0*norm/chi**(1d0/3d0)
 
 
end subroutine landau_limit 
 
 
 
 
subroutine k_abs(eb2,we,wsin,res)
!==================================================================
! Input:
! ======
! eb2=sqrt(em**2+bm**2), where em and bm are strength of electric
! and magnetic field in cgs units, respectively.
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
 ! write(*,'(A,10Es14.6)') 'k_abs lam=', lam
 ! write(*,*) lam
 
 tmin=0d0
 tmax=1d0
 call gk_adaptive(ik_abs,tmin,tmax,res,tol=1d-14)
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


! program main
 ! use eb_par_abs, only : k_abs
 ! use phys_const, only : c_light
 ! implicit none
 ! real(8) :: eb2, wsin, res
 ! real(8) :: x1, x2
 ! integer :: i, nx
 
 ! eb2=1d6
 ! x1=1d6
 ! x2=1d8
 ! nx=100
 
 ! open(1,file='lam.dat')
 ! do i=0,nx
  ! wsin=x1*(x2/x1)**(1d0*i/nx)
  ! call k_abs(eb2,wsin,res)
  ! res=c_light/res
  ! write(1,*) wsin, res
 ! end do 
 ! close(1)
 
 ! ! integer :: i, nx
 ! ! real(8) :: x, res
 ! ! real(8) :: x1, x2
 
 ! ! x1=1d-2
 ! ! x2=1d3
 ! ! nx=100
 
 ! ! open(1,file='lam.dat')
 ! ! do i=0,nx
  ! ! x=x1*(x2/x1)**(1d0*i/nx)
  ! ! call k_abs(x,res)
  ! ! write(1,*) x, res
 ! ! end do 
 ! ! close(1)

 
! end program main



module eb_gam_abs
!==================================================================
! Calculation of the total probability of the photon to annihilate 
! into electron-positron pair per second in current ref frame
!==================================================================
 use eb_par_abs, only : k_abs
 implicit none
 private

 public :: ebgam_abs,  drift_dir, drift_dir_app
 
 
 contains
 

subroutine drift_dir(ef,bf,beta,tp)
!==================================================================
! Drift trajectory direction
! tp - type of particle, is tp='elc' (electron) by default
! tp='pos' for positrons
!==================================================================
 real(8), intent(in) :: ef(3), bf(3)
 real(8), intent(out) :: beta(3)
 character(3), intent(in), optional :: tp
 character(3), parameter :: tp_def='elc'
 character(3) :: tpp
 real(8) :: ef2, bf2, ebm, ebf, q, exb(3), ss
 
 if (present(tp)) then
  tpp=tp
 else 
  tpp=tp_def
 end if  
 
 ef2=dot_product(ef,ef)
 bf2=dot_product(bf,bf)
 
 ebm=(ef2-bf2)/2
 ebf=dot_product(ef,bf)
 ss=ebf/abs(ebf)         ! relative direction
 q=sqrt(ebm**2+ebf**2)+(ef2+bf2)/2
 
 call cross_prod(ef,bf,exb)
 
 if (tpp=='elc') then
  beta=exb-(sqrt(q-bf2)*ef+sqrt(q-ef2)*ss*bf) ! electrons
 else
  beta=exb+(sqrt(q-bf2)*ef+sqrt(q-ef2)*ss*bf) ! positrons
 end if 
 beta=beta/q
 
end subroutine drift_dir

subroutine drift_dir_app(ef,bf,beta)
!==================================================================
! Drift trajectory direction (approximation if electric field
! is smaller than magnetic one)
!==================================================================
 real(8), intent(in) :: ef(3), bf(3)
 real(8), intent(out) :: beta(3)
 real(8) :: bf2, exb(3), ff
 
 bf2=dot_product(bf,bf)
 call cross_prod(ef,bf,exb)
 
 beta=exb/bf2
 
 ff=dot_product(beta,beta)
 ff=sqrt(1-ff)
 beta=beta+ff*bf/sqrt(bf2)
 
end subroutine drift_dir_app



subroutine ebgam_abs(em,bm,np,ep,res)
!==================================================================
! Input:
! ======
! em, bm - electric and magnetic fields
! np - photon direction
! ep - photon energy
! Output:
! =======
! res is total probability of the photon to annihilate into 
! electron-positron pair per second in current ref frame
!==================================================================
 real(8), intent(in) :: em(3), bm(3), np(3), ep
 real(8), intent(out) :: res 
 real(8) :: em1(3), bm1(3), nf(3), vf, gf
 real(8) :: np1(3), dopf
 real(8) :: em12, bm12, eb2, wsin, we
 
 call par_frame(em,bm,em1,bm1,nf,vf,gf)
 if (gf>1d0) then
  call n1_dir(np,nf,vf,gf,np1,dopf)
 else
  em1=em
  bm1=bm
  np1=np
  dopf=1d0
  gf=1d0  
 end if 
 
 em12=dot_product(em1,em1)
 bm12=dot_product(bm1,bm1)
 eb2=sqrt(em12+bm12)
 
 wsin=dot_product(np1,bm1)/sqrt(bm12) ! cos(theta)
 wsin=sqrt(1-wsin**2)
 we=ep*dopf
 wsin=wsin*we
 
! write(*,'(A,10Es20.12)')'we, eb2, wsin=', we, eb2, wsin
! read(*,*)
 
 call k_abs(eb2,we,wsin,res)
 res=res/gf                 ! time dilation

end subroutine ebgam_abs
 
  
subroutine n1_dir(np,nf,vf,gf,np1,dopf)
!=================================================================
! Input:
! ====== 
! np(3) - direction of velosity (unit vector) of the particle in 
! initial ref frame
! nf(3) - direction of velosity (unit vector) of the reference 
! frame in initial ref frame
! vf - velocity of the reference frame in units of speed of light
! gf - Lorentz factor factor=1/sqrt(1-vf**2)
!
! Output:
! =======
! np1(3) - direction of velosity (unit vector) of the particle in 
! moving ref frame
! dopf - 1/D, where D=1/(gamma*(1-(\vec V*\vec n))) is doppler 
! factor. The energy of photon in the new ref frame is
! E'=E*dopf, where E is photon energy in initial ref frame
!==================================================================
 real(8), intent(in) :: np(3), nf(3)
 real(8), intent(in) :: vf, gf
 real(8), intent(out) :: np1(3), dopf
 real(8) :: gv, nv, np2
 
 gv=gf*vf
 nv=dot_product(np,nf)
 
 np1=(gv*nv/(1+gf)-1)*gv*nf+np 
 np2=dot_product(np1,np1)
 np2=sqrt(np2)
 np1=np1/np2                  ! new direction of particle
 
 dopf=(1-nv)*gv+1/(gf*(1+vf)) ! 1/D, D is doppler factor 
end subroutine n1_dir

 
 
subroutine par_frame(em,bm,em1,bm1,nf,vf,gf)
!==================================================================
! Input:
! ======
! em(3) - electric field in initial ref frame
! bm(3) - magnetic field in initial ref frame
!
! Output:
! ======
! em1(3) - electric field in new ref frame, 
! where magnetic and electric field are parallel
! em1(3)|| bm1(3) and new ref frame are moving perpendicular
! to em and bm (see Landau).
! bm1(3) - magnetic field in new ref frame
! nf - is the direction of velosity (unit vector) 
! of new ref frame
! vf - velosity of a new ref frame in units of speed of light
! gf - Lorentz factor factor=1/sqrt(1-vf**2)
!==================================================================
 real(8), intent(in) :: em(3), bm(3)
 real(8), intent(out) :: em1(3), bm1(3), nf(3), vf, gf
 real(8) :: exb(3), em2, bm2, exb2, exb1, eb
 real(8) :: ebm2, rho, rho2, rs, vr, aa, aeb
 
 call cross_prod(em,bm,exb)
 exb2=dot_product(exb,exb)
 if (exb2<1d-307) then
  gf=-1d0
  return
 end if 
 em2=dot_product(em,em)
 bm2=dot_product(bm,bm)
 eb=dot_product(em,bm)
 
 exb1=sqrt(exb2)
 nf=exb/exb1           ! direction of v - velocity of moving reference frame
 
 ebm2=1/(em2+bm2)
 rho=exb1*ebm2          ! (\vec E x \vec B)/(E^2+B^2)=\vec v/(1+v^2)
 rho2=2*rho
 rs=sqrt(1-rho2**2)
 
 vr=1/(1+rs)
 vf=rho2*vr             ! module of v
 gf=sqrt((1+1/rs)/2)    ! 1/sqrt(1-v**2) - Lorentz factor 
 aa=2*vr*ebm2           ! (1+v^2)/(E^2+B^2)
 aeb=aa*eb

 em1=((1-aa*bm2)*em+aeb*bm)*gf
 bm1=((1-aa*em2)*bm+aeb*em)*gf
 

end subroutine par_frame


subroutine cross_prod(a,b,res)
! Cross product
 real(8), intent(in) :: a(3), b(3)
 real(8), intent(out) :: res(3)
 
 res(1)=a(2)*b(3)-a(3)*b(2) 
 res(2)=a(3)*b(1)-a(1)*b(3)
 res(3)=a(1)*b(2)-a(2)*b(1)
end subroutine cross_prod

end module eb_gam_abs


module eb_gam_test
!==================================================================
! Test of eb_gam_abs
!==================================================================
 use eb_gam_abs, only : ebgam_abs
 use eb_par_abs, only : landau_limit
 use phys_const, only : pi, c_light
 implicit none
 private

 public :: main_calc
 
 
 contains


subroutine main_calc
 
 call test2

end subroutine main_calc



subroutine test
 real(8) :: em(3), bm(3), np(3), ep
 real(8) :: ef, bf, res, res1, res2
 real(8) :: x1, x2
 integer :: i, nx
 
 bf=1d7
 ef=bf*0d0
 
 call unit_vec(0d0,0d0,bm)
 bm=bf*bm
 call unit_vec(0d0,0d0,em)
 em=ef*em 
 
 call unit_vec(0d0,90d0,np)
 
 
 x1=1d5
 x2=1d11
 nx=1000
 
 open(1,file='lam.dat')
 do i=0,nx
  ep=x1*(x2/x1)**(1d0*i/nx)
  call ebgam_abs(em,bm,np,ep,res)
  call landau_limit(bf,ep,res1,res2)
  write(1,*) ep*0.511d6, res, res1, res2, abs(res-res1)/res, abs(res-res2)/res
 end do 
 close(1)
 

end subroutine test


subroutine test2
 real(8) :: em(3), bm(3), np(3), ep
 real(8) :: np1(3), np2(3)
 real(8) :: em0(3), em1(3), em2(3)
 real(8) :: ef, bf, res, res1, res2
 real(8) :: x1, x2
 integer :: i, nx
 
 bf=1d11
 ef=bf
 
 call unit_vec(0d0,0d0,bm)
 bm=bf*bm
 call unit_vec(0d0,89.999d0,em)
 em0=0d0*bf*em
 em1=1d-1*bf*em
 em2=1d2*bf*em
 
 call unit_vec(0d0,10d0,np)
 call unit_vec(0d0,60d0,np1)
 call unit_vec(0d0,10d0,np2)
 
 
 x1=1d2
 x2=1d11
 nx=1000
 
 open(1,file='lam1.dat')
 do i=0,nx
  ep=x1*(x2/x1)**(1d0*i/nx)
  call ebgam_abs(em0,bm,np,ep,res)
  call ebgam_abs(em1,bm,np,ep,res1)
  call ebgam_abs(em2,bm,np,ep,res2)
  res=c_light/res
  res1=c_light/res1
  res2=c_light/res2
  write(1,*) ep*0.511d6, res, res1, res2
 end do 
 close(1)
 

end subroutine test2


! subroutine test1
 ! real(8) :: em(3), bm(3), em1(3), bm1(3)
 ! real(8) :: nf(3), vf, gf
 ! real(8) :: np(3), np1(3), dopf
 ! real(8) :: nem1(3), nbm1(3)
 ! real(8) :: eta
 
 ! eta=0.9d0
 
 ! call unit_vec(0d0,0d0,bm)
 ! call unit_vec(0d0,90d0,em)
 ! em=eta*em 
 
 ! call unit_vec(10d0,59d0,np)
 
 ! call par_frame(em,bm,em1,bm1,nf,vf,gf)
 ! call n1_dir(np,nf,vf,gf,np1,dopf)
 
 ! write(*,*) 'em1=', em1
 ! nem1= em1/sqrt(dot_product(em1,em1))
 ! write(*,*) 'nem1=', nem1, sqrt(dot_product(em1,em1))
 ! write(*,*) 'bm1=', bm1
 ! nbm1=bm1/sqrt(dot_product(bm1,bm1))
 ! write(*,*) 'nbm1=', nbm1, sqrt(dot_product(bm1,bm1))
 
 ! write(*,*) 'dot_product(bm,np)=', dot_product(bm,np)
 ! write(*,*) 'dot_product(nbm1,np1)=', dot_product(nbm1,np1)
 ! write(*,*) 'vf, gf, dopf=', vf, gf, dopf
 

! end subroutine test1


subroutine unit_vec(phi,theta,nv)
 real(8), intent(in) :: phi, theta
 real(8), intent(out) :: nv(3)
 real(8) :: phi1, theta1
 
 phi1=phi*pi/180
 theta1=theta*pi/180
  
 nv(1)=cos(phi1)*sin(theta1)
 nv(2)=sin(phi1)*sin(theta1)
 nv(3)=cos(theta1)
end subroutine unit_vec


end module eb_gam_test



! program main
 ! use eb_gam_test, only : main_calc
 
 ! call main_calc

! end program main