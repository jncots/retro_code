module pgam_photomeson
!===========================================================================
! Calculation of the production rate of secondaries according to
! Kelner, Aharonian (2008PhRvD..78c4013K).
!
! Usage:
! ======
! Initialisation:
! call init_pgam_meson(type), where 'type' is character(5) can take:
!  'gamma' - gamma
!  'posit' - positrons
!  'elect' - electrons
!  'numua' - muonic antineutrinos
!  'numuo' - muonic neutrinos
!  'nuele' - electron neutrinos
!  'nuela' - electron antineutrinos
! Calculation:
! call pgam_meson(xi,x,res),
! where xi=4*eps*Ep/(mc^2)^2, x=E/Ep is non-dimensional parameters
! eps - energy of soft photon, Ep - energy of proton, E - energy of
! secondary particle
! res is production rate in [cm^3/sec]
! Further usage:
! Let res=Phi(eta,x)
! Then integration: int f_p(E_p)*f_ph(e_ph)*Phi(eta,x)*dE_p/E_p*de_ph=dN/dE
! see Kelner, Aharonian (2008PhRvD..78c4013K) for details
!===========================================================================
 use intpol_mod, only : arr_ind_breal
 implicit none
 private
 save
 
 public :: init_pgam_meson, pgam_meson
 
 
 abstract interface
  subroutine int_rate(xi,x,res)
   real, intent(in) :: xi, x
   real, intent(out) :: res
  end subroutine int_rate
 end interface
 
 
 type pgam_int_data
   character(100) :: name ! table for 'name'
   integer :: n
   real, allocatable :: eta(:), s(:), d(:), b(:)
 end type pgam_int_data
 
 integer :: nsec, iread=0
 type(pgam_int_data), allocatable, target :: pgam_dat(:)
 type(pgam_int_data), pointer :: cdat
 procedure(int_rate), pointer :: pg_rate=>null()
 
 contains

subroutine pgam_read
! Variables
 integer :: dev, i
! Calculations 
 call get_unit(dev)
 open(dev,file='/home/anton/work/library/physlib/scr/pgam_photomeson.dat',form='unformatted')
 read(dev) nsec
 allocate(pgam_dat(nsec))
 do i=1,nsec
  read(dev) pgam_dat(i)%name
  read(dev) pgam_dat(i)%n
  allocate(pgam_dat(i)%eta(pgam_dat(i)%n))
  allocate(pgam_dat(i)%s(pgam_dat(i)%n))
  allocate(pgam_dat(i)%d(pgam_dat(i)%n))
  allocate(pgam_dat(i)%b(pgam_dat(i)%n))
  read(dev) pgam_dat(i)%eta
  read(dev) pgam_dat(i)%s
  read(dev) pgam_dat(i)%d
  read(dev) pgam_dat(i)%b
 end do 
 close(dev)
end subroutine pgam_read
 
 
 
 
subroutine init_pgam_meson(sec_type)
! Variables 
 character(5), intent(in) :: sec_type
! Calculations 
 if (iread==0) then
  call pgam_read
  iread=1
 end if
 
 select case (sec_type)
  case('gamma') ! gamma
   pg_rate=>gam
  case('posit') ! positrons
   pg_rate=>posit
  case('elect') ! electrons
   pg_rate=>elect
  case('numua') ! muonic antineutrinos
   pg_rate=>numua
  case('numuo') ! muonic neutrinos
   pg_rate=>numu
  case('nuele') ! electron neutrinos
   pg_rate=>nuel
  case('nuela') ! electron antineutrinos
   pg_rate=>nuela
  case default 
   write(*,*) 'pgam_meson_init: Incorrect type of secondary particles.'
   write(*,*) 'Choose one of following:'
   write(*,*) 'gamma, posit, elect, numuo, numua, nuele, nuela'
 end select
end subroutine init_pgam_meson


subroutine pgam_meson(xi,x,res)
!
! Conversion to real(8)
!
! Variables
 real(8), intent(in) :: xi, x
 real(8), intent(out) :: res
 real :: x1, xi1, res1
 real, parameter :: eta0=0.31344
! Calculations
 if (iread==0) then
  write(*,*) 'init_pgam_meson is not initialised'
  return
 end if

 x1=real(x)
 xi1=real(xi)/eta0 ! in the functions below xi1 is normalized to eta0
 call pg_rate(xi1,x1,res1)
 res=real(res1,8)
end subroutine pgam_meson



subroutine gam(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=1e0, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.27078, eta2=2.1332e-2
 real, parameter :: s0=7.69e-2, d0=5.44e-1, b0=6.57e-18, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 
 cdat=>pgam_dat(1)
 
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.
  return
 end if
 
 if (xi>xi_max) xi=xi_max
 
 eta=eta0*xi
 p0=sqrt((eta-eta0)*(eta+eta1))
 p0=(eta+eta2+p0)/2
 x1=eta2/p0
 x2=p0/(1+eta)
 
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 alpha=2.5+0.4*log(xi)
 if (xi<=cdat%eta(1)) then
  s=s0
  d=d0
  b=b0*(xi-1)
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 
 
 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine gam





subroutine posit(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=1e0, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.27078, eta2=2.1332e-2
 real, parameter :: s0=3.67e-1, d0=3.12e0, b0=1.86e-17, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 

 cdat=>pgam_dat(2)
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.
  return
 end if
 
 if (xi>xi_max) xi=xi_max
 
 eta=eta0*xi
 p0=sqrt((eta-eta0)*(eta+eta1))
 p0=(eta+eta2+p0)/2
 x1=0.25*eta2/p0
 x2=p0/(1+eta)
 
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 alpha=2.5+1.4*log(xi)
 if (xi<=cdat%eta(1)) then
  s=s0
  d=d0
  b=b0*(xi-1)
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 

 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine posit




subroutine elect(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=2.1361, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.6695, eta15=0.29211, eta2=2.1332e-2
 real, parameter :: s0=6.58E-1, d0=3.09E0, b0=1.983E-18, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 

 cdat=>pgam_dat(3)
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.
  return
 end if
 
 if (xi>xi_max) xi=xi_max
  
 eta=eta0*xi
 p0=sqrt((eta-eta1)*eta)
 p0=(eta-eta15+p0)/2
 x1=0.5*eta2/p0
 x2=p0/(1+eta)
 
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 
 if (xi>4.) then
  alpha=6*(1-exp(1.5*(4-xi)))
 else
  alpha=0.0
 end if
 
 if (xi<=cdat%eta(1)) then
  s=s0
  d=d0
  b=b0*(xi-2.1361)**2
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 
 
 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine elect



subroutine numua(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=1e0, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.27078, eta2=2.1332e-2
 real, parameter :: s0=3.65e-1, d0=3.09e0, b0=1.86e-17, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 
 cdat=>pgam_dat(4)
 
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.
  return
 end if
 
 if (xi>xi_max) xi=xi_max
 
 eta=eta0*xi
 p0=sqrt((eta-eta0)*(eta+eta1))
 p0=(eta+eta2+p0)/2
 x1=0.25*eta2/p0
 x2=p0/(1+eta)
 
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 alpha=2.5+1.4*log(xi)
 if (xi<=cdat%eta(1)) then
  s=s0
  d=d0
  b=b0*(xi-1)
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 
 
 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine numua




subroutine numu(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=1e0, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.27078, eta2=2.1332e-2
 real, parameter :: s0=0.0, d0=0.0, b0=2.48e-17, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 
 cdat=>pgam_dat(5)
 
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.0
  return
 end if
 
 if (xi>xi_max) xi=xi_max
 
 eta=eta0*xi
 p0=sqrt((eta-eta0)*(eta+eta1))
 p0=(eta+eta2+p0)/2
 x1=0.427*eta2/p0
 x2=p0/(1+eta)
 
 if (xi<1e1) then
  if (xi>2.14) then
   x2=(0.427+0.0729*(xi-2.14))*x2
  else
   x2=0.427*x2
  end if 
 end if 
 
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 alpha=2.5+1.4*log(xi)
 if (xi<=cdat%eta(1)) then
  s=s0
  d=d0
  b=b0*(xi-1)
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 
 
 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine numu


subroutine nuel(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=1e0, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.27078, eta2=2.1332e-2
 real, parameter :: s0=3.67e-1, d0=3.12e0, b0=2.17e-17, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 
 cdat=>pgam_dat(6)
 
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.0
  return
 end if
 
 if (xi>xi_max) xi=xi_max
 
 eta=eta0*xi
 p0=sqrt((eta-eta0)*(eta+eta1))
 p0=(eta+eta2+p0)/2
 x1=0.25*eta2/p0
 x2=p0/(1+eta)
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 alpha=2.5+1.4*log(xi)
 if (xi<=cdat%eta(1)) then
  s=s0
  d=d0
  b=b0*(xi-1)
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 
 
 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine nuel



subroutine nuela(xi0,x,res)
! Variables
 real, intent(in) :: xi0, x
 real, intent(out) :: res
 real, parameter :: xi_min=2.1361, xi_max=1e2
 real, parameter :: eta0=0.31344, eta1=0.6695, eta15=0.29211, eta2=2.1332e-2
 real, parameter :: s0=6.58e-1, d0=3.09e0, b0=1.983e-18, t0=0.69315
 integer :: n1, n2
 real :: xi, eta, p0, x1, x2, alpha, s, d, b, t1, t2, xint
! Calculations 

 cdat=>pgam_dat(7)
 xi=xi0
 
 if(xi<=xi_min) then
  res=0.
  return
 end if
 
 if (xi>xi_max) xi=xi_max
  
 eta=eta0*xi
 p0=sqrt((eta-eta1)*eta)
 p0=(eta-eta15+p0)/2
 x1=0.5*eta2/p0
 x2=p0/(1+eta)
 
 
 if(x>=x2) then
  res=0.0
  return
 end if
 
 
 if (xi>4.) then
  alpha=6*(1-exp(1.5*(4-xi)))
 else
  alpha=0.0
 end if
 
 if (xi<=cdat%eta(1)) then 
  s=s0
  d=d0
  b=b0*(xi-2.1361)**2
 else
  call arr_ind_breal(cdat%n,cdat%eta,xi,n1,n2)
  xint=(xi-cdat%eta(n1))/(cdat%eta(n2)-cdat%eta(n1))
  s=cdat%s(n1)+(cdat%s(n2)-cdat%s(n1))*xint
  d=cdat%d(n1)+(cdat%d(n2)-cdat%d(n1))*xint
  b=cdat%b(n1)+(cdat%b(n2)-cdat%b(n1))*xint
 end if 
 
 if (x>x1) then
  t1=(x-x1)/(x2-x1)
  t1=log(2/(1+t1**2))**alpha
  t2=s*log(x/x1)**d
  t2=exp(-t2) 
  res=b*t1*t2
 else
  res=b*t0**alpha	
 end if
 
end subroutine nuela



subroutine get_unit(unit)
!====================================
! Search for free 'unit' from (1..99)
!====================================
! Variables
 integer, intent(out) :: unit
 integer :: i
! Calculations
 unit=0
 do i=1,99
  if ((i/=5).and.(i/=6).and.(i/=9)) then
   open(unit=i,err=10,status='scratch')
   close(unit=i)
   unit=i
   return
  end if
10  continue
 end do
end subroutine get_unit



end module pgam_photomeson


!program main
! use pgam_int_read, only : init_pgam_meson, pgam_meson
! real(8) :: xi, x, res
! 
! call init_pgam_meson('gamma')
! xi=2.35d0
! x=0.111d0
! call pgam_meson(xi,x,res)
! write(*,*) res
!
!end program main

module pgam_phph
!=========================================================
! Convolution of production rate of secondaries from 
! p-gamma photomeson interation with soft photon field,
! and tabulation of the result for further convolution
! with proton spectrum
!=========================================================
 use pgam_photomeson, only : init_pgam_meson, pgam_meson
 use gauss_kronrod, only : gk_adaptive
 use description_line, only : descr_line, infoline
 use word_processor, only : line_break
 use timer_module, only : timer_class
 implicit none
 private
 save

 
 abstract interface
  subroutine fun1(x,res)
   real(8), intent(in) ::  x
   real(8), intent(out) :: res
  end subroutine fun1
 end interface 
 
 
 
 public :: tab_phph
 
 real(8), parameter :: mp=0.938272046d9 ! eV proton mass
 real(8), parameter :: mp2=mp**2/4
 real(8), parameter :: rpp=0.1460540999 ! m_pion/m_proton \approx 0.146
 real(8), parameter :: rpp2=rpp**2
 real(8), parameter :: etap=rpp2-2*rpp, etam=rpp2+2*rpp
 real(8), parameter :: eta0=0.31344d0
 
 
 real(8) :: xi0, xep, eph1, eph2
 procedure(fun1), pointer :: soft_ph=>null()
 real(8), allocatable :: gam_en(:), pr_en(:,:), tab_phmes(:,:)
 
 integer :: ndes
 character(len=:), allocatable :: description
 
 
 contains

 
 
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
 line='This file contains the rate of gamma-ray production in photomeson interaction &
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


end subroutine info 
 
 
 
subroutine tab_phph(ptype,nep,neg,epmin,epmax,x1,x2,ph_den,fname,inf)
 integer, intent(in) :: nep, neg ! discretisation number for proton and gamma energies
 real(8), intent(in) :: epmin, epmax ! min and max energy of protons
 real(8), intent(in) :: x1, x2 ! low and high limits of soft photon spectrum
 procedure(fun1) :: ph_den
 type(infoline) :: inf(:)
 character(500) :: fname
 character(5), intent(in) :: ptype
 integer :: i, j
 real(8) :: ep, ep1, ep2, egam, egam1, egam2, res
 real(8) :: xmm, xpp, eta
 type(timer_class) :: tmr

! Calculations
 call info(inf)
 eph1=x1
 eph2=x2
 soft_ph=>ph_den
 call init_pgam_meson(ptype)
 
 allocate(gam_en(neg+1))
 allocate(pr_en(nep+1,neg+1))
 allocate(tab_phmes(nep+1,neg+1))
 
 
 
 
 egam1=(1d0+1d-5)*rpp2*mp2/eph2
 egam2=epmax
 
 call tmr%start((nep+1)*(neg+1),1.0)
 do i=0,neg
  egam=egam1*(egam2/egam1)**(i*1d0/neg)
  gam_en(i+1)=egam
  
  eta=egam*eph2/mp2-rpp2
  ep1=(1d0+1d-3)*egam*(1+sqrt(1+4/eta))/2
  if (ep1<epmin) ep1=epmin
  ep2=3*epmax
  if (ep1>=ep2) write(*,*) 'ep1>=ep2'
  
  do j=0,nep
   ep=ep1*(ep2/ep1)**(j*1d0/nep)
   pr_en(j+1,i+1)=ep
   call phph(ep,egam,res)
   tab_phmes(j+1,i+1)=res 
   call tmr%loop
  end do
 end do
 
 
 open(1,file=fname,form='unformatted')
 write(1) ndes
 write(1) description
 write(1) size(gam_en)
 write(1) gam_en
 write(1) size(pr_en,1), size(pr_en,2)
 write(1) pr_en
 write(1) tab_phmes
 close(1)
 
 deallocate(gam_en)
 deallocate(pr_en)
 deallocate(tab_phmes)
 
end subroutine tab_phph
 
subroutine phph(ep,egam,res)
! Variables
 real(8), intent(in) :: ep, egam
 real(8), intent(out) :: res
 real(8) :: iep, tmin, tmax
! Calculations  
 iep=1/ep
 xi0=ep/mp2
 xep=egam*iep
 tmin=(rpp2/xep+xep/(1-xep))*mp2*iep
 if (tmin<eph1) tmin=eph1
 tmax=eph2
 
 if (tmax<=tmin) then
  res=0d0
  return
 end if
 
 tmin=log(tmin)
 tmax=log(tmax)
 call gk_adaptive(iphph,tmin,tmax,res)
 

end subroutine phph
 
 
 
 
subroutine iphph(x,res)
!
! Integrand
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eph, sph, xi, phi
! Calculations 
 eph=exp(x)
 call soft_ph(eph,sph)
 xi=xi0*eph
 call pgam_meson(xi,xep,phi)
 res=sph*phi*eph
end subroutine iphph



end module pgam_phph

!======================================
! Test program
!======================================
!program main
! implicit none
! 
! call tab_param
!
!end program main
!
!
!
!subroutine tab_param
! use pgam_phph, only : tab_phph
! use description_line, only : infoline
! implicit none
!! Variables
! integer :: nep, neg
! real(8) :: epmin, epmax, eph1, eph2
! character(500) :: fname
! character(5) :: ptype
! type(infoline) :: inf(3)
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
!! Information
! inf(1)%typ='sub'
! inf(1)%intro='Parameters of tabulation'
! inf(1)%sub_name='tab_param'
! inf(1)%file_name='pgam_phph.f90'
! 
! inf(2)%typ='sub'
! inf(2)%intro='Soft photon field'
! inf(2)%sub_name='nden'
! inf(2)%file_name='pgam_phph.f90'
! 
! inf(3)%typ='sub'
! inf(3)%intro='Setup of the calculations'
! inf(3)%sub_name='tab_phph'
! inf(3)%file_name='pgam_phph.f90'
!
! nep=600
! neg=600
! epmin=1d11
! epmax=1d18
! eph1=1d-10
! eph2=1d3
! fname='pgam_phph.dat'
! ptype='gamma'
! call tab_phph(ptype,nep,neg,epmin,epmax,eph1,eph2,nden,fname,inf)
!end subroutine tab_param
!
!
!subroutine nden(x,res)
!!
!! Soft photon field
!! x energy in eV
!! res number density in 1/(eV cm^3)
!!
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

module rate_mes_gam
!===========================================================
! Calculation of the production rate of secondary  
! particles (gamma, electron, positron, ...) produced in
! photomeson interactions of protons with soft photon field.
! This module convolves the rate with specific proton
! spectrum.
! 
! The information about type of particle (gamma,
! electron, ...) is in the data file which contains the
! tabulated function for the production rate of the specific
! type of particles convolved with soft photon field.
! (see information in the heading of the file).
!
! Module connection:
! ==================
! use rate_mes_gam, only : read_tab, set_prod_rate, prod_rate
!
! Reading of the data file:
! =========================
! call read_data(fname)
! where fname is character(500) is a name (or full path/name)
! of the file.
!
! Setting the proton spectrum:
! ============================
! call set_prod_rate(epmin,epmax,pr_sp)
! where epmin, epmax is minimum and maximum energy of proton
! spectrum, pr_sp(x,res) is proton spectrum function
!
! Production rate of the secondary particles:
! ===========================================
! call prod_rate(E,res) 
! E is energy of the particle in [eV]
! res is production rate in [1/(eV sec cm^3)]
!
! Example of usage is below
!===========================================================
 use intpol_mod, only : arr_ind_short, log_int, lin_int
 use tab2fun, only : tabfun
 use gauss_kronrod, only : gk_adaptive
  
 implicit none
 private
 save
 
 public :: set_prod_rate,  prod_rate, read_tab
 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 integer :: np, ng
 real(8), allocatable :: gam_en(:), pr_en(:,:), tab_phmes(:,:)
 real(8), allocatable :: g_min_max(:,:)
 real(8) :: ep1, ep2
 procedure (fun_prot), pointer :: pr_spec=>null(), phi_fun=>null()
 
 contains

 
 
subroutine read_tab(fname)
! Variables
 character(500), intent(in) :: fname
 integer, parameter :: ndev=25
 integer :: ndes
 character(len=:), allocatable :: description
! Calculations
 
 if (allocated(gam_en)) deallocate(gam_en)
 if (allocated(pr_en)) deallocate(pr_en)
 if (allocated(tab_phmes)) deallocate(tab_phmes)
 
 open(ndev,file=fname,form='unformatted')
 read(ndev) ndes
 allocate(character(len=ndes) :: description)
 read(ndev) description
 read(ndev) ng
 allocate(gam_en(ng))
 read(ndev) gam_en
 read(ndev) np, ng
 allocate(pr_en(np,ng))
 allocate(tab_phmes(np,ng))
 read(ndev) pr_en
 read(ndev) tab_phmes
 close(ndev)
end subroutine read_tab



subroutine set_prod_rate(x1,x2,prsp)
! Variables 
 real(8), intent(in) :: x1, x2
 procedure (fun_prot) :: prsp
! Calculations 
 ep1=x1 ! min and max energy of proton spectrum
 ep2=x2
 pr_spec=>prsp
end subroutine 




subroutine prod_rate(x,res)
! Variables 
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 integer :: ng1, ng2
 real(8) :: res1, res2
 type(tabfun) :: phi1, phi2
! Calculations
 if ((x<gam_en(1)).or.(x>gam_en(ng))) then
  res=0d0
  return
 end if
 
 call arr_ind_short(1,ng,gam_en,x,ng1,ng2)
 call phi1%on(pr_en(:,ng1),tab_phmes(:,ng1))
 call conv_fp(phi1,res1)
 call phi2%on(pr_en(:,ng2),tab_phmes(:,ng2))
 call conv_fp(phi2,res2)
 call log_int(gam_en(ng1),gam_en(ng2),res1,res2,x,res)
 if (isnan(res)) call lin_int(gam_en(ng1),gam_en(ng2),res1,res2,x,res)
 call phi1%off
 call phi2%off
end subroutine prod_rate


subroutine conv_fp(phic,res)
! Variables
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
 if (tmin<ep1) tmin=ep1
 if (tmax>ep2) tmax=ep2
 
 if (tmin>tmax) then
  res=0d0
  return
 end if
 
 phi_fun=>phic%fx
 tmin=log(tmin)
 tmax=log(tmax)
 call gk_adaptive(iconv_fp,tmin,tmax,res)
 phi_fun=>null()

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
 res=fp*phi
end subroutine iconv_fp



end module rate_mes_gam

!======================================
! Test program
!======================================
!program main
! use rate_mes_gam, only : set_prod_rate,  prod_rate, read_tab
! implicit none
! character(500) :: fname
! integer :: i, ng
! real(8) :: ep1, ep2, egam, res
! real(8) :: egam1, egam2
! 
! interface
!  subroutine prsp(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine prsp
! end interface
! 
! 
! 
! fname='pgam_phph.dat'
! ep1=1d12
! ep2=1d18
! call read_tab(fname)
! call set_prod_rate(ep1,ep2,prsp)
! 
! egam1=1d12
! egam2=1d18
! ng=100
! 
! do i=0,ng
!  egam=egam1*(egam2/egam1)**(1d0*i/ng)
!  call prod_rate(egam,res)
!  write(*,'(10Es14.6)') egam, res
! end do 
! 
! 
!end program main
!
!
!
!subroutine prsp(x,res)
!! Variables
! real(8), intent(in) :: x
! real(8), intent(out) :: res
!! Calculations 
! res=1/x**2
!end subroutine prsp

module pgam_mes_loss
!==================================================================
! Calculation of energy losses in photomeson interactions of 
! protons with soft photon field
! Example of usage:
! =================
! use pgam_mes_loss, only : pgmes_loss
! type(pgmes_loss) :: loss 
! call loss%sph(1d-10,1d3,sph) ! min, max and energy spectrum
! call loss%prot(l1%epmin,1d16,100) ! min, max, number of points
! call loss%solve ! calculation of losses
! do i=1,loss%np
!  write(*,*) loss%ep(i), loss%dedt(i) ! loss%dedt0(i,j) losses for 
! end do
!
! loss%dedt0(i,j) losses for each channel:
! j=[1='gamma',2='posit',3='elect',
! 4='numua',5='numuo',6='nuele',7='nuela']
! see in pgam_photomeson
!==================================================================
 use pgam_photomeson, only : init_pgam_meson, pgam_meson
 use gauss_kronrod, only : gk_adaptive
 use timer_module, only : timer_class
 implicit none
 private
 save
 
 public :: pgmes_loss
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type pgmes_loss
  real(8) :: eph1, eph2
  procedure(fun_1d), pointer, nopass :: sph=>null()
  real(8) :: epmin
  integer :: np, ntype
  real(8) :: ep1, ep2
  real(8), allocatable :: ep(:), dedt0(:,:), dedt(:)
 contains
  procedure :: set_sph, set_prot, tab_loss
 end type pgmes_loss
 
 
 
 
 real(8), parameter :: mp=0.938272046d9 ! eV proton mass
 real(8), parameter :: mp2=mp**2/4
 real(8), parameter :: rpp=0.1460540999 ! m_pion/m_proton \approx 0.146
 real(8), parameter :: rpp2=rpp**2
 real(8), parameter :: eta1=-rpp2-2*rpp, eta2=-rpp2+2*rpp
 real(8), parameter :: eta0=0.31344d0
 
 real(8) :: ep0, xep, xi0, eph10, eph20
 procedure(fun_1d), pointer :: sph0=>null()
 
 
 contains

subroutine set_sph(this,eph1,eph2,sph)
! Variables
 class(pgmes_loss) :: this
 real(8), intent(in) :: eph1, eph2
 procedure(fun_1d) :: sph
! Calculations
 this%eph1=eph1
 this%eph2=eph2
 this%sph=>sph
 this%epmin=eta0*mp2/eph2*(1d0+1d-5)
end subroutine set_sph


subroutine set_prot(this,ep1,ep2,np)
! Variables
 class(pgmes_loss) :: this
 real(8), intent(in) :: ep1, ep2
 integer, intent(in) :: np
 integer :: i
! Calculations
 if (ep1>=ep2) then
  write(*,*) 'type(pgmes_loss): ep1>=ep2'
  this%ep1=this%epmin
  this%ep2=this%epmin*1d5
 else
  this%ep1=ep1
  this%ep2=ep2
 end if
 
 if (allocated(this%ep)) deallocate(this%ep)
 allocate(this%ep(np+1))
 do i=0,np
  this%ep(i+1)=this%ep1*(this%ep2/this%ep1)**(i*1d0/np)
 end do
 this%np=np+1
end subroutine set_prot

 

subroutine tab_loss(this)
! Variables
 class(pgmes_loss) :: this
 integer :: i, j
 real(8) :: ep, res
 type(timer_class) :: tmr
 character(5), parameter :: ptype(7)=&
 ['gamma','posit','elect','numua','numuo','nuele','nuela']
! Calculations  
 this%ntype=7
 if (allocated(this%dedt0)) deallocate(this%dedt0)
 if (allocated(this%dedt)) deallocate(this%dedt)
 allocate(this%dedt0(this%np,this%ntype))
 allocate(this%dedt(this%np))
 
 eph10=this%eph1
 eph20=this%eph2
 sph0=>this%sph
 
 this%dedt=0d0
 
 call tmr%start(this%ntype*this%np)
 do j=1,this%ntype
  call init_pgam_meson(ptype(j))
  do i=1,this%np
   ep=this%ep(i)
   call loss(ep,res)
   this%dedt0(i,j)=res/ep
   call tmr%loop
  end do
  this%dedt=this%dedt+this%dedt0(:,j)
 end do
  
end subroutine tab_loss




subroutine loss(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eta, xmm, xpp, tmin, tmax
! Calculations 
 ep0=x
 eta=eph20*ep0/mp2
 if (eta<=eta0) then
  res=0d0
  return
 end if 
 xpp=(sqrt((eta+eta1)*(eta+eta2))+eta+rpp2)/2
 xmm=rpp2/xpp
 xpp=xpp/(1+eta)
 
 tmin=log(ep0*xmm)
 tmax=log(ep0*xpp)
 call gk_adaptive(phph,tmin,tmax,res)
end subroutine loss




subroutine phph(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: iep, tmin, tmax, egam
! Calculations  
 egam=exp(x)
 iep=1/ep0
 xi0=ep0/mp2
 xep=egam*iep
 tmin=(rpp2/xep+xep/(1-xep))*mp2*iep
 if (tmin<eph10) tmin=eph10
 tmax=eph20
 
 if (tmax<=tmin) then
  res=0d0
  return
 end if
 
 tmin=log(tmin)
 tmax=log(tmax)
 call gk_adaptive(iphph,tmin,tmax,res)
 res=res*egam**2
end subroutine phph 
 

subroutine iphph(x,res)
!
! Integrand
!
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eph, sph, xi, phi
! Calculations 
 eph=exp(x)
 call sph0(eph,sph)
 xi=xi0*eph
 call pgam_meson(xi,xep,phi)
 res=sph*phi*eph
end subroutine iphph
 

end module pgam_mes_loss



!==========================================
! Example of usage
!==========================================
!program main
! use pgam_mes_loss, only : pgmes_loss
! use thermal_photon_field, only : therm_phf, set_thnden, thnden
! implicit none
!! Variables
! integer :: i
! type(pgmes_loss) :: l1
! type(therm_phf) :: sph
!! Calculations 
! 
! call sph%init(4)
! sph%temp=[3d0,3d-1,6d-3,2.35d-4]
! sph%uden=[5d5,1d6,4d4,0.26d0]
! call set_thnden(sph)
! 
! call l1%set_sph(1d-10,1d3,thnden)
! call l1%set_prot(l1%epmin,1d22,200)
! call l1%tab_loss
! 
! do i=1,l1%np
!  write(*,*) l1%ep(i), l1%dedt(i)
! end do
! 
!end program main

module phmes_ratem
!===========================================================
! Calculation of the production rate of secondary  
! particles (gamma, electron, positron, ...) produced in
! photomeson interactions of protons with soft photon field.
! This module convolves the rate with specific proton
! spectrum.
! 
! The information about type of particle (gamma,
! electron, ...) is in the data file which contains the
! tabulated function for the production rate of the specific
! type of particles convolved with soft photon field.
! (see information in the heading of the file).
!
! Module connection:
! ==================
! use phmes_ratem, only : phmes_rate
! type(phmes_rate) :: f
!
! Reading of the data file:
! =========================
! f%fname='data.dat'
! call f%read
! where 'data.dat' is character(500) is a name (or full path/name)
! of the file.
!
! Setting the proton spectrum:
! ============================
! call f%prot_sp(pr_sp,epmin,epmax)
! where epmin, epmax is minimum and maximum energy of proton
! spectrum, pr_sp(x,res) is proton spectrum function
!
! Production rate of the secondary particles:
! ===========================================
! call f%rate(E,res) 
! E is energy of the particle in [eV]
! res is production rate in [1/(eV sec cm^3)]
!
! Example of usage is below
!===========================================================
 use intpol_mod, only : arr_ind_short, logl_int
 use tab2fun, only : tabfun
 use gauss_kronrod, only : gk_adaptive
 use utools_mod, only : utools
  
 implicit none
 private
 save
 
 public :: phmes_rate
 
 abstract interface
  subroutine fun1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun1d
 end interface
 
 type flim
  real(8) :: x1, x2
  procedure (fun1d), pointer, nopass :: fx=>null()
 end type flim
 
 
 type phmes_rate
  logical :: istat=.true.  
  character(500) :: fname
  integer :: ng, np ! eg(ng), ep(np,ng) energy of secondary particle and proton
  real(8), allocatable :: eg(:), ep(:,:), f(:,:) ! f(np,ng) is result of tabulation
  type(flim) :: psp ! proton spectrum
 contains
  procedure, private :: conv
  procedure :: prot_sp, rate, read=>read_tab
 end type phmes_rate
 

 procedure(fun1d), pointer :: pr_spec=>null(), phi_fun=>null() 
 
 
 contains

subroutine prot_sp(this,fx,x1,x2)
 class(phmes_rate) :: this 
 procedure(fun1d) :: fx
 real(8), intent(in) :: x1, x2

 this%psp%fx=>fx
 this%psp%x1=x1
 this%psp%x2=x2
 
end subroutine prot_sp


 
subroutine read_tab(this)
 class(phmes_rate) :: this
 integer :: ndev
 integer :: ndes
 character(len=:), allocatable :: description
 type(utools) :: ut
 
 if (allocated(this%eg)) deallocate(this%eg)
 if (allocated(this%ep)) deallocate(this%ep)
 if (allocated(this%f)) deallocate(this%f)
 
 call ut%get_unit(ndev)
 open(ndev,file=this%fname,form='unformatted')
 read(ndev) ndes
 allocate(character(len=ndes) :: description)
 read(ndev) description
 
 read(ndev) this%ng
 allocate(this%eg(this%ng))
 read(ndev) this%eg
 read(ndev) this%np, this%ng
 allocate(this%ep(this%np,this%ng))
 allocate(this%f(this%np,this%ng))
 read(ndev) this%ep
 read(ndev) this%f
 close(ndev)
 
 this%istat=.false.
end subroutine read_tab



subroutine rate(this,x,res)
 class(phmes_rate) :: this
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 integer :: ng1, ng2
 real(8) :: res1, res2
 
 
 if (this%istat) then
  write(*,*) 'phmes_rate: the data has not been read'
  res=0d0
  return
 end if
 
 if ((x<this%eg(1)).or.(x>this%eg(this%ng))) then
  res=0d0
  return
 end if
 
 call arr_ind_short(1,this%ng,this%eg,x,ng1,ng2)
 call this%conv(ng1,res1)
 call this%conv(ng2,res2)
 call logl_int(this%eg(ng1),this%eg(ng2),res1,res2,x,res)
end subroutine rate


subroutine conv(this,ng,res)
! Variables
 class(phmes_rate) :: this
 integer, intent(in) :: ng
 type(tabfun) :: phic
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations 
 call phic%on(this%ep(:,ng),this%f(:,ng))
 if (phic%zero==0) then
  res=0d0
  call phic%off
  return
 end if
 
 tmin=phic%xmin
 tmax=phic%xmax
 if (tmin<this%psp%x1) tmin=this%psp%x1
 if (tmax>this%psp%x2) tmax=this%psp%x2
 
 if (tmin>tmax) then
  res=0d0
  call phic%off 
  return
 end if
 
 pr_spec=>this%psp%fx
 phi_fun=>phic%fx
 tmin=log(tmin)
 tmax=log(tmax)
 call gk_adaptive(iconv_fp,tmin,tmax,res)
 pr_spec=>null()
 phi_fun=>null()
 call phic%off
end subroutine conv



subroutine iconv_fp(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: ep, fp, phi
! Calculations
 ep=exp(x)
 call pr_spec(ep,fp)
 call phi_fun(ep,phi)
 res=fp*phi
end subroutine iconv_fp



end module phmes_ratem

! !======================================
! ! Test program
! !======================================
! program main
 ! use phmes_ratem, only : phmes_rate
 ! use utools_mod, only : utools 
 ! implicit none
 
 ! type(phmes_rate) :: f
 ! type(utools) :: ut
 ! real(8), allocatable :: eg(:)
 ! real(8) :: ee, res
 ! integer :: i
 
 ! interface
  ! subroutine prsp(x,res)
   ! real(8), intent(in) :: x
   ! real(8), intent(out) :: res
  ! end subroutine prsp
 ! end interface
 
 
 ! f%fname='D:\project\Diffusion\SphericalDiffusion\Estimations&
 ! \xray_gc\gc_rad\xprof_dist\pgam_data\phmes\\tpmes_gam_r20pc.dat'
  
 ! write(*,*) trim(f%fname)
 ! read(*,*)
 
 ! call f%read
 ! call f%prot_sp(prsp,1d12,1d20)
 ! call ut%grid(eg,1d12,1d20,200)
 
 ! open(1,file='test_distg.dat') 
 ! do i=1,size(eg)
  ! ee=eg(i)
  
  ! call f%rate(ee,res)
  ! write(1,*) ee, res*(ee**2)
 ! end do 
 ! close(1)

! end program main



! subroutine prsp(x,res)
! ! Variables
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! ! Calculations 
! res=1/x**2
! end subroutine prsp
