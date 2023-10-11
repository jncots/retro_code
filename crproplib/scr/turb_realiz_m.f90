module turb_realiz_m
 use utools_mod, only : utools
 use vector_rot3d_m, only : vector_rot3d
 use random_gen_m, only : init_random_seed
 use phys_const, only : c_light, pi, psec

 implicit none
 private
 save

 public :: turb_realiz

 real(8), parameter :: raddeg=180/pi
 real(8), parameter :: pi2=2*pi


 
 type turb_realiz
  integer :: nmodes
  real(8) :: bms, lmin, lmax, lcor
  real(8), allocatable :: kspace(:), amp(:)
  real(8), allocatable :: cbm_n(:,:), sbm_n(:,:), kvec(:,:)
  integer, allocatable :: seed(:)
  logical :: manual_seed=.false.
 contains
  procedure :: set=>set_realiz, calc=>calc_bm
 end type turb_realiz



 contains


subroutine set_realiz(this,nmodes,lmin,lmax,bms)
!=================================================
! Precalculations of modes
!=================================================
 class(turb_realiz) :: this
 integer, intent(in) :: nmodes
 real(8), intent(in) :: lmin, lmax, bms
 integer :: i

 this%nmodes=nmodes
 this%lmin=lmin
 this%lmax=lmax
 this%bms=bms


 call set_kspace(nmodes,lmin,lmax,this%lcor,this%kspace,this%amp)

! To generate the same realization
 if (this%manual_seed) then
  
  if (allocated(this%seed)) then
!   call random_seed(put=this%seed)
    call init_random_seed(this%seed)
  else
   call init_random_seed(this%seed) ! new random realization
  end if
 
 else
  call init_random_seed(this%seed)  ! new random realization
 end if


 call generate_bmode(nmodes,this%cbm_n,this%sbm_n,this%kvec)
 
 do i=1,nmodes
  this%cbm_n(:,i)=this%cbm_n(:,i)*this%amp(i)*bms
  this%sbm_n(:,i)=this%sbm_n(:,i)*this%amp(i)*bms
  this%kvec(:,i)=this%kvec(:,i)*this%kspace(i)
 end do


end subroutine set_realiz


subroutine calc_bm(this,x,bm)
!===========================================================
! Calculation of turbulent magnetic field in the point with
! vector x
!===========================================================
 class(turb_realiz) :: this
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: bm(3)
 real(8) :: kphase, bkmode(3) 
 integer :: i


 bm=0d0
 do i=1,this%nmodes
  kphase=dot_product(this%kvec(:,i),x)
  bkmode=cos(kphase)*this%cbm_n(:,i)+sin(kphase)*this%sbm_n(:,i)
  bm=bm+bkmode
 end do
 
end subroutine calc_bm



subroutine set_kspace(nk,lmin,lmax,lcor,kspace,amp)
!===========================================================
! Setting of turbulent spectrum (Kolmogorov type)
!
! Input:
! nk - number of modes
! lmin - minimal scale of turbulence
! lmax - maximal scale of turbulence
!
! Output:
! lcor - correlation length
! kspace - module of wave-vector
! amp - amplitude of wave
!===========================================================
 integer, intent(in) :: nk
 real(8), intent(in) :: lmin, lmax
 real(8), intent(out) :: lcor
 real(8), allocatable, intent(out) :: kspace(:), amp(:)
 real(8) :: kmin, kmax
 real(8) :: tal, ki, dlk, gk_sum
 integer :: i
 type(utools) :: ut
 real(8), allocatable :: gk(:)


 tal=5d0/3d0 ! Kolmogorov
 call corr_length(tal,lmin,lmax,lcor)

 
 kmin=2*pi/lmax
 kmax=2*pi/lmin

 call ut%grid(kspace,kmin,kmax,nk)
 allocate(gk(nk))
 dlk=4*pi*log(kspace(2)/kspace(1))

 gk_sum=0d0
 do i=1,nk
  ki=kspace(i)
  gk(i)=dlk*ki**3/(1+(ki*lcor)**(11/3d0))
  gk_sum=gk_sum+gk(i)
 end do


 allocate(amp(nk))
 do i=1,nk-1
  ki=kspace(i)
  amp(i)=sqrt(gk(i)/gk_sum)
 end do

 deallocate(gk)

end subroutine set_kspace


subroutine corr_length(tal,lmin,lmax,lcor)
!===========================================================
! Calculation of correlation length for turbulent spectrum
! with min and max turbulent scales lmin and lmax and
! turbulent spectrum B^2(k)=1/k^tal
!===========================================================
 real(8), intent(in) :: tal, lmin, lmax
 real(8), intent(out) :: lcor
 real(8) :: lmm

 lmm=lmin/lmax
 lcor=(1-lmm**tal)/(1-lmm**(tal-1))
 lcor=lcor*lmax*(tal-1)/(2*tal)


end subroutine corr_length


subroutine generate_bmode(nk,cbm_n,sbm_n,kvec)
!===========================================================
! Generation of nk modes of turbulent magnetic field
!===========================================================
 integer, intent(in) :: nk
 real(8), allocatable, intent(out) :: cbm_n(:,:), sbm_n(:,:), kvec(:,:)
 real(8), allocatable :: kphi(:), ktheta(:), kalpha(:), kbeta(:)
 integer :: i


 allocate(kphi(nk))
 allocate(ktheta(nk))
 allocate(kalpha(nk))
 allocate(kbeta(nk))
 allocate(cbm_n(3,nk))
 allocate(sbm_n(3,nk))
 allocate(kvec(3,nk))
 

 call random_number(kphi)
 call random_number(ktheta)
 call random_number(kalpha)
 call random_number(kbeta)
 
 
 kphi=pi2*kphi
 kalpha=pi2*kalpha
 kbeta=pi2*kbeta
 ktheta=acos(2*ktheta-1)
 
 do i=1,nk
  call bmode(kphi(i),ktheta(i),kalpha(i),kbeta(i),cbm_n(:,i),sbm_n(:,i),kvec(:,i))
 end do


end subroutine generate_bmode




subroutine bmode(kphi,ktheta,kalpha,kbeta,cbm_n,sbm_n,kvec)
!===========================================================
! Calculation of a single mode of turbulent magnetic field
!===========================================================
 real(8), intent(in) :: kphi, ktheta, kalpha, kbeta
 real(8), intent(out) :: cbm_n(3), sbm_n(3), kvec(3)
 type(vector_rot3d) :: rot
 real(8) :: c_alpha, s_alpha, c_beta, s_beta
 real(8) :: phi, theta, cbm_n1(3), sbm_n1(3)
 real(8) :: kvec1(3)

 c_alpha=cos(kalpha)
 s_alpha=sin(kalpha)
 c_beta=cos(kbeta)
 s_beta=sin(kbeta)

 cbm_n(1)=c_alpha*c_beta
 cbm_n(2)=-s_alpha*s_beta
 cbm_n(3)=0d0

 sbm_n(1)=-c_alpha*s_beta
 sbm_n(2)=-s_alpha*c_beta
 sbm_n(3)=0d0

 phi=kphi*raddeg ! in grad
 theta=ktheta*raddeg


 call rot%set(2)   ! from k to k', where should be k'=[0,0,1]
 call rot%set_rot(1,3,-phi)
 call rot%set_rot(2,2,-theta)


 kvec1=[0d0,0d0,1d0]
 cbm_n1=cbm_n
 sbm_n1=sbm_n

 call rot%rotate_back(cbm_n1,cbm_n)
 call rot%rotate_back(sbm_n1,sbm_n)
 call rot%rotate_back(kvec1,kvec)

end subroutine bmode


end module turb_realiz_m




!module test_turb_realiz
! use turb_realiz_m, only : turb_realiz
! use phys_const, only : c_light, pi, psec
! use vector_ops_mod, only : vector_ops
! use utools_mod, only : utools
! 
! implicit none
! private
! save

! public :: main_calc
!  

! contains

!subroutine test_turb_mf
! type(turb_realiz) :: bt1
! integer :: nmodes, ns, ndev1, ndev2, i
! real(8) :: lmin, lmax, bms 
! real(8) :: x(3), bm(3), ndir1(3), ndir2(3)
! real(8), allocatable :: st(:)
! type(vector_ops) :: vop
! type(utools) :: ut
! character(500) :: fdir, fname

! nmodes=1000
! lmin=5*psec
! lmax=200*psec
! bms=1d-6



! call bt1%set(nmodes,lmin,lmax,bms)

! ns=10000
! ndir1=vop%unit_vec(10d0,60d0)
! ndir2=vop%unit_vec(3d0,30d0)
! call ut%grid(st,0d0*psec,1000*psec,ns,'lin')

! fdir='/home/anton/work/project/halo/chp_prop/traj_test/'
! fname=trim(fdir)//'turb_bm_5_1.dat'

! open(newunit=ndev1,file=fname)
! fname=trim(fdir)//'turb_bm_5_2.dat'
! open(newunit=ndev2,file=fname)
! do i=1,ns
!  x=st(i)*ndir1
!  call bt1%calc(x,bm)
!  write(ndev1,*) sqrt(dot_product(x,x))/psec, sqrt(dot_product(bm,bm))
!  x=st(i)*ndir2
!  call bt1%calc(x,bm)
!  write(ndev2,*) sqrt(dot_product(x,x))/psec, sqrt(dot_product(bm,bm))

! end do

! close(ndev1)
! close(ndev2)

!end subroutine test_turb_mf



!subroutine main_calc

! call test_turb_mf

!end subroutine main_calc


!end module test_turb_realiz



!program main
! use test_turb_realiz, only : main_calc

! call main_calc

!end program main




