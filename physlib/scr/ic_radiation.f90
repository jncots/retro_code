module inv_compton
 use gauss_kronrod, only : gk_adaptive
 use house_keeping, only : get_unit
 implicit none
 private
 
 
 
 public :: ictab_type
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 
 type ictab_type
  character(500) :: fout
  integer :: ng, ne, nel
  real(8) :: eps1, eps2, eg1, eg2, ee1, ee2
  real(8) :: eel1, eel2
  procedure(fun_1d), pointer, nopass :: sph=>null()
  real(8), allocatable :: eg(:), ee(:,:), rate(:,:)
  real(8), allocatable :: eel(:), prob(:), loss(:)
 contains
  procedure :: set_sph, set_gamr, set_elr, set_ellr, tab_rate 
 end type ictab_type
 
 
 
 
 real(8), parameter :: me=0.510998928d6
 real(8), parameter :: ime=4/me**2
! icnorm=2*Pi*re^2*c*me^2, re is the classical radius of electron, me in eV 
 real(8), parameter :: icnorm=3.90574640445437d-3
! re_norm=re^2*c, c - speed of light, re - classical radius of electrona
 real(8), parameter :: re_norm=2.38058775381309d-15 
 real(8) :: eps1, eps2
 real(8) :: egam0, ee0
 real(8) :: ee_gam, egam_ee
 real(8) :: ee_tw, ee_tl
 
 procedure(fun_1d), pointer :: sph_dist
 
 contains

 
subroutine set_sph(this,x1,x2,sph)
! Variables
 class(ictab_type) :: this
 real(8), intent(in) :: x1, x2
 procedure(fun_1d) :: sph
! Calculations 
 this%eps1=x1
 this%eps2=x2
 this%sph=>sph
! Default values
 this%eg1=1d5
 this%eg2=1d19
 this%ng=200
 
 this%ee1=1d7
 this%ee2=1d20
 this%ne=200
 
 this%eel1=1d7
 this%eel2=1d20
 this%nel=600
 
 this%fout='ic_sph_tab.dat'
 
end subroutine set_sph


subroutine set_gamr(this,eg1,eg2,ng)
! Variables
 class(ictab_type) :: this
 real(8), intent(in) :: eg1, eg2
 integer :: ng
! Calculations  
 this%eg1=eg1
 this%eg2=eg2
 this%ng=ng
end subroutine set_gamr

subroutine set_elr(this,ee1,ee2,ne)
! Variables
 class(ictab_type) :: this
 real(8), intent(in) :: ee1, ee2
 integer :: ne
! Calculations  
 this%ee1=ee1
 this%ee2=ee2
 this%ne=ne
end subroutine set_elr


subroutine set_ellr(this,ee1,ee2,ne)
! Variables
 class(ictab_type) :: this
 real(8), intent(in) :: ee1, ee2
 integer :: ne
! Calculations  
 this%eel1=ee1
 this%eel2=ee2
 this%nel=ne
end subroutine set_ellr

 
 
subroutine totw(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculation 
 ee_tw=x
 tmin=log(eps1)
 tmax=log(eps2)
 call gk_adaptive(itotw,tmin,tmax,res)
end subroutine totw
 

subroutine itotw(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eps, eta, sph, tw
! Calculations
 eps=exp(x)
 eta=eps*ee_tw*ime ! ime=4/me**2
 call sph_dist(eps,sph)
 call ic_totw(eta,tw)
 res=sph*tw*eps
end subroutine itotw
 
 
subroutine ic_totw(x,res)
!====================================================
! Analytical approximation (error within 2%) of total
! probability for IC interaction
! x=4*eps*ee/me^2
!====================================================
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x2
! Calculations
 x2=x*x
 res=log(1+2.3d-3*x2)
 res=res/(2.2d-2+2d-2*x2)
 res=1-res
 res=res*log(1.45d0+0.61d0*x)
 res=res/(1+1.78d0*x)
 res=res*22.4d0*re_norm
end subroutine ic_totw 
 

subroutine tloss(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations 
 ee_tl=x
 tmin=log(eps1)
 tmax=log(eps2)
 call gk_adaptive(itloss,tmin,tmax,res)
end subroutine tloss
 

subroutine itloss(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eps, eta, sph, tw
! Calculations
 eps=exp(x)
 eta=eps*ee_tl*ime ! ime=4/me**2
 call sph_dist(eps,sph)
 call ic_res_loss(eta,tw)
 res=sph*tw*eps*ee_tl
end subroutine itloss 
 
 
 
subroutine ic_res_loss(x,res)
!====================================================
! Analytical approximation (error within 2%) of total
! relativ energy losses for IC interaction
! x=4*eps*ee/me^2
!====================================================
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: x2
! Calculations
 x2=x*x
 res=4.6d-2*x/(1+4.9d-2*x2)
 res=1-res
 res=res*log(1+1.6d-1*x)
 res=res/(1+1.39d0*x)
 res=res*17.4d0*re_norm
end subroutine ic_res_loss


 
subroutine tab_rate(this)
! Variables
 class(ictab_type) :: this
 integer :: ndev
 integer :: i, j, ng, ne
 real(8) :: ee, egam, res
 real(8) :: egmin, eta, egmax, ee1, ee2, eemin, eemax, eg1, eg2
 character(100) :: fout
! Calculations
 
 eps1=this%eps1
 eps2=this%eps2
 sph_dist=>this%sph
 eg1=this%eg1
 eg2=this%eg2
  
 ng=this%ng
 ne=this%ne
 allocate(this%eg(ng+1))
 allocate(this%ee(ne+1,ng+1))
 allocate(this%rate(ne+1,ng+1))

 
 do i=0,ng
  egam=eg1*(eg2/eg1)**(1d0*i/ng)
  this%eg(i+1)=egam
  
  eemin=me**2/(this%eps2*egam)
  eemin=(sqrt(eemin+1)+1)*egam/2
  eemax=this%ee2
  
  do j=0,ne
   ee=eemin*(eemax/eemin)**(j*1d0/ne)
   this%ee(j+1,i+1)=ee
   call ic_sph(ee,egam,res)
   this%rate(j+1,i+1)=res
  end do
 end do
 
 
 ee1=this%eel1
 ee2=this%eel2
 ne=this%nel
 allocate(this%eel(ne+1))
 allocate(this%prob(ne+1))
 allocate(this%loss(ne+1))
 
 do i=0,ne
  ee=ee1*(ee2/ee1)**(i*1d0/ne)
  this%eel(i+1)=ee
  call totw(ee,this%prob(i+1))
  call tloss(ee,this%loss(i+1))
!  write(*,'(10Es14.6)') ee, this%prob(i+1), this%loss(i+1) 
 end do
 
 
 call get_unit(ndev)
 open(ndev,file=this%fout,form='unformatted')
 write(ndev) size(this%eg), size(this%ee,1)
 write(ndev) this%eg
 write(ndev) this%ee
 write(ndev) this%rate
 write(ndev) size(this%eel)
 write(ndev) this%eel
 write(ndev) this%prob
 write(ndev) this%loss
 close(ndev)
 
 deallocate(this%eg)
 deallocate(this%ee)
 deallocate(this%rate)
 deallocate(this%eel)
 deallocate(this%prob)
 deallocate(this%loss)
 

end subroutine tab_rate


subroutine ic_gam(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: ee
 ee=x
 call ic_sph(ee,egam_ee,res)
 res=res
end subroutine ic_gam

subroutine ic_sph(ee,egam,res)
! Variables
 real(8), intent(in) :: ee, egam
 real(8), intent(out) :: res
 real(8) :: eps_min, eps_max, tmin, tmax
! Calculations 
 ee0=ee
 egam0=egam
 
 eps_min=1/(ime*ee*(ee/egam-1))
 eps_max=ime*egam*ee**2
 
 
 if (eps_min<eps1) eps_min=eps1
 if (eps_max>eps2) eps_max=eps2
 

 tmin=log(eps_min)
 tmax=log(eps_max)
 
 if (tmax-tmin<=1d-14) then
  res=0d0
  return
 end if 
 
 call gk_adaptive(iic_sph,tmin,tmax,res)
 
 if (isnan(res).or.res<0d0) then
  write(*,'(A,10Es20.12)') 'tmin, tmax, res', tmin, tmax, res
  write(*,'(A,10Es20.12)') 'ee, egam', ee, egam
  write(*,*) eps_min, eps_max, eps_max-eps_min
  write(*,*) tmin, tmax, tmax-tmin, spacing(tmin) 
 end if 
end subroutine ic_sph

subroutine iic_sph(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eps, we, sph
! Calculations 
 eps=exp(x)
 call ic_prob(eps,ee0,egam0,we)
 call sph_dist(eps,sph)
 res=sph*we*eps
end subroutine iic_sph



subroutine ic_prob(eps,ee,egam,res)
! Variables
 real(8), intent(in) :: eps, ee, egam ! in eV
 real(8), intent(out) :: res ! in [cm^3/(s eV)]
 real(8) :: eta, q, s
! real(8) :: ime, icnorm
! Calculations
 eta=eps*ee*ime !ime=4/me**2
 q=eta*(ee/egam-1)
 q=1/q

 s=eta*q
 res=2*(1+s)
 res=s**2/res
 res=res+2*q+1
 res=res*(1-q)
 res=res+2*q*log(q)
 res=res/(eps*ee**2)
 res=icnorm*res
 
end subroutine ic_prob


end module inv_compton


!program main
! use inv_compton, only : set_sph, tab_rate
! use thermal_photon_field, only : therm_phf, set_thnden, thnden
!! Variables
! type(therm_phf) :: sf
!! Calculations 
! call sf%init(4)
! sf%temp=[3d0,3d-1,6d-3,2.35d-4]
! sf%uden=[5d4,5d4,5d4,0.26d0]
! call set_thnden(sf)
! 
! call set_sph(1d-10,1d2,thnden)
!
! 
! call tab_rate
!
!end program main


module ic_rad
 use gauss_kronrod, only : gk_adaptive
 use house_keeping, only : get_unit
 use intpol_mod, only : arr_ind_short, log_int, lin_int
 use tab2fun, only : tabfun
 implicit none
 private
 
 public :: ic_type
 
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 
 type ic_type
  character(500) :: fin
  integer :: ng, ne, nel
  real(8), allocatable :: eg(:), ee(:,:), rate(:,:)
  real(8), allocatable :: eel(:), prob(:), loss(:)
  real(8) :: ee1, ee2
  procedure(fun_1d), pointer, nopass :: spec=>null()
 contains
  procedure :: read_tab, ic_rate, set_spec, del
  procedure :: ic_loss, ic_prob
 end type ic_type
 
 type(tabfun) :: tprob
 procedure(fun_1d), pointer :: ee_spec=>null()
 
 
 
 contains
 
 
 
subroutine ic_loss(this,ee,res)
! Variables
 class(ic_type) :: this
 real(8), intent(in) :: ee
 real(8), intent(out) :: res
 type(tabfun) :: tt
! Calculations
 call tt%on(this%eel,this%loss)
 call tt%fx(ee,res)
 call tt%off
end subroutine ic_loss


subroutine ic_prob(this,ee,res)
! Variables
 class(ic_type) :: this
 real(8), intent(in) :: ee
 real(8), intent(out) :: res
 type(tabfun) :: tt
! Calculations
 call tt%on(this%eel,this%prob)
 call tt%fx(ee,res)
 call tt%off
end subroutine ic_prob

 
subroutine del(this)
! Variables
 class(ic_type) :: this
! Calculations
 
 this%ng=0
 this%ne=0
 this%nel=0
 if (allocated(this%eg)) deallocate(this%eg)
 if (allocated(this%ee)) deallocate(this%ee)
 if (allocated(this%rate)) deallocate(this%rate)
 if (allocated(this%eel)) deallocate(this%eel)
 if (allocated(this%prob)) deallocate(this%prob)
 if (allocated(this%loss)) deallocate(this%loss)
 this%ee1=0d0
 this%ee2=0d0
 this%spec=>null()
end subroutine del


subroutine set_spec(this,ee1,ee2,spec)
! Variables 
 class(ic_type) :: this
 real(8), intent(in) :: ee1, ee2
 procedure(fun_1d) :: spec
! Calculations  
 this%ee1=ee1
 this%ee2=ee2
 this%spec=>spec
end subroutine set_spec
 
 
subroutine read_tab(this)
! Variables
 class(ic_type) :: this
 integer :: ndev
! Calculations

 call this%del
 call get_unit(ndev)
 open(ndev,file=this%fin,form='unformatted')
 read(ndev) this%ng, this%ne
 allocate(this%eg(this%ng))
 allocate(this%ee(this%ne,this%ng))
 allocate(this%rate(this%ne,this%ng))
 read(ndev) this%eg
 read(ndev) this%ee
 read(ndev) this%rate
 read(ndev) this%nel
 allocate(this%eel(this%nel))
 allocate(this%prob(this%nel))
 allocate(this%loss(this%nel))
 read(ndev) this%eel
 read(ndev) this%prob
 read(ndev) this%loss
 close(ndev)
end subroutine read_tab



subroutine ic_rate(this,eg,res)
! Variables 
 class(ic_type) :: this
 real(8), intent(in) :: eg
 real(8), intent(out) :: res
 integer :: n1, n2
 real(8) :: res1, res2
! Calculations 
 
 if ((eg<this%eg(1)).or.(eg>this%eg(this%ng))) then
  res=0d0
  return
 end if 
 
 ee_spec=>this%spec
 call arr_ind_short(1,this%ng,this%eg,eg,n1,n2)
 call dist(this,n1,res1)
 call dist(this,n2,res2)
 call log_int(this%eg(n1),this%eg(n2),res1,res2,eg,res)
 if (isnan(res)) call lin_int(this%eg(n1),this%eg(n2),res1,res2,eg,res)
 

end subroutine ic_rate




subroutine dist(this,nn,res)
! Variables
 class(ic_type) :: this
 integer, intent(in) :: nn
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
! Calculations

 call tprob%on(this%ee(:,nn),this%rate(:,nn))
 
 tmin=this%ee1
 tmax=this%ee2
 if (tmin<tprob%xmin) tmin=tprob%xmin
 if (tmax>tprob%xmax) tmax=tprob%xmax

 tmin=log(tmin)
 tmax=log(tmax)
 call gk_adaptive(idist,tmin,tmax,res)
 call tprob%off

end subroutine dist



subroutine idist(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: ee, wg, fe
! Calculations
 ee=exp(x)
 call tprob%fx(ee,wg)
 call ee_spec(ee,fe)
 res=wg*fe*ee
end subroutine idist


end module ic_rad



!program main
! use ic_rad, only : ic_type
! use phys_const, only : year
! implicit none
!! Variables 
! type(ic_type) :: ff
! integer :: i, ng, ne
! real(8) :: res, egam
! real(8) :: eg1, eg2, ee1, ee2, ee
! 
! interface
!  subroutine ee_spec(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine ee_spec
! 
! end interface
! 
!! Calculations
!  ff%fin='ic_sph.dat'
!  call ff%read_tab
!  call ff%set_spec(1d9,1d20,ee_spec)
!  
!!  do i=1,size(ff%eel)
!!   write(*,'(10Es14.6)') ff%eel(i), ff%loss(i)/ff%prob(i)
!!  end do
!
! eg1=1d6
! eg2=1d17
! ng=100
! 
! open(1,file='rate_spec.dat')
! do i=0,ng
!  egam=eg1*(eg2/eg1)**(i*1d0/ng)
!  call ff%ic_rate(egam,res)
!  write(1,*) egam, res*egam**2
! end do 
! close(1)
! 
! ee1=1d8
! ee2=1d19
! ne=100
! open(1,file='loss.dat')
! do i=0,ne
!  ee=ee1*(ee2/ee1)**(i*1d0/ne)
!  call ff%ic_loss(ee,res)
!  write(1,*) ee, 1/(res*year)
! end do 
! close(1)
! 
!end program main
!
!
!subroutine ee_spec(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
!
! res=1/x**2
!end subroutine ee_spec