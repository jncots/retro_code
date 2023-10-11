module gg_elpos
!=====================================================================
! Calculation of electron-positron pair production rate and
! absorption probability.
!
! Usage:
! ======
! use gg_elpos, only : gg_elpos_rate
! type(gg_elpos_rate) :: rt ! declaration of the variable
!
! rt%fname=fname ! name of the file with data
! call rt%read ! read data from the file
! call rt%set_gam(eg1,eg2,gamsp) ! setting gamma spectrum
! call rt%wabs(egam,wabs) ! absorption probability for gamma with energy egam
!! egam in [eV], wabs in [1/cm]
! call rt%rate_epp(ee,rate) ! rate is the production rate of electron-positron
!! pairs in [1/(eV sec)], ee in [eV]
!! multiplied by 2 to account for both electrons and positrons
!
! Change history:
! ==============
! 08.04.2016 Using type(tabfun) the amount of code is reduced
! 18.08.2016 Change to object-oriented approach
!=====================================================================
 use gauss_kronrod, only : gk_adaptive
 use intpol_mod, only : arr_ind_short, logl_int
 use tab2fun, only : tabfun
 implicit none
 private
 save
 
 public :: gg_elpos_rate
 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type gg_elpos_rate
  integer :: ne, ng, nw
  real(8), allocatable :: ea(:), ga(:,:), fa(:,:)
  real(8), allocatable :: ew(:), ww(:)
  real(8) :: eg1, eg2
  procedure (fun_prot), pointer, nopass :: gam_spec=>null()
  type(tabfun) :: tww
  character(500) :: fname
 contains
  procedure :: read=>read_tab, wabs, set_gam, rate_epp, del
! read - read table, wabs - total absorption probability,
! set_gam - set gamma spectrum,
! rate_epp - calculate rate of electron-positron pairs
! del - deallocate variable
 end type gg_elpos_rate
 
 
 type(tabfun) :: rp_gam ! rate per gamma
 procedure (fun_prot), pointer :: gam_sp=>null()
 
 
 contains


subroutine del(this)
 class(gg_elpos_rate) :: this
 this%ne=0
 this%ng=0
 this%nw=0
 if (allocated(this%ea)) deallocate(this%ea)
 if (allocated(this%ga)) deallocate(this%ga)
 if (allocated(this%fa)) deallocate(this%fa)
 if (allocated(this%ew)) deallocate(this%ew)
 if (allocated(this%ww)) deallocate(this%ww)
 this%gam_spec=>null()
 call this%tww%off
end subroutine del

 
subroutine set_gam(this,eg1,eg2,gamsp)
 class(gg_elpos_rate) :: this
 real(8), intent(in) :: eg1, eg2 ! minimum and maximum arguments for gamsp(x)
 procedure (fun_prot) :: gamsp
 this%gam_spec=>gamsp
 this%eg1=eg1
 this%eg2=eg2
end subroutine set_gam

 
subroutine read_tab(this)
! Variables
 class(gg_elpos_rate) :: this
 integer, parameter :: ndev=25
 integer :: ndes, ne, ng, nw
 character(len=:), allocatable :: description
! Calculations

 call this%del
 open(ndev,file=this%fname,form='unformatted')
 read(ndev) ndes
 allocate(character(len=ndes) :: description)
 read(ndev) description
 read(ndev) ne
 allocate(this%ea(ne))
 read(ndev) this%ea
 read(ndev) ng, ne
 this%ng=ng
 this%ne=ne
 allocate(this%ga(ng,ne))
 allocate(this%fa(ng,ne))
 read(ndev) this%ga
 read(ndev) this%fa
 read(ndev) nw
 allocate(this%ew(nw))
 allocate(this%ww(nw))
 this%nw=nw
 read(ndev) this%ew
 read(ndev) this%ww
 close(ndev)
 
 
 call this%tww%on(this%ew,this%ww)
 
end subroutine read_tab


subroutine wabs(this,x,res)
!============================
! Absorption probability
!============================
! Variables
 class(gg_elpos_rate) :: this
 real(8), intent(in) :: x
 real(8), intent(out) :: res
! Calculations 
 call this%tww%fx(x,res)
end subroutine wabs



subroutine rate_epp(this,x,res)
! Variables
 class(gg_elpos_rate) :: this
 real(8), intent(in) :: x ! electron energy
 real(8), intent(out) :: res ! production rate of electron+positron
 integer :: ne1, ne2
 real(8) :: x1, x2, y1, y2
! Calculations
 
 if ((x<this%ea(1)).or.(x>this%ea(this%ne))) then
  res=0d0
  return
 end if
  
 call arr_ind_short(1,this%ne,this%ea,x,ne1,ne2)
 x1=this%ea(ne1)
 x2=this%ea(ne2)
 call rate_ee_dis(this,ne1,y1)
 call rate_ee_dis(this,ne2,y2)
 call logl_int(x1,x2,y1,y2,x,res)
 
end subroutine rate_epp




subroutine rate_ee_dis(this,nec,res)
! Variables
 class(gg_elpos_rate) :: this
 integer, intent(in) :: nec
 real(8), intent(out) :: res
 real(8) :: egam1, egam2, tmin, tmax
! Calculations 
 
 call rp_gam%on(this%ga(:,nec),this%fa(:,nec))
 if (rp_gam%zero==0) then
  res=0d0
  return
 end if
 
 egam1=rp_gam%xmin
 egam2=rp_gam%xmax 
 if (egam1<this%eg1) egam1=this%eg1
 if (egam2>this%eg2) egam2=this%eg2
 
 if (egam1>egam2) then
  res=0d0
  return
 end if 
 
 tmin=log(egam1)
 tmax=log(egam2)
 gam_sp=>this%gam_spec
 call gk_adaptive(irate_ee,tmin,tmax,res)

end subroutine rate_ee_dis



subroutine irate_ee(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: eph, rpg, fph
! Calculations 
 eph=exp(x)
 call rp_gam%fx(eph,rpg)
 call gam_sp(eph,fph)
 res=fph*rpg*eph
end subroutine irate_ee

end module gg_elpos


!============================================================
! Example of usage (test program)
!============================================================
!program main
! call test
!end program main
!
!
!
!subroutine test
! use gg_elpos, only : gg_elpos_rate
! use phys_const, only : psec, c_light
!! Variables
! type(gg_elpos_rate) :: rt
! integer :: nn, i, ne
! real(8) :: Ee1, Ee2, eel, res, egr1, egr2, egam
! real :: start, finish
! 
!  
! interface
!  subroutine gam_sp(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine
! end interface
! 
! 
!! Calculations 
! 
! rt%fname='D:\project\Diffusion\SphericalDiffusion\Estimations\xray_gc\soft_phf\soft_phf500.dat'
! call rt%read
! egr1=1d8
! egr2=1d18
! call rt%set_gam(1d8,1d18,gam_sp)
! 
! 
! Ee1=1d8
! Ee2=1d18
! nn=1000
! 
! call cpu_time(start)
! 
! open(1,file='rate_pp600.dat') 
! do i=0,nn
!  eel=Ee1*(Ee2/Ee1)**(1d0*i/nn)
!  call rt%rate_epp(eel,res)
!  write(1,*) eel, res*eel**2
! end do 
! close(1)
! 
! 
! ng=300
! open(1,file='wgam.dat')
!  do i=0,ng
!   egam=egr1*(egr2/egr1)**(1d0*i/ng)
!   call rt%wabs(egam,res)
!   write(1,*) egam, (c_light/res)/(psec)
!  end do
! close(1) 
! call cpu_time(finish)
! write(*,*) 'Elapsed time=', finish-start
! 
!end subroutine test
!
!
!
!subroutine gam_sp(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
! real(8) :: a
! 
! res=1d0/x**2
!end subroutine gam_sp
