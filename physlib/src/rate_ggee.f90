module rate_ggee
!=====================================================================
! Calculation of electron-positron pair production rate:
!
! call rate_ee(Ee,res), where
! Ee is energy of electron or positron in [eV]
! res is the production rate in [1/(eV sec)] multiplied by 2 to
! account for both electrons and positrons 
!
! Calculation of the total interaction probability of gamma rays with
! soft photon field:
!
! call wgam(Egam,res), where
! Egam is gamma-energy in [eV]
! res is the total interaction probability divided by speed of light
! (inverse mean free path) [1/cm]
!
! Information about soft photon field is contained in the data file,
! produced using module 'tab_ggee', which performs integration over
! the soft spectrum distribution. The tabulated information is read
! as following:
!
! call read_tab(fname), where
! fname is character(200) containing path and file name of data
!
! The spectrum of gamma rays is set as:
!
! call set_ggee(egmin,egmax,gam_sp), where
! egmin, egmax are minimum and maximum energy of gamma-ray spectrum
! gam_sp(x,res) is gamma-ray spectrum, where x is energy in eV
!
!
! Example of usage is in the end of the file
!
! Change history:
! ==============
! 08.04.2016 Using type(tabfun) the amount of code is reduced
!=====================================================================
 use gauss_kronrod, only : gk_adaptive
 use intpol_mod, only : log_int, lin_int, zeros_cut, arr_ind_short, arr_ind_b
 use tab2fun, only : tabfun
 implicit none
 private
 save
 
 public :: set_ggee, read_tab, rate_ee, wgam
 
 integer :: ne, ng, ng1, ng2, ngw
 real(8) :: egr1, egr2 ! minimum and maximum energies of gamma rays
 real(8), allocatable, target :: ee_arg(:), gam_arg(:,:), fres(:,:)
 real(8), allocatable :: xwgam(:), fwgam(:)

 
 abstract interface
  subroutine fun_prot(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type(tabfun) :: ww
 type(tabfun) :: rp_gam ! rate per gamma
 
 procedure (fun_prot), pointer :: gam_sp=>null()
 
 
 contains

 

subroutine set_ggee(x1,x2,gamsp)
 real(8), intent(in) :: x1, x2 ! minimum and maximum arguments for gamsp(x)
 procedure (fun_prot) :: gamsp
 gam_sp=>gamsp
 egr1=x1
 egr2=x2
end subroutine set_ggee

 
subroutine read_tab(fname)
! Variables
 character(200), intent(in) :: fname
 integer, parameter :: ndev=25
 integer :: ndes
 character(len=:), allocatable :: description
! Calculations
 
 open(ndev,file=fname,form='unformatted')
 read(ndev) ndes
 allocate(character(len=ndes) :: description)
 read(ndev) description
 read(ndev) ne
 allocate(ee_arg(ne))
 read(ndev) ee_arg
 read(ndev) ng, ne
 allocate(gam_arg(ng,ne))
 allocate(fres(ng,ne))
 read(ndev) gam_arg
 read(ndev) fres
 read(ndev) ngw
 allocate(xwgam(ngw))
 allocate(fwgam(ngw))
 read(ndev) xwgam
 read(ndev) fwgam
 close(ndev)
 
 
 call ww%on(xwgam,fwgam)
 
 
end subroutine read_tab


subroutine wgam(x,res)
!============================
! Absorbtion probability
!============================
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
! Calculations 
 call ww%fx(x,res)
end subroutine wgam




subroutine rate_ee(x,res)
! Variables
 real(8), intent(in) :: x ! electron energy
 real(8), intent(out) :: res ! production rate of electron+positron
 integer :: ne1, ne2
 real(8) :: x1, x2, y1, y2
! Calculations
 
 if ((x<ee_arg(1)).or.(x>ee_arg(ne))) then
  res=0d0
  return
 end if
  
 call arr_ind_b(ne,ee_arg,x,ne1,ne2)
 x1=ee_arg(ne1)
 x2=ee_arg(ne2)
 call rate_ee_dis(ne1,y1)
 call rate_ee_dis(ne2,y2)
 call log_int(x1,x2,y1,y2,x,res)
 if (isnan(res)) call lin_int(x1,x2,y1,y2,x,res)
 
end subroutine rate_ee




subroutine rate_ee_dis(nec,res)
! Variables
 integer, intent(in) :: nec
 real(8), intent(out) :: res
 real(8) :: egam1, egam2, tmin, tmax
! Calculations 
 
 call rp_gam%on(gam_arg(:,nec),fres(:,nec))
 if (rp_gam%zero==0) then
  res=0d0
  return
 end if
 
 egam1=rp_gam%xmin
 egam2=rp_gam%xmax 
 if (egam1<egr1) egam1=egr1
 if (egam2>egr2) egam2=egr2
 
 if (egam1>egam2) then
  res=0d0
  return
 end if 
 
 tmin=log(egam1)
 tmax=log(egam2)
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

end module rate_ggee


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
! use rate_ggee, only : set_ggee, read_tab, rate_ee, wgam
! use phys_const, only : psec, c_light
!! Variables
! character(100) :: fname
! character(100) :: fdir
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
! fdir='C:\Main\project\Diffusion\SphericalDiffusion\PrelimCalc1_1e17\ggee\\rate_ggee\'
! fname='ggee.dat'
! fname=trim(fdir)//trim(fname)
! 
! call read_tab(fname)
! egr1=1d8
! egr2=1d18
! call set_ggee(egr1,egr2,gam_sp)
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
!  call rate_ee(eel,res)
!  write(1,*) eel, res*eel**2
! end do 
! close(1)
! 
! 
! ng=300
! open(1,file='wgam.dat')
!  do i=0,ng
!   egam=egr1*(egr2/egr1)**(1d0*i/ng)
!   call wgam(egam,res)
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
