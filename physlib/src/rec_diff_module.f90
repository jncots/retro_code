module rec_diff_module
!===============================================================
! Calculation of cosmic-ray density for the source qs(t,en)
! and diffusion coefficient diff(en), where t is time since
! the start of operation, en is energy of cosmic rays
!
! Example:
! ========
! use rec_diff_module, only : rec_diff
! type(rec_diff) :: this
! call this%init(diff,qs)
! call this%cr_den(en,t,r,res) ! [t]=sec, [r]=cm, res=[1/(cm^3 eV)]
! =========
! Here:
! diff is a subroutine diff(en,res) ! in cm^2/sec
! qs is a subroutine qs(t,en,res)   ! in 1/(sec eV), time in sec
!================================================================

!================================================================
! Contains a propagator of the cosmic rays in spherically
! symmetric case, which includes both diffusion and ballistic
! regime.
! The formula is given in Aloisio et at (2009)
! (2009ApJ...693.1275A), Eq.19
!
! The module uses a function BesselK1 (dbesk1)
! from Slatec library
!================================================================

 use gauss_kronrod, only : gk_adaptive
 use phys_const, only : c_light, pi
 
 implicit none
 private
 save
 
 public :: rec_diff
 
 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 abstract interface
  subroutine fun_2d(x,y,res)
   real(8), intent(in) :: x,y
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 
 type rec_diff
  real(8), private :: en, xi, trc, tobs
  procedure (fun_1d), pointer, nopass :: diff=>null() ! diff(energy)
  procedure (fun_2d), pointer, nopass :: qs=>null() ! qsrc(time,energy)
 contains
  procedure :: init=>init_dq, cr_den 
 end type rec_diff
 
 
 type(rec_diff), pointer :: loc
 
 
 contains

 
subroutine init_dq(this,diff,qs) 
 class(rec_diff) :: this
 procedure (fun_1d) :: diff
 procedure (fun_2d) :: qs
 this%diff=>diff
 this%qs=>qs
end subroutine init_dq




subroutine cr_den(this,en,tobs,robs,res)
 class(rec_diff), target :: this
 real(8), intent(in) :: en, tobs, robs
 real(8), intent(out) :: res
 real(8), parameter :: pi4=1/(4*pi*c_light), icl=1/c_light
 real(8) :: trc, td, diff, tmin, tmax
 
 
 trc=robs*icl
 td=tobs/trc
  
 if (td<1d0) then
  res=0d0
  return
 end if 
  
 call this%diff(en,diff)
 this%en=en
 this%tobs=tobs
 this%trc=trc
 this%xi=robs*c_light/(2*diff)
 loc=>this
  
 tmin=0d0
 tmax=sqrt(td**2-1)
 call gk_adaptive(icr_den,tmin,tmax,res,tol=1d-10)
 res=res*pi4/robs**2
 
end subroutine cr_den


subroutine icr_den(x,res)
 real(8), intent(in) :: x    
 real(8), intent(out) :: res
 real(8) :: tback, tpast, pres, qres
 
 tback=sqrt(x**2+1)*loc%trc
 tpast=loc%tobs-tback
 
 call pjun0(x,loc%xi,pres)
 call loc%qs(tpast,loc%en,qres)
 res=pres*qres
end subroutine icr_den




subroutine pjun0(s,xi,res)
 real(8), intent(in) :: s, xi
 real(8), intent(out) :: res
 real(8), parameter :: pi2=sqrt(2/pi) 
 real(8), parameter :: xmin=1d-10, xmax=700d0
 real(8) :: s2, s3, ss, ss2, x
 
 interface
! dbesk1(x) is a BesselK1 from Slatec library
  real(8) function dbesk1(x)
  real(8) :: x
  end function dbesk1
 end interface
  
 
 s2=s**2
 s3=s2*s
 ss=s2+1
 ss2=sqrt(ss)
 x=xi*ss2
 
 
 if (x<xmin) then
  res=exp(-xi*(s+1/s))*x**2
  res=res/s3
  return
 end if
 
 if (x>xmax) then
  res=pi2*x*sqrt(x)
  res=res*exp(-xi*ss2/(s2+s*ss2))
  res=res/s3
  return
 end if 

 res=x/dbesk1(x)
 res=res*exp(-xi*(s+1/s))
 res=res/s3
end subroutine pjun0




end module rec_diff_module