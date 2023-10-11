module plaw_inj_mod
 use gauss_kronrod, only : gk_adaptive
 implicit none
 private
 
 public :: plaw_inj
  
 type plaw_inj
  real(8) :: power, pindex, e1, e2
  real(8) :: norm  
 contains 
  procedure :: init=>qinit, spec=>qinj
 end type plaw_inj
 
 real(8), parameter :: ev_erg=6.241509647d11 ! eV per erg
 real(8) :: e12, pindex
 
 contains

subroutine qinit(this,pow,pind,e1,e2)
 class(plaw_inj) :: this
 real(8), intent(in) :: pow, pind, e1, e2
 
 this%power=pow ! in erg/s
 this%pindex=pind
 this%e1=e1
 this%e2=e2
 call normq(this)
end subroutine qinit
 


subroutine qinj(this,x,res)
 class(plaw_inj) :: this
 real(8), intent(in) :: x
 real(8), intent(out) :: res

 if (x<this%e1) then
!  write(*,*) 'qinj: energy is out of spectral range:'
!  write(*,'(2(A,Es14.6))') 'e=', x, ' e_min=', this%e1
  res=0d0
  return
 end if
 res=exp(-x/this%e2)*x**(-this%pindex)
 res=res*this%norm         ! result in 1/(eV s)
end subroutine qinj
 
 
subroutine normq(this)
!
! Calculation of normalization factor for power-law spectrum
! with exponential cut-off
!
 class(plaw_inj) :: this
 real(8), parameter :: tmin=0d0, tmax=1d0
 real(8) :: rnorm

 e12=this%e1/this%e2
 pindex=this%pindex
 call gk_adaptive(inormq,tmin,tmax,rnorm)
 rnorm=rnorm*this%e1**(2-this%pindex)
 this%norm=this%power*ev_erg/rnorm      ! normalization in eV
end subroutine normq
 

subroutine inormq(s,res)
!
! Integrand for normq
!
 real(8), intent(in) :: s
 real(8), intent(out) :: res
 real(8), parameter :: delta=6
 real(8) :: y
! Calculation
 y=s**delta
 res=y**(pindex-2)*exp(-e12/y)*delta/s
end subroutine inormq


end module plaw_inj_mod



!program main
! use plaw_inj_mod, only : plaw_inj
! implicit none
! integer :: i, j, n
! real(8) :: e1, e2, ee, res(3)
! type(plaw_inj) :: qp(3)
! 
! call qp(1)%init(1d40,1.5d0,1d9,1d20)
! call qp(2)%init(1d40,2d0,1d9,1d18)
! call qp(3)%init(1d40,2.5d0,1d9,1d18)
! 
! e1=1d9
! e2=1d19
! n=100
! 
! open(1,file='inj_spec.dat')
! do i=0,n
!  ee=e1*(e2/e1)**(1d0*i/n)
!  do j=1,3
!   call qp(j)%spec(ee,res(j))
!  end do 
!  write(1,*) ee, (res(j),j=1,3)
! end do
! close(1)
!
!end program main