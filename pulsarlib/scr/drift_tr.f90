module drift_mod
 implicit none
 private
 save
 
 
 
 public :: drift_dir
 
 
 
 contains
 
 
!subroutine main_calc
! real(8) :: ef(3), bf(3), beta(3)
! real(8) :: bm
! 
! bm=1d0
! call unit_vec(0d0,90d0,bf)
! bf=bm*bf
! call unit_vec(0d0,40d0,ef)
! ef=1d0*bm*ef
! 
! 
! call drift_dir(ef,bf,beta,'pos')
! 
! write(*,*) dot_product(beta,beta), beta
! 
! call drift_dir1(ef,bf,beta)
! write(*,*) dot_product(beta,beta), beta
! 

!end subroutine  main_calc
 


subroutine drift_dir_old(ef,bf,beta,tp)
!
! This code works if E is not perpendicular to B
!
 real(8), intent(in) :: ef(3), bf(3)
 real(8), intent(out) :: beta(3)
 character(3), intent(in), optional :: tp
 character(3), parameter :: tp_def='elc'
 character(3) :: tpp
 real(8) :: ef2, bf2, ebm, ebf, q, exb(3), ss, qb, qe
 
 if (present(tp)) then
  tpp=tp
 else 
  tpp=tp_def
 end if  
 
 ef2=dot_product(ef,ef)
 bf2=dot_product(bf,bf)
 
 ebm=(ef2-bf2)/2
 ebf=dot_product(ef,bf)
 ss=ebf/abs(ebf)
 if (isnan(ss)) ss=ebf
 q=sqrt(ebm**2+ebf**2)+(ef2+bf2)/2
 
 call cross_prod(ef,bf,exb)
 qb=q-bf2
 qe=q-ef2
! if (qb<0d0) qb=0d0
! if (qe<0d0) qe=0d0

 
 if (tpp=='elc') then
  beta=exb-(sqrt(qb)*ef+sqrt(qe)*ss*bf)
 else
  beta=exb+(sqrt(qb)*ef+sqrt(qe)*ss*bf)
 end if

 beta=beta/q
 
end subroutine drift_dir_old




subroutine drift_dir(ef,bf,beta)
!============================================
! Approximately works when |ef|<|bf|
! There is a condition which can deal
! with a problem if |ef|>=|bf| slightly
!============================================
 real(8), intent(in) :: ef(3), bf(3)
 real(8), intent(out) :: beta(3)
 real(8) :: bf2, exb(3), ff
 
 bf2=dot_product(bf,bf)
 call cross_prod(ef,bf,exb)

 beta=exb/bf2
 ff=dot_product(beta,beta)
 if (ff<1d0) then 
   ff=sqrt(1-ff)
   beta=beta+ff*bf/sqrt(bf2)
  else
   beta=beta/sqrt(ff)
  end if
  
end subroutine drift_dir




subroutine drift_dir1(ef,bf,beta)
 real(8), intent(in) :: ef(3), bf(3)
 real(8), intent(out) :: beta(3)
 real(8) :: bf2, exb(3), ff
 
 bf2=dot_product(bf,bf)
 call cross_prod(ef,bf,exb)
 
 beta=exb/bf2
 
 ff=dot_product(beta,beta)
 if (ff<1d0) then 
  ff=sqrt(1-ff)
 else
  ff=0d0
 end if

 
 beta=beta+ff*bf/sqrt(bf2)
! beta=bf/sqrt(bf2)
 
end subroutine drift_dir1


subroutine unit_vec(phi,theta,nv)
 use phys_const, only : pi
 real(8), intent(in) :: phi, theta
 real(8), intent(out) :: nv(3)
 real(8) :: phi1, theta1
 
 phi1=phi*pi/180
 theta1=theta*pi/180
  
 nv(1)=cos(phi1)*sin(theta1)
 nv(2)=sin(phi1)*sin(theta1)
 nv(3)=cos(theta1)
end subroutine unit_vec


subroutine cross_prod(a,b,res)
! Cross product
 real(8), intent(in) :: a(3), b(3)
 real(8), intent(out) :: res(3)
 
 res(1)=a(2)*b(3)-a(3)*b(2) 
 res(2)=a(3)*b(1)-a(1)*b(3)
 res(3)=a(1)*b(2)-a(2)*b(1)
end subroutine cross_prod

end module drift_mod



!program main
! use drift_mod, only : main_calc
! 
! call main_calc

!end program main
