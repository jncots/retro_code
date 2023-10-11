module drift_mod
!============================================
! Calculation of the tangent vector to
! drift trajectory
!============================================
 use vector_ops_mod, only : vector_ops
 implicit none
 private
 save
 
 
 public :: drift_dir
! public :: test_drift_dir
 
 
 
 contains
 

subroutine drift_dir(ef,bf,npos,nneg)
!============================================
! Approximately works when |ef|<|bf|
! There is a condition which can deal
! with a problem if |ef|>=|bf| slightly
!============================================
 type(vector_ops) :: vp
 real(8), intent(in) :: ef(3), bf(3)
 real(8), intent(out) :: npos(3), nneg(3)
 real(8) :: bf2, beta(3), beta2, fb
 
 bf2=vp%vec2(bf)
 beta=vp%cross_prod(ef,bf)/bf2
 beta2=vp%vec2(beta)
 
 if (beta2<1d0) then 
  fb=sqrt((1-beta2)/bf2)
  npos=beta+fb*bf
  nneg=beta-fb*bf
 else
  npos=beta/sqrt(beta2)
  nneg=npos
 end if
  
end subroutine drift_dir



!subroutine test_drift_dir
! type(vector_ops) :: vp
! real(8) :: ef(3), bf(3), npos(3), nneg(3)
! real(8) :: bm, em
! 
! bm=100d0
! em=0.9d0*bm
! bf=bm*vp%unit_vec(0d0,90d0)
! ef=em*vp%unit_vec(85d0,85d0)
! 
! call drift_dir(ef,bf,npos,nneg)
! 
! write(*,'(A,10Es30.17)') 'n_pos=', npos
! write(*,'(A,10Es30.17)') 'n_neg=', nneg
! write(*,'(A,10Es30.17)') '|n_pos|=', 1-vp%vec_norm(npos)
! write(*,'(A,10Es30.17)') '|n_neg|=', 1-vp%vec_norm(nneg)
! 

!end subroutine  test_drift_dir


end module drift_mod



!program main
! use drift_mod, only : test_drift_dir
! 
! call test_drift_dir

!end program main
