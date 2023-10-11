module int_3d
!====================================================================
! In the version 1(v1) the ability of integration of complex
! 3d volumes is added. The limits of integration are given by
! 4 functions: fx1(y,z), fx2(y,z), fy1(z), fy2(z)
!====================================================================
 use gauss_kronrod, only : gk_adaptive
 
 implicit none
 private
 save
 
 public :: set_int_3d, calc_int_3d, volume, comp_int_3d
 
 type volume
   sequence
   real(8):: x1
   real(8):: x2
   real(8):: y1
   real(8):: y2
   real(8):: z1
   real(8):: z2
 end type volume
 
 integer :: cvol ! cvol=0, if integration is rectangular, cvol=1, if integration has complex shape
 real(8) :: x1, x2, y1, y2, z1, z2 ! limits of integration
 real(8) :: xc, yc, zc ! current value of x,y,z in integration 
 real(8) :: err_tol
 
 
 abstract interface
  subroutine fun_prototype(x,y,z,res)
   real(8), intent(in) :: x, y, z
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 abstract interface
  subroutine lim_fx(y,z,res)
   real(8), intent(in) :: y, z
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 abstract interface
  subroutine lim_fy(z,res)
   real(8), intent(in) :: z
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 procedure (fun_prototype), pointer :: fun3d=>null()
 procedure (lim_fx), pointer :: fx1=>null(), fx2=>null()
 procedure (lim_fy), pointer :: fy1=>null(), fy2=>null()
 
 
  
 contains
 
 
subroutine set_int_3d(fun,tol)
! Variables
 
 real(8), optional, intent(in) :: tol ! error tolerance in 1D integration
 real(8), parameter :: tol_def=1d-4 ! default value
 
 interface
  subroutine fun(x,y,z,res)
   real(8), intent(in) :: x, y, z
   real(8), intent(out) :: res
  end subroutine
 end interface
! Calculations 
 
 if (present(tol)) then
  err_tol=tol
 else
  err_tol=tol_def
 end if
  
 fun3d=>fun

end subroutine set_int_3d




subroutine calc_int_3d(vol,res)
!
! gives final result of 3d integral
!
! Variables 
 type(volume), intent(in) :: vol
 real(8), intent(out) :: res
! Calculations 
! Limits of integration
 x1=vol%x1
 x2=vol%x2
 y1=vol%y1
 y2=vol%y2
 z1=vol%z1
 z2=vol%z2
 cvol=0
 call gk_adaptive(iz_int,z1,z2,res,tol=err_tol)
end subroutine calc_int_3d


subroutine comp_int_3d(zz1,zz2,ffy1,ffy2,ffx1,ffx2,res)
!
! gives final result of 3d integral
!
! Variables 
 real(8), intent(in) :: zz1, zz2
 real(8), intent(out) :: res
 procedure (lim_fx) :: ffx1, ffx2
 procedure (lim_fy) :: ffy1, ffy2
! Calculations 
! Limits of integration
 z1=zz1
 z2=zz2
 fy1=>ffy1
 fy2=>ffy2
 fx1=>ffx1
 fx2=>ffx2
 cvol=1
 call gk_adaptive(iz_int,z1,z2,res,tol=err_tol)
end subroutine comp_int_3d



subroutine iz_int(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
! Calculations  
 zc=x
 if (cvol==1) then
  call fy1(zc,y1)
  call fy2(zc,y2)
 end if 
 call gk_adaptive(iy_int,y1,y2,res,tol=err_tol)
end subroutine iz_int


subroutine iy_int(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
! Calculations  
 yc=x
 if (cvol==1) then
  call fx1(yc,zc,x1)
  call fx2(yc,zc,x2)
 end if 
 call gk_adaptive(ix_int,x1,x2,res,tol=err_tol)
end subroutine iy_int


subroutine ix_int(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
! Calculations 
 xc=x
 call fun3d(xc,yc,zc,res)
end subroutine ix_int

end module int_3d