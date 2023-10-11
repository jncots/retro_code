module int_regions_mod
!==================================================
! Calculation of the integral over the region of
! сut out by the tube with radius robs in the
! the sphere lying on the axis of the tube with 
! radius rg2 minus the
! sphere radius rg1, rg1<rg2 and common centre.
! One considers that the function fxyz(x,y,z)
! is isotropic relative to the rotation around
! z axis.
! Example of usage is in the end of the file
!
! The volume cut out by the cylinder with radius d
! of the sphere with radius of R is
! V=4*Pi*R^3/3*(1-(1-(d/R)^2)^(3/2))
!==================================================
 use int_3d_mod, only : int_3dt
 implicit none
 private
 save
 
 
 public :: int_regions
 
 
 abstract interface
  subroutine fun3d(x,y,z,res)
   real(8), intent(in) :: x, y, z
   real(8), intent(out) :: res
  end subroutine fun3d
 end interface
 
 
 type int_regions
  real(8) :: rg1, rg2, robs
  real(8) :: tol1=1d-4, tol2=1d-4, tol3=1d-4
  procedure(fun3d), pointer, nopass :: fxyz=>null()
 contains
  procedure :: set_fun, set_reg, set_tol, int_res
 end type int_regions
 
 
 real(8) :: rg1, rg2, robs
  

 contains


subroutine set_fun(this,fxyz)
 class(int_regions) :: this
 procedure(fun3d) :: fxyz
 this%fxyz=>fxyz
end subroutine set_fun

subroutine set_reg(this,rg1,rg2,robs)
 class(int_regions) :: this
 real(8), intent(in) :: rg1, rg2, robs
 this%rg1=rg1
 this%rg2=rg2
 this%robs=robs
end subroutine set_reg

subroutine set_tol(this,tol1,tol2,tol3)
 class(int_regions) :: this
 real(8), intent(in) :: tol1, tol2, tol3
 this%tol1=tol1
 this%tol2=tol2
 this%tol3=tol3
end subroutine set_tol



subroutine int_res(this,res)
 class(int_regions) :: this
 real(8), intent(out) :: res
 
 rg1=this%rg1
 rg2=this%rg2
 robs=this%robs
 
 
 if (rg1>=rg2) then
  write(*,*) 'int_regions: rg2>=rg1'
  res=0d0
  return
 end if
 
 if (rg2<=robs) then
  call int_res1(this,res)
 else if ((rg1<=robs).and.(robs<rg2)) then
  call int_res2(this,res)
 else  
  call int_res3(this,res)
 end if 
 
end subroutine int_res

 
 

subroutine int_res1(this,res)
! For rg1<rg2<robs
 class(int_regions) :: this
 real(8), intent(out) :: res
 type(int_3dt) :: ical
 real(8) :: lz_min1, lz_max1, lz_min2, lz_max2
 real(8) :: lz_min3, lz_max3, lz_min4,lz_max4
 real(8) :: res1, res2, res3, res4
 
 res=0d0 
! Edges
! 1
 lz_min1=-rg2
 lz_max1=-rg1
 call ical%set(this%fxyz)
 call ical%set_l1(lz_min1,lz_max1,3,tol=this%tol1)
 call ical%set_l2(lx_min1,lx_max1,1,tol=this%tol2)
 call ical%set_l3(ly_min1,ly_max1,2,3,1,tol=this%tol3)
 call ical%calc(res1)
 res=res+res1
 
! 2 
 lz_min2=rg1
 lz_max2=rg2
 call ical%set_l1(lz_min2,lz_max2,3,tol=this%tol1) ! different from region 1
 call ical%set_l2(lx_min1,lx_max1,1,tol=this%tol2) ! the same as in region 1
 call ical%set_l3(ly_min1,ly_max1,2,3,1,tol=this%tol3)
 call ical%calc(res2)
 res=res+res2
 
 
! Centre
! 3
 lz_min3=-rg1
 lz_max3=rg1 
 call ical%set_l1(lz_min3,lz_max3,3,tol=this%tol1)
 call ical%set_l2(lx_min3,lx_max3,1,tol=this%tol2)
 call ical%set_l3(ly_min3,ly_max3,2,3,1,tol=this%tol3)
 call ical%calc(res3)
 res=res+res3

! 4 
 lz_min4=-rg1
 lz_max4=rg1 
 call ical%set_l1(lz_min4,lz_max4,3,tol=this%tol1)
 call ical%set_l2(lx_min4,lx_max4,1,tol=this%tol2)
 call ical%set_l3(ly_min4,ly_max4,2,3,1,tol=this%tol3)
 call ical%calc(res4)
 res=res+res4 

 
 res=4*res ! because we use spherical symmetry in xy plane
end subroutine int_res1 
 
 
subroutine int_res2(this,res)
! For rg1<robs<rg2
 class(int_regions) :: this
 real(8), intent(out) :: res
 type(int_3dt) :: ical
 real(8) :: lz_min21, lz_max21, lz_min22, lz_max22
 real(8) :: lz_min23, lz_max23, lz_min24, lz_max24
 real(8) :: lz_min25, lz_max25, lz_min26, lz_max26
 real(8) :: res1, res2, res3, res4, res5, res6
 
 res=0d0
 
! 1
 lz_min21=-rg2
 lz_max21=-sqrt(rg2**2-robs**2)
 call ical%set(this%fxyz)
 call ical%set_l1(lz_min21,lz_max21,3,tol=this%tol1)
 call ical%set_l2(lx_min21,lx_max21,1,tol=this%tol2)
 call ical%set_l3(ly_min21,ly_max21,2,3,1,tol=this%tol3)
 call ical%calc(res1)
 res=res+res1
 
! 2 
 lz_min22=sqrt(rg2**2-robs**2)
 lz_max22=rg2
 call ical%set_l1(lz_min22,lz_max22,3,tol=this%tol1) ! different from region 1
 call ical%set_l2(lx_min21,lx_max21,1,tol=this%tol2) ! the same as in region 1
 call ical%set_l3(ly_min21,ly_max21,2,3,1,tol=this%tol3)
 call ical%calc(res2)
 res=res+res2

! 3
 lz_min23=-sqrt(rg2**2-robs**2)
 lz_max23=-rg1
 call ical%set_l1(lz_min23,lz_max23,3,tol=this%tol1)
 call ical%set_l2(lx_min23,lx_max23,1,tol=this%tol2)
 call ical%set_l3(ly_min23,ly_max23,2,3,1,tol=this%tol3)
 call ical%calc(res3)
 res=res+res3

! 4
 lz_min24=rg1
 lz_max24=sqrt(rg2**2-robs**2)
 call ical%set_l1(lz_min24,lz_max24,3,tol=this%tol1)
 call ical%set_l2(lx_min23,lx_max23,1,tol=this%tol2) ! the same as in region 3
 call ical%set_l3(ly_min23,ly_max23,2,3,1,tol=this%tol3)
 call ical%calc(res4)
 res=res+res4 


! 5
 lz_min25=-rg1
 lz_max25=rg1
 call ical%set_l1(lz_min25,lz_max25,3,tol=this%tol1)
 call ical%set_l2(lx_min25,lx_max25,1,tol=this%tol2)
 call ical%set_l3(ly_min25,ly_max25,2,3,1,tol=this%tol3)
 call ical%calc(res5)
 res=res+res5 
 
! 6
 lz_min26=-rg1
 lz_max26=rg1
 call ical%set_l1(lz_min26,lz_max26,3,tol=this%tol1)
 call ical%set_l2(lx_min26,lx_max26,1,tol=this%tol2)
 call ical%set_l3(ly_min26,ly_max26,2,3,1,tol=this%tol3)
 call ical%calc(res6)
 res=res+res6 
 
 res=4*res ! because we use spherical symmetry in xy plane
end subroutine int_res2 
 
 
subroutine int_res3(this,res)
! For robs<rg1<rg2
 class(int_regions) :: this
 real(8), intent(out) :: res
 type(int_3dt) :: ical
 real(8) :: lz_min21, lz_max21, lz_min22, lz_max22
 real(8) :: lz_min23, lz_max23, lz_min24, lz_max24
 real(8) :: lz_min25, lz_max25, lz_min26, lz_max26
 real(8) :: res1, res2, res3, res4, res5, res6, res7, res8
 
 res=0d0
 
! 1
 lz_min21=-rg2
 lz_max21=-sqrt(rg2**2-robs**2)
 call ical%set(this%fxyz)
 call ical%set_l1(lz_min21,lz_max21,3,tol=this%tol1)
 call ical%set_l2(lx_min21,lx_max21,1,tol=this%tol2)
 call ical%set_l3(ly_min21,ly_max21,2,3,1,tol=this%tol3)
 call ical%calc(res1)
 res=res+res1
 
! 2 
 lz_min22=sqrt(rg2**2-robs**2)
 lz_max22=rg2
 call ical%set_l1(lz_min22,lz_max22,3,tol=this%tol1) ! different from region 1
 call ical%set_l2(lx_min21,lx_max21,1,tol=this%tol2) ! the same as in region 1
 call ical%set_l3(ly_min21,ly_max21,2,3,1,tol=this%tol3)
 call ical%calc(res2)
 res=res+res2

! 3
 lz_min23=-sqrt(rg2**2-robs**2)
 lz_max23=-rg1
 call ical%set_l1(lz_min23,lz_max23,3,tol=this%tol1)
 call ical%set_l2(lx_min23,lx_max23,1,tol=this%tol2)
 call ical%set_l3(ly_min23,ly_max23,2,3,1,tol=this%tol3)
 call ical%calc(res3)
 res=res+res3

! 4
 lz_min24=rg1
 lz_max24=sqrt(rg2**2-robs**2)
 call ical%set_l1(lz_min24,lz_max24,3,tol=this%tol1)
 call ical%set_l2(lx_min23,lx_max23,1,tol=this%tol2) ! the same as in region 3
 call ical%set_l3(ly_min23,ly_max23,2,3,1,tol=this%tol3)
 call ical%calc(res4)
 res=res+res4 


! 5
 lz_min25=-rg1
 lz_max25=-sqrt(rg1**2-robs**2)
 call ical%set_l1(lz_min25,lz_max25,3,tol=this%tol1)
 call ical%set_l2(lx_min25,lx_max25,1,tol=this%tol2)
 call ical%set_l3(ly_min25,ly_max25,2,3,1,tol=this%tol3)
 call ical%calc(res5)
 res=res+res5 
 
! 6
 lz_min26=-rg1
 lz_max26=-sqrt(rg1**2-robs**2)
 call ical%set_l1(lz_min26,lz_max26,3,tol=this%tol1)
 call ical%set_l2(lx_min26,lx_max26,1,tol=this%tol2)
 call ical%set_l3(ly_min26,ly_max26,2,3,1,tol=this%tol3)
 call ical%calc(res6)
 res=res+res6 
 
! 7
 lz_min25=sqrt(rg1**2-robs**2)
 lz_max25=rg1
 call ical%set_l1(lz_min25,lz_max25,3,tol=this%tol1)
 call ical%set_l2(lx_min25,lx_max25,1,tol=this%tol2)
 call ical%set_l3(ly_min25,ly_max25,2,3,1,tol=this%tol3)
 call ical%calc(res7)
 res=res+res7 
 
! 8
 lz_min26=sqrt(rg1**2-robs**2)
 lz_max26=rg1
 call ical%set_l1(lz_min26,lz_max26,3,tol=this%tol1)
 call ical%set_l2(lx_min26,lx_max26,1,tol=this%tol2)
 call ical%set_l3(ly_min26,ly_max26,2,3,1,tol=this%tol3)
 call ical%calc(res8)
 res=res+res8  
 
 res=4*res ! because we use spherical symmetry in xy plane
end subroutine int_res3 
 

! For rg1<rg2<robs 
! function describing limits for 
! region 1 and 2
subroutine lx_min1(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=0d0
end subroutine lx_min1

subroutine lx_max1(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2)
end subroutine lx_max1

subroutine ly_min1(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=0d0
end subroutine ly_min1

subroutine ly_max1(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2-x**2)
end subroutine ly_max1

! region 3
subroutine lx_min3(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=0d0
end subroutine lx_min3

subroutine lx_max3(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg1**2-z**2)
end subroutine lx_max3

subroutine ly_min3(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(rg1**2-z**2-x**2)
end subroutine ly_min3

subroutine ly_max3(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2-x**2)
end subroutine ly_max3


! region 4
subroutine lx_min4(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg1**2-z**2)
end subroutine lx_min4

subroutine lx_max4(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2)
end subroutine lx_max4

subroutine ly_min4(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=0d0
end subroutine ly_min4

subroutine ly_max4(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2-x**2)
end subroutine ly_max4

! For rg1<robs<rg2

subroutine lx_min21(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=0d0
end subroutine lx_min21

subroutine lx_max21(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2)
end subroutine lx_max21


subroutine ly_min21(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=0d0
end subroutine ly_min21

subroutine ly_max21(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(rg2**2-z**2-x**2)
end subroutine ly_max21


 
subroutine lx_min23(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=0d0
end subroutine lx_min23

subroutine lx_max23(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=robs
end subroutine lx_max23


subroutine ly_min23(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=0d0
end subroutine ly_min23

subroutine ly_max23(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(robs**2-x**2)
end subroutine ly_max23

 
subroutine lx_min25(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=0d0
end subroutine lx_min25

subroutine lx_max25(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg1**2-z**2)
end subroutine lx_max25


subroutine ly_min25(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(rg1**2-z**2-x**2)
end subroutine ly_min25

subroutine ly_max25(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(robs**2-x**2)
end subroutine ly_max25

 
subroutine lx_min26(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=sqrt(rg1**2-z**2)
end subroutine lx_min26

subroutine lx_max26(z,res)
 real(8), intent(in) :: z
 real(8), intent(out) :: res
 res=robs
end subroutine lx_max26


subroutine ly_min26(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=0d0
end subroutine ly_min26

subroutine ly_max26(z,x,res)
 real(8), intent(in) :: z, x
 real(8), intent(out) :: res
 res=sqrt(robs**2-x**2)
end subroutine ly_max26 
 


end module int_regions_mod


! module test_int_regions
 ! use int_regions_mod, only : int_regions
 ! implicit none
 ! private
 ! save
 
 ! public :: main_calc
 
 ! contains
 
! subroutine main_calc
 ! type(int_regions) :: fg
 ! real(8) :: res

 ! call fg%set_fun(fun)
 ! call fg%set_reg(50d0,100d0,4d0)
 ! call fg%set_tol(1d-3,1d-3,1d-3)
 ! call fg%int_res(res)
 
 ! write(*,*) res

! end subroutine main_calc


! subroutine fun(x,y,z,res) 
 ! real(8), intent(in) :: x, y, z
 ! real(8), intent(out) :: res
 ! res=1d0
! end subroutine fun


! end module test_int_regions



! program main
 ! use test_int_regions, only : main_calc
 ! call main_calc
! end program main
