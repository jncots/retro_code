module frame_rot3d_m
!=======================================================
! Calculation of the coordinates of the vector in the 
! new frame obtained by sequence of rotations around 
! determined axes
!=======================================================
 implicit none
 private
 save
 
 public :: rot2d, frame_rot3d
 
 real(8), parameter :: pi=3.14159265358979324d0 ! pi number
 real(8), parameter :: degrad=pi/180
 real(8), parameter :: imat(3,3)= reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/), (/3,3/))


 type rot2d
  integer :: nax=3      ! by default rotation around z
  real(8) :: phi=0d0    ! by default no rotation
  real(8) :: rm(3,3)=imat, rm_inv(3,3)=imat ! by default no rotation
 contains
  procedure :: calc=>get_rm,  clean=>clean_rm
 end type rot2d
 

 type frame_rot3d
  integer :: nrot
  type(rot2d), allocatable :: tran(:) ! transformations
  real(8) :: trm(3,3)=imat, trm_inv(3,3)=imat ! total rotation matrix and its inverse
 contains
  procedure :: set=>set_nrot, del=>del_nrot, set_rot, calc_rm
  procedure :: to_new_frame, from_new_frame
 end type frame_rot3d



 contains

subroutine del_nrot(this)
 class(frame_rot3d) :: this
 if (allocated(this%tran)) deallocate(this%tran)
 this%nrot=0
 this%trm=imat
 this%trm_inv=imat
end subroutine del_nrot

subroutine set_nrot(this,nrot)
 class(frame_rot3d) :: this
 integer, intent(in) :: nrot
 integer :: i

 if (this%nrot==nrot) then
  do i=1,this%nrot
   call this%tran(i)%clean
  end do
  this%trm=imat
  this%trm_inv=imat
 else  
  call this%del
  this%nrot=nrot
  allocate(this%tran(nrot))
 end if

end subroutine set_nrot


subroutine set_rot(this,rotn,nax,phi)
 class(frame_rot3d) :: this
 integer, intent(in) :: rotn, nax
 real(8), intent(in) :: phi
 
 if (rotn>this%nrot) then
  write(*,*) 'Rotation number is greater than number of allocated array'
  return
 end if

 this%tran(rotn)%nax=nax
 this%tran(rotn)%phi=phi
 call this%tran(rotn)%calc


end subroutine set_rot


subroutine calc_rm(this)
 class(frame_rot3d) :: this
 integer :: i

 this%trm=this%tran(1)%rm
 do i=2,this%nrot
  this%trm=matmul(this%trm,this%tran(i)%rm) 
! we use transposed matrices, so we multiply in opposite order
 end do
 
 this%trm_inv=this%tran(this%nrot)%rm_inv
 do i=this%nrot-1,1,-1
  this%trm_inv=matmul(this%trm_inv,this%tran(i)%rm_inv)
 end do
end subroutine calc_rm

subroutine to_new_frame(this,x,xnew)
 class(frame_rot3d) :: this
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: xnew(3)
! Change coordinates to new reference frame:  
 xnew(1)=dot_product(this%trm(:,1),x)
 xnew(2)=dot_product(this%trm(:,2),x)
 xnew(3)=dot_product(this%trm(:,3),x)
end subroutine to_new_frame

subroutine from_new_frame(this,xnew,x)
 class(frame_rot3d) :: this
 real(8), intent(in) :: xnew(3)
 real(8), intent(out) :: x(3)
! Change coordinates to new reference frame:  
 x(1)=dot_product(this%trm_inv(:,1),xnew)
 x(2)=dot_product(this%trm_inv(:,2),xnew)
 x(3)=dot_product(this%trm_inv(:,3),xnew)
end subroutine from_new_frame



subroutine clean_rm(this)
 class(rot2d) :: this

 this%nax=3
 this%phi=0d0
 this%rm=imat
 this%rm_inv=imat

end subroutine clean_rm


subroutine get_rm(this)
 class(rot2d) :: this
 real(8) :: dphi, cp, sp
 
 dphi=degrad*this%phi   ! from degrees to radians
 cp=cos(dphi)
 sp=sin(dphi)

 select case (this%nax)
  case(1)
! Rotation around x from y to z
  this%rm(:,1)=[1d0,0d0,0d0]
  this%rm(:,2)=[0d0,cp,sp]
  this%rm(:,3)=[0d0,-sp,cp]

  this%rm_inv(:,1)=[1d0,0d0,0d0]
  this%rm_inv(:,2)=[0d0,cp,-sp]
  this%rm_inv(:,3)=[0d0,sp,cp]

  case(2)
! Rotation around y from z to x
  this%rm(:,1)=[cp,0d0,-sp]
  this%rm(:,2)=[0d0,1d0,0d0]
  this%rm(:,3)=[sp,0d0,cp]

  this%rm_inv(:,1)=[cp,0d0,sp]
  this%rm_inv(:,2)=[0d0,1d0,0d0]
  this%rm_inv(:,3)=[-sp,0d0,cp]


  case(3)
! Rotation around z from x to y
  this%rm(:,1)=[cp,sp,0d0]
  this%rm(:,2)=[-sp,cp,0d0]
  this%rm(:,3)=[0d0,0d0,1d0]

  this%rm_inv(:,1)=[cp,-sp,0d0]
  this%rm_inv(:,2)=[sp,cp,0d0]
  this%rm_inv(:,3)=[0d0,0d0,1d0]


 case default
  write(*,*) 'Axis number is out of possible values=1,2,3'
  write(*,*) 'Axis number is set to default nax=3'
! Rotation around z from x to y
  this%rm(:,1)=[cp,sp,0d0]
  this%rm(:,2)=[-sp,cp,0d0]
  this%rm(:,3)=[0d0,0d0,1d0]

  this%rm_inv(:,1)=[cp,-sp,0d0]
  this%rm_inv(:,2)=[sp,cp,0d0]
  this%rm_inv(:,3)=[0d0,0d0,1d0]

 end select


end subroutine get_rm



end module frame_rot3d_m


!program main
! use frame_rot3d_m, only : rot2d, frame_rot3d
! type(frame_rot3d) :: rt3d
! real(8) :: x(3), xnew(3)
!  
! x=[0d0,1d0,0d0]

! call rt3d%set(4)
! call rt3d%set_rot(1,2,90d0)
! call rt3d%set_rot(2,3,19d0)
! call rt3d%set_rot(3,2,42d0)
! call rt3d%set_rot(1,1,-50d0)
! call rt3d%calc_rm


! call rt3d%to_new_frame(x,xnew)
! write(*,'(A,3Es14.6)') 'xnew=', xnew
! call rt3d%from_new_frame(xnew,x)
! write(*,'(A,3Es14.6)') 'xold=', x

!end program main
