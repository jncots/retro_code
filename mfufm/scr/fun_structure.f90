module fun_structure
!
! Module contains the derived data types for storage and processing arrays
! for tabulated function of 1,2,3 or 4 arguments
!
! Example of usage is below
!
 implicit none
 private
 save
 
 public :: fun_data1d, fun_data2d, fun_data3d, fun_data4d
 
 type fun_data1d
  integer :: n1
  real(8), allocatable :: x1(:)
  real(8), allocatable :: fun(:)
 contains
  procedure :: initialize => init1d, read => read_1d, write => write_1d
  procedure :: terminate => term1d
 end type fun_data1d
 
 
 type fun_data2d
  integer :: n1, n2
  real(8), allocatable :: x1(:), x2(:)
  real(8), allocatable :: fun(:,:)
 contains
  procedure :: initialize => init2d, read => read_2d, write => write_2d
  procedure :: terminate => term2d
 end type fun_data2d
 
 
 type fun_data3d
  integer :: n1, n2, n3
  real(8), allocatable :: x1(:), x2(:), x3(:)
  real(8), allocatable :: fun(:,:,:)
 contains
  procedure :: initialize => init3d, read=>read_3d, write => write_3d
  procedure :: terminate => term3d
 end type fun_data3d
 
 type fun_data4d
  integer :: n1, n2, n3, n4
  real(8), allocatable :: x1(:), x2(:), x3(:), x4(:)
  real(8), allocatable :: fun(:,:,:,:)
 contains
  procedure :: initialize => init4d, read=>read_4d, write => write_4d
  procedure :: terminate => term4d
 end type fun_data4d
 
 
 
 contains


subroutine init1d(fdata,n1)
! Allocation of fun_data1d type variable 
 class(fun_data1d) :: fdata
 integer, intent(in) :: n1
! Calculation
 fdata%n1=n1
 if (.not. allocated(fdata%x1)) allocate(fdata%x1(n1))
 if (.not. allocated(fdata%fun)) allocate(fdata%fun(n1))
end subroutine init1d


subroutine term1d(fdata)
! Allocation of fun_data1d type variable 
 class(fun_data1d) :: fdata
! Calculation
 fdata%n1=0
 if (allocated(fdata%x1)) deallocate(fdata%x1)
 if (allocated(fdata%fun)) deallocate(fdata%fun)
end subroutine term1d



subroutine read_1d(fdata,ndev)
! Read of fdata from file connected to ndev
 class(fun_data1d) :: fdata
 integer, intent(in) :: ndev
 integer :: n1
! Calculations 
 read(ndev) n1
 call init1d(fdata,n1)
 read(ndev) fdata%x1
 read(ndev) fdata%fun
end subroutine read_1d

subroutine write_1d(fdata,ndev)
! Write of fdata from file connected to ndev
 class(fun_data1d) :: fdata
 integer, intent(in) :: ndev
! Calculations 
 write(ndev) fdata%n1
 write(ndev) fdata%x1
 write(ndev) fdata%fun
end subroutine write_1d


subroutine init2d(fdata,n1,n2)
! Allocation of fun_data2d type variable 
 class(fun_data2d) :: fdata
 integer, intent(in) :: n1, n2 
! Calculation
 fdata%n1=n1
 fdata%n2=n2
 if (.not. allocated(fdata%x1)) allocate(fdata%x1(n1))
 if (.not. allocated(fdata%x2)) allocate(fdata%x2(n2))
 if (.not. allocated(fdata%fun)) allocate(fdata%fun(n1,n2))
end subroutine init2d


subroutine term2d(fdata)
! Allocation of fun_data1d type variable 
 class(fun_data2d) :: fdata
! Calculation
 fdata%n1=0
 fdata%n2=0
 if (allocated(fdata%x1)) deallocate(fdata%x1)
 if (allocated(fdata%x2)) deallocate(fdata%x2)
 if (allocated(fdata%fun)) deallocate(fdata%fun)
end subroutine term2d




subroutine read_2d(fdata,ndev)
! Read of fdata from file connected to ndev
 class(fun_data2d) :: fdata
 integer, intent(in) :: ndev
 integer :: n1, n2
! Calculations 
 read(ndev) n1
 read(ndev) n2
 call init2d(fdata,n1,n2)
 read(ndev) fdata%x1
 read(ndev) fdata%x2
 read(ndev) fdata%fun
end subroutine read_2d


subroutine write_2d(fdata,ndev)
! Write of fdata from file connected to ndev
 class(fun_data2d) :: fdata
 integer, intent(in) :: ndev
! Calculations 
 write(ndev) fdata%n1
 write(ndev) fdata%n2
 write(ndev) fdata%x1
 write(ndev) fdata%x2
 write(ndev) fdata%fun
end subroutine write_2d





subroutine init3d(fdata,n1,n2,n3)
! Allocation of fun_data3d type variable 
 class(fun_data3d) :: fdata
 integer, intent(in) :: n1, n2, n3
! Calculation
 fdata%n1=n1
 fdata%n2=n2
 fdata%n3=n3
 if (.not. allocated(fdata%x1)) allocate(fdata%x1(n1))
 if (.not. allocated(fdata%x2)) allocate(fdata%x2(n2))
 if (.not. allocated(fdata%x3)) allocate(fdata%x3(n3))
 if (.not. allocated(fdata%fun)) allocate(fdata%fun(n1,n2,n3))
end subroutine init3d



subroutine term3d(fdata)
! Allocation of fun_data1d type variable 
 class(fun_data3d) :: fdata
! Calculation
 fdata%n1=0
 fdata%n2=0
 fdata%n3=0
 if (allocated(fdata%x1)) deallocate(fdata%x1)
 if (allocated(fdata%x2)) deallocate(fdata%x2)
 if (allocated(fdata%x3)) deallocate(fdata%x3)
 if (allocated(fdata%fun)) deallocate(fdata%fun)
end subroutine term3d



subroutine read_3d(fdata,ndev)
! Read of fdata from file connected to ndev
 class(fun_data3d) :: fdata
 integer, intent(in) :: ndev
 integer :: n1, n2, n3
! Calculations 
 read(ndev) n1
 read(ndev) n2
 read(ndev) n3
 call init3d(fdata,n1,n2,n3)
 read(ndev) fdata%x1
 read(ndev) fdata%x2
 read(ndev) fdata%x3
 read(ndev) fdata%fun
end subroutine read_3d


subroutine write_3d(fdata,ndev)
! Write of fdata from file connected to ndev
 class(fun_data3d) :: fdata
 integer, intent(in) :: ndev
! Calculations 
 write(ndev) fdata%n1
 write(ndev) fdata%n2
 write(ndev) fdata%n3
 write(ndev) fdata%x1
 write(ndev) fdata%x2
 write(ndev) fdata%x3
 write(ndev) fdata%fun
end subroutine write_3d



subroutine init4d(fdata,n1,n2,n3,n4)
! Allocation of fun_data4d type variable 
 class(fun_data4d) :: fdata
 integer, intent(in) :: n1, n2, n3, n4
! Calculation
 fdata%n1=n1
 fdata%n2=n2
 fdata%n3=n3
 fdata%n4=n4
 
 if (.not. allocated(fdata%x1)) allocate(fdata%x1(n1))
 if (.not. allocated(fdata%x2)) allocate(fdata%x2(n2))
 if (.not. allocated(fdata%x3)) allocate(fdata%x3(n3))
 if (.not. allocated(fdata%x4)) allocate(fdata%x4(n4))
 if (.not. allocated(fdata%fun)) allocate(fdata%fun(n1,n2,n3,n4))
end subroutine init4d


subroutine term4d(fdata)
! Allocation of fun_data1d type variable 
 class(fun_data4d) :: fdata
! Calculation
 fdata%n1=0
 fdata%n2=0
 fdata%n3=0
 fdata%n4=0
 if (allocated(fdata%x1)) deallocate(fdata%x1)
 if (allocated(fdata%x2)) deallocate(fdata%x2)
 if (allocated(fdata%x3)) deallocate(fdata%x3)
 if (allocated(fdata%x4)) deallocate(fdata%x4)
 if (allocated(fdata%fun)) deallocate(fdata%fun)
end subroutine term4d


subroutine read_4d(fdata,ndev)
! Read of fdata from file connected to ndev
 class(fun_data4d) :: fdata
 integer, intent(in) :: ndev
 integer :: n1, n2, n3, n4
! Calculations 
 read(ndev) n1
 read(ndev) n2
 read(ndev) n3
 read(ndev) n4
 call init4d(fdata,n1,n2,n3,n4)
 read(ndev) fdata%x1
 read(ndev) fdata%x2
 read(ndev) fdata%x3
 read(ndev) fdata%x4
 read(ndev) fdata%fun
end subroutine read_4d



subroutine write_4d(fdata,ndev)
! Write of fdata from file connected to ndev
 class(fun_data4d) :: fdata
 integer, intent(in) :: ndev
! Calculations 
 write(ndev) fdata%n1
 write(ndev) fdata%n2
 write(ndev) fdata%n3
 write(ndev) fdata%n4
 write(ndev) fdata%x1
 write(ndev) fdata%x2
 write(ndev) fdata%x3
 write(ndev) fdata%x4
 write(ndev) fdata%fun
end subroutine write_4d

end module fun_structure



!
!program example
! use fun_structure, only : fun_data4d
! integer :: ndev=15
! integer :: n1, n2, n3, n4
! type(fun_data4d) :: fd
! 
! n1=10
! n2=20
! n3=30
! n4=20
! call fd%initialize(n1,n2,n3,n4)
! fd%x4=10
! fd%fun=40
! 
! write(*,*) fd%n4
! open(ndev,file='fdtest.dat',form='unformatted')
! call fd%write(ndev)
! close(ndev)
! 
! open(ndev,file='fdtest.dat',form='unformatted')
! call fd%read(ndev)
! close(ndev)
! 
! write(*,*) fd%fun(2,3,3,2), fd%x4(2)
! call fd%terminate()
!
!
!end program example