module type_info
!================================================================
! This module contains description of type 'info', which is
! just character variable, which saves long line in array of
! lines with length=100 characters. Addition of the text 'text'
! to the line 'ff' is just:
!   call ff%add(len(text),text)
! deallocation, writing or reading to the device with unit='unit':
!   call ff%del
!   call ff%write(unit)
!   call ff%read(unit)
! The addition is without trailing spaces
!================================================================
 implicit none
 private
 
 public :: info
 
 type info
  integer :: nmes
  character(100), allocatable :: message(:)
 contains
  procedure :: del=>del_info , write=>write_info, read=>read_info, add
 end type info
 
 
 contains

subroutine del_info(ff)
! Variables
 class(info) :: ff
 ff%nmes=0
 if (allocated(ff%message)) deallocate(ff%message)
end subroutine del_info

subroutine write_info(ff,ndev)
! Variables
 class(info) :: ff
 integer, intent(in) :: ndev
 integer :: i
! Calculations 
 if (allocated(ff%message)) then
  write(ndev) ff%nmes
  write(ndev) ff%message
 else
  ff%nmes=0
  write(ndev) ff%nmes
 end if
end subroutine write_info



subroutine read_info(ff,ndev)
! Variables
 class(info) :: ff
 integer, intent(in) :: ndev
 integer :: i
! Calculations 
 read(ndev) ff%nmes
 if (ff%nmes>0) then
  if (allocated(ff%message)) deallocate(ff%message)
  allocate(ff%message(ff%nmes))
  read(ndev) ff%message
 end if 
end subroutine read_info


subroutine add(ff,nchar,line)
!=====================================
! Add text of line to the variable ff
!=====================================
! Variables
 class(info) :: ff
 integer, intent(in) :: nchar ! lenth of line
 character(nchar), intent(in) :: line
 integer, parameter :: lmes=100 ! length of character in type info
 character(lmes), allocatable :: mes1(:)
 integer :: i, nmes, nfill, nemp, nline, nadd, ntot, n1, n2
! Calculations

! If there is not yet filled last line 
 if (allocated(ff%message)) then
  nmes=ff%nmes
  nfill=len_trim(ff%message(nmes))
  if (nfill<lmes) then
   nemp=lmes-nfill-1
   ff%message(nmes)(nfill+2:lmes)=line(1:nemp)
  end if
 else
  nmes=0
  nemp=0
 end if
 
 nline=len_trim(line) 
 nadd=ceiling((nline-nemp)*1d0/lmes)
 ntot=nmes+nadd
 allocate(mes1(ntot))
 
 if (nmes>0) then
  do i=1,nmes
   mes1(i)=ff%message(i)
  end do
 end if
 
 n1=1+nemp
 n2=n1+lmes-1
 do i=nmes+1,ntot
  if (i==1) then
   n2=n2-1
   mes1(i)='\n'//line(n1:n2) ! division of the line
  else
   mes1(i)=line(n1:n2)
  end if 
  n1=n2+1
  n2=n1+lmes-1
  if (n2>nline) n2=nline
 end do
 
 if (allocated(ff%message)) deallocate(ff%message)
 allocate(ff%message(ntot))
 ff%nmes=ntot
 ff%message=mes1
 deallocate(mes1) 
 

end subroutine add

end module type_info








module type_fun_data
!====================================================================
!
! Module contains the type 'fun_data', which describes
! tabulated function of 1,2,3, or 4 dimentions
! Each dimension characterised by data (x or f array) and
! metadata such as name, short name, units of measurement and other
! information. Each variable of 'fun_data' type can be deallocated
! read, or written to file='fname' with len(fname)=lf as follows
!  call f%del
!  call f%read(lf,fname)
!  call f%write(lf,fname)
! For initialization the array of sizes in each dimention should be
! given, e.g. n(3)=[n1,n2,n3], then
!  call f%init(n)
!
!====================================================================
 use type_info, only : info
 implicit none
 private
 save
 
 public :: fun_data
 
 type finfo
  character(100) :: name      ! name of variable, e.g. 'Energy'
  character(100) :: short     ! short name, e.g. 'E'
  character(100) :: unitm     ! units, e.g. 'eV', 'GeV', 'TeV'
  type(info) :: extra         ! additional information   
 end type finfo
 
 
 type farg
  integer :: nx               ! array size
  real(8), allocatable :: x(:)! data array
  character(100) :: name      ! name of variable, e.g. 'Energy'
  character(100) :: short     ! short name, e.g. 'E'
  character(100) :: unitm     ! units, e.g. 'eV', 'GeV', 'TeV'
  type(info) :: extra         ! additional information 
 end type farg 
 
 type fun_data 
  integer :: ndim                  ! dimention of the function
  type(finfo), allocatable :: info ! information about function
  type(farg), allocatable :: a(:)  ! argument variable
  real(8), allocatable :: f1(:), f2(:,:), f3(:,:,:), f4(:,:,:,:)
 contains
  procedure :: init, del, write=>write_data, read=>read_data, check
 end type fun_data
 
 
 contains
 

subroutine check(fdat,res)
! Variables
 class(fun_data) :: fdat
 character(100) :: res
 integer :: i
! Calculations
 
 do i=1,fdat%ndim
  if (fdat%a(i)%nx==size(fdat%a(i)%x)) then
   res='OK'
  else
   write(res,'(A,I5)') 'No correspondence between x and nx along i=', i
   return
  end if
 end do
 
 select case (fdat%ndim)
  case(1)
   do i=1,fdat%ndim
    if (fdat%a(i)%nx==size(fdat%f1,i)) then
     res='OK'
    else
     write(res,'(A,I5)') 'No correspondence between f and nx along i=', i
     return
    end if
   end do
  case(2) 
   do i=1,fdat%ndim
    if (fdat%a(i)%nx==size(fdat%f2,i)) then
     res='OK'
    else
     write(res,'(A,I5)') 'No correspondence between f and nx along i=', i
     return
    end if
   end do
  case(3)
   do i=1,fdat%ndim
    if (fdat%a(i)%nx==size(fdat%f3,i)) then
     res='OK'
    else
     write(res,'(A,I5)') 'No correspondence between f and nx along i=', i
     return
    end if
   end do
  case(4)
   do i=1,fdat%ndim
    if (fdat%a(i)%nx==size(fdat%f4,i)) then
     res='OK'
    else
     write(res,'(A,I5)') 'No correspondence between f and nx along i=', i
     return
    end if
   end do
  case default 
   write(*,*) 'type(fun_data): max of ndim=4'
  end select 
 
end subroutine check
 

subroutine del(fdat)
! Variables
 class(fun_data) :: fdat
! Calculations 
 fdat%ndim=0 
 if (allocated(fdat%info)) deallocate(fdat%info)
 if (allocated(fdat%a)) deallocate(fdat%a)
 if (allocated(fdat%f1)) deallocate(fdat%f1)
 if (allocated(fdat%f2)) deallocate(fdat%f2)
 if (allocated(fdat%f3)) deallocate(fdat%f3)
 if (allocated(fdat%f4)) deallocate(fdat%f4) 
end subroutine del


 
 
subroutine init(fdat,n)
! Variables
 class(fun_data) :: fdat
 integer, intent(in) :: n(:)
 integer :: ndim, i
! Calculations 
 call del(fdat)
 ndim=size(n)
 fdat%ndim=ndim
 
 
 allocate(fdat%info) 
 fdat%info%name='none'
 fdat%info%short='none'
 fdat%info%unitm='none'
 fdat%info%extra%nmes=0
 
 
 allocate(fdat%a(ndim)) 
 do i=1,ndim
  fdat%a(i)%nx=n(i)
  if (allocated(fdat%a(i)%x)) deallocate(fdat%a(i)%x)
  allocate(fdat%a(i)%x(n(i)))
  fdat%a(i)%name='none'
  fdat%a(i)%short='none'
  fdat%a(i)%unitm='none'
  fdat%a(i)%extra%nmes=0
 end do
 
 select case (ndim)
  case(1)
   allocate(fdat%f1(n(1)))
  case(2) 
   allocate(fdat%f2(n(1),n(2)))
  case(3)
   allocate(fdat%f3(n(1),n(2),n(3)))
  case(4)
   allocate(fdat%f4(n(1),n(2),n(3),n(4)))
  case default 
   write(*,*) 'type(fun_data): max of ndim=4'
  end select
   
end subroutine init



subroutine write_data(fdat,lf,fname)
! Read of fdata from file connected to ndev
 class(fun_data) :: fdat
 integer :: lf
 character(lf) :: fname
 integer :: ndev, i, j
! Calculations
 call get_unit(ndev)
 open(ndev,file=fname,form='unformatted')
 write(ndev) fdat%ndim
 call fdat%info%extra%write(ndev)
 write(ndev) fdat%info%name
 write(ndev) fdat%info%short
 write(ndev) fdat%info%unitm
 
 do i=1,fdat%ndim
  write(ndev) fdat%a(i)%nx
  write(ndev) fdat%a(i)%x
  write(ndev) fdat%a(i)%name
  write(ndev) fdat%a(i)%short
  write(ndev) fdat%a(i)%unitm
  call fdat%a(i)%extra%write(ndev)
 end do
 
 
 select case (fdat%ndim)
  case(1)
   write(ndev) fdat%f1
  case(2) 
   write(ndev) fdat%f2
  case(3)
   write(ndev) fdat%f3
  case(4)
   write(ndev) fdat%f4
  case default 
   write(*,*) 'type(fun_data): max of ndim=4'
 end select
 close(ndev)
end subroutine write_data



subroutine read_data(fdat,lf,fname)
! Read of fdata from file connected to ndev
 class(fun_data) :: fdat
 integer :: lf
 character(lf) :: fname
 integer :: ndev, i, j
! Calculations
 call get_unit(ndev)
 call del(fdat)
 open(ndev,file=fname,form='unformatted')
 read(ndev) fdat%ndim
 allocate(fdat%info)
 call fdat%info%extra%read(ndev)
 read(ndev) fdat%info%name
 read(ndev) fdat%info%short
 read(ndev) fdat%info%unitm
 
 allocate(fdat%a(fdat%ndim))
 do i=1,fdat%ndim
  read(ndev) fdat%a(i)%nx
  allocate(fdat%a(i)%x(fdat%a(i)%nx))
  read(ndev) fdat%a(i)%x
  read(ndev) fdat%a(i)%name
  read(ndev) fdat%a(i)%short
  read(ndev) fdat%a(i)%unitm
  call fdat%a(i)%extra%read(ndev)
 end do
 
 
 select case (fdat%ndim)
  case(1)
   allocate(fdat%f1(fdat%a(1)%nx))
   read(ndev) fdat%f1
  case(2) 
   allocate(fdat%f2(fdat%a(1)%nx,fdat%a(2)%nx))
   read(ndev) fdat%f2
  case(3)
   allocate(fdat%f3(fdat%a(1)%nx,fdat%a(2)%nx,fdat%a(3)%nx))
   read(ndev) fdat%f3
  case(4)
   allocate(fdat%f4(fdat%a(1)%nx,fdat%a(2)%nx,fdat%a(3)%nx,fdat%a(4)%nx))
   read(ndev) fdat%f4
  case default 
   write(*,*) 'type(fun_data): max of ndim=4'
 end select
 close(ndev)
end subroutine read_data


subroutine get_unit(unit)
!====================================
! Search for free 'unit' from (1..99)
!====================================
! Variables
 integer, intent(out) :: unit
 integer :: i
! Calculations
 unit=0
 do i=1,99
  if ((i/=5).and.(i/=6).and.(i/=9)) then
   open(unit=i,err=10,status='scratch')
   close(unit=i)
   unit=i
   return
  end if
10  continue
 end do
end subroutine get_unit
 

end module type_fun_data

! Tests
!program main
! use type_fun_data, only : test
! 
! call test
!
!end program main
!
!
!
! 
!subroutine test
! type(fun_data) :: fdat, fdat1
! integer :: n(2), i
! character(100) :: fname
! character(:), allocatable :: line
! character(2000) :: ttt
! 
! n=[10,10]
! 
! call fdat%init(n)
! 
! fdat%a(1)%x=5
! fdat%a(1)%name='energy'
! fdat%a(2)%name='distance'
! fdat%f2(1,2)=4
! 
! ttt='In the past decade, mobile computing           has gone from a niche market for well-heeled &
! enterprises with large field organisations to the fastest growing, and often most popular, &
! way for employees of organisations of all sizes to do business computing. &
! The near-universal adoption of mobile devices by consumers-who are also &
! employees-has forced one of the most major shifts that corporate IT has ever seen.'
! line=ttt
! call line_break(70,line)
! ttt=line
! 
!! write(*,*) len(ttt)
!! write(*,*) trim(ttt)
! 
! if (allocated(fdat%info%extra%message)) then
!  write(*,*) 'yes'
! else
!  write(*,*) 'no'
! end if 
! 
! call fdat%info%extra%add(len(ttt),ttt)
! !fdat%info%extra%message(1)='\n This \n is \n something'
! write(*,*) fdat%info%extra%nmes
! write(*,*) fdat%info%extra%message
! read(*,*)
! 
! fname='test.dat'
! call fdat%write(fname)
! call fdat1%read(fname)
! 
! write(*,'(Es14.6)') fdat1%a(1)%x
! do i=1,2
! write(*,*) trim(fdat1%a(i)%name)
! end do
! 
! write(*,'(Es14.6)') fdat1%f2(1,2)
! 
! write(*,*) fdat1%info%extra%message
! write(*,'(A)') trim(fdat1%a(1)%short)
! write(*,'(A)') trim(fdat1%a(1)%unitm)
! 
!end subroutine test
