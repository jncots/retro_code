!  =========== DOCUMENTATION ===========
!
!  Definition:
!  ===========
!    descr_line (subroutine) - creates a string with information
!    infoline (type) - data structure used by descr_line
!
!  Syntax:
!  =======
!    call descr_line(nitem,note,ndes,description)
!
!  Arguments:
!  ==========
!    nitem(input,integer)  
!
!  Example:
!  ========
!   subroutine info
!   !
!   ! Write information of calculation into description(ndes) variable
!   !
!   ! Variables
!    integer, parameter :: nitem=3
!    type(infoline) :: note(nitem)
!    integer :: ndes
!    character(len=:), allocatable :: description
!   ! Calculations
!    note(1)%typ='des'
!    note(1)%intro='Introduction'
!    note(1)%info='This file contains radial profile of surface brightness\n of gamma radiation.'
! 
!    note(2)%typ='sub'
!    note(2)%intro='Setup of the calculations'
!    note(2)%sub_name='set_profile'
!    note(2)%file_name='rad_profile.f90'
! 
!    note(3)%typ='sub'
!    note(3)%intro='Description of the data structure of the current file'
!    note(3)%sub_name='out_profile'
!    note(3)%file_name='rad_profile.f90'
! 
!    call descr_line(nitem,note,ndes,description)
! 
!   end subroutine info
!
!  Result:
!  =======
!   ---------------------------------------------------------------------
!   
!                 Introduction                                                                                        
!
!   ---------------------------------------------------------------------
!   This file contains radial profile of surface brightness
!    of gamma radiation.
!   ---------------------------------------------------------------------
!   
!   ---------------------------------------------------------------------
!   
!                 Setup of the calculations                                                                           
!   
!   ---------------------------------------------------------------------
!   subroutine set_profile
!    ...(Content of the subroutine)...
!   end subroutine set_profile
!   ---------------------------------------------------------------------
!   
!______________________________________________________________________________
!   
!                 Description of the data structure of the current file                                               
!______________________________________________________________________________
!
!   subroutine out_profile
!    ...(Content of the subroutine)...
!   end subroutine out_profile
!_______________________________________________________________________________
!
!
!
!_______________________________________________________________________________
! call descr_line(nitem,note,ndes,description)
!-------------------------------------------------------------------------------
! Situation: 
! 1. Need for writing down information in the form:
!      -------------------------
!            Heading
!      -------------------------
!           Information ....
!           ................
!      -------------------------
!     
!
!
!
!
module description_line
 implicit none
 private
 save
 
 public :: descr_line, infoline
 
 type infoline
   sequence
   character(10) :: typ
   character(100) :: intro
   character(100) :: sub_name
   character(500) :: file_name
   character(2000) :: info
 end type infoline
 
 
 contains
 
 
subroutine descr_line(nitem,note,ndes,description)
! Variables
 integer, intent(in) :: nitem
 type(infoline), intent(in) :: note(nitem)
 integer, intent(out) :: ndes
 character(len=:), allocatable, intent(out) :: description
 character(len=:), allocatable :: des
 integer, parameter :: nres=10000
 character(nres) :: res
 character(100) :: sub_name, file_name, intro
 character(2000) :: info
 integer :: i
! Calculations

 do i=1,nitem
  if (note(i)%typ=='sub') then
   intro=note(i)%intro
   sub_name=note(i)%sub_name
   file_name=note(i)%file_name
   call print_subroutine(nres,sub_name,file_name,intro,res)
  end if
  if (note(i)%typ=='des') then
   intro=note(i)%intro
   info=note(i)%info
   call print_info(nres,intro,info,res)
  end if
  
  if (i==1) then
   ndes=len_trim(trim(res))
   allocate(character(len=ndes) :: description)
   description=trim(res)
  else
   ndes=len_trim(trim(description)//trim(res))
   allocate(character(len=ndes) :: des)
   des=trim(description)//trim(res)
   deallocate(description)
   allocate(character(len=ndes) :: description)
   description=trim(des)
   deallocate(des)
  end if 
 end do
end subroutine descr_line


subroutine print_info(nres,intro,info,res)
! Variables
 integer, intent(in) :: nres
 character(100), intent(in) :: intro
 character(2000), intent(in) :: info
 character(nres), intent(out) :: res
 
! Write result variable
 res='\n---------------------------------------------------------------------\n'
 res=trim(res)//'\n              '//intro//'\n'
 res=trim(res)//'\n---------------------------------------------------------------------\n'
 res=trim(res)//trim(info)//'\n'
 res=trim(res)//'---------------------------------------------------------------------\n'
end subroutine print_info



subroutine print_subroutine(nres,sub_name,file_name,intro,res)
! Variables
 integer, intent(in) :: nres
 character(100), intent(in) :: sub_name, file_name, intro
 character(nres), intent(out) :: res
 
 integer, parameter :: ndev=10
 integer :: n, istat, i, ns
 integer :: nstart, nend
 character(1000) :: str
 character(1000), allocatable :: fstr(:)
! Calculations

! Read file 
 open(ndev,file=file_name)
 
 n=0
 do
  read(ndev,'(A)',iostat=istat) str
  if (istat<0) exit
  n=n+1
 end do

 rewind(ndev)
 allocate(fstr(n))
 do i=1,n
  read(ndev,'(A)') fstr(i)
 end do 
 
 close(ndev)

 nstart=0
 nend=0
 
! Search the line 
 do i=1,n
  ns=index(trim(fstr(i)),trim(sub_name))
  if (ns/=0) then
   ns=index(trim(fstr(i)),'subroutine')
   if (ns/=0) then
    ns=index(trim(fstr(i)),'end')
    
    if (ns/=0) then
     nend=i
    else
     nstart=i
    end if
       
   end if 
  end if 
 end do
 
 if ((nstart==0).or.(nend==0)) then
  write(*,'(A)') 'print_subroutine: '//trim(sub_name)//' is absent in the file '//trim(file_name)
  res=''
  return
 end if 
 
! Write result variable
 res='\n---------------------------------------------------------------------\n'
 res=trim(res)//'\n              '//intro//'\n'
 res=trim(res)//'\n---------------------------------------------------------------------\n'
 
 do i=nstart,nend
  res=trim(res)//trim(fstr(i))//'\n'
 end do
 
 res=trim(res)//'---------------------------------------------------------------------\n'
  
 deallocate(fstr)
end subroutine print_subroutine

end module description_line
