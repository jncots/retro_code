module txt_file_m
 implicit none
 private
 save

 public :: txt_file !, main_calc


 type txt_file
  character(500) :: fname='file.dat'
  integer :: nlines=0
  character(1000), allocatable :: line(:) 
 contains
  procedure :: read=>read_file, fscript, del
 end type txt_file


 contains

subroutine del(this)
 class(txt_file) :: this
 
 if (allocated(this%line)) deallocate(this%line)
 this%nlines=0
 this%fname='file.dat'

end subroutine del


subroutine read_file(this)
!===========================================
! Read file to array of strings,
! each string is line
!===========================================
 class(txt_file) :: this
 integer :: ndev, n, istat, i
 character(1000) :: str

 open(newunit=ndev,file=this%fname)
 
 n=0
 do
  read(ndev,'(A)',iostat=istat) str
  if (istat<0) exit
  n=n+1
 end do
 rewind(ndev)
 
 this%nlines=n
 if (allocated(this%line)) deallocate(this%line)
 allocate(this%line(n))
 do i=1,n
  read(ndev,'(A)') this%line(i)
 end do 
 
 close(ndev)
end subroutine read_file


subroutine fscript(this,fname)
!===========================================
! Write the peace of code for executable
! in bash (sh) variable
!===========================================
 class(txt_file) :: this
 character(500), intent(in) :: fname
 integer :: ndev, i

 open(newunit=ndev,file=fname)
 
 do i=1,this%nlines
  write(ndev,'(A)')  'comm=trim(comm)//'//"'"//trim(this%line(i))//"\\n"//"'"
 end do 
 close(ndev)
end subroutine fscript




!subroutine main_calc
! character(500) :: fdir, scr, fort
! type(txt_file) :: tf

! fdir='/home/anton/work/project/halo/chp_prop/'
! scr=trim(fdir)//'search_data.sh'

! tf%fname=scr
! call tf%read

! fort=trim(fdir)//'search_data2.f90'

! call tf%fscript(fort)
! 
!end subroutine main_calc


end module txt_file_m

!program main
! use txt_file_m, only : main_calc
! 
! call main_calc

!end program main
