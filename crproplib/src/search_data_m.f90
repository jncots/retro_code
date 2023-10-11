module search_data_m
 use txt_file_m, only : txt_file

 implicit none
 private
 save

 public :: search_data !, main_calc


 type search_data
  character(500) :: search_scr='C:/work/library/crproplib/scr/search_data.sh'
  character(500) :: fout='data_data.dat'
  character(500) :: search_dir, prefix
  integer :: nfiles
  character(500), allocatable :: dfile(:)
 contains
  procedure :: search, del
 end type search_data
 

 contains


subroutine del(this)
 class(search_data) :: this
 
 if (allocated(this%dfile)) deallocate(this%dfile)
 this%nfiles=0

end subroutine del


subroutine search(this)
 class(search_data) :: this
 character(500) :: cwd, fout
 character(500) :: sdir, spref, sout
 type(txt_file) :: tf
 integer :: i
 character(1000) :: comm

 call getcwd(cwd)
 fout=trim(cwd)//'/'//trim(this%fout)
 sdir=' "'//trim(this%search_dir)//'"'
 spref=' "'//trim(this%prefix)//'"'
 sout=' "'//trim(fout)//'"'

 comm='bash '//trim(this%search_scr)//trim(sdir)//trim(spref)//trim(sout)

! write(*,*) 'Command=', trim(comm)
 call execute_command_line(trim(comm))

 tf%fname=fout
 call tf%read
 comm='rm '//trim(fout)
 call execute_command_line(trim(comm))

 
 call this%del
 this%nfiles=tf%nlines
 allocate(this%dfile(this%nfiles))
 this%dfile=tf%line

 call tf%del

end subroutine search


!subroutine main_calc
! type(search_data) :: sd
! integer :: i


! sd%search_dir='/home/anton/work/project/halo/chp_prop/traj/traj_17/'
! sd%prefix='traj_part' 

! call sd%search

! do i=1,sd%nfiles
!  write(*,*) trim(sd%dfile(i))
! end do

!end subroutine main_calc



end module search_data_m

!program main
! use search_data_m, only : main_calc

! call main_calc

!end program main
