module file_manager_m
!=====================================================
! Module contains functions for manipulations
! with files in the system
!=====================================================
 implicit none
 private
 save

 public :: file_manager !, test_calc

 type file_manager
 contains
  procedure, nopass :: numbered_fname
  procedure, nopass :: save_old_data
  procedure, nopass :: copy_to, create_dir
 end type file_manager



 contains

function numbered_fname(dir,prefix,num)
 character(500), intent(in) :: dir, prefix
 integer, intent(in) :: num
 character(500) :: fname, numbered_fname

 write(fname,'(I20)') num
 fname=trim(prefix)//trim(adjustl(fname))//'.dat'
 fname=trim(dir)//trim(fname)
 numbered_fname=fname

end function numbered_fname



subroutine save_old_data(ftypes,dir,sdir)
!===========================================================
! Copy old data, if exists, with extensions 
! given in ftypes (e.g. ftypes='*.f90 *.dat *.sh')
! from the directory dir to numberated subdirectory sdir 
!===========================================================
 character(500), intent(in) :: ftypes, dir, sdir
 character(500) :: dir1
 character(3000) :: comm

 dir1=dir(:len_trim(dir)-1)
 comm="#! /bin/bash\n"
 comm=trim(comm)//"dir_check='"//trim(dir1)//"'\n"
 comm=trim(comm)//"dir_try='"//trim(sdir)//"'\n"
 comm=trim(comm)//'if [ -d "$dir_check" ]; then\n'
 comm=trim(comm)//'if ls "$dir_check"/'//trim(ftypes)//' 1> /dev/null 2>&1; then # check if any *.dat files exist\n'
 comm=trim(comm)//'cdir0=$(pwd)\n'
 comm=trim(comm)//'cdir="$dir_check"\n'
 comm=trim(comm)//'cd $cdir\n'
 comm=trim(comm)//'isw=0\n'
 comm=trim(comm)//'i=1\n'
 comm=trim(comm)//'while [ "$isw" -lt "1" ]; do\n'
 comm=trim(comm)//'ndir="$dir_try$i"\n'
 comm=trim(comm)//'if [ ! -d "$ndir" ]; then\n'
 comm=trim(comm)//'mkdir $ndir\n'
 comm=trim(comm)//"mv  "//trim(ftypes)//"  $ndir 2>/dev/null\n"
 comm=trim(comm)//'isw=1\n'
 comm=trim(comm)//'fi\n'
 comm=trim(comm)//'i=$(( i + 1))\n'
 comm=trim(comm)//'done\n'
 comm=trim(comm)//'cd $cdir0\n'
 comm=trim(comm)//'fi\n'
 comm=trim(comm)//'fi\n'

 call execute_command_line(trim(comm))

end subroutine save_old_data


subroutine create_dir(dir)
 character(500), intent(in) :: dir
 call execute_command_line('mkdir -p '//trim(dir))
end subroutine create_dir

subroutine copy_to(files,dir)
 character(500), intent(in) :: dir, files
 call execute_command_line('mkdir -p '//trim(dir))
 call execute_command_line('cp '//trim(files)//' '//trim(dir))
end subroutine copy_to


subroutine test_name
 type(file_manager) :: fm
 character(500) :: dir, prefix, ftypes, sdir, dir1, files
 integer :: num, i

 dir='/home/anton/work/project/halo/chp_prop/in_jf_field/reg_cent/parall_launch/node_test/'
 prefix='file_'
 
 do i=10000,10005
  write(*,*) trim(fm%numbered_fname(dir,prefix,i))
 end do

! ftypes='*.svg'
! sdir='save'
! call fm%save_old_data(ftypes,dir,sdir)

! dir1=trim(dir)//'other'
! files='*.f90 *.sh'
! call fm%copy_to(files,dir1)
 

end subroutine test_name


subroutine test_calc

 call test_name

end subroutine test_calc



end module file_manager_m


!program main
! use file_manager_m, only : test_calc

! call test_calc

!end program main
