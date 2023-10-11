module format_file_tools
 implicit none
 private
 save
 
 public :: format_file_data
 
 
 type format_file_data
  character(500) :: fin
  integer :: nrow, ncol
  real(8), allocatable :: fdata(:,:) ! fdata(nrow,ncol)
 contains
  procedure :: read=>data_from_file, del
 end type format_file_data
 
 
 
 contains


subroutine del(this)
! Variables
 class(format_file_data) :: this
! Calculations
 deallocate(this%fdata)
 this%nrow=0
 this%ncol=0
 this%fin=''
end subroutine del


subroutine data_from_file(this)
!
! Reading numerical data fdat(nrow,ncol) from
! the file fname. The first  lines can contain
! comments starting from '#'.
!
! Variables
 class(format_file_data) :: this
 integer, parameter :: ndev=14
 integer :: nskip, i, j
 character(1000) :: str
! Calculations

 
 call file_row_col(this%fin,this%nrow,this%ncol,nskip)
 
 if (allocated(this%fdata)) deallocate(this%fdata)
 allocate(this%fdata(this%nrow,this%ncol))
 
 open(ndev,file=this%fin)
 
 if (nskip/=0) then
 do i=1,nskip
  read(ndev,'(A)') str
 end do
 end if
 
 do i=1+nskip,this%nrow+nskip
  read(ndev,*) (this%fdata(i-nskip,j),j=1,this%ncol)
 end do
 
 close(ndev)
 

end subroutine data_from_file



subroutine file_row_col(fname,nrow,ncol,nskip)
!
! Calculations of number of rows ('nrow') and
! columns ('ncol') in the file 'fname' containing
! numerical data. The first 'nskip' lines can contain
! comments starting from '#'.
!
! Variables
 character(500), intent(in) :: fname
 integer, intent(out) :: nrow, ncol, nskip
 integer, parameter :: ndev=14
 integer :: istat, lstr, nchar, i
 character(1000) :: str, str1
 character(15), parameter :: numstr='0123456789Ee.-+'
! Calculations

! Number of rows 
 nrow=0
 nskip=0
 open(ndev,file=fname)
 do
  read(ndev,'(A)',iostat=istat) str
  if (istat<0) exit
  str1=adjustl(str)
  if (str1(1:1)=='#') then
   nskip=nskip+1
  else 
   nrow=nrow+1
  end if 
 end do
 close(ndev)
 
! Number of columns 
 lstr=len_trim(str)
 nchar=0
 ncol=0
 
 do i=1,lstr
  if ((index(numstr,str(i:i))/=0).and.(i<lstr)) then
   nchar=1
  else
   if (nchar==1) then
    nchar=0
    ncol=ncol+1
   end if
  end if
 end do

end subroutine file_row_col



end module format_file_tools


!================================================================
! Test program
!================================================================
!program main
! use format_file_tools, only : data_from_file
! use intpol_mod, only : plaw_index
! implicit none
! 
! integer, parameter :: ndev=15 
! integer :: i, j
! integer :: nrow, ncol, icol
! character(100) :: fname, fname_out
! real(8), allocatable :: f_se_pp(:,:), f_se_gam(:,:), f_se_tot(:,:)
!! Calculations 
! 
! fname='syn_totsp3_plt.dat'
! call data_from_file(fname,f_se_pp,nrow,ncol)
! 
! 
! fname='syn_gam_tot3_plt.dat'
! call data_from_file(fname,f_se_gam,nrow,ncol)
! write(*,*) nrow, ncol
! 
! call data_from_file(fname,f_se_tot,nrow,ncol)
! 
! do i=2,ncol
!  f_se_tot(:,i)=f_se_pp(:,i)+f_se_gam(:,i)
! end do
! 
! fname_out='syn_pg_tot3.dat'
! open(ndev,file=fname_out)
! do i=1,nrow
!  write(ndev,'(10Es14.6)') (f_se_tot(i,j),j=1,ncol)
! end do
! close(ndev)
!
!
!
!
!end program main