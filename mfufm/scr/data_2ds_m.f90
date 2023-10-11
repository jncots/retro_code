module data_2ds_m
!====================================================
! Simple data structure for storage of function with
! 2 arguments
!====================================================
 use utools_mod, only : utools
 implicit none
 private
 save
 
 public :: data_2ds
 
 type data_2ds
  integer :: nx, ny
  real(8), allocatable :: xa(:), ya(:), fa(:,:)
  character(500) :: fname='data.dat'
 contains 
  procedure :: del, set, write=>write_tab, read=>read_tab
 end type data_2ds 
 
 
 type(utools) :: ut

 contains 
 

subroutine del(this)
 class(data_2ds) :: this
 if (allocated(this%xa)) deallocate(this%xa)
 if (allocated(this%ya)) deallocate(this%ya)
 if (allocated(this%fa)) deallocate(this%fa)
 this%nx=0
 this%ny=0
end subroutine del

subroutine set(this,nx,ny)
 class(data_2ds) :: this
 integer, intent(in) :: nx, ny
 call this%del
 allocate(this%xa(nx))
 allocate(this%ya(ny))
 allocate(this%fa(nx,ny))
 this%nx=nx
 this%ny=ny
end subroutine set


subroutine write_tab(this)
 class(data_2ds) :: this
 integer :: ndev
 
 if (this%fname=='data.dat') then
  write(*,*) 'data_2ds%write: data_2ds%fname=?, choose the name'
  return
 end if 
 call ut%get_unit(ndev)
 open(ndev,file=this%fname,form='unformatted')
 write(ndev) this%nx
 write(ndev) this%ny
 write(ndev) this%xa
 write(ndev) this%ya
 write(ndev) this%fa
 close(ndev)

end subroutine write_tab
 

subroutine read_tab(this)
 class(data_2ds) :: this
 integer :: ndev
 integer :: nx, ny
 
 if (this%fname=='data.dat') then
  write(*,*) 'data_2ds%read: data_2ds%fname=?, choose the name'
  return
 end if 
 
 call ut%get_unit(ndev)
 open(ndev,file=this%fname,form='unformatted')
 read(ndev) nx
 read(ndev) ny
 call this%set(nx,ny)
 read(ndev) this%xa
 read(ndev) this%ya
 read(ndev) this%fa
 close(ndev)
end subroutine read_tab 
 

end module data_2ds_m