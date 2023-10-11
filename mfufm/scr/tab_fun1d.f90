module tab_fun1d_m
 implicit none
 private
 save
 
 
 public :: tab_fun1d
 
 type tab_fun1d
  integer :: nx=0
  real(8), allocatable :: xa(:), fa(:)  
  character(500) :: fname='data.dat'
 contains 
  procedure :: set, del, set_xa, read=>read_tab, write=>write_tab
 end type tab_fun1d
 
 
 contains


subroutine del(this)
 class(tab_fun1d) :: this
 
 if (allocated(this%xa)) deallocate(this%xa)
 if (allocated(this%fa)) deallocate(this%fa)
 this%nx=0
end subroutine del 

subroutine set(this,nx)
 class(tab_fun1d) :: this
 integer, intent(in) :: nx
 
 call this%del
 allocate(this%xa(nx))
 allocate(this%fa(nx))
 this%nx=nx
end subroutine set

subroutine set_xa(this,x1,x2)
 class(tab_fun1d) :: this
 real(8), intent(in) :: x1, x2

 if (this%nx==0) then
  write(*,'(A)') 'tab_fun1d: set_xa: array is not allocated'
  return
 end if
 
 call set_arr_log(x1,x2,this%nx,this%xa)
end subroutine set_xa


subroutine set_arr_log(a1,a2,na,aa)
 real(8), intent(in) :: a1, a2
 integer, intent(in) :: na
 real(8) :: aa(na)
 integer :: i
 real(8) :: da
 
 da=(a2/a1)**(1d0/(na-1))
 
 do i=0,na-2
  aa(i+1)=a1*da**i
 end do
 aa(na)=a2
 
end subroutine set_arr_log


subroutine write_tab(this)
 class(tab_fun1d) :: this
 integer, parameter :: ndev=15
 
 if (this%fname=='data.dat') then
  write(*,*) 'write: choose the name for the object of "tab_fun1d" type'
  return
 end if 
 
 open(ndev,file=this%fname,form='unformatted')
 write(ndev) this%nx
 write(ndev) this%xa
 write(ndev) this%fa
 close(ndev)

end subroutine write_tab



subroutine read_tab(this)
 class(tab_fun1d) :: this
 integer, parameter :: ndev=15
 integer :: nx
 
 if (this%fname=='data.dat') then
  write(*,*) 'read: choose the name for the object of "tab_fun1d" type'
  return
 end if 
 
 open(ndev,file=this%fname,form='unformatted')
 read(ndev) nx
 call this%set(nx)
 read(ndev) this%xa
 read(ndev) this%fa
 close(ndev)
 
end subroutine read_tab


end module tab_fun1d_m



! program main
 ! use tab_fun1d_m, only : tab_fun1d
 ! type(tab_fun1d) :: pp
 ! integer :: i
 
 ! pp%fname='test.dat'
 ! call pp%set(100)
 ! call pp%set_xa(1d0,100d0)
 ! pp%fa=2*pp%xa
 ! call pp%write
 ! call pp%del
 
 ! write(*,*) pp%nx, trim(pp%fname)
 ! read(*,*)
 ! call pp%read
 
 ! do i=1,pp%nx
  ! write(*,'(I5,10Es14.6)') i, pp%xa(i), pp%fa(i)
 ! end do
 
 
! end program main