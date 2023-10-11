module fifo_rec_m
!=================================================
! Realisation of FIFO (first in, first out)
! data structure
!=================================================
 implicit none
 private
 save


 public :: fifo_rec
! public :: test_calc


 type fifo_rec
  integer :: nsys
  integer :: ntot, n1, n2
  logical :: empty=.true.
  integer :: def_ntot=1000
  real(8), allocatable :: y(:,:)
 contains
  procedure :: del=>del_fifo, set=>set_fifo
  procedure :: set_ntot, add=>add_fifo, arrange
 end type fifo_rec


 contains

subroutine del_fifo(this)
 class(fifo_rec) :: this
 if (allocated(this%y)) deallocate(this%y)
 this%nsys=0
 this%ntot=0
 this%n1=0
 this%n2=0
 this%empty=.true.
end subroutine del_fifo


subroutine set_ntot(this,ntot)
 class(fifo_rec) :: this
 integer, intent(in) :: ntot
 
 this%def_ntot=ntot

end subroutine set_ntot


subroutine set_fifo(this,nsys,ntot)
 class(fifo_rec) :: this
 integer, intent(in) :: nsys, ntot

 call this%del
 allocate(this%y(0:nsys,ntot))
 this%nsys=nsys
 this%ntot=ntot
 this%n1=1
 this%n2=1
 this%empty=.false.

end subroutine set_fifo


subroutine add_fifo(this,xc,yc)
 class(fifo_rec) :: this
 real(8), intent(in) :: xc, yc(:)

 if (this%empty) then
  call this%set(size(yc),this%def_ntot)
  this%y(0,this%n2)=xc
  this%y(1:,this%n2)=yc
  return
 end if
 
 call move_pt(this%ntot,this%n1,this%n2)

 this%y(0,this%n2)=xc
 this%y(1:,this%n2)=yc

end subroutine add_fifo


subroutine move_pt(ntot,n1,n2)
 integer, intent(in) :: ntot
 integer :: n1, n2

 n2=n2+1
 if (n2>ntot) n2=1
 if (n1==n2) n1=n1+1
 if (n1>ntot) n1=1

end subroutine move_pt

subroutine arrange(this)
 class(fifo_rec) :: this
 type(fifo_rec) :: that
 integer :: nc
 nc=this%n1

 call that%set_ntot(this%def_ntot)
 do
  call that%add(this%y(0,nc),this%y(1:,nc))
  if (nc==this%n2) exit
  nc=nc+1
  if (nc>this%ntot) nc=1 
 end do

 call this%del
 call this%set_ntot(that%n2)
 
 do nc=1,that%n2
  call this%add(that%y(0,nc),that%y(1:,nc))
 end do 

 call that%del

end subroutine arrange


!subroutine test_calc
! type(fifo_rec) :: arr
! real(8) :: x, y(1)
! integer :: i


! call arr%set_ntot(10)


! do i=1,100
!  x=i
!  y=x
!  call arr%add(x,y)
! end do

! call arr%arrange
! 
! do i=1,arr%ntot
!  write(*,*) i, arr%y(:,i)
! end do

!end subroutine test_calc




end module fifo_rec_m


!program main
! use fifo_rec_m, only : test_calc

! call test_calc

!end program main

