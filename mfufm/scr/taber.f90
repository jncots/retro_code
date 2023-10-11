module taber
 use intpol_mod, only : arr_ind_short, logl_int
 implicit none
 private
 save
 
 
 public :: tab_er
 
 
 type tab_er
  integer :: ne=0, nr=0
  real(8), allocatable :: ea(:), ra(:), fa(:,:)
  character(500) :: fname='data.dat'
 contains
  procedure :: set, del, set_ea, set_ra, read=>read_tab, write=>write_tab
  procedure :: copy_from, add
 end type tab_er

 contains

 
subroutine add(this,that)
  class(tab_er) :: this
  type(tab_er), intent(in) :: that
  integer :: ir, ie, nr1, nr2, ne1, ne2
  real(8) :: res1, res2, res
  
  do ir=1,this%nr
   if ((this%ra(ir)>=that%ra(1)).and.(this%ra(ir)<=that%ra(that%nr))) then
    call arr_ind_short(1,that%nr,that%ra,this%ra(ir),nr1,nr2)
    do ie=1,this%ne
     if ((this%ea(ie)>=that%ea(1)).and.(this%ea(ie)<=that%ea(that%ne))) then 
      call arr_ind_short(1,that%ne,that%ea,this%ea(ie),ne1,ne2)
	  call logl_int(that%ea(ne1),that%ea(ne2),&
	  that%fa(ne1,nr1),that%fa(ne2,nr1),this%ea(ie),res1)
	  call logl_int(that%ea(ne1),that%ea(ne2),&
	  that%fa(ne1,nr2),that%fa(ne2,nr2),this%ea(ie),res2)
	  call logl_int(that%ra(nr1),that%ra(nr2),&
	  res1,res2,this%ra(ir),res)
	  this%fa(ie,ir)=this%fa(ie,ir)+res
     end if
    end do
   end if
  end do
end subroutine add
 
 
 
subroutine copy_from(this,that) 
  class(tab_er) :: this
  type(tab_er), intent(in) :: that 
  
  call this%set(that%ne,that%nr)
  this%ea=that%ea
  this%ra=that%ra
  this%fa=that%fa
  this%fname=that%fname
end subroutine copy_from 

subroutine read_tab(this)
 class(tab_er) :: this
 integer, parameter :: ndev=15
 integer :: ne, nr
 
 if (this%fname=='data.dat') then
  write(*,*) 'read_tab: choose the name for the object of "tab_er" type'
  return
 end if 
 
 open(ndev,file=this%fname,form='unformatted')
 read(ndev) ne
 read(ndev) nr
 call set(this,ne,nr)
 read(ndev) this%ea
 read(ndev) this%ra
 read(ndev) this%fa
 close(ndev)

end subroutine read_tab


subroutine write_tab(this)
 class(tab_er) :: this
 integer, parameter :: ndev=15
 
 if (this%fname=='data.dat') then
  write(*,*) 'write_tab: choose the name for the object of "tab_er" type'
  return
 end if 
 
 open(ndev,file=this%fname,form='unformatted')
 write(ndev) this%ne
 write(ndev) this%nr
 write(ndev) this%ea
 write(ndev) this%ra
 write(ndev) this%fa
 close(ndev)

end subroutine write_tab



subroutine set(this,ne,nr)
 class(tab_er) :: this
 integer, intent(in) :: ne, nr
 
 if (allocated(this%ea)) deallocate(this%ea)
 if (allocated(this%ra)) deallocate(this%ra)
 if (allocated(this%fa)) deallocate(this%fa)
 
 this%ne=ne
 this%nr=nr
 
 allocate(this%ea(ne))
 allocate(this%ra(nr))
 allocate(this%fa(ne,nr))
end subroutine set
 
subroutine del(this)
 class(tab_er) :: this
 
 if (allocated(this%ea)) deallocate(this%ea)
 if (allocated(this%ra)) deallocate(this%ra)
 if (allocated(this%fa)) deallocate(this%fa)
 
 this%ne=0
 this%nr=0
end subroutine del

subroutine set_ea(this,e1,e2)
 class(tab_er) :: this
 real(8), intent(in) :: e1, e2
 integer :: i
 real(8) :: e21, ne

 if (this%ne==0) then
  write(*,'(A)') 'set_ea: array is not allocated'
  return
 end if
 
 e21=e2/e1
 ne=1d0/(this%ne-1)
 
 do i=0,this%ne-1
  this%ea(i+1)=e1*e21**(i*ne)
 end do
 
end subroutine set_ea

subroutine set_ra(this,r1,r2)
 class(tab_er) :: this
 real(8), intent(in) :: r1, r2
 integer :: i
 real(8) :: r21, nr

 if (this%nr==0) then
  write(*,'(A)') 'set_ra: array is not allocated'
  return
 end if
 
 r21=r2/r1
 nr=1d0/(this%nr-1)
 
 do i=0,this%nr-1
  this%ra(i+1)=r1*r21**(i*nr)
 end do
 
end subroutine set_ra
 
 

end module taber


!program main
! use taber, only : tab_er
! type(tab_er) :: gden
! 
! gden%fname='gden.dat'
! call gden%set(100,100)
! call gden%set_ea(1d8,1d16)
! call gden%set_ra(1d-2,1d2)
! 
! do i=1,gden%nr
!  do j=1,gden%ne
!   gden%fa(j,i)=gden%ea(j)*gden%ra(i)
!  end do
! end do
! 
! call gden%write
! 
! call gden%del
! 
! call gden%read
! 
! do i=1,gden%nr
!  do j=1,gden%ne
!   write(*,'(2I5,10Es14.6)') i, j, gden%fa(j,i), gden%ea(j), gden%ra(i)
!  end do
!  read(*,*)
! end do
!
!end program main



