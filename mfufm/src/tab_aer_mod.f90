module arr1d_ipol
 use intpol_mod, only : arr_ind_short
 implicit none
 private
 save
 
 public :: arr_1d
 
 type arr_1d
  integer :: nx, nx0 ! size of x and x0
  integer :: i=0, n1=0, n2=0
  integer :: er=0 ! er=1 is OK, er=0 is not OK
  real(8) :: v1, v2, v 
  real(8), pointer :: x(:), x0(:) ! x0 is old array, in which we are searching place (n1,n2) for x(i)
 contains
  procedure :: set, find=>wfind
 end type arr_1d
 
 contains

subroutine set(this,x,x0)
 class(arr_1d) :: this
 real(8), intent(in), target :: x(:), x0(:)
 this%nx=size(x)
 this%nx0=size(x0)
 this%x=>x
 this%x0=>x0
 this%er=0
 this%i=0
 this%v1=0d0
 this%v2=0d0
 this%v=0d0
end subroutine set

subroutine wfind(this,i)
 class(arr_1d) :: this
 integer, intent(in) :: i
 real(8) :: x
 
 if (this%i/=i) then
  this%i=i
  this%er=1  
  if ((this%i<1).or.(this%i>this%nx)) then
   this%er=0
   return
  end if
  
  x=this%x(this%i)
  if ((x<this%x0(1)).or.(x>this%x0(this%nx0))) then
   this%er=0
   return
  end if
  call arr_ind_short(1,this%nx0,this%x0,x,this%n1,this%n2)
 end if 
end subroutine wfind


end module arr1d_ipol



! program main
 ! use arr1d_ipol, only : arr_1d
 ! type(arr_1d) :: aa
 ! real(8), allocatable :: a1(:), a2(:)
 ! real(8) :: x1, x2, x
 ! integer :: i, nx
 
 ! x1=1d0
 ! x2=1d2
 ! nx=100
 ! allocate(a1(nx))
 
 ! do i=0,nx-1
  ! a1(i+1)=x1*(x2/x1)**(i*1d0/(nx-1))
 ! end do
 
 ! x1=1d1
 ! x2=1d2
 ! nx=50
 ! allocate(a2(nx))
 
 ! do i=0,nx-1
  ! a2(i+1)=x1*(x2/x1)**(i*1d0/(nx-1))
 ! end do
 
 ! call aa%set(a2,a1)
 ! do i=0,60
  ! call aa%find(i)
  ! write(*,*) aa%er, aa%n1, aa%n2, a2(i), a1(aa%n1), a1(aa%n2)
 ! end do

 
! end program main



module tab_aer_mod
!==========================================================================
! Describes object 'tab_aer' of 3 dimentional data
! see example program below
!==========================================================================
 use arr1d_ipol, only : arr_1d
 use intpol_mod, only : logl_int
 implicit none
 private
 save
 
 
 public :: tab_aer
 
 
 type tab_aer
  integer :: na=0, ne=0, nr=0
  real(8), allocatable :: aa(:), ea(:), ra(:), fa(:,:,:)
  character(500) :: fname='data.dat'
 contains
  procedure :: del, set, set_ea, set_ra, set_aa, read=>read_tab, write=>write_tab
  procedure :: copy_from, map=>map_that
 end type tab_aer
 
 
 contains
 
subroutine del(this)
 class(tab_aer) :: this
 
 if (allocated(this%aa)) deallocate(this%aa)
 if (allocated(this%ea)) deallocate(this%ea)
 if (allocated(this%ra)) deallocate(this%ra)
 if (allocated(this%fa)) deallocate(this%fa)
 
 this%na=0
 this%ne=0
 this%nr=0
 
end subroutine del


subroutine set(this,na,ne,nr)
 class(tab_aer) :: this
 integer, intent(in) :: na, ne, nr
 
 call this%del
 
 allocate(this%aa(na))
 allocate(this%ea(ne))
 allocate(this%ra(nr))
 allocate(this%fa(na,ne,nr))
 this%na=na
 this%ne=ne
 this%nr=nr

end subroutine set

subroutine set_from(this,that)
 class(tab_aer) :: this
 type(tab_aer), intent(in) :: that
 
 call this%set(that%na,that%ne,that%nr)

end subroutine set_from


subroutine copy_from(this,that)
 class(tab_aer) :: this
 type(tab_aer), intent(in) :: that
 
 call this%set(that%na,that%ne,that%nr)
 this%aa=that%aa
 this%ea=that%ea
 this%ra=that%ra
 this%fa=that%fa
end subroutine copy_from



subroutine set_aa(this,x1,x2)
 class(tab_aer) :: this 
 real(8), intent(in) :: x1, x2 ! is angles in degrees
 real(8), allocatable :: ang(:)
 real(8), parameter :: pi=3.14159265358979324d0 ! pi number 
 integer :: i
 
 if (this%na==0) then
  write(*,'(A)') 'set_aa: array is not allocated'
  return
 end if
 
 allocate(ang(this%na))
 call set_arr_log(x1,x2,this%na,ang)
 do i=1,this%na
  this%aa(i)=1-cos(ang(i)*pi/180)
 end do
 deallocate(ang)
 
end subroutine set_aa

subroutine set_ea(this,x1,x2)
 class(tab_aer) :: this 
 real(8), intent(in) :: x1, x2
 
 if (this%ne==0) then
  write(*,'(A)') 'set_ea: array is not allocated'
  return
 end if
 
 call set_arr_log(x1,x2,this%ne,this%ea)
end subroutine set_ea

subroutine set_ra(this,x1,x2)
 class(tab_aer) :: this
 real(8), intent(in) :: x1, x2

 if (this%nr==0) then
  write(*,'(A)') 'set_ra: array is not allocated'
  return
 end if
 
 call set_arr_log(x1,x2,this%nr,this%ra)
end subroutine set_ra


subroutine set_arr_lin(a1,a2,na,aa)
 real(8), intent(in) :: a1, a2
 integer, intent(in) :: na
 real(8) :: aa(na)
 integer :: i
 real(8) :: da
 
 da=(a2-a1)/(na-1)
 do i=0,na-1
  aa(i+1)=a1+da*i
 end do
 
end subroutine set_arr_lin

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
 class(tab_aer) :: this
 integer, parameter :: ndev=15
 
 if (this%fname=='data.dat') then
  write(*,*) 'write_tab: choose the name for the object of "tab_aer" type'
  return
 end if 
 
 open(ndev,file=this%fname,form='unformatted')
 write(ndev) this%na
 write(ndev) this%ne
 write(ndev) this%nr
 write(ndev) this%aa
 write(ndev) this%ea
 write(ndev) this%ra
 write(ndev) this%fa
 close(ndev)

end subroutine write_tab



subroutine read_tab(this)
 class(tab_aer) :: this
 integer, parameter :: ndev=15
 integer :: na, ne, nr
 
 if (this%fname=='data.dat') then
  write(*,*) 'read_tab: choose the name for the object of "tab_aer" type'
  return
 end if 
 
 open(ndev,file=this%fname,form='unformatted')
 read(ndev) na
 read(ndev) ne
 read(ndev) nr
 call set(this,na,ne,nr)
 read(ndev) this%aa
 read(ndev) this%ea
 read(ndev) this%ra
 read(ndev) this%fa
 close(ndev)
 
end subroutine read_tab




subroutine map_that(this,that,order)
!=================================================================
! mapping that to this data using logarithmic interpolation
! (logarithmic everywhere, except points with value=0)
! order is so, that the interpolation starts with the first index
!=================================================================
 class(tab_aer), target :: this
 type(tab_aer), intent(in), target :: that
 integer, optional, intent(in) :: order(3)
 type(arr_1d) :: sa(3)
! order is so, that the interpolation starts with the first index
 integer, parameter :: order_def(3)=[1,2,3]
 real(8), pointer :: af(:)
 integer :: i, i1, i2, i3, j2, j3, k1, k2
 integer :: ii(3), nt(3), ti(3), r1(3), r2(3)
 integer :: q(3,2)
 real(8) :: res2(2,2), res1(2), res
 real(8), pointer :: fr(:,:,:), f0(:,:,:)  
  
! Setting data  
 call sa(1)%set(this%aa,that%aa)
 call sa(2)%set(this%ea,that%ea)
 call sa(3)%set(this%ra,that%ra)
 nt(1)=this%na
 nt(2)=this%ne
 nt(3)=this%nr  
 fr=>this%fa
 f0=>that%fa
  
! Checking order 
 if (present(order)) then
  do i=1,3
   if ((order(i)<1).or.(order(i)>3)) then
    write(*,*) 'Order is given in a wrong format:'
	write(*,*) 'order(',i,')=',order(i)
	ii=order_def
	exit
   end if
   i1=i+1
   if (i>2) i1=1
   if (order(i)==order(i1)) then
    write(*,*) 'Order is given in a wrong format:'
	write(*,*) 'order(',i,'=order(',i,')=',order(i)
    ii=order_def
	exit
   end if
  end do
  ii=order
 else
  ii=order_def
 end if

! Change order, thus the first one is the innermost one.
 i=ii(1)
 ii(1)=ii(3)
 ii(3)=i 
  
  do i1=1,nt(ii(1))
   call sa(ii(1))%find(i1)
   if (sa(ii(1))%er==0) then
	ti(ii(1))=i1
	do j2=1,nt(ii(2))
	 ti(ii(2))=j2
	 do j3=1,nt(ii(3))
	  ti(ii(3))=j3
	  fr(ti(1),ti(2),ti(3))=0d0
	 end do
    end do
   	cycle 
   end if
   q(ii(1),1)=sa(ii(1))%n1
   q(ii(1),2)=sa(ii(1))%n2
   
   do i2=1,nt(ii(2))
    call sa(ii(2))%find(i2)	
	if (sa(ii(2))%er==0) then
	 ti(ii(1))=i1
	 ti(ii(2))=i2
	 do j3=1,nt(ii(3))
	  ti(ii(3))=j3
	  fr(ti(1),ti(2),ti(3))=0d0
	 end do
     cycle	 
    end if
    q(ii(2),1)=sa(ii(2))%n1
    q(ii(2),2)=sa(ii(2))%n2
	
	
	
	do i3=1,nt(ii(3))
	 call sa(ii(3))%find(i3)	 
	 ti(ii(1))=i1
	 ti(ii(2))=i2
	 ti(ii(3))=i3
	  
	 if (sa(ii(3))%er==0) then
	  fr(ti(1),ti(2),ti(3))=0d0
	  cycle
	 end if 
	 q(ii(3),1)=sa(ii(3))%n1
	 q(ii(3),2)=sa(ii(3))%n2
	 
     do k1=1,2
      do k2=1,2	 
	   r1(ii(1))=q(ii(1),k1)
	   r2(ii(1))=q(ii(1),k1)
	   r1(ii(2))=q(ii(2),k2)
	   r2(ii(2))=q(ii(2),k2)
	   r1(ii(3))=q(ii(3),1)
	   r2(ii(3))=q(ii(3),2)
	   
	   call logl_int(sa(ii(3))%v1,sa(ii(3))%v2,&
	   f0(r1(1),r1(2),r1(3)),f0(r2(1),r2(2),r2(3)),sa(ii(3))%v,res2(k1,k2))
	  end do
	 end do 

     do k1=1,2
      call logl_int(sa(ii(2))%v1,sa(ii(2))%v2,&
	  res2(k1,1),res2(k1,2),sa(ii(2))%v,res1(k1))
	 end do 
	 
	 call logl_int(sa(ii(1))%v1,sa(ii(1))%v2,&
	  res1(1),res1(2),sa(ii(1))%v,res) 
		 
	 fr(ti(1),ti(2),ti(3))=res
	  
	end do
   end do
  end do   
  
  fr=>null()
  f0=>null()
  
end subroutine map_that


end module tab_aer_mod


! program main
 ! use tab_aer_mod, only : tab_aer
 ! type(tab_aer) :: f1, f2
 ! integer :: n1, n2, n3
 
 ! n1=10
 ! n2=20
 ! n3=30
 ! call f1%set(n1,n2,n3)
 ! call f1%set_aa(1d1,1d2)
 ! call f1%set_ea(1d8,1d12)
 ! call f1%set_ra(1d-2,1d3)
 ! f1%fa=10d0
 
 
 ! n1=5
 ! n2=3
 ! n3=15
 ! call f2%set(n1,n2,n3)
 ! call f2%set_aa(1d1,1d2)
 ! call f2%set_ea(1d2,1d12)
 ! call f2%set_ra(1d-2,1d3)
 ! call f2%map(f1)
 
 ! write(*,*) f2%fa(1,2,3)
 
 
 
! end program main

