module plot_fun_m
!=================================================
! Calculation of the function on a cartesian
! grid
!=================================================
 use utools_mod, only : utools
 implicit none
 private
 save


 public :: plot_fun
 public :: test_calc


 abstract interface
  subroutine fun_to_plot(x,res)
   real(8), intent(in) :: x(:)
   real(8), intent(out) :: res(:)
  end subroutine
 end interface


 type one_dim
  real(8), allocatable :: x(:)
 contains 
  procedure :: set=>set_one_dim, del=>del_one_dim
 end type one_dim


 type plot_fun
  character(500) :: fname
  integer :: ndout
  type(one_dim), allocatable :: ax(:), fun(:)
  real(8), allocatable :: units(:)
  procedure(fun_to_plot), pointer, nopass :: funx=>null()
 contains
  procedure :: add_dim, calc=>calc_on_grid
  procedure :: ndim=>set_ax, del=>del_ax
 end type plot_fun



 contains


subroutine del_one_dim(this)
 class(one_dim) :: this
 if (allocated(this%x)) deallocate(this%x)
end subroutine del_one_dim


subroutine set_one_dim(this,nx)
 class(one_dim) :: this 
 integer, intent(in) :: nx
 call this%del
 allocate(this%x(nx))
end subroutine set_one_dim


subroutine del_ax(this)
 class(plot_fun) :: this
 if (allocated(this%ax)) deallocate(this%ax)
 if (allocated(this%units)) deallocate(this%units)
end subroutine del_ax


subroutine set_ax(this,nax,ndout)
 class(plot_fun) :: this
 integer, intent(in) :: nax, ndout
 call this%del
 allocate(this%ax(nax))
 allocate(this%units(nax))
 this%ndout=ndout
end subroutine set_ax



subroutine add_dim(this,num,xmin,xmax,nx,scl,units)
 class(plot_fun) :: this
 integer, intent(in) :: num, nx
 real(8), intent(in) :: xmin, xmax
 character(3), intent(in), optional :: scl
 real(8), intent(in), optional :: units
 character(3) :: scl0
 type(utools) :: ut

 if (present(scl)) then
  scl0=scl
 else
  scl0='lin'
 end if

 if (present(units)) then
  this%units(num)=units
 else 
  this%units(num)=1d0
 end if
 
 call ut%grid(this%ax(num)%x,xmin,xmax,nx,scl0)
end subroutine add_dim


subroutine calc_on_grid(this,funx)
 class(plot_fun) :: this
 procedure(fun_to_plot) :: funx
 integer, allocatable :: nx(:), ix(:), nxx(:) 
 integer :: ntot, num, i, i1, i2, j
 real(8), allocatable :: x(:), res(:)
 integer :: nunit

 num=size(this%ax)
 if (allocated(nx)) deallocate(nx)
 allocate(nx(num))

 if (allocated(nxx)) deallocate(nxx)
 allocate(nxx(num))

 if (allocated(ix)) deallocate(ix)
 allocate(ix(num))

 if (allocated(x)) deallocate(x)
 allocate(x(num))

 if (allocated(res)) deallocate(res)
 allocate(res(this%ndout))

 
 ntot=1
 do i=1,num
  nx(i)=size(this%ax(i)%x)
  nxx(num-i+1)=ntot
  ntot=ntot*nx(i)
 end do


 if (allocated(this%fun)) deallocate(this%fun)
 allocate(this%fun(ntot))

 do i=1,ntot
  i1=i
  
  do j=1,num
   i2=mod(i1,nxx(j))
   if (i2==0) then
    ix(num-j+1)=i1/nxx(j)
    i2=nxx(j)
   else
    ix(num-j+1)=i1/nxx(j)+1
   end if
   i1=i2
  end do
 
  do j=1,num
   x(j)=this%ax(j)%x(ix(j)) 
  end do
  call funx(x,res)
  call this%fun(i)%set(this%ndout)
  this%fun(i)%x=res
 end do

 open(newunit=nunit,file=this%fname,form='unformatted',access='stream',status='replace')
 write(nunit) num, this%ndout
 do i=1,num
  write(nunit) nx(i)
 end do 
 write(nunit) ntot
 do i=1,num
  write(nunit) real(this%ax(i)%x/this%units(i)) 
 end do

 do i=1,ntot
  write(nunit) real(this%fun(i)%x)
 end do

end subroutine calc_on_grid


subroutine test_calc
 type(plot_fun) :: pf
 integer :: i


 call pf%ndim(3,3)
 call pf%add_dim(1,1d0,10d0,100,'lin')
 call pf%add_dim(2,1d0,10d0,100,'lin')
 call pf%add_dim(3,1d0,10d0,100,'lin')
 pf%fname='test.dat'
 call pf%calc(mff)

 
end subroutine test_calc


subroutine mff(x,res)
 real(8), intent(in) :: x(:)
 real(8), intent(out) :: res(:)
 
 res=2*x

end subroutine mff



end module plot_fun_m


!program main
! use plot_fun_m, only : test_calc

! call test_calc

!end program main
