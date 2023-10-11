module cont_fun1d_m
 implicit none
 private
 save
 
 public :: cont_fun1d
 
 abstract interface
  subroutine fun1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type cont_fun1d
  real(8) :: x1, x2
  procedure(fun1d), pointer, nopass :: fx=>null()
 contains
  procedure :: set, del 
 end type cont_fun1d

 contains
 
subroutine set(this,fx,x1,x2)
 class(cont_fun1d) :: this
 procedure(fun1d) :: fx
 real(8), intent(in) :: x1, x2
 this%fx=>fx
 this%x1=x1
 this%x2=x2
end subroutine set

subroutine del(this)
 class(cont_fun1d) :: this
 this%fx=>null()
 this%x1=0d0
 this%x2=0d0
end subroutine del

end module cont_fun1d_m