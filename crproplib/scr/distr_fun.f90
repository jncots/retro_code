module distr_fun_m
 use sort_methods_m, only : sort_methods
 use utools_mod, only : utools
 implicit none
 private
 save

 public :: distr_fun
! public :: test_calc

 type distr_fun
  integer :: narr=0, np=0
  real(8), allocatable :: dat(:)
  character(500) :: fname='noname'
 contains
  procedure :: del=>del_df, set=>set_df
  procedure :: double=>double_df, add=>add_df
  procedure :: sort=>sort_df
  procedure :: hist=>hist_df, norm_hist=>nhist_df
  procedure :: aver_val, var_val
  procedure :: dfun=>dfun_df
  procedure :: print=>print_res
 end type distr_fun


 contains

subroutine del_df(this)
 class(distr_fun) :: this
 this%narr=0
 this%np=0
 if (allocated(this%dat)) deallocate(this%dat)
 this%fname='noname'
end subroutine del_df

subroutine set_df(this,narr)
 class(distr_fun) :: this
 integer, intent(in), optional :: narr
 integer :: narr0

 if (present(narr)) then
  narr0=narr
 else
  narr0=100
 end if
 
 call this%del
 allocate(this%dat(narr0))
 this%narr=narr0

end subroutine set_df

subroutine double_df(this)
 class(distr_fun) :: this
 type(distr_fun) :: that
 integer :: i, narr2

 narr2=2*this%narr

 call that%set(this%narr)
 that%np=this%np

 do i=1,this%np
  that%dat(i)=this%dat(i)
 end do

 call this%set(narr2)
 this%np=that%np

 do i=1,that%np
  this%dat(i)=that%dat(i)
 end do
 
 call that%del

end subroutine double_df

subroutine add_df(this,pt)
 class(distr_fun) :: this
 real(8) :: pt
 integer :: nss

 if (this%narr==0) call this%set
 nss=this%np+1
 if (nss>this%narr) call this%double

 this%np=nss
 this%dat(nss)=pt

end subroutine add_df


subroutine sort_df(this)
 class(distr_fun) :: this
 type(sort_methods) :: sm

 call sm%merge_sort(this%dat(1:this%np))

end subroutine sort_df


subroutine hist_df(this,nbin,dist,scl)
 class(distr_fun) :: this
 integer, intent(in) :: nbin
 real(8), allocatable :: dist(:,:)
 character(3), intent(in) :: scl
 type(utools) :: ut
 real(8), allocatable :: x(:)
 integer :: i, j

 call this%sort
 if (allocated(dist)) deallocate(dist)
 allocate(dist(2,nbin+1))

 call ut%grid(x,this%dat(1),this%dat(this%np),nbin+1,scl)
 
 dist(1,:)=x
 dist(2,:)=0d0

 j=1
 do i=2,nbin+1  
  
  do
   if ((this%dat(j)<=x(i)).and.(j<=this%np)) then
    dist(2,i)=dist(2,i)+1
    j=j+1
   else
    exit
   end if

  end do
 end do
 
 deallocate(x)

end subroutine hist_df

subroutine nhist_df(this,nbin,dist,scl)
 class(distr_fun) :: this
 integer, intent(in) :: nbin
 real(8), allocatable :: dist(:,:)
 character(3), intent(in) :: scl


 call this%hist(nbin,dist,scl)
 dist(2,:)=dist(2,:)/this%np

end subroutine nhist_df


subroutine dfun_df(this,nbin,dist,scl)
 class(distr_fun) :: this
 integer, intent(in) :: nbin
 real(8), allocatable :: dist(:,:)
 character(3), intent(in) :: scl
 integer :: i

 call this%norm_hist(nbin,dist,scl)

 do i=2,nbin+1
  dist(2,i)=dist(2,i)/(dist(1,i)-dist(1,i-1))
 end do

end subroutine dfun_df

function aver_val(this)
 class(distr_fun) :: this
 real(8) :: aver_val
 integer :: i 

 aver_val=0d0

 do i=1,this%np
  aver_val=aver_val+this%dat(i)
 end do

 aver_val=aver_val/this%np

end function aver_val

function var_val(this)
 class(distr_fun) :: this
 real(8) :: var_val, av
 integer :: i

 av=this%aver_val()
 var_val=0d0

 do i=1,this%np
  var_val=var_val+(this%dat(i)-av)**2
 end do 

 var_val=sqrt(var_val/this%np)

end function var_val


subroutine print_res(this,nbin,scl,zero,icase,average)
 class(distr_fun) :: this
 integer, intent(in) :: nbin
 character(3), intent(in), optional :: scl
 real(8), intent(in), optional :: zero
 integer, intent(in), optional :: icase
 logical, intent(in), optional :: average
 real(8), allocatable :: dist1(:,:), dist2(:,:), dist3(:,:)
 integer :: nunit, i
 character(3) :: scl0  
 real(8) :: zero0, min2, min3
 integer:: icase0
 logical :: average0

 if (this%fname=='noname') then
  write(*,*) 'distr_fun: print_res: The name of the output file is not given'
  write(*,*) 'distr_fun: print_res: Please set it as this%fname=fname'
  return
 end if

 if (present(scl)) then
  scl0=scl
 else
  scl0='lin'
 end if

 if (present(zero)) then
  zero0=zero
 else
  zero0=0d0
 end if

 if (present(icase)) then
  icase0=icase
 else
  icase0=3
 end if

 if (present(average)) then
  average0=average
 else
  average0=.true.
 end if


 open(newunit=nunit,file=this%fname)

 select case (icase0)
  case(1)
   call this%hist(nbin,dist1,scl0)
   do i=1,size(dist1,2)
    if (dist1(2,i)<1d0) dist1(2,i)=zero0
    write(nunit,*) dist1(1,i), dist1(2,i)
   end do

  case(2)
   call this%hist(nbin,dist1,scl0)
   call this%dfun(nbin,dist2,scl0)
   min2=min_gt0(dist2(2,:))

   do i=1,size(dist1,2)
    if (dist1(2,i)<1d0) dist1(2,i)=zero0
    if (dist2(2,i)<min2) dist2(2,i)=zero0*min2
    write(nunit,*) dist1(1,i), dist1(2,i), dist2(2,i)
   end do
  case(3)
   call this%hist(nbin,dist1,scl0)
   call this%dfun(nbin,dist2,scl0)
   call this%norm_hist(nbin,dist3,scl0)
   min2=min_gt0(dist2(2,:))
   min3=min_gt0(dist3(2,:))
   do i=1,size(dist1,2)
    if (dist1(2,i)<1d0) dist1(2,i)=zero0
    if (dist2(2,i)<min2) dist2(2,i)=zero0*min2
    if (dist3(2,i)<min3) dist3(2,i)=zero0*min3
    write(nunit,*) dist1(1,i), dist1(2,i), dist2(2,i), dist3(2,i)
   end do
  case default
   write(*,*) 'distr_fun: print_res: There is only 3 levels of output'
   write(*,*) 'distr_fun: print_res: You gave icase=', icase
   write(*,*) 'distr_fun: print_res: Output is for icase=3'
   call this%hist(nbin,dist1,scl0)
   call this%dfun(nbin,dist2,scl0)
   call this%norm_hist(nbin,dist3,scl0)
   min2=min_gt0(dist2(2,:))
   min3=min_gt0(dist3(2,:))

   do i=1,size(dist1,2)
    if (dist1(2,i)<1d0) dist1(2,i)=zero0
    if (dist2(2,i)<min2) dist2(2,i)=zero0*min2
    if (dist3(2,i)<min3) dist3(2,i)=zero0*min3
    write(nunit,*) dist1(1,i), dist1(2,i), dist2(2,i), dist3(2,i)
   end do
 end select


 if (average0) then
  write(nunit,'(A,10I10)')    '# npoints=', this%np
  write(nunit,'(A,10Es14.6)') '# average=', this%aver_val() 
  write(nunit,'(A,10Es14.6)') '# sigma=  ', this%var_val()
 end if
 close(nunit)

 if (allocated(dist1)) deallocate(dist1)
 if (allocated(dist2)) deallocate(dist2)
 if (allocated(dist3)) deallocate(dist3)

end subroutine print_res



function min_gt0(arr)
 real(8), intent(in) :: arr(:)
 real(8) :: min_gt0, res, val 
 integer :: i

 res=maxval(arr)
 
 do i=1,size(arr)
  val=arr(i)
  if ((val>0d0).and.(val<res)) res=val
 end do
 min_gt0=res
end function min_gt0


!subroutine test_calc
! use random_gen_m, only : init_random_seed
! type(sort_methods) :: sm
! type(distr_fun) :: df
! integer :: i, nt
! real(8) :: a, dsum
! real(8), allocatable :: dist(:,:)
! 
! call init_random_seed()

! 
! nt=21570
! do i=1,nt
!  call random_number(a)
!!  a=2*a
!  call df%add(a)
! end do


! df%fname='print_hist.dat'
! call df%print(100,zero=0.5d0,scl='lin',icase=3)

!! call df%norm_hist(100,dist,'lin')

! call df%dfun(4,dist,'lin')
! open(1,file='hist_df.dat')

! dsum=0d0
! do i=1,size(dist,2)

!  write(1,'(10Es14.6)') dist(:,i)
!  dsum=dsum+int(dist(2,i))
!!  write(1,'(I10,10Es14.6)') i, df%dat(i)
! end do
!  write(*,*) 'dsum=', int(dsum)
!  write(*,*) df%aver_val(), df%var_val(), 1d0/sqrt(12d0)
! close(1)




!! call sm%merge_sort(df%dat(1:df%np),ist=50,ifin=1500)
! 
!! open(1,file='test_df.dat')

!! do i=1,df%np

!!!  write(*,'(I5,10Es14.6)') i, df%dat(i)
!!  write(1,'(I10,10Es14.6)') i, df%dat(i)
!! end do
!! 
!! close(1)


!end subroutine test_calc



end module distr_fun_m


!program main
! use distr_fun_m, only : test_calc 

! call test_calc

!end program main
