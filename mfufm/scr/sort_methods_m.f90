module sort_methods_m
!=======================================
! Contains fast sorting algorithms
!=======================================
 implicit none
 private
 save
 
 public :: sort_methods
! public :: test_calc

 type sort_methods
 contains
   procedure, nopass :: merge_sort=>bup_merge_sort
 end type sort_methods


 contains

subroutine bup_merge_sort(ar,ast,ist,ifin)
!=======================================
! Implimentation of bottom-up
! merge sort algorithm
!=======================================
 real(8) :: ar(:)
 integer, intent(in), optional :: ast, ist, ifin
 real(8), allocatable :: br(:)
 integer :: st, wd, wd2, mid, fin
 integer :: ast0, afin0, ist0, ifin0, wd0

 if (present(ast)) then
  ast0=ast
 else
  ast0=1
 end if

 if (present(ist)) then
  ist0=ist
 else
  ist0=ast0
 end if

 afin0=size(ar)+ast0-1

 if (present(ifin)) then
  if (ifin<afin0) then
   ifin0=ifin
  else
   ifin0=afin0
  end if
 else
  ifin0=afin0
 end if


 allocate(br(ast0:afin0))

 wd0=ifin0-ist0+1
 wd=1 ! width
 do
  if (wd>wd0) exit
  st=ist0
  do
   if (st>ifin0) exit

   mid=st+wd
   if (mid>ifin0) mid=ifin0
   
   wd2=2*wd
   fin=st+wd2-1
   if (fin>ifin0) fin=ifin0
   
   call bottom_up_merge(st,mid,fin,ar,br)
   st=fin+1
 
  end do
  ar(ist0:ifin0)=br(ist0:ifin0)
  wd=wd2
 end do
 
 deallocate(br)

end subroutine bup_merge_sort



subroutine bottom_up_merge(start,middle,finish,ahalf,afull)
 integer, intent(in) :: start, middle, finish 
 real(8), intent(in) :: ahalf(:)
 real(8), intent(out) :: afull(:)
 integer :: i, j, k

 i=start
 j=middle

 do k=start,finish
  
  if ((i<middle).and.((j>finish).or.(ahalf(i)<=ahalf(j)))) then
    afull(k)=ahalf(i)
    i=i+1
  else
    afull(k)=ahalf(j)
    j=j+1
  end if

 end do

end subroutine bottom_up_merge


!subroutine test_calc
! use random_gen_m, only : init_random_seed
! type(sort_methods) :: sm
! integer, parameter :: nt=17
! real(8) :: a(nt), b(nt)
! integer :: i

! call init_random_seed()
! call random_number(a)

! a=a*1d3-5d2
!! a=[3d0,-2d0,6d0,35d0,8d0,1d0,17d0,33d0,-12d0,11d0,-33d0]

!! do i=1,nt
!!  write(*,'(10I5)') i, int(a(i))
!! end do

! b=a
! call sm%merge_sort(a)

! do i=1,nt
!  write(*,'(10I5)') i, int(a(i)), int(b(i))
! end do

!end subroutine test_calc



end module sort_methods_m



!program main
! use sort_methods_m, only : test_calc
! call test_calc

!end program main
