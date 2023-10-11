module house_keeping
 implicit none
 private
 save
 
 public :: get_unit
 
 
 contains

subroutine get_unit(unit)
!====================================
! Search for free 'unit' from (1..99)
!====================================
! Variables
 integer, intent(out) :: unit
 integer :: i
! Calculations
 unit=0
 do i=1,99
  if ((i/=5).and.(i/=6).and.(i/=9)) then
   open(unit=i,err=10,status='scratch')
   close(unit=i)
   unit=i
   return
  end if
10  continue
 end do
end subroutine get_unit


end module house_keeping