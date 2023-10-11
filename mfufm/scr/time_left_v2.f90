module time_left
!============================================
! Changes:
! 31.03.2016 add optional argument del to
! subroutine timer, which is delay in seconds
!
! 08.04.2016 move present(del) condition in
! str=='start' condition  and delete variable
! sw
!============================================
 implicit none
 private
 save
 
 public :: timer
 
 
 contains

subroutine timer(str,nt,del)
! Variables
 character(5), intent(in) :: str
 integer, intent(in) :: nt
 real, optional, intent(in) :: del
 integer :: ntot, ncount
 real :: start, finish, finish1, elap1, elap, left, delay
 real, parameter :: delay_def=5.0 ! default value 
 character(100) :: s1, s2
 save
! Calculations 
 
 if (str=='start') then
  if (present(del)) then
   delay=del
  else
   delay=delay_def
  end if
  
  call cpu_time(start)
  finish1=start
  ntot=nt
  return
 end if

 if (str=='left') then
  call cpu_time(finish)
  elap1=finish-finish1
  if (elap1<delay) then
   continue
  else
  finish1=finish
  ncount=nt
  elap=finish-start
  left=elap*(ntot-ncount)/ncount
  call sec2watch(elap,s1)
  call sec2watch(left,s2)
  write(*,'(A)') 'Elapsed time = '//trim(s1)//'   Time left = '//trim(s2)
 end if 
 end if

end subroutine timer



subroutine sec2watch(times,stime)
! Variables 
 real, intent(in) :: times
 character(100), intent(out) :: stime
 real, parameter :: minute=60
 real, parameter :: hour=60*minute
 real, parameter :: day=24*hour
 real, parameter :: year=365*day
 real :: time
 integer :: ny, nd, nh, nm
 character(100) :: str
 
 time=times
 ny=0
 nd=0
 nh=0
 nm=0
 
 if (time>=year) then
  ny=int(time/year)
  time=time-ny*year
 end if
 
 if (time>=day) then
  nd=int(time/day)
  time=time-nd*day
 end if
 
 if (time>=hour) then
  nh=int(time/hour)
  time=time-nh*hour
 end if
 
 if (time>=minute) then
  nm=int(time/minute)
  time=time-nm*minute
 end if
 
 write(stime,'(f6.3,A)') time,'s'
 
 if ((nm>0).or.(nh>0).or.(nd>0).or.(ny>0)) then
  write(str,'(I2,A)') nm,'m:'
  stime=trim(adjustl(str))//trim(adjustl(stime))
 end if
 
 if ((nh>0).or.(nd>0).or.(ny>0)) then
  write(str,'(I2,A)') nh,'h:'
  stime=trim(adjustl(str))//trim(adjustl(stime))
 end if
 
 if ((nd>0).or.(ny>0)) then
  write(str,'(I3,A)') nd,'d:'
  stime=trim(adjustl(str))//trim(adjustl(stime))
 end if
 
 if (ny>0) then
  write(str,'(I6,A)') ny,'yr:'
  stime=trim(adjustl(str))//trim(adjustl(stime))
 end if

end subroutine sec2watch

end module time_left