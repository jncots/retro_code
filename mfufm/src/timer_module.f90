module timer_module
!===================================================
! Modification of timer for object-oriented approach
! Usage:
! =====
! use timer_module, only : timer_class
! type(timer_class) :: tmr
! call tmr%start(ni*nj,1.0) ! 1.0 is interval in sec
! do i=1,ni
!  do j=1,nj
!   ......
!   call tmr%loop          ! add number to counter and
! ! show information
!  end do
! end do
!===================================================
 implicit none
 private
 
 public :: timer_class 
 
 type timer_class
  character(500) :: proc_name='none'
  integer, private :: ntot, ncount, ncount1
  real, private :: delay, time0, rec_time, time1
 contains
   procedure :: start, loop
 end type timer_class
 
 
 contains



subroutine start(this,nt,del,pname)
! Variables
 class(timer_class) :: this
 integer, intent(in) :: nt
 real, optional, intent(in) :: del
 character(500), optional, intent(in) :: pname
 real, parameter :: delay_def=5.0 ! default value 
! Calculations
 if (present(del)) then
  this%delay=del
 else
  this%delay=delay_def
 end if
 
 if (present(pname)) then
  this%proc_name=pname
 end if
 
 call cpu_time(this%time0)
 this%rec_time=this%time0
 this%ntot=nt
 this%ncount=0
 this%ncount1=0
end subroutine


subroutine show(this)
! Variables
 class(timer_class) :: this
 integer :: a
! Calculations
 a=this%ntot
end subroutine show


subroutine loop(this,ostr)
! Variables
 class(timer_class) :: this
 character(500), optional, intent(out) :: ostr
 real :: now, elap, left, left1
 character(100) :: s1, s2, s3, s4, str
 integer :: procent
! Calculations
 this%ncount=this%ncount+1
 call cpu_time(now)
 if (now-this%rec_time<this%delay) then
  if (present(ostr)) ostr=''
  return
 end if 

 this%time1=this%rec_time
 this%rec_time=now
 elap=this%rec_time-this%time0
! left=elap*(this%ntot-this%ncount)/this%ncount
 left=(this%ntot-this%ncount)*(this%rec_time-this%time1)/(this%ncount-this%ncount1)
 procent=floor(this%ncount*1d2/this%ntot)
 call sec2watch(elap,s1)
 call sec2watch(left,s2)
 write(s3,'(I3,A)') procent,'%'
 s3=trim(adjustl(s3))
 call str_time(left,s4)
 
 if (this%proc_name=='none') then
  str=' ['//trim(s3)//', '//trim(s1)//']'//'  <'//trim(s2)
  str=trim(str)//'>  ['//trim(s4)//']'
 else
  str=' ['//trim(this%proc_name)//']'
  str=trim(str)//' ['//trim(s3)//', '//trim(s1)//']'//'  <'//trim(s2)
  str=trim(str)//'>  ['//trim(s4)//']'
 end if 
 
 if (present(ostr)) then
  ostr=trim(str)
 else
  write(*,'(A)') trim(str)
 end if
 this%ncount1=this%ncount
end subroutine loop




subroutine str_time(left,res)
! Variables
 real, intent(in) :: left
 character(100), intent(out) :: res
 character(100) :: sn, sl, day
 character(2) :: sndig(3), sldig(3), sday
 integer :: i, now(9), later(9), snow, slater
 character(9), parameter :: week_day(0:6)=['Sunday   ','Monday   ',&
 'Tuesday  ','Wednesday','Thursday ','Friday   ','Saturday ']
! Calculations
   
 snow=time()
 slater=snow+ceiling(left)
 call ltime(snow,now)
 call ltime(slater,later)


 write(day,'(A,I4)') '/',later(6)+1900
 call dig2str(later(5)+1,sday)
 day='/'//sday//trim(day)
 call dig2str(later(4),sday)
 day=sday//trim(day)

 do i=1,3
  call dig2str(now(i),sndig(i))
  call dig2str(later(i),sldig(i))
 end do
 sn=sndig(3)//':'//sndig(2)//':'//sndig(1)
 sl=sldig(3)//':'//sldig(2)//':'//sldig(1)
 res=trim(sn)//' - '//trim(sl)//', '//trim(week_day(later(7)))
 res=trim(res)//', '//trim(day)
end subroutine str_time

subroutine dig2str(dig,str)
! Variables
 integer, intent(in) :: dig
 character(2), intent(out) :: str
! Calculations 
 if (dig<10) then
  write(str,'(A,I1)') '0',dig
 else
  write(str,'(I2)') dig
 end if
end subroutine


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
! Calculations 
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
 
 
 write(stime,'(I2,A)') floor(time),'s'
 stime=trim(adjustl(stime))
 
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



end module timer_module
