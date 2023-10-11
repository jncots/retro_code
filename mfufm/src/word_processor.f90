!  =========== DOCUMENTATION ===========
!
!  Module:
!  =======
!
!    use word_processor, only : del_break, del_2space, line_break
!
!  Definition:
!  ===========
!
!    (subroutine) del_break(line) - delete all line breaks "\n" in "line" creating single line
!
!    (subroutine) del_2space(line) - delete more than 1 space in a row leaving only 1 space
!
!    (subroutine) line_break(nleng,line) - break a string in sequence of line not longer "nleng"
!    deleting excess spaces
!
!  Syntax:
!  =======
!
!    call del_break(line)
!    call del_2space(line)
!    call line_break(nleng,line)
!
!  Arguments:
!  ==========
!
!    LINE:     character(:), allocatable :: line
!    NLENG:    integer :: nleng
!
!  Example:
!  ========
!
!   integer :: nleng
!   character(:), allocatable :: line
!   character(1000) :: str
!
!   str='In the past decade, mobile computing           has gone from a niche market for well-heeled &
!   enterprises with large field organisations to the fastest growing, and often most popular, &
!   way for employees of organisations of all sizes to do business computing. &
!   The near-universal adoption of mobile devices by consumers—who are also &
!   employees—has forced one of the most major shifts that corporate IT has ever seen.'  
!
!   line=str
!   nleng=50
!   call line_break(nleng,line)
!   write(*,*) line
!
!  Result:
!  =======
!
!    In the past decade, mobile computing has gone 
!    from a niche market for well-heeled enterprises 
!    with large field organisations to the fastest 
!    growing, and often most popular, way for 
!    employees of organisations of all sizes to do 
!    business computing. The near-universal adoption 
!    of mobile devices by consumers—who are also 
!    employees—has forced one of the most major shifts 
!    that corporate IT has ever seen.
!
!
!

module word_processor
 implicit none
 private
 save
 
 public :: del_break, del_2space, line_break
 
 
 contains
 
 
subroutine del_break(line)
! Variables
 character(:), allocatable :: line
 integer :: nline, nb
 character(len=:), allocatable :: str
! Calculations
 
 nb=index(line,'\n') ! in case of absence of line break
 if (nb==0) return
 
 nline=len_trim(line)
 allocate(character(len=nline) :: str)
 str=line
 line=''
 
 do
  nb=index(trim(str),'\n')
  if (nb==0) then
   line=line//trim(str)
   exit
  end if
  line=(line)//str(1:nb-1)
  if (nb>=len_trim(str)) exit
  str=trim(str(nb+1:))
  if (trim(str)=='') exit
 end do
 deallocate(str)
end subroutine del_break





subroutine del_2space(line)
! Variables
 character(:), allocatable :: line
 integer :: nn, nnn
! Calculations
 
 do 
  nnn=len_trim(line)
  nn=index(line,'  ')
  if (nn==0) exit
  line=line(1:nn)//line(nn+2:nnn)
 end do

end subroutine del_2space




subroutine line_break(nleng,line)
! Variables  
 integer, intent(in) :: nleng
 character(:), allocatable :: line
 character(len=:), allocatable :: str
 logical, parameter :: back=.true.
 integer :: nb, nline
! Calculations
 
 call del_break(line)
 call del_2space(line)
 
 nline=len_trim(line)
 allocate(character(len=nline) :: str)
 str=trim(line)
 
 line=''
 
 
 do
  if (len_trim(str)<=nleng) then
   line=line//trim(str(1:))
   exit
  else 
   nb=index(str(1:nleng),' ',back)
   if (nb==0) nb=nleng
   line=line//str(1:nb)//'\n'
   str=trim(str(nb+1:))
   if (trim(str)=='') exit
  end if 
 end do
 
 line=adjustl(line)//'\n'
 deallocate(str)
 
end subroutine line_break



end module word_processor