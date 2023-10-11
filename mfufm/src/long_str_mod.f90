module long_str_mod
!===========================================================
! Module contains class long_str, which is describes
! dinamical string.
! Usage:
! ======
! call s%add('Some text') ! add to the string
! call s%del ! delete string
! call s%copy(3,15,word) ! copy string content
!                        ! in position from 3 to 15
!                        ! to word(string)
! call s%scopy(3,15,s1)  ! the same but s1 is long_str
! call s%sadd(s1) ! add to string s1(long_str)
! call s%find('hello',n) ! gives position of the first
!                        ! occurrence of the substring
!                        ! 'hello' in the string
! call s%find_back('hello',n) ! last occurrence
! call s%sbreak(50) ! break the string in lines not longer
!                   ! the the length 50
! call s%del_break ! delete all break line symbols '\n'
! call s%del_2space ! delete more than 2 spaces in a row
! call s%salign(50) ! formatting the string in the nice form
! call ff%sbeg(word,50) ! formatting string adding 'word'
!                       ! in the beginning of the each line
!===========================================================
 implicit none
 private
 
 public :: long_str
 
 type long_str
  integer :: slen=0
  character(100), allocatable :: s(:)
 contains
  procedure :: del, add, sadd
  procedure :: find, find_back, copy, scopy
  procedure :: sbreak, del_break, del_2space
  procedure :: salign, sbeg
 end type long_str

 
 contains

 
 
 
subroutine copy(this,n1,n2,word)
! Variables
 class(long_str) :: this
 integer, intent(in) :: n1, n2 
 character(*) :: word
 integer :: i, j, ia, is, n22
! Calculations
 
 if (n2-n1+1>len(word)) then
  word=''
  write(*,*) 'copy: The string is too small to copy in'
  return
 end if
 word=''
 n22=n2
 if (n2>this%slen) n22=this%slen
 do i=n1,n22
  call ind_conv(i,ia,is)
  j=i+1-n1
  word(j:j)=this%s(ia)(is:is)
 end do

end subroutine copy
 

subroutine scopy(this,n1,n2,that) 
! Variables
 class(long_str) :: this
 integer, intent(in) :: n1, n2
 type(long_str) :: that
 integer :: n22, i, j, ia, is, ia1, is1
! Calculations 
 call that%del
 n22=n2
 if (n2>this%slen) n22=this%slen
 that%slen=n22+1-n1
 call ind_conv(that%slen,ia,is)
 allocate(that%s(ia))
 that%s(ia)=''
 do i=n1,n22
  call ind_conv(i,ia,is)
  j=i+1-n1
  call ind_conv(j,ia1,is1)
  that%s(ia1)(is1:is1)=this%s(ia)(is:is)
 end do
end subroutine scopy 
 
 

subroutine find_back(this,word,res)
! Variables
 class(long_str) :: this
 character(*) :: word
 integer, intent(out) :: res
 integer :: nw, ns, nsw, i, j, ia, is, iw, ic 
! Calculations
 nw=len(word)
 ns=this%slen
 
 nsw=ns-nw+1 
 if (nsw<1) then
  res=0
  return
 end if 
  
 do i=1,nsw
  j=nsw-i+1
  call ind_conv(j,ia,is)
   if (this%s(ia)(is:is)==word(1:1)) then
   iw=1
   ic=j
   do
    
    if (iw==nw) then
     res=j
     return
    end if
    
    iw=iw+1
    ic=ic+1
    call ind_conv(ic,ia,is)
    if (this%s(ia)(is:is)/=word(iw:iw)) exit
   
   end do 
  end if
 end do
 
 
end subroutine find_back

 

subroutine find(this,word,res)
! Variables
 class(long_str) :: this
 character(*) :: word
 integer, intent(out) :: res
 integer :: nw, ns, i, ia, is, iw, ic
! Calculations
 
 nw=len(word)
 ns=this%slen
 
 do i=1,ns
  call ind_conv(i,ia,is)
  if (this%s(ia)(is:is)==word(1:1)) then
   iw=1
   ic=i
   do
    
    if (iw==nw) then
     res=i
     return
    end if
    
    if (ic==ns) then
     res=0
     return
    end if 
    
    iw=iw+1
    ic=ic+1
    call ind_conv(ic,ia,is)
    if (this%s(ia)(is:is)/=word(iw:iw)) exit
   
   end do 
  end if
 end do
 
end subroutine find



subroutine ind_conv(ic,ia,is)
! Variables 
 integer, intent(in) :: ic
 integer, intent(out) :: ia, is
 integer, parameter :: n0=100
! Calculations 
 ia=ic/n0
 is=ic-ia*n0
 if (is==0) then
  is=n0
 else
  ia=ia+1
 end if
end subroutine ind_conv
 

subroutine del(this)
! Variables
 class(long_str) :: this
! Calculations
 if (allocated(this%s)) deallocate(this%s)
 this%slen=0
end subroutine del
 
subroutine add(this,word)
! Variables
 class(long_str) :: this
 character(100), allocatable :: sloc(:)
 character(*) :: word
 integer, parameter :: n0=100
 integer :: ns, nf, nws, nw, np, i
! Calculations

 if (allocated(this%s)) then
  ns=size(this%s)
  nf=len_trim(this%s(ns))
  nw=len(word)
  if (nf<n0) then
   nws=n0-nf
   if (nw<=nws) then
    this%s(ns)(nf+1:nf+nw)=word(1:nw)
    return
   else 
    this%s(ns)(nf+1:n0)=word(1:n0-nf)
   end if 
   nws=nws+1
  else
   nws=1
  end if
  
  nw=len(word(nws:))
  np=ceiling(nw*1d0/n0)
  allocate(sloc(ns))
  sloc=this%s
  deallocate(this%s)
  allocate(this%s(ns+np))
  
  do i=1,ns
   this%s(i)=sloc(i)
  end do
  deallocate(sloc)
  do i=ns+1,ns+np-1
   this%s(i)=word(nws:nws+n0-1)
   nws=nws+n0
  end do
  this%s(ns+np)=word(nws:)
 
 else
  
  nw=len(word)
  np=ceiling(nw*1d0/n0)
  allocate(this%s(np))
  nws=1
  do i=1,np-1
   this%s(i)=word(nws:nws+n0-1)
   nws=nws+n0
  end do 
  this%s(np)=word(nws:nw)
 
 end if
 
 call sleng(this)
end subroutine add

subroutine sleng(this)
! Variables
 class(long_str) :: this
 integer :: ns
! Calculations  
 ns=size(this%s)
 this%slen=100*(ns-1)+len_trim(this%s(ns))
end subroutine sleng



subroutine sadd(this,that)
! Variables
 class(long_str) :: this
 type(long_str), intent(in) :: that
 character(100), allocatable :: sloc(:)
 integer :: tot, ia, is, i, j, ia1, is1
! Calculations
 
 if (allocated(this%s)) then
  tot=this%slen+that%slen
  allocate(sloc(size(this%s)))
  sloc=this%s
  deallocate(this%s)
  call ind_conv(tot,ia,is)
  allocate(this%s(ia))
  this%s=''
  do i=1,tot
   call ind_conv(i,ia,is)
   if (i<=this%slen) then
    this%s(ia)(is:is)=sloc(ia)(is:is)
   else
    j=i-this%slen
    call ind_conv(j,ia1,is1)
    this%s(ia)(is:is)=that%s(ia1)(is1:is1)
   end if 
  end do
  this%slen=tot
 else
  allocate(this%s(size(that%s)))
  this%slen=that%slen
  do i=1,this%slen
   call ind_conv(i,ia,is)
   this%s(ia)(is:is)=that%s(ia)(is:is)
  end do
 end if
 call sleng(this)
end subroutine sadd


subroutine sbreak(this,nlen)
! Variables
 class(long_str) :: this
 integer, intent(in) :: nlen
 integer :: ns, i, ia, is, nc, nl, nnn
! Calculations 
 
 ns=1
 nnn=this%slen
 call ind_conv(nnn,ia,is)
 if (ia<=size(this%s)) nnn=nnn+1
 
 do i=1,nnn
  call ind_conv(i,ia,is)
  if (this%s(ia)(is:is)==' ') then
   nc=i-ns+1
   if (nc>nlen) then
    ns=nl+1
    call ind_conv(nl,ia,is)
    this%s(ia)(is:is)='\n'
   end if
   nl=i
  end if
 end do
 call sleng(this)
end subroutine sbreak


subroutine del_break(this)
! Variables
 class(long_str) :: this
 integer :: i, ia, is
! Calculations 
 
 do i=1,this%slen
  call ind_conv(i,ia,is)
  if (this%s(ia)(is:is)=='\n') this%s(ia)(is:is)=' '
 end do
 call sleng(this)
end subroutine del_break

subroutine del_2space(this)
! Variables
 class(long_str) :: this
 integer :: i, j, ia, is, ia1, is1, nb
! Calculations
 i=1
 j=1
 nb=0
 do
  call ind_conv(i,ia,is)
  call ind_conv(j,ia1,is1)
  this%s(ia)(is:is)=this%s(ia1)(is1:is1)
  if (this%s(ia1)(is1:is1)==' ') then
   nb=nb+1
  else
   nb=0
  end if
  
  j=j+1
  if (j>this%slen) exit
  if (nb<2) i=i+1
 end do
 
 i=i+1
 do
  call ind_conv(i,ia,is)
  this%s(ia)(is:is)=' '
  i=i+1
  if (i>this%slen) exit
 end do 
 call sleng(this) 
end subroutine del_2space


subroutine salign(this,nlen)
! Variables
 class(long_str) :: this
 integer, intent(in) :: nlen
! Calculations 
 call this%del_break
 call this%del_2space
 call this%sbreak(nlen)
 call sleng(this)
end subroutine salign



subroutine sbeg(this,word,nlen)
! Variables
 class(long_str) :: this
 character(*), intent(in) :: word
 integer, intent(in) :: nlen
 type(long_str) :: that
 integer :: nw, nn, i, j, k, ia, is, ia1, is1, tot
! Calculations
 nw=len(word)
 if (nw>=nlen) then
  write(*,*) 'sbeg: the line is too short'
  call this%salign(nlen)
  return
 end if
 call this%salign(nlen-nw)
 
 nn=0
 do i=1,this%slen
  call ind_conv(i,ia,is)
  if (this%s(ia)(is:is)=='\n') nn=nn+1
 end do
 
 nn=nn+1
 tot=nn*nw+this%slen
 call that%sadd(this)
 call this%del
 call ind_conv(tot,ia,is)
 allocate(this%s(ia))
 this%s=''
 
 i=0
 j=1
 
 do k=1,nw
  i=i+1
  call ind_conv(i,ia,is)
  this%s(ia)(is:is)=word(k:k)
 end do
 
 do
  call ind_conv(i,ia,is)
  call ind_conv(j,ia1,is1)
  this%s(ia)(is:is)=that%s(ia1)(is1:is1)
  if (that%s(ia1)(is1:is1)=='\n') then
  do k=1,nw
   i=i+1
   call ind_conv(i,ia,is)
   this%s(ia)(is:is)=word(k:k)
  end do
  end if
  j=j+1
  if (j>that%slen) exit
  i=i+1
 end do
 
 call sleng(this)
 call that%del
end subroutine sbeg



end module long_str_mod



!program main
! use long_str_mod, only : long_str
! implicit none
! integer :: i, ns1, ns2
! character(1000) :: word
! type(long_str) :: ff, ff1
! call ff%add('The announcement was met with')
! call ff%add(' the speech and many still held feelings. ')
! call ff%add(' The announcement was met with polite applause.&
! However, Ars informally polled some audience members &
! after the speech and many still held feelings of &
! resentment toward the Obama administration &
! apparent lack of action. The          steps being taken by Huerta &
! come at the 11th hour for the administration, and &
! some audience members commented to Ars that nothing &
! would get done unless the next administration takes the same &
! sorts of steps that other countries have taken to accelerate innovation around drones.')
! write(*,*) ff%s
! call ff%find('  ',ns1)
! call ff%find_back('Ars',ns2)
! call ff%scopy(ns1,ns1+19,ff1)
! write(*,*) ff%slen
! call ff%sadd(ff1)
! call ff%sadd(ff1)
! call ff%sadd(ff1)
! write(*,*) '70'
! call ff%sbreak(70)
! call ff%add(' something elsessssss')
! write(*,*) ff%s
! call ff%del_break
! call ff%del_2space
! call ff%sbreak(70)
! write(*,*) ff%s
! call ff%salign(60)
! call ff%add(' something elsessssss')
! write(*,*) ff%s
! call ff%sbeg('***!!      ',50)
! write(*,*) ff%s
! 
!end program main