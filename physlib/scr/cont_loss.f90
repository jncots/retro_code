module cont_loss
 use gauss_kronrod, only : gk_adaptive
 use format_file_tools, only : format_file_data
 use tab2fun, only : tabfun
 use house_keeping, only : get_unit
 use intpol_mod, only : arr_ind_short, logl_int
 implicit none
 private
 
 public :: closs

 abstract interface
  subroutine fun_1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine fun_1d
 end interface 
  
  
 type closs
  character(500) :: floss='losses.dat', ftab='tab_losses.dat'
  integer :: nn, n1=600, n2=600
  real(8), allocatable :: en(:), be(:)
  real(8), allocatable :: e1(:), e2(:,:), tl(:,:)
  real(8) :: eq1, eq2
  procedure(fun_1d), pointer, nopass :: qsp=>null()
  integer :: imp=0        ! continous injection, imp=1 is impulsive
  real(8) :: time=0d0     ! stationary problem, time dependent if this%time>0
  integer :: tabon=0      ! tabulated data is not read, 1 - is read
  integer :: eff_err
  real(8) :: max_time
  type(tabfun) :: befun
 contains
  procedure :: eff_en, del, tab, set_inj, fdist
 
 end type closs
 
 
 type(tabfun) :: loss
 procedure (fun_1d), pointer :: qspec=>null()
 
 contains

 
 
subroutine del(this) 
 class(closs) :: this
 call this%befun%off
 if (allocated(this%en)) deallocate(this%en)
 if (allocated(this%be)) deallocate(this%be)
 if (allocated(this%e1)) deallocate(this%e1)
 if (allocated(this%e2)) deallocate(this%e2)
 if (allocated(this%tl)) deallocate(this%tl)
 this%n1=600
 this%n2=600
 this%imp=0
 this%time=0d0
 this%tabon=0
 this%floss='losses.dat'
 this%ftab='tab_losses.dat'
 this%qsp=>null()
end subroutine del
 
 
 

subroutine tab(this)
 class(closs) :: this
 type(format_file_data) :: ldat
 integer :: i, j, ndev
 integer :: n1, n2, na1, na2
 real(8) :: e1, e2, eh1, eh2, en1, en2, res
 
 ldat%fin=this%floss
 call ldat%read
 call loss%on(ldat%fdata(:,1),ldat%fdata(:,2))
 e1=loss%xmin
 e2=loss%xmax
 n1=this%n1
 n2=this%n2
 na1=n1+1
 na2=n2+1
 
 allocate(this%e1(na1))
 allocate(this%e2(na2,na1))
 allocate(this%tl(na2,na1))
 do i=0,n1
  en1=e1*(e2/e1)**(1d0*i/n1)
  if (i==n1) en1=e2
  this%e1(i+1)=en1
  eh1=en1
  eh2=e2
  do j=0,n2
   en2=eh1*(eh2/eh1)**(1d0*j/n2)
   if (j==n2) en2=e2
   this%e2(j+1,i+1)=en2
   call tbe(en1,en2,res)
   this%tl(j+1,i+1)=res
   if (isnan(res).or.res<0d0) write(*,*) res
  end do
 end do

 call get_unit(ndev)
 open(ndev,file=this%ftab,form='unformatted')
 write(ndev) size(this%e1)
 write(ndev) size(this%e2,2)
 write(ndev) this%e1
 write(ndev) this%e2
 write(ndev) this%tl
 write(ndev) size(ldat%fdata(:,1))
 write(ndev) ldat%fdata(:,1)
 write(ndev) ldat%fdata(:,2)
 close(ndev)
 
 call loss%off
 deallocate(this%e1)
 deallocate(this%e2)
 deallocate(this%tl)
 
end subroutine tab


subroutine tbe(e1,e2,res)
 real(8), intent(in) :: e1, e2
 real(8), intent(out) :: res
 real(8) :: tmin, tmax
 
 if (abs(e2-e1)<1d-307) then
  res=0d0
  return
 end if 
 
 tmin=log(e1)
 tmax=log(e2)
 call gk_adaptive(ibe,tmin,tmax,res)
end subroutine tbe


subroutine ibe(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: ee
 
 ee=exp(x)
 call loss%fx(ee,res) 
 res=ee/res
end subroutine ibe



subroutine read_tab(this) 
 class(closs) :: this
 integer :: ndev
 
 call get_unit(ndev)
 open(ndev,file=this%ftab,form='unformatted')
 
 read(ndev) this%n1
 read(ndev) this%n2
 if (allocated(this%e1)) deallocate(this%e1)
 if (allocated(this%e2)) deallocate(this%e2)
 if (allocated(this%tl)) deallocate(this%tl)
 allocate(this%e1(this%n1))
 allocate(this%e2(this%n2,this%n1))
 allocate(this%tl(this%n2,this%n1))
 read(ndev) this%e1
 read(ndev) this%e2
 read(ndev) this%tl
 read(ndev) this%nn
 if (allocated(this%en)) deallocate(this%en)
 if (allocated(this%be)) deallocate(this%be)
 allocate(this%en(this%nn))
 allocate(this%be(this%nn))
 read(ndev) this%en
 read(ndev) this%be 
 
 close(ndev)
 call this%befun%on(this%en,this%be) ! initiation of the function on losses
end subroutine read_tab



subroutine eff_en(this,e0,t0,res)
! Variables
 class(closs) :: this
 real(8), intent(in) :: e0, t0
 real(8), intent(out) :: res
 integer :: n1, n2, ne1, ne2, nes, nts
 real(8) :: res1, res2
! Calculations

 this%eff_err=0
 
 
 nes=size(this%e1)
 if ((e0<this%e1(1)).or.(e0>this%e1(nes))) then
  res=e0
  this%eff_err=1
  return
 end if
 
 call arr_ind_short(1,nes,this%e1,e0,ne1,ne2)
 nts=size(this%tl(:,ne1))
 this%max_time=this%tl(nts,ne1)
 if (t0>this%max_time) then
  res=this%e2(size(this%e2(:,ne1)),ne1)
  this%eff_err=2
  return
 else
  call arr_ind_short(1,nts,this%tl(:,ne1),t0,n1,n2)
  call logl_int(this%tl(n1,ne1),this%tl(n2,ne1),this%e2(n1,ne1),this%e2(n2,ne1),t0,res1)
 end if
 
 nts=size(this%tl(:,ne2))
 this%max_time=this%tl(nts,ne2)
 if (t0>this%max_time) then
  res=this%e2(size(this%e2(:,ne2)),ne2)
  this%eff_err=2
  return
 else
  call arr_ind_short(1,nts,this%tl(:,ne2),t0,n1,n2)
  call logl_int(this%tl(n1,ne2),this%tl(n2,ne2),this%e2(n1,ne2),this%e2(n2,ne2),t0,res2)
 end if
 
 call logl_int(this%e1(ne1),this%e1(ne2),res1,res2,e0,res)

end subroutine eff_en


subroutine set_inj(this,eq1,eq2,qsp)
 class(closs) :: this
 real(8), intent(in) :: eq1, eq2
 procedure (fun_1d) :: qsp
 integer :: ns
 
 if (this%tabon==0) then
  call read_tab(this)
  this%tabon=1
 end if
 
 this%eq1=eq1
 this%eq2=eq2
 this%qsp=>qsp 
 
 ns=size(this%e1)
 if (this%eq2>this%e1(ns)) then
  write(*,'(A)')'=================Calculations could be WRONG !!!============================='
  write(*,'(A)') 'set_inj: emax_inj > emax_loss:'
  write(*,'(A,Es14.6,A)') 'maximum energy of the injection spectrum: emax_inj=', this%eq2,' eV'
  write(*,'(A,Es14.6,A)') 'maximum energy of the tabulated losses: emax_loss=', this%e1(ns),' eV'
  write(*,'(A)') 'it should be emax_inj <= emax_loss'
  write(*,'(A)')'============================================================================='
 end if
 
end subroutine set_inj



subroutine fdist(this,en,res)
 class(closs), target :: this
 real(8), intent(in) :: en
 real(8), intent(out) :: res
 
 if (en<this%en(1)) then
  call zero_loss(this,en,res)
  write(*,'(A)')'=================Calculations could be WRONG !!!============================='
  write(*,'(A)') 'fdist: input energy is out of range of the tabulated losses'
  write(*,'(A,Es14.6,A)') 'input energy: en=', en, ' eV'
  write(*,'(A,Es14.6,A)') 'minimum energy of the tabulated losses: emin_loss=', this%en(1), ' eV'
  write(*,'(A)')'============================================================================='
  return
 end if
 
 if (this%imp==0) then
  call qinj_cont(this,en,res)
 else
  call qinj_imp(this,en,res)
 end if 

end subroutine fdist



subroutine zero_loss(this,en,res)
 class(closs), target :: this
 real(8), intent(in) :: en
 real(8), intent(out) :: res
 
 call this%qsp(en,res)
 if (this%imp==0) res=res*this%time  
end subroutine zero_loss




subroutine qinj_imp(this,en,res)
 class(closs) :: this
 real(8), intent(in) :: en
 real(8), intent(out) :: res
 real(8) :: eff, beff, be
  
 if (this%time>0d0) then
  call this%eff_en(en,this%time,eff)
 else
  eff=this%en(size(this%en))
 end if
 
! write(*,'(10Es14.6)') en, eff, this%time, this%max_time
! read(*,*)
 
 call this%befun%fx(eff,beff)
 call this%befun%fx(en,be)
 if ((eff<this%eq1).or.(eff>this%eq2)) then
  res=0d0
  return
 end if 
 call this%qsp(eff,res)
 res=res*beff/be
end subroutine qinj_imp



subroutine qinj_cont(this,en,res)
 class(closs), target :: this
 real(8), intent(in) :: en
 real(8), intent(out) :: res
 real(8) :: eff, emax, emin, e1, e2
 real(8) :: tmin, tmax, be
 
 qspec=>this%qsp
 if (this%time>0d0) then
  call this%eff_en(en,this%time,eff)
  if (this%eff_err>0) eff=this%eq2
 else
  eff=this%eq2
 end if 
 
 
 if (en<this%eq1) then
  e1=this%eq1
 else
  e1=en
 end if
 
 if (eff>this%eq2) then
  e2=this%eq2
 else
  e2=eff
 end if
 
 if (e1>=e2) then
  res=0d0
  return
 end if 
 
 tmin=log(e1)
 tmax=log(e2)
 call gk_adaptive(iqinj,tmin,tmax,res)
 call this%befun%fx(en,be)
 res=res/be
end subroutine qinj_cont
 
 
 
subroutine iqinj(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: en
 en=exp(x)
 call qspec(en,res)
 res=res*en
end subroutine iqinj



end module cont_loss



!program main
! use cont_loss, only : closs
! use phys_const, only : year
! implicit none
! type(closs) :: ff
! integer :: i, ne
! real(8) :: res, e1, e2, ee
! 
! interface
!  subroutine qspec(x,res)
!   real(8), intent(in) :: x
!   real(8), intent(out) :: res
!  end subroutine qspec
! end interface
! 
! 
! ff%floss='sum_lossh.dat'
!! ff%ftab='tab_lossh.dat'
!! call ff%tab
!! call ff%read_tab
!! call ff%eff_en(1d12,4.03991036d8,res)
!! write(*,'(10Es20.12)') res, ff%max_time
! 
! e1=1d9
! e2=1d17
! ne=1000
! 
! 
! call ff%set_inj(1d14,1d17,qspec)
!! ff%imp=1
! ff%time=1d1*year
! open(1,file='dist.dat')
! do i=0,ne
!  ee=e1*(e2/e1)**(1d0*i/ne)
!  call ff%fdist(ee,res)
!  write(1,*) ee, res*ee**2
! end do
! close(1)
!end program main
!
!
!subroutine qspec(x,res)
! real(8), intent(in) :: x
! real(8), intent(out) :: res
!
! res=1/x**2
!
!end subroutine qspec