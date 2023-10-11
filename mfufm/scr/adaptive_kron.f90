! usage line:
! use gauss_kronrod, only : gk_adaptive
module gauss_kronrod
!====================================================
! Description:
! ============
! Calculation of integrals by adaptive method using
! gauss-kronrod rule with 15 points
!
! Version history:
! ================
! 23.03.2016:
! The problem for negative valued functions is
! fixed.
! delta_tot=abs(eps_tot/res_tot) now
!  instead of
!  delta_tot=eps_tot/res_tot was
! The calculation before probably was not always
! accurate. 
! 11.04.2016
! Add the condition:
! if (isnan(delta_tot)) fin=1
!====================================================
 implicit none
 private
 public :: gk_adaptive
  
 type segment
  real(8) :: a
  real(8) :: b
  real(8) :: res
  real(8) :: eps
 end type segment
  
  
 contains


recursive subroutine gk_adaptive(fun,a,b,res,tol,eps,nfun)
!
! Main program for adaptive integration
!
 implicit none
! Arguments
 integer, optional, intent(out) :: nfun
 real(8), intent(in) :: a, b
 real(8), optional, intent(in) :: tol
 real(8), optional, intent(out) :: eps
 real(8), intent(out) :: res
! Parameters
 integer, parameter :: nfmax_def=10000 ! default value
 real(8), parameter :: tol_def=1d-4 ! default value 
! Variables
 type(segment) :: arr(nfmax_def)
 integer :: nfmx,nfn, fin, nmax
 real(8) :: delta_tot
 real(8) :: res_tot, eps_tot, toler
 type(segment) :: seg, seg1, seg2
! Interface 
 interface
   subroutine fun(x,res)
    real(8), intent(in) :: x
    real(8), intent(out) :: res
   end subroutine
 end interface   
! Calculations
 
 
 if (present(tol)) then
  toler=tol
 else
  toler=tol_def
 end if
 
 nfmx=nfmax_def
 
 nfn=0
 nmax=0
 fin=0
 
 res_tot=0d0
 eps_tot=0d0
 
 call gauss_kron15(fun,a,b,seg)
 
 nfn=15
 res_tot=seg%res
 eps_tot=seg%eps
 delta_tot=abs(eps_tot/res_tot)

 
 if (delta_tot<=toler) then
   res=res_tot
   if (present(eps)) eps=delta_tot
   if (present(nfun)) nfun=nfn
   return
 end if

 call add_to_array(seg,arr,nmax)

 
 do  
  seg=arr(nmax)
  call subdiv(fun,seg,seg1,seg2)
  call int_val(seg,seg1,seg2,toler,res_tot,eps_tot,fin)
  if (fin==1) exit
  nmax=nmax-1
  call add_to_array(seg1,arr,nmax)
  call add_to_array(seg2,arr,nmax)
  nfn=nfn+30
  if (nfn>nfmx) then
!   write(*,*) 'The limit for number of function evaluations is attained: nfmax=', nfmx
!   write(*,*) abs(eps_tot/res_tot), res_tot
!   read(*,*)
   exit
  end if
  
 end do
 
 res=res_tot
 
 if (present(eps)) eps=abs(eps_tot/res_tot)
 if (present(nfun)) nfun=nfn


end subroutine gk_adaptive



recursive subroutine int_val(seg,seg1,seg2,toler,res_tot,eps_tot,fin)
!
! Calculation of sum at all segments and check for the accuracy
!
 implicit none
! Arguments
 type(segment), intent(in) :: seg, seg1, seg2
 real(8), intent(in) :: toler
 real(8) :: res_tot, eps_tot
 integer :: fin
! Variables
 real(8) :: delta_tot
 real(8) :: res_plus, eps_plus
! Calculations 
 
 res_plus=seg1%res+seg2%res-seg%res
 res_tot=res_tot+res_plus
 
 eps_plus=seg1%eps+seg2%eps-seg%eps
 eps_tot=eps_tot+eps_plus
 
 delta_tot=abs(eps_tot/res_tot)
 
 if (delta_tot<=toler) fin=1
 if (isnan(delta_tot)) fin=1
 
end subroutine int_val



recursive subroutine add_to_array(x,arr,nmax)
!
! Addition of a new element to array
!
 implicit none
! Arguments
 type(segment), intent(in) :: x
 type(segment) :: arr(10000)
 integer :: nmax
! Variables
 integer :: i, ii
! Calculations 
 
 arr(nmax+1)=x
 
 if (nmax==0) then
  nmax=1
  return
 end if 
 
 do i=0,nmax-1
  ii=nmax+1-i 
  if (arr(ii)%eps>=arr(ii-1)%eps) then
   exit
  else
   arr(ii)=arr(ii-1)
   arr(ii-1)=x
  end if  
 end do
 nmax=nmax+1
end subroutine add_to_array
 

recursive subroutine subdiv(fun,seg,seg1,seg2)
!
! Division of a segment
!
 implicit none
! Arguments
 type(segment), intent(in) :: seg
 type(segment), intent(out) :: seg1, seg2
! Variables
 real(8) :: a, b, ab, rr
 
 interface
   subroutine fun(x,res)
    real(8), intent(in) :: x
    real(8), intent(out) :: res
   end subroutine fun
 end interface
! Calculation
 a=seg%a
 b=seg%b
 ab=(a+b)/2
 call gauss_kron15(fun,a,ab,seg1)
 call gauss_kron15(fun,ab,b,seg2)
end subroutine subdiv




recursive subroutine gauss_kron15(f,a,b,seg)
!
! Calculation of integral according gauss-kronrod rule with 15 points
!
 implicit none
! Arguments
 real(8), intent(in) :: a, b
 type(segment), intent(out) :: seg
! Parameters 
 real(8), parameter :: xw(0:7)=&
 [0d0,0.2077849550078984676007d0,0.405845151377397166907d0,0.5860872354676911302941d0,&
  0.741531185599394439864d0,0.86486442335976907279d0,0.9491079123427585245262d0,&
  0.991455371120812639207d0]
 real(8), parameter :: wk(0:7)=&
 [0.209482141084727828013d0,0.204432940075298892414d0,0.1903505780647854099133d0,&
  0.169004726639267902827d0,0.140653259715525918745d0,0.1047900103222501838399d0,&
  0.063092092629978553291d0,0.0229353220105292249637d0]
 real(8), parameter :: wg(0:3)=&
 [0.4179591836734693877551d0,0.3818300505051189449504d0,0.2797053914892766679015d0,&
  0.129484966168869693271d0]
! Variables
 integer :: ig, i
 real(8) :: h, x0, res0, xl, xr, resl, resr, sk, sg, resg, res
! Interface
 interface
   subroutine f(x,res)
    real(8), intent(in) :: x
    real(8), intent(out) :: res
   end subroutine
 end interface 
! Calculations 
 ig=0 
 h=0.5d0*(b-a)
 x0=a+h
 call f(x0,res0)
 sk=wk(0)*res0
 sg=wg(0)*res0
 
 do i=1,7
  xl=x0-h*xw(i)
  xr=x0+h*xw(i)
  
  call f(xl,resl)
  call f(xr,resr)
  sk=sk+wk(i)*(resl+resr)
   
  if (mod(i,2)==0) then
   ig=ig+1
   sg=sg+wg(ig)*(resl+resr)
  end if
 end do
 
 resg=sg*h
 res=sk*h
 
 seg%a=a
 seg%b=b
 seg%res=res
 seg%eps=abs(resg-res)
 
end subroutine gauss_kron15


end module gauss_kronrod