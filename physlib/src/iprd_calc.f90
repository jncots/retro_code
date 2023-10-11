module iprd_calc_m
!
! Calculation of the cosmic-ray distributions
!
 use plaw_inj_mod, only : plaw_inj
 use rec_diff_module, only : rec_diff
 use timer_module, only : timer_class
 use taber, only : tab_er
 use phys_const, only : psec, year

 implicit none
 private
 save
 
 public :: iprd_calc

 type iprd_calc
  real(8) :: diff_ind=0.3d0
  real(8) :: diff_d0=1d28
  real(8) :: diff_en0=1d9
 
  real(8) :: qsp_lum=1d40 ! erg/s
  real(8) :: qsp_ind=2d0
  real(8) :: qsp_elow=1d9
  real(8) :: qsp_ecut=1d18
  real(8) :: qsp_tobs=1d5*year
 
  integer :: res_ne=250
  real(8) :: res_e1=1d9, res_e2=1d19 
  
  integer :: res_nr=100
  real(8) :: res_r1=1d-2*psec, res_r2=300d0*psec
  
  character(500) :: cname='prot_spec'
 
 contains 
  procedure :: calc, show
 end type iprd_calc
 
 

 type(plaw_inj) :: qp
 type(rec_diff) :: rdiff
 type(timer_class) :: tmr
 real(8) :: d0, dind, den0
 
 
 
 contains
 

subroutine show(this)
 class(iprd_calc) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'Cosmic-ray density calculation for:'
 write(*,*) 'Diffusion coefficient:'
 write(*,'(A,f14.6)') ' Index=', this%diff_ind
 write(*,'(A,Es14.6)') ' D0=', this%diff_d0
 write(*,'(A,Es14.6)') ' E0=', this%diff_en0
 write(*,*)
 write(*,*) 'Spectrum:'
 write(*,'(A,Es14.6)') ' Luminosity=', this%qsp_lum
 write(*,'(A,Es14.6)') ' Power-law index=', this%qsp_ind
 write(*,'(A,Es14.6)') ' E0=', this%qsp_elow
 write(*,'(A,Es14.6)') ' Ecut=', this%qsp_ecut
 write(*,'(A,Es14.6)') ' tobs=', this%qsp_tobs/year
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' ne=', this%res_ne
 write(*,'(2(A,Es14.6))') ' E1=', this%res_e1, ' E2=', this%res_e2 
 write(*,'(A,I5)') ' nr=', this%res_nr
 write(*,'(2(A,Es14.6))') ' r1=', this%res_r1/psec, ' r2=', this%res_r2/psec
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show
 
 
 
subroutine calc(this,prden)
 class(iprd_calc) :: this
 type(tab_er) :: prden
 integer :: ir, ie
 real(8) :: lum, ind, elow, ecut, tobs
 real(8) :: r0, e0, res
 
 lum=this%qsp_lum
 ind=this%qsp_ind
 elow=this%qsp_elow
 ecut=this%qsp_ecut
 tobs=this%qsp_tobs
 
 d0=this%diff_d0
 dind=this%diff_ind
 den0=this%diff_en0
 

 call qp%init(lum,ind,elow,ecut)
 call rdiff%init(diff,qs)
 
 call prden%set(this%res_ne,this%res_nr)
 call prden%set_ea(this%res_e1,this%res_e2)
 call prden%set_ra(this%res_r1,this%res_r2)
 
 tmr%proc_name=trim(this%cname)
 call tmr%start(prden%ne*prden%nr,1.0)
 do ir=1,prden%nr
  r0=prden%ra(ir)
  do ie=1,prden%ne
   e0=prden%ea(ie)
   call rdiff%cr_den(e0,tobs,r0,res)
   prden%fa(ie,ir)=res
   call tmr%loop
  end do
 end do
 
 
end subroutine calc

 
subroutine qs(t,en,res)
 real(8), intent(in) :: t, en    
 real(8), intent(out) :: res
 real(8) :: a
 a=t
 call qp%spec(en,res)
end subroutine qs


subroutine diff(en,res)
 real(8), intent(in) ::  en    
 real(8), intent(out) :: res
 res=d0*(en/den0)**dind
end subroutine diff


end module iprd_calc_m


module ipp_rate_m
!
! Calculation of pp production rate for 
! gamma-rays and electrons
!
 implicit none
 private
 save
 
 
 public :: ipp_rate
 
 
 abstract interface
  subroutine fun1d(x,res)
   real(8), intent(in) :: x
   real(8), intent(out) :: res
  end subroutine
 end interface
 
 
 type ipp_rate
  character(10) :: ptype='gam' ! elec
  real(8) :: gas_den=100d0 ! gas number density  
  integer :: ne=250
  real(8) :: e1=1d7, e2=1d19
  character(500) :: cname='ipp_rate'
 contains 
  procedure :: calc, show
 end type ipp_rate
 
 

 contains


 
subroutine show(this)
 class(ipp_rate) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'pp calculation for:'
 write(*,'(A,A)') ' Particle=', this%ptype
 write(*,'(A,Es14.6)') ' Gas density=', this%gas_den
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' ne=', this%ne
 write(*,'(2(A,Es14.6))') ' E1=', this%e1, ' E2=', this%e2 
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show
 


subroutine calc(this,prden,grate)
 use timer_module, only : timer_class
 use taber, only : tab_er
 use tab2fun, only : tabfun
 use pp_sec_prod, only : set_protsp, set_protsp_min_max, pp_gam_dist, pp_el_dist
 
 class(ipp_rate) :: this
 type(tab_er), intent(in) :: prden
 type(tab_er) :: grate
 integer :: ir, ie
 real(8) :: pp_norm, res
 type(tabfun) :: tpden
 type(timer_class) :: tmr
 procedure (fun1d), pointer :: pp_ge=>null()

 
 if (this%ptype=='gam') then
  pp_ge=>pp_gam_dist
 else
  pp_ge=>pp_el_dist
 end if
 
 call grate%set(this%ne,prden%nr)
 call grate%set_ea(this%e1,this%e2)
 grate%ra=prden%ra
 
 tmr%proc_name=trim(this%cname)//':'//trim(this%ptype)
 call tmr%start(grate%ne*grate%nr,1.0)
 
 do ir=1,grate%nr
  call tpden%on(prden%ea,prden%fa(:,ir))
  call set_protsp(this%gas_den,tpden%fx,pp_norm)
  call set_protsp_min_max(tpden%xmin,tpden%xmax)
  
  do ie=1,grate%ne
   call pp_ge(grate%ea(ie),res)
   grate%fa(ie,ir)=pp_norm*res
   call tmr%loop
  end do
 end do
 
 pp_ge=>null()
 call tpden%off

end subroutine calc
 
 

end module ipp_rate_m


module ibh_rate_m
!
! Calculation of the cosmic-ray distributions
!
 implicit none
 private
 save
 
 
 public :: ibh_rate
 
 type ibh_rate
  character(500) :: bhtab='none'  
  integer :: ne=250
  real(8) :: e1=1d7, e2=1d19
  character(500) :: cname='ibh_rate'
 contains 
  procedure :: calc, show
 end type ibh_rate
 
 contains


 
subroutine show(this)
 class(ibh_rate) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'Bethe-Heitler rate calculation for:'
 write(*,'(A,A)') ' Table=', trim(this%bhtab)
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' ne=', this%ne
 write(*,'(2(A,Es14.6))') ' E1=', this%e1, ' E2=', this%e2 
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show
 

subroutine calc(this,prden,erate)
 use timer_module, only : timer_class
 use taber, only : tab_er
 use tab2fun, only : tabfun
 use beheit_module, only : beheit
 
 class(ibh_rate) :: this
 type(tab_er), intent(in) :: prden
 type(tab_er) :: erate
 type(timer_class) :: tmr
 type(beheit),save :: bhvar
 type(tabfun) :: tden
 integer :: ir, ie
 
 if (this%bhtab=='none') then
  write(*,*) 'ibh_rate: this%bhtab="none"'
  write(*,*) 'Give to this%bhtab correct path to the table file'
  return
 end if 
 
 if (bhvar%tab_file/=this%bhtab) then
  bhvar%tab_file=this%bhtab
  call bhvar%read_tab
 end if 
 
 call erate%set(this%ne,prden%nr)
 call erate%set_ea(this%e1,this%e2)
 erate%ra=prden%ra
 
 tmr%proc_name=this%cname
 call tmr%start(erate%nr*erate%ne,1.0)
 
 do ir=1,erate%nr
  call tden%on(prden%ea,prden%fa(:,ir))
  call bhvar%set_protsp(tden%xmin,tden%xmax,tden%fx)
  
  do ie=1,erate%ne
   call bhvar%bh_rate(erate%ea(ie),erate%fa(ie,ir))
   call tmr%loop
  end do
 
 end do
 call tden%off
 
end subroutine calc

 
end module ibh_rate_m




module ipg_rate_m
!
! Calculation of the photo-meson production rate
!
 implicit none
 private
 save
 
 
 public :: ipg_rate
 
 type ipg_rate
  character(500) :: tab='none'  
  integer :: ne=250
  real(8) :: e1=1d7, e2=1d19
  character(500) :: cname='ipg_rate'
 contains 
  procedure :: calc, show
 end type ipg_rate
 
 contains


 
subroutine show(this)
 class(ipg_rate) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'Photo-meson rate calculation for:'
 write(*,'(A,A)') ' Table=', trim(this%tab)
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' ne=', this%ne
 write(*,'(2(A,Es14.6))') ' E1=', this%e1, ' E2=', this%e2 
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show
 

subroutine calc(this,prden,erate)
 use rate_mes_gam, only : set_prod_rate, prod_rate, read_tab
 use timer_module, only : timer_class
 use taber, only : tab_er
 use tab2fun, only : tabfun
 
 class(ipg_rate) :: this
 type(tab_er), intent(in) :: prden
 type(tab_er) :: erate
 type(timer_class) :: tmr
 type(tabfun) :: tden
 integer :: ir, ie
 
 if (this%tab=='none') then
  write(*,*) 'ibh_rate: this%bhtab="none"'
  write(*,*) 'Give to this%bhtab correct path to the table file'
  return
 end if 
 
 call read_tab(this%tab) 
 call erate%set(this%ne,prden%nr)
 call erate%set_ea(this%e1,this%e2)
 erate%ra=prden%ra
 
 tmr%proc_name=this%cname
 call tmr%start(erate%nr*erate%ne,1.0)
 
 do ir=1,erate%nr
  call tden%on(prden%ea,prden%fa(:,ir))
  call set_prod_rate(tden%xmin,tden%xmax,tden%fx)
  
  do ie=1,erate%ne
   call prod_rate(erate%ea(ie),erate%fa(ie,ir))
   call tmr%loop
  end do
 
 end do
 call tden%off
 
end subroutine calc

 
end module ipg_rate_m


module elcool_den_m
!
! Cooling of electrons
!
 implicit none
 private
 save
 
 
 public :: elcool_den
 
 

 type elcool_den
  real(8) :: time=0d0
  character(500) :: tab='none' 
  character(500) :: tabout='tabout.dat'  
  character(500) :: cname='elcool'
 contains 
  procedure :: calc, show
 end type elcool_den
 
 contains


 
subroutine show(this)
 class(elcool_den) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'Electron cooling calculation for:'
 write(*,'(A,A)') ' Table=', trim(this%tab)
 write(*,*)
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show
 

subroutine calc(this,erate,eden)
 use timer_module, only : timer_class
 use taber, only : tab_er
 use tab2fun, only : tabfun
 use cont_loss, only : closs
 
 class(elcool_den) :: this
 type(tab_er), intent(in) :: erate
 type(tab_er) :: eden
 type(closs), save :: cl
 type(tabfun) :: trate
 type(timer_class) :: tmr
 integer :: ir, ie
 
 
 if (this%tab=='none') then
  write(*,*) 'elcool_den: this%tab="none"'
  write(*,*) 'Give to this%tab correct path to the table file'
  return
 end if 
 
 
 if (cl%floss/=this%tab) then
  call cl%del
  cl%floss=this%tab
  cl%ftab=this%tabout
  call cl%tab
 end if 
 cl%time=this%time
 
 call eden%set(erate%ne,erate%nr)
 eden%ea=erate%ea
 eden%ra=erate%ra
 
 
 tmr%proc_name=trim(this%cname)
 call tmr%start(eden%nr*eden%ne,1.0)
 do ir=1,eden%nr
  call trate%on(erate%ea,erate%fa(:,ir))
  call cl%set_inj(trate%xmin,trate%xmax,trate%fx)
  
  do ie=1,eden%ne
   call cl%fdist(eden%ea(ie),eden%fa(ie,ir))
   call tmr%loop
  end do
 
 end do
 call trate%off
 
end subroutine calc
 
 

end module elcool_den_m


module iggep_rate_m
!========================================================================================
! Calculation of density of electron-positron pairs 
! from gamma-rays
! Example:
! ========
! subroutine calc_eppair
! use iggep_rate_m, only : iggep_rate
! type(iggep_rate) :: cep
!  
! cep%tab='D:\project\Diffusion\SphericalDiffusion\Estimations\xray_gc\soft_phf\soft_phf500.dat'
! cep%cname='eppair'
! cep%rmin=0d0
! cep%rmax=10*psec
! 
! cep%nr=100
! cep%r1=1d-2*psec
! cep%r2=10d0*psec
! 
! cep%ne=250
! cep%e1=1d7
! cep%e2=1d18
! 
! call cep%show
! call cep%calc(pggam_rate,pggam_ep_rate)
!
! end subroutine calc_eppair
!========================================================================================
 use tab2fun, only : tabfun
 use phys_const, only : psec
 
 implicit none
 private
 save
 
 
 public :: iggep_rate
 
 

 type iggep_rate
  integer :: nr ! number of points in space
  integer :: ne ! number of points for el-pos pair spectrum
  real(8) :: rmin=0d0, rmax=1d2*psec ! size of the gamma-ray region
  real(8) :: r1, r2 ! calculation of the density in this region
  real(8) :: e1, e2 ! spectrum range for el-pos spectrum
  character(500) :: tab='none' ! table for gg->ep production
  character(500) :: cname='ggep_rate'
 contains 
  procedure :: calc, show
 end type iggep_rate
 
 
 type(tabfun), pointer :: t0fun
 
 
 contains


 
subroutine show(this)
 class(iggep_rate) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'Calculation of electron-positron pair production:'
 write(*,'(A,A)') ' Data table=', trim(this%tab)
 write(*,*)
 write(*,*) 'Take gamma-ray radiation in the region'
 write(*,'(2(A,Es14.6))') ' rmin=', this%rmin/psec, ' rmax=', this%rmax/psec
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' nr=', this%nr
 write(*,'(2(A,Es14.6))') ' r1=', this%r1/psec, ' r2=', this%r2/psec
 write(*,'(A,I5)') ' ne=', this%ne
 write(*,'(2(A,Es14.6))') ' E1=', this%e1, ' E2=', this%e2
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show
 


subroutine calc(this,ggrate,eprate)
!
! Calculation of electron-positron rate from gamma rate
!
 use taber, only : tab_er
 use gg_elpos, only : gg_elpos_rate
 use gam_den, only : gamden, set_gamden
 use timer_module, only : timer_class
 
 class(iggep_rate) :: this
 type(tab_er), intent(in) :: ggrate
 type(tab_er) :: eprate
 
 real(8) :: lambda, ee, res
 integer :: ie, ir
 type(tab_er) :: ggden
 type(gg_elpos_rate), save :: rtep
 type(tabfun), target :: tgrate
 type(tabfun) :: tgden
 type(timer_class) :: tmr

! Calculations

 if (this%tab=='none') then
  write(*,*) 'iggep_rate: this%tab="none"'
  write(*,*) 'Give to this%tab correct path to the table file'
  return
 end if 
 
 
 if (rtep%fname/=this%tab) then
  rtep%fname=this%tab
  call rtep%read
 end if
 
 
! Calculation of the density of gamma-rays
 call ggden%set(ggrate%ne,this%nr)
 call ggden%set_ra(this%r1,this%r2)
 ggden%ea=ggrate%ea
  
 tmr%proc_name=trim(this%cname)//':ggden'
 call tmr%start(ggden%ne*ggden%nr,1.0)
 
 do ie=1,ggden%ne
  call tgrate%on(ggrate%ra,ggrate%fa(ie,:))
  call rtep%wabs(ggrate%ea(ie),lambda)
  t0fun=>tgrate
  call set_gamden(this%rmin,this%rmax,t0ext)
  do ir=1,ggden%nr
   call gamden(lambda,ggden%ra(ir),res)
   ggden%fa(ie,ir)=res
   call tmr%loop
  end do
 end do
 
 t0fun=>null()
 call tgrate%off
 
 
! Calculation of the rate of electron-positron pairs 
 call eprate%set(this%ne,ggden%nr)
 call eprate%set_ea(this%e1,this%e2)
 eprate%ra=ggden%ra
 
 tmr%proc_name=trim(this%cname)//':eerate'
 call tmr%start(eprate%nr*eprate%ne,1.0)
 
 do ir=1,eprate%nr
  call tgden%on(ggden%ea,ggden%fa(:,ir))
  call rtep%set_gam(tgden%xmin,tgden%xmax,tgden%fx)
  
  do ie=1,eprate%ne
   ee=eprate%ea(ie)
   call rtep%rate_epp(ee,res)
   eprate%fa(ie,ir)=res
   call tmr%loop
  end do
 
 end do
 
 call tgden%off
 call ggden%del
 
end subroutine calc


subroutine t0ext(x,res)
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 
 if (x<t0fun%xmin) then
  call t0fun%fx(t0fun%xmin,res)
  return
 else
  call t0fun%fx(x,res)
 end if

end subroutine t0ext

 
 
end module iggep_rate_m


module syn_env_m
!
! Cooling of electrons
!
 implicit none
 private
 save
 
 
 public :: syn_env
 
 

 type syn_env
  integer :: ne
  real(8) :: bm
  real(8) :: e1, e2 
  character(500) :: cname='synrad'
 contains 
  procedure :: calc, show
 end type syn_env
 
 contains


 
subroutine show(this)
 class(syn_env) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'Synchrotron radiation calculation:'
 write(*,'(A,Es14.6)') ' For magnetic field=', this%bm
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' ne=', this%ne
 write(*,'(2(A,Es14.6))') ' E1=', this%e1, ' E2=', this%e2
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show





subroutine calc(this,eden,grate)
 use timer_module, only : timer_class
 use taber, only : tab_er
 use tab2fun, only : tabfun
 use synch_module, only : synch_class
 
 class(syn_env) :: this
 type(tab_er), intent(in) :: eden
 type(tab_er) :: grate
 type(synch_class), save :: synvar
 type(tabfun) :: tden
 type(timer_class) :: tmr
 integer :: ir, ie
 
 
 call synvar%set_bm(this%bm)
 
 call grate%set(this%ne,eden%nr)
 grate%ra=eden%ra 
 call grate%set_ea(this%e1,this%e2)
 
 tmr%proc_name=trim(this%cname)
 call tmr%start(grate%nr*grate%ne,1.0)
 
 do ir=1,grate%nr
  call tden%on(eden%ea,eden%fa(:,ir))
  call synvar%set_elsp(tden%xmin,tden%xmax,tden%fx)
  do ie=1,grate%ne
   call synvar%syn_rad(grate%ea(ie),grate%fa(ie,ir))
   call tmr%loop
  end do
 end do
 
 call tden%off
 
end subroutine calc
 
 
end module syn_env_m


module ic_env_m
!
! Cooling of electrons
!
 implicit none
 private
 save
 
 
 public :: ic_env
 
 
 type ic_env
  integer :: ne
  real(8) :: e1, e2
  character(500) :: tab='none'
  character(500) :: cname='icrad'
 contains 
  procedure :: calc, show
 end type ic_env
 
 contains


 
subroutine show(this)
 class(ic_env) :: this
 write(*,*) '---------------------------------------'
 write(*,*) 'IC radiation calculation:'
 write(*,'(A,A)') ' Data table=', trim(this%tab)
 write(*,*)
 write(*,*) 'Results on the grid:'
 write(*,'(A,I5)') ' ne=', this%ne
 write(*,'(2(A,Es14.6))') ' E1=', this%e1, ' E2=', this%e2
 write(*,*) 'Name of calculation = ', trim(this%cname)
 write(*,*) '---------------------------------------'
end subroutine show



subroutine calc(this,eden,grate)
 use timer_module, only : timer_class
 use taber, only : tab_er
 use tab2fun, only : tabfun
 use ic_rad, only : ic_type
 
 class(ic_env) :: this
 type(tab_er), intent(in) :: eden
 type(tab_er) :: grate
 type(ic_type), save :: icvar
 type(tabfun) :: tden
 type(timer_class) :: tmr
 integer :: ir, ie
 
 
 if (this%tab=='none') then
  write(*,*) 'ic_env: this%tab="none"'
  write(*,*) 'Give to this%tab correct path to the table file'
  return
 end if 
 
 if (icvar%fin/=this%tab) then
  icvar%fin=this%tab
  call icvar%read_tab
 end if 
 
 
 call grate%set(this%ne,eden%nr)
 grate%ra=eden%ra 
 call grate%set_ea(this%e1,this%e2)
 
 tmr%proc_name=trim(this%cname)
 call tmr%start(grate%nr*grate%ne,1.0)
 
 do ir=1,grate%nr
  call tden%on(eden%ea,eden%fa(:,ir))
  call icvar%set_spec(tden%xmin,tden%xmax,tden%fx)
  do ie=1,grate%ne
   call icvar%ic_rate(grate%ea(ie),grate%fa(ie,ir))
   call tmr%loop
  end do
 end do
 
 call tden%off
end subroutine calc 

end module ic_env_m









! program main
 ! use iprd_calc_m, only : iprd_calc
 ! use taber, only : tab_er
 ! use phys_const, only : psec, year
 
 ! implicit none
 ! type(iprd_calc) :: prc
 ! type(tab_er) :: pden
 
 
 ! prc%diff_ind=0.5d0
 ! prc%cname='pden05'
 ! call prc%show
 ! call prc%calc(pden)
 
 
! end program main