module thermal_photon_field
!==================================================================
! Description:
! ============ 
! Calculation of number density of photon field consisting of the
! set of thermal photon fields with given termperatures and
! energy density.
!
! Usage:
! ======
! use thermal_photon_field, only : therm_phf, set_thnden, thnden
!
! Set photon field:
! =================
! call field%init(nf) ! number of thermal fields, let nf=4, then
! field%temp=[3d0,3d-1,6d-3,2.35d-4] ! temperature in eV
! field%uden=[5d3,4d4,4d4,0.26d0] ! energy density in eV/cm^3,
! where 'field' is an object of type(therm_phf)
!
! Photon field number density:
! ============================
! call field%nden(x,res) ! x is energy in eV, res is number density
! in 1/(eV cm^3)
! or, if simple function of two arguments is needed:
!
! call set_thnden(field)
! call thnden(x,res) ! calculates number density for 'field'
!
!==================================================================
 implicit none
 private
 save

 public :: therm_phf, set_thnden, thnden
 
 type therm_phf
  integer :: nf ! number of thermal fields
  real(8), allocatable :: temp(:), uden(:) ! temperature (in [eV]) and energy density (in [eV/cm^3])
 contains
  procedure :: init=>set_thermal_nden, nden=>thermal_nden
 end type therm_phf

 
 
 type(therm_phf), pointer :: fieldc
 
 
 contains
 
 
subroutine set_thnden(field) 
! Variables
 type(therm_phf), intent(in), target :: field
! Calculations  
 fieldc=>field
end subroutine set_thnden


subroutine thnden(x,res)
! Variables
 real(8), intent(in) :: x
 real(8), intent(out) :: res
! Calculations
 call fieldc%nden(x,res)
end subroutine thnden

subroutine set_thermal_nden(field,nfield)
! Variables
 class(therm_phf) :: field
 integer, intent(in) :: nfield
! Calculations
 field%nf=nfield
 if (allocated(field%temp)) deallocate(field%temp)
 allocate(field%temp(field%nf))
 if (allocated(field%uden)) deallocate(field%uden)
 allocate(field%uden(field%nf))
end subroutine set_thermal_nden



subroutine thermal_nden(field,x,res)
!=========================================
! Soft photon field
! x energy in eV
! res number density in 1/(eV cm^3)
!=========================================
! Variables
 class(therm_phf), intent(in) :: field
 real(8), intent(in) :: x
 real(8), intent(out) :: res
 real(8) :: nd
 integer :: i
! Calculations 
 res=0d0
 do i=1,field%nf
  call therm_den(field%uden(i),field%temp(i),x,nd)
  res=res+nd
 end do 
 
end subroutine thermal_nden

 
 
subroutine therm_den(uden,temp,x,res)
! Variables
! energy density [eV/cm^3], temperature [eV], photon energy [eV]
 real(8), intent(in) :: uden, temp, x
 real(8), intent(out) :: res ! number density [1/cm^3 eV] in all directions
! inorm=1/a, where a=int(x^3/(exp(x)-1),x=0..infinity)=Pi^4/15
 real(8), parameter :: inorm=0.153989733820265028d0 
! Calculations
 res=x**2/(exp(x/temp)-1)
 res=res*uden*inorm/temp**4
end subroutine therm_den


end module thermal_photon_field

!!======================================
!! Test program
!!======================================
!
!program main
! use thermal_photon_field, only : therm_phf, set_thnden, thnden
! implicit none
! integer :: nf=4
! type(therm_phf) :: sph1
! real(8) :: x, res
! 
! nf=4
! call sph1%init(nf)
! sph1%temp=[3d0,3d-1,6d-3,2.35d-4]
! sph1%uden=[5d3,4d4,4d4,0.26d0]
! 
! call set_thnden(sph1)
! x=1d0
! call thnden(x,res)
! write(*,*) x, res
!
!end program main



