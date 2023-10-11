module phys_const
!
! Physical and mathematical constants
!
 implicit none
 private
 save
 
 
 real(8), parameter, public :: psec=3.08567758d18 ! cm, parsec
 real(8), parameter, public :: year=3.155695200d7 ! s,  Gregorian year in seconds 3600*24*365.2425
 real(8), parameter, public :: pi=3.14159265358979324d0 ! pi number 
 real(8), parameter, public :: millibarn=1d-27 ! [cm^2/millibarn]
 real(8), parameter, public :: c_light=2.99792458d10 ! cm/s
 real(8), parameter, public :: erg_eV=1.602176487d-12
 real(8), parameter, public :: me_ev=0.510998928d6 ! eV, electron mass
 real(8), parameter, public :: e_charge=4.80320425d-10 ! charge of electron in CGS system
 real(8), parameter, public :: mprot=0.938272046d9 ! eV, proton mass
 real(8), parameter, public :: hbar=1.054571800d-27 ! erg*s, plank constant
 real(8), parameter, public :: alpha_fs=0.0072973525664d0 ! e^2/(hbar*c) fine-structure constant
 real(8), parameter, public :: kbol_evk=8.6173303d-5 ! Boltzmann constant eV/Kelvin
 real(8), parameter, public :: kbol_ergk=1.38064852d-16 ! Boltzmann constant erg/Kelvin
 real(8), parameter, public :: solar_lum=3.828d33 ! solar luminosity in erg/s
 
end module phys_const
