module surfrad

!This code is based on the paper:
!X. Yin (1998) Temporally-aggregated atmospheric optical properties as a function of common climatic information:
!Systems development and application, Meteorol. Atmos. Phys. 68, 99-113
!Jed Kaplan, EPFL, 2008
!Adapted for ARVE JM July 28 08

!This needs to be called once a day before surfaceheatflux. The integrated SW value from this is then used in
!surfaceheatflux for the short timestep calculations. The PET that is used is the previous day's value.

implicit none

public :: surface_rad
public :: simplesurfrad

contains

subroutine surface_rad(j,i)

use arveparams, only : dp,pi,Pstd,daysec
use iovariables, only : ibuf
use metvarsmod, only : met,atm
use statevars, only : sv,dtime
use soilstate_vars, only : surf
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j        !gridcell index
integer, intent(in) :: i        !timestep

!pointers
real(dp), pointer :: Patm       !atmospheric pressure (Pa)
real(dp), pointer :: Pjj        !precipitation equitability index
real(dp), pointer :: mbar       !mean daily air mass
real(dp), pointer :: mc         !airmass at medium cosine Z angle
real(dp), pointer :: ml         !airmass at bottom-quarter cos Z angle
real(dp), pointer :: mo         !airmass at max cosine zenith angle
real(dp), pointer :: dsw_t      !top-of-atmosphere insolation (kJ m-2 d-1)
real(dp), pointer :: prec       !precipitation mm/day
real(dp), pointer :: cldf       !24 hour mean cloud cover fraction
real(dp), pointer :: direct     !direct-beam downwelling shortwave (kJ m-2 d-1)
real(dp), pointer :: diffuse    !diffuse downwelling shortwave (kJ m-2 d-1)
real(dp), pointer :: lastpet    !potential evapotranspiration from last day (mm/timestep)
logical, pointer :: polar       !polar night if true.
real(dp), pointer :: dswb       !24 hour total downwelling shortwave (KJ m-2)
real(dp), pointer, dimension(:) :: tmp  !monthly mean temp (C)

!variables
real(dp) :: tau   !direct insolation atmospheric turbidity factor
real(dp) :: zeta0 !diffuse insolation atmospheric turbidity factor
real(dp) :: x     !tropics indicator (tropical = 1, north (lowest T below 10 C) = 0, otherwise btwn 0 & 1)
real(dp) :: fm    !atmospheric transmittance function
real(dp) :: p     !relative atmospheric pressure (1=sea level)
real(dp) :: tcm   !minimum temperature in the year (used as tropics indicator)
real(dp) :: sun   !bright sunshine duration fraction, n/N (fraction)

!parameters
real(dp), parameter :: kp  = 0.500_dp !links absorption coeff. to trans. coeff.
real(dp), parameter :: kag = 3.300_dp
real(dp), parameter :: kan = 2.320_dp
real(dp), parameter :: kn  = 0.686_dp !cloud parameter
real(dp), parameter :: ag = 0.17_dp   !Surface shortwave albedo (average=0.17)

!point pointers
mbar            => atm%mbar
mc              => atm%mc
ml              => atm%ml
mo              => atm%mo
lastpet         => sv(j)%lastpet
dsw_t           => surf%dsw_t
Patm            => sv(j)%Patm
cldf            => met(j,0)%cldf
direct          => surf%ra_dp(1)
diffuse         => surf%ra_dp(2)
Pjj             => atm%Pjj
tmp             => ibuf(j)%tmp
prec            => met(j,0)%prec
polar           => sv(j)%polar
dswb            => met(j,0)%dswb

!---------------

if (polar .or. i == 1) then  !daytime only calculation

  p = min(1._dp,Patm / Pstd)

  tcm = minval(tmp) !find min temp in the year (used as a tropics indicator)

  if (tcm < 10.0_dp) then
    x = 0.0_dp
  else if (tcm > 20.0_dp) then !tropics
    x = 1.0_dp
  else
    x = sin(pi / 2.0_dp * (tcm / 10.0_dp - 1.0_dp))
  end if

  !--
  !sun is simply 1 - cloud cover
  sun = 1._dp - cldf

  !--
  !find direct insolation atmospheric turbidity factor

  tau = exp(-0.115_dp * p * ((2.15_dp - 0.713_dp * x + exp(-6.74_dp / (prec + 1._dp))) &
                    * exp(0.0971_dp * lastpet) - 0.65_dp * (1._dp - x) * Pjj))  !Eqn. 4.1

  !find atmospheric transmittance function
  fm = 0.01452_dp * (mbar + ml) * exp(1.403_dp * tau) - 0.1528_dp * mo + mc + 0.48700_dp * (mc - ml) + 0.2323_dp   !Eqn. 2.4 2nd term

  !direct beam downwelling shortwave (kJ m-2 d-1)
  direct = sun * tau**kp * dsw_t * tau**fm   !Eqn. 2.4

  !find diffuse insolation atmospheric turbidity factor
  zeta0 = 0.503_dp * exp(-1.2_dp * p * exp(-0.633_dp / (prec + 1._dp) - 0.226_dp * lastpet)) &
           * 3.3_dp**ag * 2.32_dp**(1._dp - sun) * (1._dp - 0.686_dp * (1._dp - sun))  !Eqn. 4.2

  !diffuse downwelling shortwave (kJ m-2 d-1)
  diffuse = zeta0 * kag**ag * kan**(1.0_dp - sun) * (1 - kn * (1.0_dp - sun)) * (tau**kp * dsw_t - direct)   !Eqn. 2.5

  !for gridded data, assign the total radiation (in kJ m-2 d-1) to met (used by phenology)
      dswb = diffuse + direct  !FLAG maybe not anymore, possibly extraneous JM 14.12.2010

  !convert to W m-2 per timestep
  if (polar) then
          direct = direct * 1000. / daysec     !polar day spreads radiation over whole 24 hr period
          diffuse = diffuse * 1000. / daysec
  else !normal day
          direct = direct * 1000. / dtime      !convert kJ m-2 day to W m-2 timestep-1
          diffuse = diffuse * 1000. / dtime
  end if

  !reset lastpet to 0
  lastpet = 0._dp

else  !night time set radiation to 0

        diffuse = 0._dp                       !no radiation when night
        direct = 0._dp
        dswb = 0._dp

end if

end subroutine surface_rad

!----------------------------------------

subroutine simplesurfrad(j,i)

!ARVE-Point
!This is a simple way to split the downwelling solar radiation at the surface into direct and diffuse segments. This is used when
!we have station data that only records the total downwelling shortwave radiation at the surface. The reference is: Wang et al. 2002
!Ecological Modelling, 00, 1-14. JM Oct 2 2008

!NOTE: for the Amazon dataset (it is in Wm-2) you will need to make this so that it does not convert the expected input
!in kJ m-2 d-1 into W m-2.

use arveparams, only : dp,solarc,daysec
use statevars, only : sv,dtime
use metvarsmod, only : met
use soilstate_vars, only : surf

implicit none

!arguments
integer, intent(in) :: j        !gridcell index
integer, intent(in) :: i        !timestep

!variables
real(dp) :: Kt          !parameter giving estimate of 'cloudiness'
real(dp) :: Sday        !temporary variable for the direct downwelling shortwave (kJ m-2 d-1)
real(dp) :: epsilon     !fraction of diffusive radiation in total

!pointers
real(dp), pointer :: direct     !direct-beam downwelling shortwave (kJ m-2 d-1)
real(dp), pointer :: diffuse    !diffuse downwelling shortwave (kJ m-2 d-1)
real(dp), pointer :: avgZ       !cosine of the zenith angle (rads)
real(dp), pointer :: dsw_b      !downwelling shortwave (KJ m-2 day-1) (1 is direct, 2 is diffuse)
logical, pointer :: polar      !polar night if true.

!point pointers
direct          => surf%ra_dp(1)
diffuse         => surf%ra_dp(2)
avgZ            => surf%avg_cZ
dsw_b           => met(j,0)%dswb
polar           => sv(j)%polar

!---------------
!begin calculations

if (polar .or. (i == 1 .and. dsw_b > 0._dp)) then

     !when the downwelling shortwave at the surface is read in from the station data,
     !it is put into the direct compartment of the dswb.

 if (polar) then
         direct = dsw_b * 1000. / daysec     !polar day spreads radiation over whole 24 hr period
 else
         direct = dsw_b * 1000. / dtime      !convert kJ m-2 day to W m-2 timestep-1
 end if


Sday = direct

Kt = Sday / (solarc * avgZ)

if (Kt <= 0.1_dp) then
   epsilon = 0.98_dp
else
   epsilon = max((0.91_dp + 1.154_dp * Kt - 4.936_dp * Kt**2 + 2.848_dp * Kt**3), 0.15_dp)
end if

!find the diffuse component of the total
diffuse = epsilon * Sday

!the remainder is the direct downwelling shortwave
direct = Sday - diffuse

else !night

direct = 0._dp                       !no radiation when night
diffuse = 0._dp

end if

end subroutine simplesurfrad

end module surfrad
