module airmass

use arveparams, only : dp,pi,d2r
use iovariables, only : slope

implicit none

public  :: get_airmass
private :: m
private :: F
private :: elev_corr

contains

!-------

subroutine get_airmass(j)

!Based on X. Yin (1997) Optical air mass: Daily integration and its applications, Meteorol. Atmos. Phys. 63, 227-233
!Jed Kaplan, EPFL, 2008

use statevars, only : sv,dayl
use soilstate_vars, only : surf
use metvarsmod, only : atm

implicit none

!pointers and arguments
integer , intent(in) :: j       !gridcell index

real(dp), pointer :: lat          !latitude (degrees)
real(dp), pointer :: delta        !solar declination (degrees)
real(dp), pointer :: elevation    !elevation above sea level (m)
real(dp), pointer :: mbar        !daytime mean optical air mass (unitless, 1 at equatorial noon)
real(dp), pointer :: mo          !air mass at cosine zenith angle maximum
real(dp), pointer :: mc          !air mass at cosine zenith angle medium
real(dp), pointer :: ml          !air mass at cosine zenith angle bottom quarter range point
real(dp), pointer :: rZ0         !zenith angle at solar noon (radians)

!parameters
real(dp), parameter :: m0  =  1.0_dp  !air mass at 0 degree solar zenith angle
real(dp), parameter :: m80 =  5.6_dp  !air mass at 80 degree solar zenith angle
real(dp), parameter :: m90 = 39.7_dp  !air mass at 90 degree solar zenith angle
real(dp), parameter :: cos80 = 0.173648177666930442_dp !cos(80) (degrees)
real(dp), dimension(3), target :: c00 = [ 0.008307_dp, 1.0213_dp, -0.01288_dp ] !air mass coefficients for solar zenith angle <=80 degrees
real(dp), dimension(3), target :: c80 = [ 0.03716_dp,  1.538_dp, -1.6973_dp   ] !air mass coefficients for solar zenith angle  >80 degrees
real(dp), parameter :: w   = 15.0_dp  !solar angular velocity (degrees hr-1)

!local variables
real(dp) :: rlat        !latitude (radians)
real(dp) :: rdelta      !solar declination (radians)
real(dp) :: t1          !number of hours between sunrise/sunset and solar noon (hr)
real(dp) :: t80         !solar hour corresponding to the 80 degree zenith angle
real(dp) :: Z           !solar zenith angle (degrees)
real(dp) :: Zn          !lesser of solar zenith angle at sunset or at midnight (degrees)
real(dp) :: Z0          !zenith angle at solar noon (degrees)
real(dp) :: cosZ        !cosine solar zenith angle (fraction), used in calculation of instantaneous air mass
real(dp) :: sinlat
real(dp) :: coslat
real(dp) :: sindel
real(dp) :: cosdel
real(dp)                        :: a   !values in equation 2.6b
real(dp)                        :: b
real(dp), pointer, dimension(:) :: c
real(dp) :: tmp1
real(dp) :: tmp2
real(dp) :: tmp3
real(dp) :: tinv
real(dp) :: rZn

!point pointers
lat             => sv(j)%lat
delta           => surf%sdecl
elevation       => slope(j)%elevation
mbar            => atm%mbar
mo              => atm%mc
mc              => atm%ml
ml              => atm%mo
rZ0             => surf%rZ0

!----------------
!These are parameters, below shows how the values are calculated

!cos80 = cos(80._dp * d2r)

!c00(1) = 0.008307_dp
!c00(2) = (m0 - m80) * (c00(1) + 1.0_dp) * (c00(1) + cos80) / (cos80 - 1.0_dp)
!c00(3) = m0 - c00(2) / (c00(1) + 1.0_dp)

!c80(1) = 0.037160_dp
!c80(2) = (m90 - m80) * c80(1) * (c80(1) + cos80) / cos80
!c80(3) = m90 - c80(2) / c80(1)

!----------------

!calculate daily mean air mass (mbar)

if (dayl(1) <= 3600.0_dp) then  !NOTE: this was changed to one hour since our smallest timestep is one hour. JM Dec 15 08
  mbar = m90
  mc   = m90
  ml   = m90

else

  !basic setup
  rlat   = d2r * lat
  rdelta = d2r * delta

  sinlat = sin(rlat)
  sindel = sin(rdelta)
  coslat = cos(rlat)
  cosdel = cos(rdelta)

  !------

  !Eqn. 2.5
  if (abs(lat - delta) < 90.0_dp .and. abs(lat + delta) >= 90.0_dp) then
   t1 = 12.0_dp
  else
   t1 = (12.0_dp / pi) * acos(-tan(rlat) * tan(rdelta))
  end if

  tinv = 1._dp / t1

  !Eqn. 2.9
  if (abs(lat + delta) >= 90.0_dp) then
    Zn = acos(sin(rlat) * sin(rdelta) - cos(rlat) * cos(rdelta)) / d2r
  else
    Zn = 90.0_dp
  end if

  !Eqn. 2.10
  if (abs(lat - delta) >= 90.0_dp) then
    Z0 = 90.0_dp
  else
    Z0 = abs(lat - delta)
  end if

  rZ0 = Z0 * d2r  !convert to radians
  rZn = Zn * d2r

  !-----
  b = coslat * cosdel

  if (t1 == 0._dp) then

    mbar = m90

  else if (abs(Zn) <= 80._dp) then

    c => c00
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else if (abs(Z0) >= 80._dp) then

    c => c80
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else

    t80 = 1._dp / w * acos((cos80 - sinlat * sindel) / (coslat * cosdel)) / d2r  !Eqn. 2.8

    c => c00
    a = c(1) + sinlat * sindel
    tmp1 = F(t80,a,b,c)

    c => c80
    a = c(1) + sinlat * sindel
    tmp2 = F(t1,a,b,c)

    c => c80
    a = c(1) + sinlat * sindel
    tmp3 = F(t80,a,b,c)

    mbar = tinv * (tmp1 + tmp2 - tmp3)

  end if

  !---------
  !calculate instantaneous air mass at max, mid, and bottom quarter solar zenith angle (m0, mc, ml)

  Z = Z0
  cosZ = cos(Z * d2r)

  if (Z <= 80._dp) then
    c => c00
  else
    c => c80
  end if

  mo = m(cosZ,c)

  !--

  Z = (Z0 + Zn) / 2._dp
  cosz = (cos(rZ0) + cos(rZn)) / 2._dp

  if (Z <= 80._dp) then
    c => c00
  else
    c => c80
  end if

  mc = m(cosZ,c)

  !--

  Z = (Z0 + 3._dp * Zn) / 4._dp
  cosz = (cos(rZ0) + 3._dp * cos(rZn)) / 4._dp

  if (Z <= 80._dp) then
    c => c00
  else
    c => c80
  end if

  ml = m(cosZ,c)

end if

!correct calculated air mass for elevation

        mbar = elev_corr(mbar,elevation)
        mo = elev_corr(mo,elevation)
        mc = elev_corr(mc,elevation)
        ml = elev_corr(ml,elevation)

end subroutine get_airmass

!-----------------------------------

real(dp)function m(cosZ,c)

!Instantaneous air mass m, equation 2.1 in Yin, 1997

implicit none

real(dp),               intent(in) :: cosZ
real(dp), dimension(:), intent(in) :: c

m = c(2) / (c(1) + cosZ) + c(3)

end function m

!---------------------------------

real(dp) function F(t1,a,b,c)

!integral air mass function F, equation 2.6b in Yin, 1997
!section inside curly braces only - multiply result by 1/t1 to get mbar

implicit none

real(dp),               intent(in) :: t1
real(dp),               intent(in) :: a
real(dp),               intent(in) :: b
real(dp), dimension(:), intent(in) :: c

real(dp) :: wt1
real(dp) :: wpi
real(dp) :: e1
real(dp) :: e2

real(dp), parameter :: w   = 15.0_dp  !solar angular velocity (degrees hr-1)
real(dp), parameter :: rw  = d2r * w !solar angular velocity (radians hr-1)

wpi  = 180._dp / (pi * w)
wt1  = rw * t1

if (a > b) then

  F = wpi * c(2) / sqrt(a**2 - b**2) * acos((b + a * cos(wt1)) / (a + b * cos(wt1))) + c(3) * t1

else if (a < b) then

  e1 = sqrt((b + a) * (1._dp + cos(wt1))) + sqrt((b - a) * (1._dp - cos(wt1)))
  e2 = sqrt((b + a) * (1._dp + cos(wt1))) - sqrt((b - a) * (1._dp - cos(wt1)))

  F = wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2) + c(3) * t1

else

  F = wpi * c(2) / a * tan(wt1 / 2._dp) + c(3) * t1

end if

end function F

!------------------------------

real(dp) function elev_corr(m,elevation)

implicit none

real(dp), intent(in) :: m
real(dp), intent(in) :: elevation

elev_corr = m * exp(-elevation / 8000.0_dp)

end function elev_corr

!------------------------------

end module airmass
