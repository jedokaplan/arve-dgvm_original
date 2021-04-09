module coszen

! Compute cosine solar zenith angle using day value where a round day (such as 213.0)
! refers to 0z at Greenwich longitude. This is adapted from CSM 3.0 - JM Oct 6 08
! Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
! Processes in Meterology and Climatology, Elsevier Scientific
! Publishing Company, New York  p. 57, p. 62,63.

implicit none

public :: avgZen
public :: disZen

contains

!----

subroutine avgZen(j)

!use for diurnally averaged values.

use arveparams, only : dp,dayspy,pi,d2r
use statevars, only : sv,doy
use soilstate_vars, only : surf

implicit none

!arguments
integer, intent(in) :: j        !gridcell index

!pointers
real(dp), pointer :: lat        !latitude (degrees)
real(dp), pointer :: avg_cZ     !diurnally averaged cos of the zenith angle (rad) of the incident beam
real(dp), pointer :: delta      !declination (deg)

!variables
real(dp) :: rdelta              !Solar declination angle  (radians)
real(dp) :: arg
real(dp) :: calday              !calendar day
real(dp) :: frac                !Daylight fraction
real(dp) :: tsun                !temporary term in diurnal averaging
real(dp) :: rlat                !latitude (radians)

!point pointers
lat             => sv(j)%lat
avg_cZ          => surf%avg_cZ
delta           => surf%sdecl

!-------
!begin calculations
rlat = lat * d2r
calday = real(doy)

!Solar declination in radians (precalculated in either ARVE-Point (daylength) or ARVE-Grid (dayins):
  rdelta = delta * d2r

!for diurnal averaging, compute the average local cosine solar
! zenith angle using formulas from paltridge and platt 1976  p. 57, p. 62,63.

arg = -(sin(rlat) / cos(rlat) * (sin(rdelta) / cos(rdelta)))

if (arg < -1._dp) then
   frac = 1._dp
else if (arg > 1._dp) then
   frac = 0._dp
else
   frac = (1._dp / pi) * acos(arg)
end if

tsun = pi * frac

if (tsun > 0._dp) then
   avg_cZ =  sin(rlat) * sin(rdelta) + (cos(rlat) * cos(rdelta) * sin(tsun)) / tsun
else
   avg_cZ = 0._dp
end if


end subroutine avgZen

!----------------------------------

subroutine disZen(j)

!use for discrete time periods in the day.

use arveparams, only : dp,dayspy,pi,d2r,npft
use statevars, only : sv,counter,dt,counter_lim,doy,dtime
use pftparametersmod, only : prm
use soilstate_vars, only : surf

implicit none

!arguments
integer, intent(in) :: j        !gridcell index

!pointers
real(dp), pointer :: lat        !latitiude (degrees)
real(dp), pointer :: zen        !zenith angle for this timestep (rads)
real(dp), pointer :: delta     !declination (deg)

!variables
!real(dp) :: theta               !Earth orbit seasonal angle (radians)
real(dp) :: rdelta               !Solar declination angle  (radians)
real(dp) :: calday              !calendar day
real(dp) :: rlat                !latitude (radians)
real(dp) :: hdl                 !half day length (hrs)
real(dp) :: sunrise             !time of sunrise in hrs
real(dp) :: zed                 !time of day in hours
real(dp) :: phi                 !adjusted calendar day

!point pointers
lat             => sv(j)%lat
zen             => surf%zen
delta           => surf%sdecl

!-------
!begin calculations
rlat = lat * d2r

!find the what the fraction of day we are in this timestep.
hdl = 0.5 * dtime / 3600.  !hrs
sunrise = 12. - hdl
zed = sunrise + (dt * (0.5 + real(counter - 1)) / 3600.)

calday = real(doy) + zed / 24._dp !calendar day plus the fraction of the present day.

!Solar declination in radians:
rdelta = delta * d2r  !NOTE since I changed this to the Berger delta, we lose the ability to have a changing delta as the day goes on (yet gain
                      !the ability to accurately do paeloruns). This is likely a very small change through the day.JM Oct 30 08

  !below is a different calculation for delta. However this calculation is only valid for modern simulations.
  !find the Earth orbit seasonal angle in radians
  !theta = 2._dp * pi * calday / dayspy
  !delta = 0.006918_dp - 0.399912_dp * cos(theta) + 0.070257_dp * sin(theta) - &
  !        0.006758_dp * cos(2._dp * theta) + 0.000907_dp * sin(2._dp * theta) - &
  !        0.002697_dp * cos(3._dp * theta) + 0.001480_dp * sin(3._dp * theta)

! Calday is the calender day for Greenwich, including fraction
! of day; the fraction of the day represents a local time at
! Greenwich; to adjust this to produce a true instantaneous time
! For other longitudes, we must correct for the local time change:
! local time based on the longitude and day of year
! then compute the local cosine solar zenith angle
! HOWEVER, we do not treat individual timezones discretely In ARVE, we keep only day/night timesteps
! thus day is day at any point on the globe (not moving as is the case in reality). We do not need to
! correct for the longitude as a result.

phi = calday !+ (real(lon(i)-1)/real(plon))

zen = sin(rlat) * sin(rdelta) - cos(rlat) * cos(rdelta) * cos(2._dp * pi * phi)

end subroutine disZen

end module coszen
