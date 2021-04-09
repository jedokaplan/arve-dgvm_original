module diurnaltempmod

use arveparams, only : dp,pi,d2r,daysec

implicit none

public :: diurnaltemp
public :: instantemp


real(dp) :: sunrise     !time of (relative to solar noon (rad))
real(dp) :: sundown        !time of (relative to solar noon (rad))
real(dp) :: peakt       !time of peak temperature
real(dp) :: t0                !temperature at the sunset hour
real(dp) :: r                !difference between max temp and temp at sunset (C)
real(dp) :: a           !amplitude of temp change (max - min) in that day (C)
real(dp) :: b                !eqn 10 from Cesaraccio paper (ref below)
real(dp) :: hni         !length of the night (hrs)
real(dp) :: tfollow     !the next days minimum temperature (assigned based on grid or point mode)

!-------
contains

! Formulas used herein are adapted from Cesaraccio, C. Spano, D., Pierpaolo, D., Snyder R. Int. J. Biometerol.
! 2001 (45) 161-169
!NOTE there is an error in the paper,  7a should have a sin after the alpha.
!--

subroutine diurnaltemp(j)

use statevars,   only : sv,gridded,dayl
use metvarsmod,  only : met,dm,atm

!arguments
integer, intent(in) :: j          !gridcell index

!pointers
real(dp), pointer :: tmin                       !min temperature (C)
real(dp), pointer :: tmax                       !max temperature (C)
real(dp), pointer :: tnextmin                      !minimum temperature of the following day
real(dp), pointer :: tmean                      !derived mean temperature for the 24 hour period (C)
real(dp), pointer, dimension(:) :: temp         !outputs the night/day temperature to here!
real(dp), pointer :: tnext                       !minimum temperature of the following day.
real(dp), pointer, dimension(:) :: runtemp      !holds one month of temp data for mean monthly calc.

!parameter
real(dp), parameter :: tfpk = 1./6.                 !delay between solar noon and peak temperature (fraction)

!local variables
real(dp) :: hdl,hdlnext         !half day length (sec) (this day and next)
real(dp) :: ti1                 !midnight till sunup
real(dp) :: ti                          !sundown till midnight
real(dp) :: tam                 !daytime till noon
real(dp) :: tpm                 !daytime post noon till sundown
real(dp) :: morn                !sunrise - peakt
real(dp) :: tday                !integrated mean daytime temperature
real(dp) :: tnight              !integrated mean nighttime temperature
real(dp) :: sunrise_next        !relative to solar noon (rad)

!point pointers to global data
tmin    => met(j,0)%tmin
tmax    => met(j,0)%tmax
tmean   => met(j,0)%temp24
tnextmin   => met(j,1)%tmin
temp    => met(j,0)%temp
tnext   => dm(j)%tnext
runtemp => atm%runtemp

!------
!initial assignments

!ARVE grid looks for the next days min temp in a different area than ARVE point
if (gridded) then
  tfollow = tnextmin
else
  tfollow = tnext
end if

!these are prep for other calculations below
!tfollow = minimum temperature of the following day
t0 = tmax - 0.39 * (tmax - tfollow)
a = tmax - tmin
r = tmax - t0

!find sunrise and sundown
hdl = 0.5 * dayl(1) / 3600.  !hrs
sunrise = 12. - hdl  !hrs, fixes time from noon till sun up by using half the day length
sundown = 12. + hdl

!find next days sunrise
hdlnext = 0.5 * dayl(2) / 3600.
sunrise_next = 12 - hdlnext

!this gets the time of peak temperature
peakt = 12. + 2. * hdl * tfpk

if (dayl(1) <= 3600.) then
  !FLAG this could be a bit of an ugly way of doing this... it could bias the weather to be mostly
  !hot since the day (which is the longest) will be set to the max temp of the day. The main reason for
  !this problem is that for gridded data we do not have the mean daily temp calced out in weathergen
  !like we do for the max and min. Will think of better ways to do this -JM Oct 29 08
      !FLAG check on this again after we fix weathergen JM Jan 10 09
  tday = tmax
  tnight = tmin

else if (dayl(1) < daysec) then !daylength is less than 24 hours

        !has a night and a day, calculate night first

          !find the length of the night
          hni = (43200. - 0.5 * dayl(1) + 43200. - 0.5 * dayl(2))  / 3600. !hrs

          b = (tfollow - t0) / sqrt(hni) !eqn 10 from Cesaraccio paper

          ti  = t0 * sundown                                            !sundown
          ti1 = t0 * (sundown + hni) + 2./3. * b * (hni**(3./2.))   !sunrise (next morn)

          tnight = (ti1 - ti) / hni

          if (dayl(1) > 0.) then  !regular night and day

                    !morning integral (ti is at sunrise, ti1 is at temperature peak time)
                    morn = sunrise - peakt

                    ti  = (tmin * sunrise) + (1. / pi) * (2. * a * morn)
                    ti1 = (tmin * peakt)   + (1. / pi) * (2. * a * morn * cos(pi/2. * (sunrise - peakt) / morn))

                    tam = (ti1 - ti) / (-morn)

                    !afternoon integral (ti is at temperature peak time, ti1 is at sundown)
                    ti  = t0 * peakt   - (1. / pi) * 8. * r * cos(pi / 8. * (-4.))
                    ti1 = t0 * sundown - (1. / pi) * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))

                    tpm = (ti1 - ti) / (sundown - peakt)

                    tday = (tam + tpm) / 2.

          else

                tday = tnight      !only night, day = night

          end if

else !no night, only day

          !morning integral (ti is at sunrise, ti1 is at temperature peak time)
          morn = sunrise - peakt

          ti  = tmin * sunrise + 1. / pi * (2. * a * morn)
          ti1 = tmin * peakt   + 1. / pi * (2. * a * morn * cos(pi / 2. * (sunrise - peakt) / morn))

          tam = (ti1 - ti) / (-morn)

          !afternoon integral (t10 is at temperature peak time, ti1 is at sundown)
          ti  = t0 * peakt   - 1. / pi * 8. * r * cos(pi / 8. * (-4.))
          ti1 = t0 * sundown - 1. / pi * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))


          tpm = (ti1 - ti) / (sundown - peakt)

          tday = (tam + tpm) / 2.

          tnight = tday

end if

temp(1) = tday
temp(2) = tnight

if (gridded) then !gridded data has no mean temperature for the day so find that here.
        tmean = temp(1) * (dayl(1) / daysec) + temp(2) * ((daysec - dayl(1)) / daysec)

end if

!write(*,'(a,f12.4,a,2f12.4)')'================tday',temp(1),'night              ',temp(2),tfollow

!set the value for the running mean
runtemp(1) = runtemp(1) + tday
runtemp(2) = runtemp(2) + tnight

end subroutine diurnaltemp

!---------------------------------------------------------------------------

function instantemp(i,j,count)
  !Added May 08 to match the temperature variation to the radiation variation for the short timestep
  !Based on the same paper as diurnal temp variation. Coded by Joe Melton

  !Added in wind speed calculation. This wind speed is a diurnal variation closely following the temperature.
  !Obviously not a rigorous approach but workable. Added Oct 08 08 JM

use metvarsmod, only : met,dm,atm
use statevars, only : dt,sv,doy,dayl

implicit none

!arguments
integer, intent(in) :: i        !timestep
integer, intent(in) :: j        !grid-cell index
integer, intent(in) :: count

!pointers
real(dp), pointer :: tmax               !max temperature (C)
real(dp), pointer :: tmin               !min temperature (C)
real(dp), pointer :: tmean              !mean temperature (C)
real(dp), pointer :: winds              !parameterized wind speed (m s-1)
real(dp), pointer :: tmean_next          !mean temp for the following day

!local variables
real(dp) :: instantemp
real(dp) :: hr
real(dp) :: w0
real(dp) :: bw
real(dp) :: rw
real(dp) :: c                        !like the b variable but for the fix I made (JM)
real(dp), save :: tsave                !save the sunset temperature to use as the starting point for the night
real(dp), save :: wsave                !winds as for temp above

!parameter
real(dp), parameter :: aw = 0.6_dp             !amplitude of diurnal wind speed variation (m s-1) (Dai & Deser JGR 1999)
real(dp), parameter :: wmean = 3.28           !global 'mean' wind speed (m s-1) (Archer & Jacobson JGR 2005)
real(dp), parameter :: wmin = wmean - aw      !global 'min' wind speed (m s-1)

!point pointers
tmin    => met(j,0)%tmin
tmax    => met(j,0)%tmax
tmean   => met(j,0)%temp24
winds   => atm%winds
tmean_next => dm(j)%tmean_next

!-------
!begin calculations

w0 = (wmean + aw) - 0.39_dp * (wmean + aw - wmin)
rw = wmean + aw - w0

if (dayl(1) <= 3600._dp) then  !polar night

  sunrise = 11.5 !(because the daylength is only 1 hour)
  hr = sunrise + (dt * (0.5_dp + real(count - 1)) / 3600._dp)

  instantemp = tmean * (1.5_dp - hr / 24._dp) + tmean_next * (1._dp - (1.5_dp - hr / 24._dp))
  winds = wmean * (1.5_dp - hr / 24._dp) + wmean * (1._dp - (1.5_dp - hr / 24._dp))

else   !not polar night

  if (i == 1) then !normal day

    hr = sunrise + (dt * (0.5_dp + real(count - 1)) / 3600._dp)

  else

    if (dayl(1) <82800._dp) then !normal night
      hr = sundown + (dt * (0.5_dp + real(count - 1)) / 3600._dp)
    else  !night of polar day (model defined as 1hr)
      hni = 0.5_dp
      sundown = 23._dp
      hr = 23._dp +  (dt * (0.5_dp + real(count - 1)) / 3600._dp)
    end if

  end if

  if (i == 1) then !normal day

      if (hr < peakt) then  !between sunrise and peak temperature time

        instantemp = tmin + a * sin(((hr - sunrise) / (peakt - sunrise)) * pi * 0.5)
        winds = wmin + aw * sin(((hr - sunrise) / (peakt - sunrise)) * pi * 0.5)

      else if (hr > peakt) then !between peak temperature time and sunset

        instantemp = t0 + r * sin(pi * 0.5 + (hr - peakt) * 0.25 * pi * 0.5)
        winds = w0 + rw * sin(pi * 0.5 + (hr - peakt) * 0.25 * pi * 0.5)
        tsave = instantemp
        wsave = winds

      else

        instantemp = tmax
        winds = wmean + aw

      end if

  else !night

      !NOTE: This was not working well, it was creating fairly large jumps between the last day
      !and first timestep of the night. I have changed it to look at the last temperature at sunset
      ! (tsave) and to take the slope of the function using that. This seems to allow for much smoother transitions.
      !JM 10.06.2010

      !OLD WAY: between sunset and sunrise of the next day.
      !instantemp = t0 + b * (hr - sundown)**0.5   !eqn 10 from Cesaraccio paper

      !NEW WAY
      c = (tfollow - tsave) / sqrt(hni)
      instantemp = tsave + c * (hr - sundown)**0.5

      !also done for wind
      bw = (wmin - wsave) / sqrt(hni)

      winds = wsave + bw * (hr - sundown)**0.5
  end if

end if

! write(*,'(a8,3i5,2f12.4)')'wind',doy,i,count,winds,instantemp

end function instantemp

end module diurnaltempmod
