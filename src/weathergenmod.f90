module weathergenmod

!This module includes subroutines to calculate daily maximum and minimum temperature and cloud cover fraction
!based on an annual timeseries of monthly values of these variables
!The weather generator is based on the WGEN model (Richardson, 1981) with extension to use monthly summary
!data from Geng et al., 1986, and Geng and Auburn 1986.
!Additional statistical relationships for both temperature and cloudiness have been produced
!by J.O. Kaplan using global weather station datasets (GSOD and global synoptic cloud reports).

!Coded in 2007-2009 by Jed Kaplan and Joe Melton, ARVE Group, EPFL/UVic, jed.kaplan@epfl.ch

use arveparams, only : dp

implicit none

public  :: weathergen
public  :: rmsmooth
private :: meancv

!-------------------------------

!local module variables

real(dp) :: tmax_mn     !maximum temperature monthly mean (K)
real(dp) :: tmin_mn     !minimum temperature mothly mean (K)
real(dp) :: cldf_mn     !mean cloud fraction (fraction)

real(dp) :: tmax_cv     !coefficients of variation of above variables
real(dp) :: tmin_cv     !
real(dp) :: cldf_cv     !

contains

!---------------------------------------------------------------------------------------------------

subroutine weathergen(j,m,d)

use statevars,     only : sv
use arveparams,    only : mdays,tfreeze
use randomdistmod, only : ran_uni,ran_normal,ran_gamma
use iovariables,   only : ibuf,doweather
use metvarsmod,    only : met,dm

implicit none

!arguments
integer, intent(in) :: j                          !the pixel being calculated
integer, intent(in) :: d                          !day of the year
integer, intent(in) :: m                          !month number

!pointers to input variables
!raw input: monthly means
real(dp), pointer :: pre                            !monthly total precipitation amount (mm)
real(dp), pointer :: wet                            !number of days in month with precipitation (days)
!interpolated, smoothed pseudo-daily values
real(dp), pointer :: tmx                            !maximum temperture (C)
real(dp), pointer :: tmn                            !minumum temperture (C)
real(dp), pointer :: cld                            !cloud fraction (0=clear sky, 1=overcast) (fraction)
real(dp), pointer, dimension(:) :: met_res        !weather residuals

!local variables
real(dp), dimension(3) :: lag1             !previous day's weather residuals
!simulated actual daily output values (written to met at end)
logical   :: pday_local                    !precipitation status
real(dp)  :: tmin_local                    !24 hour minimum temperature (C)
real(dp)  :: tmax_local                    !24 hour maximum temperature (C)
real(dp)  :: prec_local                    !24 hour total precipitation (mm)
real(dp)  :: cldf_local                    !24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
real(dp), dimension(3) :: unorm            !vector of uniformly distributed random numbers (0-1)
real(dp) :: pbar                             !mean amount of precipitation per wet day (mm)
real(dp) :: pwd                              !transition probability of a wet day following a dry day (fraction)
real(dp) :: pww                              !transition probability of a wet day following a wet day (fraction)
real(dp) :: alpha                            !shape parameter for the precipitation amount function
real(dp) :: beta                             !shape parameter for the precipitation amount function
real(dp) :: u                                !uniformly distributed random number (0-1)
real(dp) :: wetf                             !fraction of wet days in month (fraction)
integer :: i

!parameters
real(dp), parameter :: pmin = 2.16_dp / 0.83_dp !minimum value for pbar when using Geng linear relationship,&
                                              !below this value we use a 1:1 line, see below
real(dp), parameter :: small = 5.e-5

real(dp), dimension(9), parameter :: corva =   [ 0.567_dp, 0.086_dp,-0.002_dp, &   !lag day-1 correlation coefficients
                                                    0.253_dp, 0.504_dp,-0.050_dp, &
                                                   -0.006_dp,-0.039_dp, 0.244_dp ]

real(dp), dimension(3,3), parameter :: cor_a = reshape(corva,[3,3])

real(dp), dimension(9), parameter :: corvb =   [ 0.781_dp, 0.000_dp, 0.000_dp, &   !current day correlation coefficients
                                                    0.328_dp, 0.637_dp, 0.000_dp, &
                                                    0.238_dp,-0.341_dp, 0.873_dp ]

real(dp), dimension(3,3), parameter :: cor_b = reshape(corvb,[3,3])

!Point pointers
pre => ibuf(j)%pre(m)
wet => ibuf(j)%wet(m)
tmn => dm(j)%tmn(d)
tmx => dm(j)%tmx(d)
cld => dm(j)%cld(d)
met_res => met(j,0)%met_res

!--------------------
if (doweather) then     !perform the weather generator

!Assign previous days rain state
!FLAG met position 1 is the next days weather! not previous (met position -1), what is intended here? JM 
! I changed it to be -1. JM 16.04.2011
pday_local  =  met(j,-1)%pday

!Begin calculations

!1) Precipitation occurrence

!if there is precipitation this month, calculate the precipitation state for today

if (wet > 0._dp .and. pre > 0._dp) then

  !calculate fraction of wetdays in the month
  wetf = wet / mdays(m)

  !calculate transitional probabilities for dry to wet and wet to wet days
  !Relationships from Geng & Auburn, 1986, Weather simulation models based on summaries of long-term data
  pwd = 0.75_dp * wetf
  pww = 0.25_dp + pwd

  !determine the precipitation state of the current day using the Markov chain approach
  call ran_uni(u)

  !the precip status of the current day is conditioned on the status of the previous day
  if (pday_local) then   !previous day's precip saved from last call to this subroutine

    if (u - pww > 0._dp) then
      pday_local = .false.
    else
      pday_local = .true.
    end if

  else if (u - pwd > 0._dp) then
    pday_local = .false.
  else
    pday_local = .true.
  end if

!-----

!2) precipitation amount

  if (pday_local) then  !today is a wet day, calculate the rain amount

    !calculate parameters for the distribution function of precipitation amount
    pbar = pre / wet

    if (pbar > pmin) then
      beta = -2.16_dp + 1.83_dp * pbar
    else
      beta = pbar
    end if

    beta  = max(beta,small)   !put here to avoid infinity values of alpha
    alpha = pbar / beta

    prec_local = ran_gamma(alpha,beta,.true.)

  else

    prec_local = 0._dp

  end if

else

  pday_local = .false.
  prec_local = 0._dp

end if

!-----

!3) temperature min and max, cloud fraction

!calculate a baseline mean and cv for today's weather dependent on precip status
!want temp in K here so convert from C.
call meancv(pday_local,tmn+tfreeze,tmx+tfreeze,cld)

!use random number generator for the normal distribution
do i = 1,3
  call ran_normal(unorm(i))
end do

!load yesterday's residuals into lag1
lag1 = met_res

!calculate today's residuals for weather variables
met_res = matmul(cor_a,lag1) + matmul(cor_b,unorm)  !Richardson 1981, eqn 5; WGEN tech report eqn. 3

tmax_local = tmax_mn * (met_res(1) * tmax_cv + 1._dp)      !WGEN tech report eqn. 13
tmin_local = tmin_mn * (met_res(2) * tmin_cv + 1._dp)
cldf_local = cldf_mn * (met_res(3) * cldf_cv + 1._dp)

!correct tmax if it is less than tmin and vice versa
tmax_local = max(tmax_local,tmin_local)
tmin_local = min(tmax_local,tmin_local)  !FLAG added fix by JM Mar 10 09

!set cldf to be inside its range
cldf_local = max(cldf_local,0._dp)
cldf_local = min(cldf_local,1._dp)

! Put the weather variables into the met array
met(j,:)%pday = eoshift(met(j,:)%pday,1,pday_local)
met(j,:)%tmin = eoshift(met(j,:)%tmin,1,tmin_local-Tfreeze)  !make tmin in deg. C
met(j,:)%tmax = eoshift(met(j,:)%tmax,1,tmax_local-Tfreeze)  !make tmax in deg C.
met(j,:)%prec = eoshift(met(j,:)%prec,1,prec_local)
met(j,:)%cldf = eoshift(met(j,:)%cldf,1,cldf_local)

else

  !no weather generated!
  !use pseudo-daily smoothed values from rmsmooth, flag is set in joboptions file.

met(j,:)%tmin = eoshift(met(j,:)%tmin,1,tmn)  !tmin in deg. C
met(j,:)%tmax = eoshift(met(j,:)%tmax,1,tmx)  !tmax in deg C.
met(j,:)%prec = eoshift(met(j,:)%prec,1,pre/mdays(m))  !put the rain evenly throughout the month
met(j,:)%cldf = eoshift(met(j,:)%cldf,1,cld)

end if


end subroutine weathergen

!-------------------------------------------------------------------

subroutine meancv(pday,tmn,tmx,cld)
!calculate the mean and CV for a single day value of tmax, tmin, and cloud fraction
!requires temperatures in K

implicit none

!arguments
logical,  intent(in) :: pday  !precipitation status
real(dp), intent(in) :: tmn   !temperature (K)
real(dp), intent(in) :: tmx
real(dp), intent(in) :: cld   !fraction (0-1)

!parameters
!coefficients for the temperature wet day:dry day mean split (based on temperatures in K)
real(dp), dimension(3), parameter :: tmnc = [ -7.424e+01,  4.573e-01, -6.944e-04 ]
real(dp), dimension(3), parameter :: tmxc = [ -3.794e+02,  2.604e+00, -4.440e-03 ]

!coefficients for the temperature mean:CV
real(dp), dimension(3), parameter :: cvtmnc = [ 2.699e-01, -1.362e-03, 1.596e-06 ]
real(dp), dimension(3), parameter :: cvtmxc = [ 2.399e-01, -1.220e-03, 1.517e-06 ]

real(dp), parameter :: a = 2.414532_dp    !coefficient for the cloud dry day mean split
real(dp), parameter :: b = 2.436341_dp    !coefficient for the cloud wet day mean split
real(dp), parameter :: cc = -1.0479262_dp  !slope of the regression line through the cloud wet day mean:CV relationship

!coefficients for the regression line through the cloud dry day:CV relationship
real(dp), dimension(3), parameter :: cd = [ 4.084e+00, 4.165e-03, 1.694e-01 ]

!local variables
real(dp) :: tmaxdiff
real(dp) :: tmindiff

!----------------------------------

tmaxdiff = tmxc(1) + tmxc(2) * tmx + tmxc(3) * tmx**2
tmindiff = tmnc(1) + tmnc(2) * tmn + tmnc(3) * tmn**2

if (pday) then   !wet day

  !---tmax---
  tmax_mn = tmx - tmaxdiff                                              !mean
  tmax_cv = cvtmxc(1) + cvtmxc(2) * tmx + cvtmxc(3) * tmx**2            !CV

  !---tmin---
  tmin_mn = tmn - tmindiff                                              !mean
  tmin_cv = cvtmnc(1) + cvtmnc(2) * tmn + cvtmnc(3) * tmn**2            !CV

  !---cloud---
  cldf_mn = 1._dp + (cld - 1._dp) / (b * cld + 1._dp)                      !mean
  cldf_cv = cc * cld                                                    !CV

else   !dry day

  !---tmax---
  tmax_mn = tmx + tmaxdiff                                              !mean
  tmax_cv = cvtmxc(1) + cvtmxc(2) * tmx + cvtmxc(3) * tmx**2            !CV

  !---tmin---
  tmin_mn = tmn + tmindiff                                              !mean
  tmin_cv = cvtmnc(1) + cvtmnc(2) * tmn + cvtmnc(3) * tmn**2            !CV

  !---cloud---
  cldf_mn = - (a * cld) / (cld - a - 1._dp)                          !mean
  cldf_cv = 1._dp /(cd(1) * (cld + cd(2))) - cd(3)                   !CV

end if

!write(0,'(a,l4,4f9.3)')'meancv',pday,tmn,tmin_mn,tmin_cv,tmin_mn*tmin_cv
!write(0,'(a,l4,4f9.3)')'meancv',pday,tmx,tmax_mn,tmax_cv,tmax_mn*tmax_cv

end subroutine meancv

!---------------------------------------------------------------------------------------

subroutine rmsmooth(m,dmonth,bc,r)

!Iterative, mean preserving method to smoothly interpolate mean data to pseudo-sub-timestep values
!From Rymes, M.D. and D.R. Myers, 2001. Solar Energy (71) 4, 225-231

implicit none

!arguments
real(dp), dimension(:), intent(in)  :: m          !vector of mean values at super-time step (e.g., monthly), minimum three values
integer,  dimension(:), intent(in)  :: dmonth         !vector of number of intervals for the time step (e.g., days per month)
real(dp), dimension(2), intent(in)  :: bc         !boundary conditions for the result vector (1=left side, 2=right side)
real(dp), dimension(:), intent(out) :: r          !result vector of values at chosen time step

!parameters
real(dp), parameter :: ot = 1._dp / 3

!local variables
integer :: n
integer :: ni
integer :: a
integer :: b
integer :: i
integer :: j
integer :: k
integer :: l
integer, dimension(size(r)) :: g
real(dp) :: ck

!----------

n  = size(m)
ni = size(r)

!initialize the result vector
i = 1
do a = 1,n
  j = i
  do b = 1,dmonth(a)
    r(i) = m(a)
    g(i) = j
    i = i + 1
  end do
end do

!iteratively smooth and correct the result to preserve the mean

!iteration loop
do i = 1,ni

  do j = 2,ni-1
    r(j) = ot * (r(j-1) + r(j) + r(j+1))   !Eqn. 1
  end do

  r(1)  = ot * (bc(1)   + r(1)  +  r(2))    !Eqns. 2
  r(ni) = ot * (r(ni-1) + r(ni) + bc(2))

  j = 1
  do k = 1,n                               !calculate one correction factor per super-timestep

    a = g(j)                                       !index of the first timestep value of the super-timestep
    b = g(j) + dmonth(k) - 1                           !index of the last timestep value of the super-timestep

    ck = sum(m(k) - r(a:b)) / ni           !Eqn. 4

    do l = 1,dmonth(k)                         !apply the correction to all timestep values in the super-timestep
      r(j) = r(j) + ck
      j = j + 1
    end do

  end do
end do

end subroutine rmsmooth

!-------------------------------

end module weathergenmod
