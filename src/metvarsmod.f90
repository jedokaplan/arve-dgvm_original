module metvarsmod

! This module contains the MET, DM, and ATM types
! Contains the read-in (ARVE-Point) and derived meteorological variables (ARVE-Grid)
! called by modules below the main driver hierarchically
! (calcdaily, calcannual, etc.)
! Also contains some average and summed meteorological and state of the atmosphere variables

!-------------------------------------------------------------
! ARVE-Point:
! read-only variables that are externally prescribed
! these are daily values from weather station data
! read directly into metvars (met)

! ARVE-GRID:
! monthly data for arve_grid is read into iovariables (ibuf), this is then smoothed to
! daily values using rmsmooth (in weathergenmod) and put in dmetvars(dm). 'Weather' is added and
! The weather generator puts the daily values here in metvars (met)

!=======================================

use arveparams, only : dp

implicit none

! DM
! Meteorological variables that are derived (in the program) from the ibuf variables

type dmetvars

  !vectors for a year of annual, pseudo-daily, smoothed values
  real(dp), dimension(365) :: tmx          ! minimum temperature deg C)
  real(dp), dimension(365) :: tmn          ! maximum temperature (deg C)
  real(dp), dimension(365) :: cld          ! cloud cover (fraction)
  real(dp) :: p_wet_m                      ! precipitation of the wettest month (mm)
  real(dp) :: p_dry_m                      ! precipitation of the driest month (mm)
  real(dp) :: p_mo_avg                     ! mean monthly precip (mm)
  logical :: init                          ! initialization flag for the weather generator
  real(dp) :: prec_ra                      ! precipitation running sum over 30 days (mm)
  real(dp), dimension(30) :: mprec         ! quasi-monthly precipitation total (mm)
  real(dp) :: tnext                        ! minimum temperature of the following day.
  real(dp) :: tmean_next                   ! mean temp for the following day

end type dmetvars

type(dmetvars), allocatable, target, dimension(:) :: dm

!=======================================
! MET
! Met has three positions, -1 = previous day, 0 = current day, 1 = next day
type metvars

  logical  :: pday         !precipitation status
  real(dp) :: tmin         !24 hour minimum temperature (C)
  real(dp) :: tmax         !24 hour maximum temperature (C)
  real(dp) :: temp24       !24 hour mean temperature (C)
  real(dp) :: prec         !24 hour total precipitation (mm)
  real(dp) :: cldf         !24 hour mean cloud cover fraction
  real(dp) :: dswb         !24 hour integral downwelling shortwave (KJ m-2)

  real(dp), dimension(3) :: met_res        !residual values (see weathergen)

  real(dp), dimension(2)   :: temp          !vector for mean temperature of the day and night timesteps: 1=day, 2=night (deg C)

end type metvars

type(metvars), allocatable, target, dimension(:,:) :: met

!=======================================
! ATM
! State of atmosphere varibles
! these values are NOT retained year over year!

type atm_variables

        real(dp)                  :: Pjj              !precipitation equitability index
        real(dp)                  :: mbar             !mean daily air mass
        real(dp)                  :: mc               !airmass at medium cosine Z angle
        real(dp)                  :: ml               !airmass at bottom-quarter cos Z angle
        real(dp)                  :: mo               !airmass at max cosine zenith angle
        real(dp)                  :: winds            !instantaneous 'wind speed' diurnally varying (m s-1)
        real(dp)                  :: Tdew             !dewpoint temperature (K)
        real(dp), dimension(365)  :: daytempv         !vector of daytime averaged temperatures
        real(dp)                  :: twm              !maximum monthly temperature for the year being simulated
        real(dp)                  :: tminm            !minimum monthly temperature for the year being simulated
        real(dp), dimension(2)    :: runtemp          !holds one month of temp data for mean monthly calc.
        real(dp)                  :: aprec            !annual precipitation (mm)
        real(dp)                  :: MAT              !mean annual temperature (C)
        real(dp)                  :: MAT_d            !stores running total of MAT (C)

end type atm_variables

type (atm_variables), save, target :: atm

end module metvarsmod
