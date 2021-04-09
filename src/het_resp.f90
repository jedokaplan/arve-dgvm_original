!SUBROUTINE het_resp
!Litter and soil decomposition
!Incorporates analytical solution for soil pool sizes once litter inputs
!are (assumed to be) at equilibrium, reducing spin-up time for carbon
!fluxes due to soil respiration.
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code and converted from monthly to daily call.

subroutine het_resp(j,year)

use arveparams,only : dp,npft,ncvar,daysec
use statevars, only : sv,ov,dtime
use pft_state_variables, only : veg
use soilstate_vars, only : surf

implicit none

!ARGUMENTS
integer, intent(in)   :: j        !gridcell]
integer, intent(in)   :: year

!pointers
real(dp), pointer, dimension(:,:) :: litter_ag_fast        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg                !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: cpool_fast         !fast-decomposing soil C pool (gC/m2)
real(dp), pointer, dimension(:)   :: cpool_slow                !slow-decomposing soil C pool (gC/m2)
real(dp), pointer, dimension(:)   :: drh                !daily heterotrophic respiration (gC/m2 timestep-1)
real(dp), pointer, dimension(:)   :: Tliq               !soil volumetric water content (fraction)
real(dp), pointer, dimension(:)   :: Tsat               !soil volumetric pore space (fraction)
real(dp), pointer, dimension(:)   :: temp_soil                !top-layer soil temperature
logical, pointer, dimension(:)    :: present            !true if PFT is present
real(dp), pointer, dimension(:) :: clay    !soil clay content (mass percent)
real(dp), pointer, dimension(:) :: bulk    !soil bulk density (kg m-3)
real(dp), pointer, dimension(:) :: dz      !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: litter_decom_ave
real(dp), pointer :: k_fast_ave
real(dp), pointer :: k_slow_ave

!LOCAL VARIABLES
integer  :: pft,c                      !counters
real(dp) :: k_litter_fast               !fast litter decomposition rate [timestep)
real(dp) :: k_litter_slow               !slow litter decomposition rate [timestep)
real(dp) :: k_soil_fast                      !fast pool decomposition rate [timestep)
real(dp) :: k_soil_slow                      !slow pool decomposition rate [timestep)
real(dp) :: temp_resp                        !temperature response of decomposition
real(dp) :: moist_resp                        !moisture response of decomposition
real(dp) :: cflux_litter_soil          !litter decomposition flux to soil
real(dp) :: cflux_litter_atmos         !litter decomposition flux to atmosphere
real(dp) :: cflux_fast_atmos               !soil fast pool decomposition flux to atmosphere
real(dp) :: cflux_slow_atmos               !soil slow pool decomposition flux to atmosphere
real(dp) :: litterdag_fast               !above-ground component of fast litter decomp
real(dp) :: litterdag_slow               !above-ground component of slow litter decomp
real(dp) :: litter_decom_bg               !below-ground component of litter decomposition
real(dp) :: temp                       !temporary variable
real(dp) :: Tmoist                       !temp variable for mean soil moisture over top 3 layers
real(dp) :: Soiltemp                       !temp variable for mean soil temperature over top 3 layers
real(dp), dimension(ncvar) :: litter_decom !litter decomposition
real(dp) :: fastfrac  !fraction of litter entering fast soil decomposition pool (opposed to slowC)
real(dp) :: slowfrac  !fraction of litter entering slow soil decomposition pool
real(dp) :: ek_lf
real(dp) :: ek_ls
real(dp) :: ek_sf
real(dp) :: ek_ss
real(dp) :: claytot   !total mass of clay in the soil column (to 3m, kg)

!PARAMETERS
real(dp), parameter :: s_per_yr        = 365._dp * 24._dp * 60._dp * 60._dp  !number of seconds in one year
real(dp), parameter :: i_s_per_yr      = 1._dp / s_per_yr                         !
real(dp), parameter :: k_litter_fast10 = 1./ (2.*s_per_yr)      !fast litter decomposition rate at 10 deg C (s-1)(2. /year)
real(dp), parameter :: k_litter_slow10 = 1./ (20.*s_per_yr)     !slow litter decomposition rate at 10 deg C (s-1)(20 /year)
real(dp), parameter :: k_soil_fast10   = 1./ (20.*s_per_yr)     !fast pool decomposition rate at 10 deg C (s-1)(20 /year)
real(dp), parameter :: k_soil_slow10   = 1./ (1000.*s_per_yr)   !slow pool decomposition rate at 10 deg C (s-1)(1000 /year)
real(dp), parameter :: atmfrac         = 0.65           !fraction of litter decomposition going directly into atmosphere
real(dp), parameter :: soilfrac = 1._dp - atmfrac   !fraction of litter decomposition going to soil C pools
real(dp), parameter :: sb  = 0.18  !parameters for the moisture effect on respiration equation below
real(dp), parameter :: sfc = 0.7   !units relative soil wetness S
real(dp), parameter :: t0  = 0.15  !maximum SOM partitioning coefficient between fast and slow soil pools
integer,  parameter :: soil_equil_year = 750      !number of years until pool sizes for soil decomposition solved analytically
real(dp), parameter :: seyr = real(soil_equil_year) * s_per_yr

!point pointers
temp_soil        => sv(j)%Tsoil
Tliq             => sv(j)%Tliq
Tsat             => sv(j)%Tsat
litter_ag_fast   => sv(j)%litter_ag_fast
litter_ag_slow   => sv(j)%litter_ag_slow
litter_bg        => sv(j)%litter_bg
cpool_fast       => sv(j)%cpool_fast
cpool_slow       => sv(j)%cpool_slow
drh              => veg%drh
present          => sv(j)%present
dz               => sv(j)%dz
clay             => sv(j)%clay
bulk             => sv(j)%bulk
litter_decom_ave => sv(j)%litter_decom_ave
k_fast_ave       => sv(j)%k_fast_ave
k_slow_ave       => sv(j)%k_slow_ave

!-----------------------------
!begin calculations
temp = 0._dp
litter_decom = 0._dp

!Temperature response function is a modified Q10 relationship
!Lloyd & Taylor 1994)

!this uses the temperature of the very top soil layer
!I changed it to be for the top 3 soil layers, similar depth to LPJ old soil model -JM Jan 26 08

Soiltemp = sum(temp_soil(1:3)) / 3

!if (Soiltemp < -46._dp) then !avoid division by zero
!to keep in line with LPJ2 changed to -25 degC (JM 16.07.2010)
if (Soiltemp < 248.15) then !avoid division by zero
  temp_resp = 0._dp
else
  temp_resp = exp(308.56_dp * ((1 / 56.02_dp) - (1 /(Soiltemp - 227.13_dp)))) !Lloyd & Taylor 1994, eqn 11
end if

!Moisture response based on soil layer 1 moisture content (Foley 1995)
!I changed it to be for the top 3 soil layers, again similar depth to old LPJ soil model -JM Jan 26 08

Tmoist = sum(Tliq(1:3) / Tsat(1:3)) / 3

!moisture response based on soil layer 1 moisture content (eqn. 7, Manzoni & Porporato, Soil Biol. Biochem. 39, 2007)
!this does not work in lpj because mw1 represents only the range between wilting point and field capacity, not total soil wetness. This works in ARVE.

!NOTE this is different than LPJ (JM 16.07.2010)
if (Tmoist <= sb) then
  moist_resp = 0.
else if (Tmoist <= sfc) then
  moist_resp = (Tmoist - sb) / (sfc - sb)
else
  moist_resp = sfc / Tmoist
end if

!original LPJ function after Foley et al. (1995)
!moist_resp = (0.25 + (0.75 * Tmoist))

!Calculate decomposition rates (k, timestep-1) as a function of temperature and moisture
!k = k_10 * temp_resp * moist_resp

k_litter_fast = k_litter_fast10 * temp_resp * moist_resp * dtime
k_litter_slow = k_litter_slow10 * temp_resp * moist_resp * dtime
k_soil_fast = k_soil_fast10 * temp_resp * moist_resp * dtime
k_soil_slow = k_soil_slow10 * dtime            !following LPJ2 remove temp response on this one.

ek_lf = 1. - exp(-k_litter_fast)
ek_ls = 1. - exp(-k_litter_slow)
ek_sf = 1. - exp(-k_soil_fast)
ek_ss = 1. - exp(-k_soil_slow)

!calculate partitioning fractions for fast and slow SOM based on soil clay content

!find the total mass of clay in the soil column (to 3m)

claytot  = sum(dz(1:10) * 0.01 * clay(1:10)* bulk(1:10))  !final units kg m-2

!the partitioning between fast and slow SOM pools is a function of total soil clay mass
!after Jobbagy & Jackson (Ecol. App. 10, 2000)

slowfrac = t0 * claytot / 300._dp !300 kg m-2 is ca. the maximum value in the ISRIC WISE data, extrapolating layer 5 to 3m
fastfrac = 1. - slowfrac

!Calculate daily litter decomposition using equation
! (1) dc/dt = -kc        where c=pool size, t=time, k=decomposition rate from (1),
! (2) c = c0*exp(-kt) where c0=initial pool size from (2), decomposition in any month given by
! (3) delta_c = c0 - c0*exp(-k) from (4)
! (4) delta_c = c0*(1.0-exp(-k))

!-------
!litter decomposition

do pft = 1,npft
  if (present(pft)) then

    !eqn 4;
    litterdag_fast  = litter_ag_fast(pft,1) * ek_lf
    litterdag_slow  = litter_ag_slow(pft,1) * ek_ls
    litter_decom_bg = litter_bg(pft,1) *      ek_lf                !belowground litter has the same turnover rate as fast litter
    litter_decom(1) = litter_decom(1) + litterdag_fast + litterdag_slow + litter_decom_bg

    do c = 2,ncvar
      litter_decom(c) = litter_decom(c) + litter_ag_fast(pft,c) * litterdag_fast  &
                      + litter_ag_slow(pft,c) * litterdag_slow + litter_bg(pft,c) &
                      * litter_decom_bg
    end do

    !Update the litter pools
    litter_ag_fast(pft,1) = max(litter_ag_fast(pft,1) - litterdag_fast,0._dp)
    litter_ag_slow(pft,1) = max(litter_ag_slow(pft,1) - litterdag_slow,0._dp)
    litter_bg(pft,1)      = max(litter_bg(pft,1) -      litter_decom_bg,0._dp)

  end if
end do  !PFT loop

!-------
!soil organic matter decomposition

if (litter_decom(1) > 0._dp) then
  litter_decom(2) = litter_decom(2) / litter_decom(1)
  litter_decom(3) = litter_decom(3) / litter_decom(1)
end if

!Calculate carbon flux to atmosphere and soil

cflux_litter_atmos = atmfrac  * litter_decom(1)
cflux_litter_soil  = soilfrac * litter_decom(1)

!Further subdivide soil fraction between fast and slow soil pools

temp = cpool_fast(1)
cpool_fast(1) = cpool_fast(1) + fastfrac * cflux_litter_soil

if (cpool_fast(1) > 0._dp) then
  do c = 2,ncvar
    cpool_fast(c) = (cpool_fast(c) * temp + fastfrac * litter_decom(c) * cflux_litter_soil) / cpool_fast(1)
  end do
end if

temp = cpool_slow(1)
cpool_slow(1) = cpool_slow(1) + slowfrac * cflux_litter_soil

if (cpool_slow(1) > 0._dp) then
  do c = 2,ncvar
    cpool_slow(c) = (cpool_slow(c) * temp + (1._dp - fastfrac) * cflux_litter_soil * litter_decom(c)) / cpool_slow(1)
  end do
end if

!Calculate timestep soil decomposition to the atmosphere
        cflux_fast_atmos = cpool_fast(1) * ek_sf  !eqn 4
        cflux_slow_atmos = cpool_slow(1) * ek_ss  !eqn 4

!Update the soil pools
        cpool_fast(1) = max(cpool_fast(1) - cflux_fast_atmos,0._dp)
        cpool_slow(1) = max(cpool_slow(1) - cflux_slow_atmos,0._dp)

!-------
!total microbial decomposition
        drh(1) = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos

!write(0,'(5es12.4)')cflux_litter_atmos,cflux_fast_atmos,cflux_slow_atmos,drh(1)
!read(*,*)

!13C and 14C
if (drh(1) > 0._dp) then
  do c = 2,ncvar
    drh(c) = (cflux_litter_atmos * litter_decom(c) + cflux_fast_atmos * cpool_fast(c) + cflux_slow_atmos * cpool_slow(c)) / drh(1)
  end do
end if

!write(0,'(a,7f10.5)')'drh',drh(1),cflux_litter_atmos,cflux_fast_atmos,cflux_slow_atmos,cpool_fast(1),cpool_slow(1),litter_decom(1)

!Empty soil pools below a minimum threshold

!if (cpool_fast(1) < 1.e-5) then
!  do c = 1,ncvar
!    cpool_fast(c) = 0._dp
!  end do
!end if

!if (cpool_slow(1) < 1.e-5) then
!  do c = 1,ncvar
!    cpool_slow(c) = 0._dp
!  end do
!end if

!SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
!Analytical solution of differential flux equations for fast and slow
!soil carbon pools.  Implemented after (soil_equil_year) simulation
!years, when annual litter inputs should be close to equilibrium.  Assumes
!average climate (temperature and soil moisture) from all years up to
!soil_equil_year.

if (year == soil_equil_year + 1 .and. k_fast_ave > 0. .and. k_slow_ave > 0.) then

  !Analytically calculate pool sizes

  !Rate of change of soil pool size = litter input - decomposition
  !  (5) dc/dt = litter_decom - kc

  !At equilibrium,
  !  (6) dc/dt = 0

  !From (5) & (6),
  !  (7) c = litter_decom / k

  cpool_fast(1) = (soilfrac * fastfrac * litter_decom_ave(1)) / k_fast_ave   !eqn 7

  cpool_slow(1) = (soilfrac * slowfrac * litter_decom_ave(1)) / k_slow_ave   !eqn 7

  !FLAG what about the isotopes?? JM 31.08.2010


else if (year <= soil_equil_year) then

  !Update running average respiration rates and litter input

  k_fast_ave = k_fast_ave + k_soil_fast * (seyr / dtime)
  k_slow_ave = k_slow_ave + k_soil_slow * (seyr / dtime)

  litter_decom_ave(1) = litter_decom_ave(1) + litter_decom(1) * (seyr / dtime)

end if

end subroutine het_resp
