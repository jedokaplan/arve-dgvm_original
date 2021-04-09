module soilstate_vars

use arveparams, only : dp,ns,nl,band,npft

implicit none

public :: init_surf_state_vars

! FLAG I want to put these into the surf type below! JM 29.03.2011

! Soil instantaneous state variables
! These are NOT retained year on year.

!------------------------------------------------
!common variables used in the soil physics and surface energy balance routines
  real(dp) :: lam                                !latent heat of phase change
  real(dp) :: fsnow_mo                          !monthly sum of gridcell covered by snow
  real(dp), dimension(2,band) :: asurf_g         !overal surface albedo (direct/diffuse:vis/nir) (fraction)
  real(dp), dimension(2) :: albedo_mo_tot        !overall surface albedo monthly cumulative diffuse and direct (1=vis, 2=nir)
  real(dp), dimension(npft) :: Eveg                   !surface emissivity (fraction) vegetation
  real(dp) :: Egrnd                                   !ground emissivity(fraction)

  !instantaneous state of the atmosphere calcs
  real(dp) :: TairK                                  !surface (2m) air temperature (K)
  real(dp) :: tair                                    !temperature in C degrees
  real(dp) :: eatm                                     !atmospheric vapor pressure (Pa)
  real(dp) :: emb                                      !atmospheric vapor pressure, millibar
  real(dp) :: ratm                                   !atmospheric density (kg m-3)
  real(dp) :: theta                                  !atmospheric potential temperature

  !surface heat flux variables

  !Shortwave
  real(dp) :: Swg_b                             !shortwave rad absorbed by ground bare ground (W m-2)
  real(dp), dimension(npft) :: Swg_v            !shortwave rad absorbed by ground beneath vegetation (W m-2)
  real(dp), dimension(npft) :: Swv              !shortwave rad absorbed by vegetation (W m-2)

  !Longwave
  real(dp) :: Latm                              !downwelling atmospheric longwave radiation (W m-2)
  real(dp), dimension(npft) :: Lwg_v            !canopy net longwave rad flux for the ground (W m-2)
  real(dp) :: Lwg_b
  real(dp) :: dLgdT_b                             !derivative of the longwave radiation flux into the surface (W m-2) w.r.t. Temp
  real(dp), dimension(npft) :: dLgdt_v                !derivative of the longwave radiation flux into the surface (W m-2) w.r.t. veg temp
  real(dp), dimension(npft) :: Lwvd                !downward longwave below vegetation (W m-2)
  real(dp), dimension(npft) :: Lwv              !net longwave rad flux for vegetation (W m-2)
  real(dp), dimension(npft) :: dLvdT            !vegetation derivative of the longwave radiation flux into surface (W m-2) w,r,t Temp

  !Sensible
  real(dp) :: Hg                                     !total sensible heat flux into the surface (W m-2)
  real(dp) :: dHgdT                                  !total derivative of the sensible heat flux into the surface (W m-2) w.r.t. temp
  real(dp) :: Hg_b                                   !sensible heat flux into the surface (W m-2) bare surfaces
  real(dp) :: dHgdt_b                                  !derivative of the sensible heat flux into the surface (W m-2) w.r.t temp (bare surf)
  real(dp), dimension(npft) :: Hg_v                   !sensible heat flux into the surface (W m-2) vegetated surfaces
  real(dp), dimension(npft) :: dHgdt_v          !derivative of the sensible heat flux into surface (W m-2) w.r.t temp (vegetated surf)

  !Latent
  real(dp) :: Eg                                     !latent heat flux into the surface (W m-2)
  real(dp) :: dEgdT                                  !derivative of the latent heat flux into the surface (W m-2) with respect to Temperature
  real(dp) :: Eg_b                                   !latent heat flux into the surface (W m-2)(bare surfaces)
  real(dp) :: dEgdt_b                                 !derivative of the latent heat flux into surface (W m-2) w.r.t temp (bare surf)
  real(dp), dimension(npft) :: Eg_v                    !latent heat flux into the surface (W m-2)(vegetated surfaces)
  real(dp), dimension(npft) :: dEgdt_v                 !derivative of the latent heat flux into the surface (W m-2) w.r.t temp (veg surf)
  real(dp) :: Beta_soil                                !CLM 4.0 Eqn 5.68 represent molecular diffusion processes in evaporation from soil
  real(dp) :: dqsath_dT                                !saturated water vapour specific humidity (CLM 4 Eqn 5.143) derivative w.r.t grnd temp.
  real(dp) :: dqgdTg                                !derivative of the specific humidity of the soil surface (CLM 4.0 Eqn 5.74)
  real(dp) :: qg                                !ground specific humidity
  real(dp) :: qs                                !canopy specific humidity

  !Total
  real(dp) :: hs                                     !net energy flux into the surface (W m-2)
  real(dp) :: dhsdT                                  !derivative of hs with respect to temperature

  ! soil thermal conductivity and heat capacity
  real(dp), dimension(ns:nl) :: Kl      !soil thermal conductivity across layer boundary (W m-1 K-1)
  real(dp), dimension(ns:nl) :: Kh      !soil thermal conductivity at layer midpoint (W m-1 K-1)
  real(dp), dimension(ns:nl) :: Cl      !instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)

  real(dp), dimension(ns:nl) :: Ku      !soil water instantaneous (unsaturated) conductivity across layer boundary (m s-1)

  !used in soilphasechg routine
  real(dp), dimension(ns:nl) :: fact    !factor used in computing tridiagonal coefficients
  real(dp), dimension(ns:nl) :: Fhti    !heat flux across soil layer boundary (W m-2)
  integer,  dimension(ns:nl) :: ithaw
  real(dp) :: qsnomelt                  !snow melt (kg m-2 sec -1)
  
  
!-------------  

  type surface_variables
  
    real(dp), dimension(2)    :: ra_dp            !instantaneous downwelling shortwave radiation (W m-2)(1 is direct, 2 is diffuse)
    real(dp)                  :: dsw_t            !integrated calculated downwelling shortwave radiation top of atmos(KJ m-2 day)
    real(dp)                  :: absrad_annual    !annual total of absorbed radiation (J/m2)
    real(dp)                  :: sdecl            !solar declination (deg)
    real(dp)                  :: rZ0              !zenith angle at solar noon (radians)
    real(dp)                  :: zen              !cosine of the zenith angle (rads)
    real(dp)                  :: avg_cZ           !diurnally averaged cos of the zenith angle (rad) of the incident beam

    !resistance terms
    real(dp)                  :: raw_b            !aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
    real(dp)                  :: rah_b            !aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)
    real(dp), dimension(npft) :: rawp             !aerodynamic resistance for latent heat transfer btwn canopy air & ground (s m-1)
    real(dp), dimension(npft) :: raw              !aerodynamic resistance for latent heat transfer btwn canopy air & atmosphere (s m-1)
    real(dp), dimension(npft) :: rah              !aerodynamic resistance for sensible heat transfer btwn canopy air & atmosphere (s m-1)
    real(dp), dimension(npft) :: rahp             !aerodynamic resistance for sensible heat transfer btwn canopy air & ground (s m-1)
    real(dp), dimension(npft) :: rb               !leaf boundary layer resistance (s m-1)
    real(dp), dimension(npft) :: rdbp             !fraction of potential evapotranspiration
    real(dp), dimension(npft) :: rdp_dry          !term for contributions by wet leaves and transpiration, limited by avail H2O and pot. evapo.
    real(dp)                  :: ustar_b          !friction velocity estimated from scalar wind speed (m s-1)
  
    real(dp), dimension(365)   :: ald             !active layer depth (m)
    real(dp), dimension(365)   :: Aliq            !daily soil water content for layer 1
    
    real(dp) :: zsno_mo                           !monthly cumulative snow depth (m) 
    real(dp) :: realsoilT                         !used when weather station data is used. Can read in the real soil T measured at the site.

    real(dp) :: fsnow                             !fraction of the gridcell covered by snow (fraction)
  
    real(dp) :: mo_runoff_tot                     !monthly total overland runoff (mm)
    real(dp) :: qseva_tot                         !timestep soil evaporation (mm)
    real(dp) :: qsubl_tot                         !timestep soil sublimation (mm)
    real(dp) :: qsdew_tot                         !timestep soil dew  (mm)
    real(dp) :: qfrost_tot                        !timestep soil frost (mm)

    real(dp) :: qseva                             !soil evaporation flux (mm s-1)     
    real(dp) :: qsubl                             !soil sublimation flux (mm s-1)     
    real(dp) :: qsdew                             !soil surface dew flux (mm s-1)          
    real(dp) :: qfrost                            !soil surface frost flux (mm s-1)        

    real(dp) :: qliq                              !liquid precipitation (kg m-2 sec-1)                                      
    real(dp) :: qover                             !liquid water surface runoff (kg m-2 sec-1)                               
    real(dp) :: qsno                              !snowfall (kg m-2 sec -1)                                                 
    real(dp), dimension(npft) :: qgrnd_l          !total rate of liquid precip reaching the ground under canopy(mm s-1)     
    real(dp), dimension(npft) :: qgrnd_s          !total rate of solid precip reaching the ground under canopy(mm s-1)      

  end type surface_variables

type (surface_variables), save, target :: surf

!------------------------------------------------------

contains

subroutine init_surf_state_vars

surf%Aliq(:) = 0._dp
surf%rdbp = 1._dp 
surf%mo_runoff_tot = 0._dp

end subroutine init_surf_state_vars

end module soilstate_vars
