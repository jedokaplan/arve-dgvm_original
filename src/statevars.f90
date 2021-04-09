module statevars

!this model contains the sv and ov types as well as some select other variables

! All of the state variables that are needed year upon year are contained in this one structure.
! This structure is then declared in the main program.
! In the subroutines, pointers point to variables contained in this
! structure.

use arveparams, only : dp,sp,npft,ncvar,climbuf,ns,nl,band,numtimestep,nts

implicit none

public :: initstatevars
public :: area

   !--------------------------------------------------------

type state_variables

  ! Physical variables
  logical                       :: valid_cell                   !cell is valid to calculate                                                   
  real(dp)                      :: lat                          !latitude of gridcell (degrees)                                               
  real(dp)                      :: lon                          !longitude of gridcell (degrees)                                              
  real(dp)                      :: grid_area                    !gridcell area in m2                                                          
  logical                       :: polar                        !Cell is located in polar region. polar night if true.                        
  real(dp)                      :: Patm                         !atmospheric pressure from elevation (Pa)                                     
  real(dp)                      :: Patm30                       !atmospheric pressure at reference height (30m) from elevation (Pa)           

  real(dp)                      :: dsw_tn                       !the next day's dsw_t                                                                     
  real(dp)                      :: maxdl                        !gridcell's longest daylength (s)                                                         
  real(dp)                      :: sdecl_n                      !the next day's solar decl (deg)
  real(dp), dimension(climbuf)  :: mtemp_min_buf                !buffer to store 'climbuf' years of coldest month temperatures
  real(dp), dimension(climbuf)  :: mtemp_max_buf                !buffer to store 'climbuf' years of warmest month temperatures

  ! Plant tissue compartments (per individual)
  real(dp), dimension(npft,ncvar) :: lm_ind                     !individual leaf mass (gC)
  real(dp), dimension(npft,ncvar) :: lm_max                     !maximal leaf mass for the plant based on stem and roots (gC)
  real(dp), dimension(npft,ncvar) :: rm_ind                     !individual fine root mass (gC)
  real(dp), dimension(npft,ncvar) :: rm_max                     !maximal individual fine root mass (gC)
  real(dp), dimension(npft,ncvar) :: sm_ind                     !individual sapwood mass (gC)
  real(dp), dimension(npft,ncvar) :: hm_ind                     !individual heartwood mass (gC)
  real(dp), dimension(npft,ncvar) :: NSC_storage                !individual plant non-structural carbon stores (gC)

  ! Litter pools and decomposition
  real(dp), dimension(npft,ncvar)     :: litter_ag_fast         !gridcell above-ground litter (gC/m2)
  real(dp), dimension(npft,ncvar)     :: litter_ag_slow         !gridcell above-ground litter (gC/m2)
  real(dp), dimension(npft,ncvar)     :: litter_bg              !gridcell below-ground litter (gC/m2)
  real(dp), dimension(ncvar)          :: cpool_fast             !fast-decomposing soil C pool (gC/m2)     
  real(dp), dimension(ncvar)          :: cpool_slow             !slow-decomposing soil C pool (gC/m2)     
  real(dp)                            :: k_fast_ave             !running average fast litter respiration rates
  real(dp)                            :: k_slow_ave             !running average slow litter respiration rates
  real(dp), dimension(ncvar)          :: litter_decom_ave       !running average litter input

  ! Vegetation physical structure
  logical, dimension(npft)  :: present                          !whether PFT present in gridcell
  logical, dimension(npft)  :: virtual_leaf                     !true if the leaf being calculated is a virtual one to determine phenology.
  real(dp), dimension(npft) :: crownarea                        !crown area (m2)   
  real(dp), dimension(npft) :: fpc_gmax                         !maximal gridcell foliar projective cover (FPC) per pft for gridcell (m2 m-2)
  real(dp), dimension(npft) :: stema                            !stem area index (m2 m-2)
  real(dp), dimension(npft) :: stem_sno                         !stema with snow burial
  real(dp), dimension(npft) :: height                           !tree height (m)
  real(dp), dimension(npft) :: lai_max                          !maximal individual leaf area index
  real(dp), dimension(npft) :: lai_ind                          !instantaneous individual leaf area index based on present lm_ind and phenological state(m2 m-2)
  real(dp), dimension(npft) :: lai_sno                          !lai_ind with snow burial
  real(dp), dimension(npft) :: lai_vir                          !virtual LAI for phenology decisions with snow burial
  real(dp), dimension(npft) :: lai_vir_no_burial                !virtual LAI for phenology decisions with no snow burial
  real(dp), dimension(npft) :: nind                             !gridcell individual density (indiv/m2)
  real(dp), dimension(npft) :: needed_sap                       !the minimum amount of sapwood for the plant (at start of dormant season) gC indiv

  real(dp), dimension(nl,npft) :: rootfracl                     !root fraction in that soil layer
  real(dp), dimension(npft)    :: rootdepth                     !root depth (m)
  real(dp), dimension(npft)    :: betaT                         !soil moisture function limiting transpiration (fraction)

  real(dp), dimension(npft)     :: Wcan_old                     !amount of water in the canopy in the previous timestep (mm)                                       
  integer, dimension(npft)      :: yrpres                       !number of years that the pft was present                                                          
  real(dp), dimension(2,npft)   :: temp_veg                     !temperature of the vegetation (K) per pft (for the previous (1) and present (2) timesteps)        
  real(dp), dimension(2,npft)   :: Tveg_ave                     !timestep average vegetation temperature (K) (1) previous timestep, (2) present timestep           

  ! Phenology
  real(dp), dimension(npft)     :: aleafdays                    !number of days with foliage
  integer,  dimension(npft)     :: pstate                       !instantaneous phenological state for deciduous plants 
                                                                ! (0=no leaves, 1=full leaf out, 2=leafing out, 3=senescing)
  integer, dimension(npft)      :: leafonday                    !day of year of leaf onset
  integer, dimension(npft)      :: leafoffday                   !day of year of completed senescence
  real(dp), dimension(nts)      :: meanT_ra                     !daily mean temperature (C) 
  real(dp), dimension(npft,nts) :: soilT_ra                     !daily soil temperature in root zone (C) per pft 
  real(dp), dimension(nts)      :: fsnow_ra                     !daily fractional snow cover array 

  ! Phenology (annual C allocation version)
  real(dp), dimension(npft) :: phdays                           !number of days the pft has been in the given phenological state
  real(dp), dimension(nts)  :: dsw_ra                           !dsw running average over nts timesteps
  real(dp), dimension(npft,nts) :: smp_ra                       !soil matric potential in the root zone running average over nts timesteps
  real(dp), dimension(nts) :: minT_ra                           !daily min temperature (C) running average over nts timesteps
  real(dp) :: chilld                                            !winter chilling days (# of days with daily temp <= 5C) reset Aug 1 (N. Hemi)
  real(dp), dimension(npft)     :: phfrac                       !instantaneous phenological state fraction (0=no leaves, 1=full leaf out)
  integer, dimension(npft)      :: leafout                      !status of leaf out (1) able to leaf out or presently leafed out, (2) leaf shed due to cold/daylight
                                                                ! (3) leaf shed due to moisture stress, (4) leaf shed due to leaf longevity
  real(dp), dimension(npft)     :: gddsum                       !growing degree-day sum (reset at senesence)
  real(dp), dimension(npft)     :: gddleaf                      !growing degree-day value at leaf out initiation
  real(dp), dimension(npft,nts) :: npp_ra                         !FLAG   
  real(dp), dimension(npft) :: NPP7     !FLAG

  ! Spatial references
  integer                      :: snl                           !index value of top snow layer (negative is more layers)
  integer                      :: gnl                           !index value of lowest 'soil' layer
  real(dp)                     :: soildepth                     !depth to bedrock (m)
  real(dp), dimension(ns:nl)   :: dz                            !thickness of the soil or snow layers (m)
  real(dp), dimension(ns:nl)   :: zpos                          !z coordinate of the middle of the soil layer relative to soil surface (m), pos down
  real(dp), dimension(ns-1:nl) :: zipos                         !z coordinate of the interface between soil layers relative to soil surface (m), pos down
  real(dp), dimension(1:nl)    :: dzmm          !in mm
  real(dp), dimension(1:nl+1)  :: zposmm        !in mm

  ! Soil physical components
  real(dp), dimension(nl)    :: sand                            !soil sand content (mass percent)
  real(dp), dimension(nl)    :: silt                            !soil silt content (mass percent)
  real(dp), dimension(nl)    :: clay                            !soil clay content (mass percent)
  real(dp), dimension(nl)    :: sorg                            !soil organic matter content (mass percent)
  real(dp), dimension(nl)    :: rock                            !soil rock content (mass FRACTION)
  real(dp), dimension(ns:nl) :: bulk                            !soil bulk density (kg m-3)
  real(dp), dimension(nl)    :: Vcf                             !FRACTION of coarse fragments (rocks) by volume
  real(dp), dimension(nl)    :: Vom                             !volumetric fraction
  real(dp), dimension(nl)    :: Vsand                           !volumetric fraction
  real(dp), dimension(nl)    :: Vclay                           !volumetric fraction
  real(dp), dimension(nl)    :: Vsilt                           !volumetric fraction
  logical                    :: peat                            !peatland flag

  ! Soil temperature variables
  real(dp), dimension(1:nl)  :: Tsoil_ave                       !avg. soil temp across night/day (K)
  real(dp), dimension(ns:nl) :: Tsoil                           !soil temperature (K)
  real(dp), dimension(ns:nl) :: Tsoiln                          !soil temperature (K) at previous timestep
  real(dp), dimension(ns:nl) :: Tnsoi                           !updated soil temperature at timestep n+1 (used in soiltemperature)

  ! Soil hydrology
  ! the dimensions of these arrays in practice depend on the soil depth and are allocated in subroutine 'soillayers'.
  real(dp), dimension(ns:nl) :: Tliq                             !soil liquid water content (fraction)
  real(dp), dimension(ns:nl) :: Tice                             !soil ice content (fraction)
  real(dp), dimension(ns:nl) :: Tpor                             !soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
  real(dp)                   :: theta_fc                         !the volumetric water content at field capacity (fraction)
  real(dp), dimension(ns:nl) :: Wliq                             !soil liquid water content at layer midpoint (mm)
  real(dp), dimension(ns:nl) :: Wice                             !soil ice content at layer midpoint (mm)
  real(dp), dimension(ns:nl) :: Psi                              !soil water potential at layer midpoint (mm)
  real(dp), dimension(1:nl+1)  :: Psi_eq                         !equilibrium soil matric potential (mm)
  real(dp), dimension(1:nl)  :: rPsi                             !Psi calc used in root extraction calculations
  real(dp), dimension(ns:nl) :: Ku                               !soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
  real(dp), dimension(ns:nl) :: fice0                            !layer ice fraction, previous timestep
  real(dp), dimension(ns:nl) :: Ksat                             !soil water saturated conductivity at layer midpoint (mm s-1)
  real(dp), dimension(ns:nl) :: Tsat                             !soil water volumetric water content at saturation (fraction)
  real(dp), dimension(ns:nl) :: Bexp                             !soil water B exponent used in the Brooks & Corey Pedotransfer Functions
  real(dp), dimension(ns:nl) :: Psat                             !soil water matric potential at saturation (mm)
  real(dp), dimension(nl)    :: Kdry                             !soil thermal conductivity of dry natural soil (W m-1 K-1)
  real(dp), dimension(nl)    :: Ksolid                           !soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
  real(dp), dimension(nl)    :: Csolid                           !soil solids volumetric heat capacity at layer midpoint (J m-3 K-1)
  real(dp), dimension(npft)  :: wetsum                           !Plant wilting factor (CLM 4.0 eq 8.18) weighted by root zone
  real(dp)                   :: lastpet                          !potential evapotranspiration from last day (mm/timestep)
  
  ! Snow variables
  real(dp)                :: Wsno                                !snow water equivalent of the snowpack (mm)                               
  real(dp)                :: Wsno_old                            !snow water equivalent of the snowpack (mm) from previous timestep        
  real(dp)                :: zsno                                !snow total depth (m)                                                     
  real(dp)                :: tausno                              !snow age (non-dimensional)                                               
  real(dp)                :: tausno_old                          !snow age of previous timestep                                            
  real(dp), dimension(nl) :: Ffrz                                !fractional impermeable area as a function of soil ice content at a layer 
  
  ! Aquifer
  integer  :: jint                                               !index of the soil layer directly above the water table                   
  real(dp) :: Waquif_a                                           !water in unconfined aquifer below soil (mm)                              
  real(dp) :: Waquif_t                                           !water in aquifer within soil column (mm)                                 
  real(dp) :: zw                                                 !groundwater depth (m)                                                    
  
  !soil moisture averages used in data assimilation mode (none-base ARVE code)   
  !real(dp), dimension(12)  :: meansoilwater  !long term running mean soil moisture
  !real(dp), dimension(365) :: dmsoilm        !smooth daily interpolation of mean soil moisture above variable
  
end type state_variables

type(state_variables), save, allocatable, target, dimension(:) :: sv

!------------------------------------------

type output_variables

  real(sp), dimension(12)    :: mw1                     !monthly root zone soil water content fraction
  real(sp), dimension(12)    :: mts                     !monthly root zone soil temperature
  real(sp), dimension(12)    :: runoff                  !monthly overland runoff (mm /day)
  real(sp), dimension(12,3)  :: mrhc                    !het resp carbon
  real(sp), dimension(12,2)  :: albedo                  !monthly mean surface albedo (fraction)(1 = vis, 2 = nir)
  real(sp), dimension(12)    :: snow_frac               !snow fractional ground cover mean monthly fraction
  real(sp), dimension(12)    :: snow_depth              !mean monthly snow depth (m)
  real(sp)                   :: annual_hetresp          ! BENMODS
  real(sp)                   :: soilc                   !soil carbon
  real(sp)                   :: soili                   !soil isotope ratio
  real(sp), dimension(npft)  :: leafondays              !days with leaves on the plants (days/yr)
  real(sp)                   :: plant_carbon            !total plant structural carbon (gC/m2)
  real(sp), dimension(ncvar) :: grid_npp                !gridcell total npp (gC/m2)
  real(sp), dimension(ncvar) :: grid_gpp                !gridcell total gpp (gC/m2)
  real(sp), dimension(ncvar) :: grid_autoresp           !gridcell total autotrophic respiration (gC/m2)
  real(sp)                   :: abs_rad_annual          !annual absorbed total radiation (GJ m-2)
  real(sp)                   :: grndwat_depth           !depth to groundwater (m)
  real(sp), dimension(npft)  :: aaet_out                !annual actual evapotranspiration (mm/year)
  real(sp)                   :: max_ald                 !deepest active layer depth in the year (m)

  real(sp), dimension(ncvar) :: acflux_fire             !C flux to atmosphere due to fire (gC/m2)
  real(sp), dimension(ncvar) :: CH4flux_fire            !CH4 flux to atmosphere due to fire (gCH4/m2)
  real(sp)                   :: afire_frac              !fraction of gridcell burnt this year
  real(sp)                   :: COfire_flux             !CO produced by fire (g/m2)
  real(sp)                   :: VOCfire_flux            !VOC produced by fire (g/m2)
  real(sp)                   :: TPMfire_flux            !TPM produced by fire (g/m2)
  real(sp)                   :: NOxfire_flux            !NOx produced by fire (g/m2)

end type output_variables

type(output_variables), save, allocatable, target, dimension(:) :: ov

!--------------------------------------------

!Variables that are the same regardless of gridcell 
!or denote time and are calculated daily

real(dp), dimension(ncvar)      :: co2                  !Carbon dioxide (concentration in ppm)
logical                         :: gridded              !flag for type of simulation
integer                         :: counter              !
integer                         :: counter_lim          !
real(dp)                        :: dt                   !short time step (sec)
integer                         :: doy                  !day of the year
real(dp)                        :: dtime                !current timestep (sec)
real(dp), dimension(2)          :: dayl                 !daylength (sec) (1 is present timestep, 2 is next)

!--------------------------------------------

contains

subroutine initstatevars()

implicit none

integer :: i,k,n
integer :: j
integer :: l
integer :: m

sv%valid_cell = .true.
sv%peat = .false.
sv%tausno_old   = 0._dp
sv%tausno       = 0._dp
sv%Wsno          = 0._dp
sv%Wsno_old     = 0._dp
sv%k_fast_ave = 0._dp        
sv%k_slow_ave = 0._dp    
sv%lastpet = 3.e-3   !FLAG ok value to start? JM Apr 8 09   

forall (l=1:nl)
  sv%Tsoil_ave(l) = 0._dp
end forall

forall (m=1:climbuf)
  sv%mtemp_min_buf(m) = 0._dp
  sv%mtemp_max_buf(m) = 0._dp
end forall

forall (i=1:npft)
  sv%virtual_leaf(i) = .false.
  sv%rootdepth(i) = 0._dp
  sv%Wcan_old(i) = 0._dp
  sv%crownarea(i) = 0._dp
  sv%height(i)    = 1._dp
  sv%lai_ind(i)   = 0._dp
  sv%nind(i)      = 0._dp
  sv%gddsum(i)    = 0._dp
  sv%aleafdays(i) = 0._dp
  sv%yrpres(i)    = 0
  sv%leafout(i) = 1
  sv%present(i) = .false.
  sv%stem_sno(i) = 0._dp
  sv%lai_sno(i) = 0._dp
  sv%lai_vir(i) = 0._dp
  sv%lai_vir_no_burial(i) = 0._dp
  sv%stema(i) = 0._dp          
  sv%stem_sno(i) = 0._dp      
  sv%NPP7(i) = 0._dp  
  sv%leafonday(i) = 0
  sv%leafoffday(i) = 0
        
  forall (l=1:nl)
    sv%rootfracl(l,i) = 0._dp
  end forall

  forall (j=1:ncvar)
    sv%hm_ind(i,j) = 0._dp
    sv%lm_ind(i,j) = 0._dp
    sv%sm_ind(i,j) = 0._dp
    sv%rm_ind(i,j) = 0._dp
    sv%litter_ag_fast(i,j) = 0._dp 
    sv%litter_ag_slow(i,j) = 0._dp
    sv%litter_bg(i,j) = 0._dp    
    sv%lm_max(i,j) = 0._dp 
    sv%NSC_storage(i,j) = 0._dp
  end forall

  forall (k=1:nts)
   sv%soilT_ra(i,k) = 0._dp
   sv%smp_ra(i,k) = -5.e4
   sv%npp_ra(i,k) = 0._dp
  end forall

end forall

forall (j=1:ncvar)
 sv%cpool_slow(j) = 0._dp
 sv%cpool_fast(j) = 0._dp
 sv%litter_decom_ave(j) = 0._dp  
end forall

forall (n=1:nts)
   sv%dsw_ra(n) = 0._dp
   sv%meanT_ra(n) = 0._dp
   sv%minT_ra(n) = 0._dp
   sv%fsnow_ra(n) = 0._dp
end forall

!OV variables

forall (j=1:ncvar)
  ov%acflux_fire(j) = 0.
  ov%CH4flux_fire(j) = 0.
end forall

ov%afire_frac = 0.
ov%COfire_flux = 0.
ov%VOCfire_flux = 0.
ov%TPMfire_flux = 0.
ov%NOxfire_flux = 0.

end subroutine initstatevars

!-----------------------------------------

!this function returns the size of a regular grid cell in square meters.
real function area(lat,dlon,dlat)

use arveparams, only : pi,dp,d2r,radius

implicit none

real(dp), intent(in) :: lat
integer, intent(in) :: dlon
integer, intent(in) :: dlat

real(dp) :: cellarea
!real(dp) :: resolution
real(dp) :: deltalat
real(dp) :: deltalon
real(dp) :: elevation
real(dp) :: dlon_real
real(dp) :: dlat_real
!-------
dlon_real = 360. / real(dlon)
dlat_real = 180. / real(dlat)

elevation = d2r * (lat + (dlat_real / 2.0_dp))

deltalon = d2r * dlon_real
deltalat = d2r * dlat_real

cellarea = (2.0_dp * radius**2 * deltalon * cos(elevation) * sin((deltalat / 2.0_dp)))

area = cellarea * 1.0e6

end function area

end module statevars
