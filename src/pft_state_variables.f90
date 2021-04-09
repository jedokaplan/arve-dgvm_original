module pft_state_variables

! this module contains the type veg

use arveparams, only : dp,npft,ncvar,band

implicit none

public :: initpftvars

!--------------------------

! This is the basic data structure that contains the state variables
! for the Plant Functional Type (PFT).These variables are not retained beyond present year!

type veg_struct

  real(dp), dimension(npft) :: APAR_sun               !visible only shortwave rad absorbed by vegetation averaged for leaf area (W m-2) sunlight
  real(dp), dimension(npft) :: APAR_sha               !visible only shortwave rad absorbed by vegetation averaged for leaf area (W m-2) shaded
  real(dp), dimension(npft) :: ratio                  !the lambda ratio of stomatal opening (fraction)
  real(dp), dimension(npft) :: CO2comp                !CO2 compensation point in the absence of mitochondrial respiration
  real(dp), dimension(npft) :: gt_sun                 !instantaneous canopy conductance - sunlight leaves
  real(dp), dimension(npft) :: gt_sha                 !instantaneous canopy conductance - shaded leaves
  real(dp), dimension(npft)  :: fsun                  !sunlit fraction of canopy
  real(dp), dimension(npft)  :: fsha                  !shaded fraction of canopy  
  real(dp), dimension(npft)  :: f_wet                 !fraction of the canopy that is wet
  real(dp), dimension(npft)  :: f_dry                 !fraction of the canopy that is dry
 
  real(dp), dimension(npft)  :: gdd_bioclim           !growing degree days in one year (bioclimatic)

  real(dp), dimension(npft)  :: Wcan                  !amount of water in the canopy (mm)
  real(dp), dimension(npft)  :: Wcanmax               !maximum quantity of water the canopy can hold (mm)
  real(dp), dimension(npft)  :: can_evap              !canopy evaporation (mm s-1)
  real(dp), dimension(npft)  :: xs_can_drip           !excess canopy water (mm)

  logical, dimension(npft) :: pft_estab               !whether PFT within bioclimatic limits for establishment
  logical, dimension(npft) :: pft_survive             !whether PFT within bioclimatic limits for survival

  real(dp), dimension(npft)         :: pet            !potential evapotranspiration for canopy surface (mm s-1) per pft
  real(dp), dimension(npft)         :: red_pet        !potential evapotranspiration for canopy surface reduced by canopy evaporation (mm s-1) per pft
  real(dp), dimension(npft)         :: underpet       !potential evapotranspiration for ground beneath canopy (mm s-1) per pft
  real(dp), dimension(npft)         :: dpet           !potential evapotranspiration (mm/ day) per pft
  real(dp)                          :: gpet           !PET for bare ground

  real(dp), dimension(npft)         :: supply         !water supply to plants (mm s-1)
  real(dp), dimension(npft)         :: demand         !plant water demand (mm s-1)
  real(dp), dimension(npft)         :: supply_ts      !water supply to plants  timestep average
  real(dp), dimension(npft)         :: demand_ts      !plant water demand  timestep average
  real(dp), dimension(npft)         :: awscal         !annual average water scalar
    
  !daily
  real(dp), dimension(npft)  :: aaet                  !annual actual evapotranspiration (mm/year)
  real(dp), dimension(npft)  :: daet                  !calculated daily actual evapotranspiration (mm s-1)
  real(dp), dimension(npft)  :: dreprod
  real(dp), dimension(npft,ncvar) :: dresp            !daily maintenance respiration (g C m-2 timestep-1)
  real(dp), dimension(npft,ncvar) :: gresp            !daily growth resp. (g C m-2 timestep-1)
  real(dp), dimension(npft,ncvar) :: dgpp             !daily gross primary productivity (g C m-2 timestep-1)
  real(dp), dimension(npft,ncvar) :: dnpp             !daily net primary productivity (gC m-2 timestep-1)
  real(dp), dimension(npft)  :: dturnover_ind         !daily total turnover of living biomass per individual (gC)
  !short timestep
  real(dp), dimension(npft)  :: aet                   !short timestep calculated actual evapotranspiration (mm s-1)
 real(dp), dimension(npft,ncvar) :: sgpp              !short timestep gross primary productivity (g C m-2 timestep-1)
  real(dp), dimension(ncvar) :: drh                   !daily heterotrophic respiration (gC/m2 timestep-1)
  real(dp), dimension(ncvar) :: arh                   !annual heterotrophic respiration(gC/m2 yr-1)

  !albedo on plants
  real(dp), dimension(band,npft) :: ftid              !downward diffuse fluxes per unit incident direct beam(surface albedos)
  real(dp), dimension(band,npft) :: ftii              !downward diffuse fluxes per unit incident diffuse beam (surface albedos)
  real(dp), dimension(band,npft) :: scatcf_tot        !fraction of intercepted radiation that is scattered (0 to 1)
  real(dp), dimension(band,npft) :: scatcf_tot_vir    !fraction of intercepted radiation that is scattered (0 to 1) (virtual leaves)
  real(dp), dimension(band,npft) :: ftdd              !down direct flux below veg per unit dir flx
  real(dp), dimension(band,npft) :: fabd              !flux absorbed by vegetation
  real(dp), dimension(band,npft) :: fabd_vir          !flux absorbed by vegetation (virtual leaves)
  real(dp), dimension(band,npft) :: fabi              !flux absorbed by veg per unit diffuse flux
  real(dp), dimension(band,npft) :: fabi_vir          !flux absorbed by veg per unit diffuse flux (virtual leaves)
  real(dp), dimension(npft)      :: bigK              !the optical depth of direct beam per unit leaf and stem area

  real(dp), dimension(npft,ncvar) :: lm_sapl          !initial (sapling) leaf mass (gC)
  real(dp), dimension(npft,ncvar) :: rm_sapl          !initial (sapling) fine root mass (gC)
  real(dp), dimension(npft,ncvar) :: sm_sapl          !initial (sapling) sapwood mass (gC)
  real(dp), dimension(npft,ncvar) :: hm_sapl          !initial (sapling) heartwood mass (gC)
  real(dp), dimension(npft,ncvar) :: NSC_sapl         !initial (sapling) plant non-strucutral carbon store (gC)
  real(dp), dimension(npft) :: sla                    !invariant PFT specific leaf area (m2/gC)
  real(dp), dimension(npft) :: sla_sun                !specific leaf area sunlight (m2 gC-1) (variant)
  real(dp), dimension(npft) :: sla_sha                !specific leaf area shaded (m2 gC-1) (variant)
  
  real(dp), dimension(npft)       :: lm_tot_max       !running max leaf mass for the present year (gC)
  real(dp), dimension(npft)       :: rm_tot_max       !running max individual fine root mass for the present year (gC) 
  
  real(dp), dimension(npft)       :: frac             !fraction of sapwood and heartwood that is below ground (in coarse roots)(gC)
  real(dp), dimension(npft)       :: NSC_limit        !allometric plant non-structural carbon stores size (gC/indiv)      
                                                
  real(dp), dimension(npft,ncvar) :: acflux_estab     !annual biomass increment due to establishment (gC/m2)
  real(dp), dimension(npft)       :: turnover_ind     !total turnover of living biomass per individual (gC / yr)                                             
  real(dp), dimension(npft)       :: heatstress       !reduction in individual density (& establishment) due to heat induced mortality  (indiv/m2)  
  real(dp), dimension(npft)       :: reprod           !reproduction cost (gC/m2)    
  real(dp), dimension(npft,ncvar) :: bm_inc           !annual biomass increment (gC/m2) minus reproduction cost
  real(dp), dimension(npft,ncvar) :: annual_npp       !the total annual NPP accumulated by this PFT (gC/m2)
  real(dp), dimension(npft,ncvar) :: annual_gpp       !the total annual GPP for this PFT (gC/m2)
  real(dp), dimension(npft,ncvar) :: annual_resp      !annual gridcell leaf respiration (gC/m2) (annual autotropic + growth total respiration)-JM Apr 9 08
  real(dp), dimension(npft) :: lai_sun                !instantaneous sunlit leaf area index from lai_sno (m2 m-2)                    
  real(dp), dimension(npft) :: lai_sha                !instantaneous shaded leaf area index from lai_sno (m2 m-2)                    
  real(dp)                  :: fpc_sum                !instantaneous cumulative fpc_grid over all pfts for gridcell                  
  real(dp), dimension(npft) :: fpc_ind                !foliar protective cover per individual (m2 m-2)                               
  real(dp), dimension(npft) :: fpc_grid               !instantaneous gridcell foliar projective cover (FPC) for gridcell (m2 m-2)    
  real(dp), dimension(npft) :: wscal                  !mean daily water scalar (among leaf-on days) (0-1 scale)                      

end type veg_struct

type (veg_struct), save, target :: veg

!-------------------------
contains

subroutine initpftvars()

implicit none

integer :: k

veg%ratio         = 0._dp
veg%gt_sun        = 0._dp
veg%gt_sha        = 0._dp
veg%aet           = 0._dp
veg%daet          = 0._dp
veg%dreprod = 0._dp
veg%gdd_bioclim  = 0._dp
veg%awscal = 1._dp
veg%sla_sun = 0._dp          
veg%sla_sha = 0._dp 
veg%lai_sun = 0._dp        
veg%lai_sha = 0._dp   
veg%fpc_sum = 0._dp  
veg%fpc_grid  = 0._dp             
veg%fpc_ind = 0._dp    
veg%wscal = 1._dp    

forall (k = 1:ncvar)

veg%drh(k)          = 0._dp
veg%dresp(:,k)      = 0._dp
veg%gresp(:,k)      = 0._dp
veg%dgpp(:,k)       = 0._dp
veg%dnpp(:,k)       = 0._dp
veg%sgpp(:,k)       = 0._dp
veg%annual_gpp(:,k) = 0._dp
veg%annual_npp(:,k) = 0._dp
veg%annual_resp(:,k) = 0._dp

end forall

end subroutine initpftvars

end module pft_state_variables
