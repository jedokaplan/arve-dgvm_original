module gppmod

! Contains gpp, photosynthesis and Vcmax subroutines.

use arveparams,          only : dp,npft,pliq,g_m,daysec,Tfreeze,Rgas,umol2g,c12mass,minplant
use statevars,           only : sv,co2,dt,dtime,counter,counter_lim  !counters just brought in for testing. JM 13.01.2011
use pft_state_variables, only : veg
use pftparametersmod,    only : prm

implicit none

public  :: gpp
private :: optgpp
public  :: photosynthesis
public  :: calc_Vcmax
private :: ft

! module level parameters
real(dp), parameter :: drCO2H2O  = 1.6_dp    !volumetric diffusivity ratio of H2O / CO2
real(dp), parameter :: ftbase    = 2._dp     !base for use in the ft function
real(dp), parameter :: alpV_Ha_m = 1.391_dp  !maximum Priestley-Taylor coefficient (Gerten et al.2004)

!--------------------------------------

contains

subroutine gpp(j,i)

implicit none

! arguments
integer, intent(in)  :: j  !grid-cell index
integer, intent(in)  :: i  !subdaily timestep (1=day, 2=night)

! pointers to global variables
real(dp), pointer :: Patm                             !atmospheric pressure (Pa)
logical,  pointer, dimension(:) :: c4                 !true if pft is c4
real(dp), pointer, dimension(:) :: APAR_sun           !visible only shortwave rad absorbed by vegetation 
                                                      ! averaged for leaf area (W m-2) sunlight 
real(dp), pointer, dimension(:) :: APAR_sha           !visible only shortwave rad absorbed by vegetation 
                                                      ! averaged for leaf area (W m-2) shaded
real(dp), pointer, dimension(:) :: aet                !actual evapotranspiration (mm s-1)
real(dp), pointer, dimension(:) :: betaT              !soil moisture function limiting transpiration (fraction)
real(dp), pointer, dimension(:) :: demand             !plant water demand (mm s-1)
real(dp), pointer, dimension(:) :: gmin               !canopy conductance component (mm/s) not associated with photosynthesis
real(dp), pointer, dimension(:,:) :: gpp_shorttimestep  !Total timestep net photosynthesis at the gridscale, gC/m2/timestep
real(dp), pointer, dimension(:) :: gt_sun             !stomatal conductance (mm/s) sunlight
real(dp), pointer, dimension(:) :: gt_sha             !stomatal conductance (mm/s) shaded
real(dp), pointer, dimension(:) :: lai_sno            !instantaneous leaf area index per pft (m2 m-2)
real(dp), pointer, dimension(:) :: lai_vir            !VIRTUAL leaf area index per pft (m2 m-2) (for phenology)
real(dp), pointer, dimension(:) :: ratio              !stomatal conductance (gc) resulting from the maximum lambda value (ci/ca ratio)
real(dp), pointer, dimension(:) :: red_pet            !potential evapotranspiration for canopy reduced by canopy evaporation (mm s-1)
real(dp), pointer, dimension(:) :: supply             !water supply to the plants (mm s-1)
real(dp), pointer, dimension(:,:) :: leaftemp         !timestep average vegetation temperature (K)
real(dp), pointer, dimension(:) :: lai_sun            !instantaneous sunlight leaf area index per pft (m2 m-2)
real(dp), pointer, dimension(:) :: lai_sha            !instantaneous shaded leaf area index per pft (m2 m-2)
logical, pointer, dimension(:)  :: virtual_leaf       !true if the leaf being calculated is a virtual one to determine phenology.
logical, pointer, dimension(:)  :: present            !true if pft is present
real(dp), pointer, dimension(:) :: fpc_grid           !instantaneous gridcell foliar projective cover (FPC) per pft for gridcell (m2 m-2)
logical, pointer, dimension(:)  :: tree               !true if pft is tree
real(dp), pointer, dimension(:) :: fsun               !sunlit fraction of canopy
real(dp), pointer, dimension(:) :: fsha               !shaded fraction of canopy
real(dp), pointer, dimension(:) :: crownarea          !crown area (m2)
real(dp), pointer, dimension(:) :: nind               !individual density (indiv m-2)
                                                        
! parameters
real(dp), parameter :: emax   = 5._dp / 86400._dp  !maximum daily transpiration rate (mm/day converted to mm/s)

! local variables
real(dp) :: A_gross                             !gross photosynthesis !(umol CO2 m-2 s-1)
real(dp) :: Vcmax                               !maximum rate of Rubisco activity (umol m-2 s-1)
real(dp) :: adtmm                               !volume of gas exchanged during photosynthesis (mm timestep-1)
real(dp) :: Anet                                !net photosynthesis (umol CO2 m-2 s-1)
real(dp) :: Rd                                  !leaf respiration (umol m-2 s-1)
real(dp) :: max_gt                              !maximum stomatal conductance (mm/s)
real(dp) :: upper_bound                         !
real(dp) :: lower_bound                         !
real(dp) :: rtemp                               !leaf temperature for respiration calculation
integer :: pft,n                                !plant functional type
real(dp) :: lai_tot_temp                        !temporary sunlight and shaded LAI variable
real(dp) :: lai_suntemp                         !temporary sunlight LAI variable
real(dp) :: lai_shatemp                         !temporary shaded LAI variable
real(dp) :: lai                                 !temp lai variable
real(dp), dimension(2) :: gpp_temp              !
real(dp), dimension(2) :: supply_temp           !
real(dp), dimension(2) :: demand_temp           !
real(dp), dimension(2) :: aet_temp              !
real(dp), dimension(2) :: gt_temp               !
real(dp) :: APAR                                !temp APAR for sunlight or shaded leaves (W m-2)
real(dp) :: fpc_temp                            !temp fpc_grid (m2 veg m-2 grnd)

! point pointers
APAR_sun          => veg%APAR_sun
APAR_sha          => veg%APAR_sha
Patm              => sv(j)%Patm
aet               => veg%aet
betaT             => sv(j)%betaT
c4                => prm%c4
demand            => veg%demand
gmin              => prm%gmin
gt_sun            => veg%gt_sun
gt_sha            => veg%gt_sha
lai_sno           => sv(j)%lai_sno
lai_vir           => sv(j)%lai_vir
leaftemp          => sv(j)%temp_veg
ratio             => veg%ratio
red_pet           => veg%red_pet
supply            => veg%supply
gpp_shorttimestep => veg%sgpp
lai_sun           => veg%lai_sun
lai_sha           => veg%lai_sha  
virtual_leaf      => sv(j)%virtual_leaf 
present           => sv(j)%present    
fpc_grid          => veg%fpc_grid 
tree              => prm%tree
fsun              => veg%fsun
fsha              => veg%fsha
crownarea         => sv(j)%crownarea
nind              => sv(j)%nind

!-----------------------------------------------

!begin calculations

do pft = 1,npft
        
  if (present(pft)) then

      !set lai's for real or virtual leaves
       if (.not. virtual_leaf(pft)) then
         lai_tot_temp = lai_sno(pft)
         lai_suntemp = lai_sun(pft)
         lai_shatemp = lai_sha(pft)
         fpc_temp = fpc_grid(pft)
       else
         lai_tot_temp = lai_vir(pft)
         lai_suntemp = lai_tot_temp
         lai_shatemp = 0._dp    !virtual leaves have no shaded leaf area
         fpc_temp = 1.!crownarea(pft) * nind(pft) * (1._dp - exp(-0.5_dp * lai_tot_temp)) !FLAG try? JM 17.03.2011
       end if
  
   !loop over the sunlight and shaded leaves
   do n = 1,2  !1 = sunlight leaves, 2 = shaded leaves

         !Assign APAR the sunlight or shaded APAR value
         if (n == 1) then 
           APAR = APAR_sun(pft)
           lai = lai_suntemp
         else 
           APAR = APAR_sha(pft)
           lai = lai_shatemp
         end if 
   
    ! Calculate gross leaf level photosynthesis for the theoretical maximum value of ratio   
    if (APAR > 0._dp) then     ! photosynthesis only occurs when there is downwelling solar

         if (.not. c4(pft)) then
           ratio(pft) = 0.8_dp !initial value for c3
         else
           ratio(pft) = 0.4_dp !initial value for c4
         end if

         call photosynthesis(j,pft,APAR,n,Vcmax,A_gross)

    else ! no incoming solar, so no photosynthesis, but still dark respiration

         ! Perform Vcmax calculation for mitochondrial (dark) respiration calculation (umol m2 s-1)
         ! set it to the sunlight calculation (1) !FLAG ok? JM 06.03.2011
         call calc_Vcmax(j,pft,1,Vcmax)

         ! Set the gross leaf level photosynthesis to 0
         A_gross = 0._dp

         ! prevent unconstrained water loss during night by setting red_pet to zero
         red_pet(pft) = 0._dp 

    end if

    !leaftemp in deg C
    rtemp = leaftemp(2,pft) - Tfreeze

    ! Calc rate of mitochondrial respiration (c3) using CLM 3.0 formulation
    if (.not. c4(pft)) then        
      Rd = 0.015 * Vcmax * ft(ftbase,rtemp)   !(umol m-2 s-1) 
    else
      Rd = 0.02  * Vcmax * ft(ftbase,rtemp)   !(umol m-2 s-1)
    end if

    ! Net photosynthesis at LEAF LEVEL (umol CO2 leaf m-2 s-1)  
    Anet = A_gross - Rd

    ! change units of net photosynthesis (g C leaf m-2 timestep-1)
    gpp_temp(n) = Anet * umol2g * dt 
                  
    ! Net photosynthesis in terms of gas loss (mm s-1) at the LEAF level
    ! this is the absolute amount of gas that passes through the stomata.
    adtmm = ABS(Anet * umol2g / c12mass * Rgas * leaftemp(2,pft) / Patm)

    ! Relate photosynthetic uptake to transpiration (gt = stomatal conductance) (mm H2O s-1)
    ! converts [CO2] from ppm to mole fraction.This follows Gerten et al. 2004 formulation.

    if (lai > 0._dp) then
     
      ! In CLM photosynthesis canopy conductance is 1/rs * lai, so multiply by lai here. JM 08.03.2011
      gt_temp(n) = (gmin(pft) + drCO2H2O * adtmm / (co2(1) * 1.e-6 * (1._dp - ratio(pft)))) * lai 

      ! Calculate demand for water (Gerten et al. 2004 J. Hydrology)(already reduced demand
      ! by the canopy wetness in latentsens) at leaf level
      demand_temp(n) =  red_pet(pft) * alpV_Ha_m / (1._dp + g_m / gt_temp(n))  !mm s-1

    else

      gt_temp(n) = 0._dp
      demand_temp(n) = 0._dp

    end if
    
    ! Water supply is a function of the soil moisture in the root zone (mm s-1)
    supply_temp(n) = betaT(pft) * emax

    ! Determine actual evapotranspiration (AET) by comparing supply and demand, if demand is greater
    ! than supply reduce the gpp by reducing ratio until the demand equation is satisfied

    if (demand_temp(n) > supply_temp(n)) then 

          ! Supply-limited canopy conductance (mm s-1) (Haxeltine & Prentice 1996 Eqn. 25)
          max_gt = -g_m * log(1._dp - supply_temp(n) / (red_pet(pft) * alpV_Ha_m))

          ! Actual evapotranspiration at leaf level 
          aet_temp(n) = supply_temp(n)

          upper_bound = 1._dp
          lower_bound = 0.02_dp

          if (max_gt > 0._dp) then

            ! Enter optgpp to determine the optimal ratio of stomatal opening
            call optgpp(j,pft,n,max_gt,APAR,upper_bound,lower_bound,gt_temp(n),gpp_temp(n),lai)

          else

            ! this means that essentially supply is 0
            gt_temp(n)             = 0._dp
            gpp_temp(n)            = max(0._dp,gpp_temp(n))  !NOTE I made it so that gpp can be negative on days of 0 supply.

          end if !mag_gt

        else !demand is less than supply

          ! actual evapotranspiration at leaf level 
          aet_temp(n) = demand(pft)

    end if !demand/supply
       
    if (i == 2) then !night
    
      ! GPP is negative at night, assign total lai to lai_sun since calculation should
      ! exit before doing the shaded leaves (at night there is no shaded/sunlight). 
      lai_suntemp = lai_tot_temp
      lai_shatemp = 0._dp  
      gt_temp(1:2) = 0._dp  !no stomatal conductance under the assumption of complete stomatal closure at night.
   
      exit
      
    end if
      
  end do !sunlight/shaded
        
    !Now add the sunlight and shaded fractions as well as scale to the actual ground surface vegetated.    

    ! Find the total GPP (gC / m2 of vegetated surface)
    if (tree(pft)) then

    gpp_shorttimestep(pft,1) = (gpp_temp(1) * lai_suntemp + gpp_temp(2) * lai_shatemp) * fpc_temp

    ! total canopy level actual evapotranspiration
    aet(pft) = (aet_temp(1) * lai_suntemp + aet_temp(2) * lai_shatemp) * fpc_temp
   
    else  !grass
    
    ! Grass are treated as a 'green carpet'. The light interception of grass accounts for the scaling
    ! of how much of the ground surface is grass. Therefore further scaling here, like in the case
    ! of trees is not required. We scale using fpc_ind. JM 03.03.2011 
   
    gpp_shorttimestep(pft,1) = (gpp_temp(1) * (1. - exp(-0.5 * lai_suntemp)) + gpp_temp(2) * (1. - exp(-0.5 * lai_shatemp)))
       
    ! total canopy level actual evapotranspiration
    aet(pft) = (aet_temp(1) * lai_suntemp + aet_temp(2) * lai_shatemp) 
    
    end if  

!write(*,'(a5,2i4,2es12.4)')'gpp',sv(j)%d,pft,gpp_shorttimestep(pft),fpc_temp   

    ! supply and demand are compared to each other to determine wscal. The exact scaling is not 
    ! important thus fsun and fsha are sufficient to scale appropriately.
    
    ! total canopy level supply 
    supply(pft) = supply_temp(1) * fsun(pft) + supply_temp(2) * fsha(pft)
    
    ! total canopy level demand
    demand(pft) = demand_temp(1) * fsun(pft) + demand_temp(2) * fsha(pft) 
  
    ! assign the canopy level stomatal conductances sunlight/shaded
    gt_sun(pft) = gt_temp(1) * lai_suntemp
    gt_sha(pft) = gt_temp(2) * lai_shatemp
 
    ! reset values
    gpp_temp(:) = 0._dp
    supply_temp(:) = 0._dp
    demand_temp(:) = 0._dp
    aet_temp(:) = 0._dp
              
   ! If this was for a virtual leaf keep only canopy level gpp, the rest, reset to 0.
   if (virtual_leaf(pft)) then 
        gt_sun(pft)             = 0._dp
        gt_sha(pft)             = 0._dp
        aet(pft)                = 0._dp
        demand(pft)             = 0._dp
        supply(pft)             = 0._dp
   end if

 end if  !lai
end do !pft loop

end subroutine gpp

!---------------------------------------------------------------------------
!---Comments---(Old)
!from old gpp.f for when decay.f is used.
!assign DELTA 14C value from atmospheric time series to plant production
!    mgpp(:,pft,3)=co2(3)
!    agpp(pft,3)=co2(3)
!  else
!    mgpp(:,pft,2:3)=0.0
!    agpp(pft,2:3)=0.0
!---------------------------------------------------------------------------

subroutine optgpp(j,pft,n,max_gt,APAR,upper_bound,lower_bound,gt_temp,gpp_temp,lai)

implicit none

!arguments
integer,  intent(in)  :: j                              !grid index
integer,  intent(in)  :: pft                            !pft
real(dp), intent(in)  :: max_gt                         !
integer,  intent(in)  :: n                              !1 = sunlight calc. 2 = shaded
real(dp), intent(in)  :: APAR                           !visible only shortwave rad absorbed by vegetation averaged 
                                                        ! for leaf area (W m-2) sunlight or shaded
real(dp), intent(out) :: gt_temp                        !
real(dp), intent(out) :: gpp_temp                       !
real(dp), intent(in)  :: lai                            !lai, either sunlight or shaded

!pointers
real(dp), pointer                       :: Patm         !atmospheric pressure (Pa)
logical,  pointer, dimension(:)         :: c4           !true if pft is c4
real(dp), pointer, dimension(:)         :: ratio        !stomatal conductance (gc) resulting from the maximum lambda value (ci/ca ratio)
real(dp), pointer, dimension(:,:)       :: leaftemp     !timestep average vegetation temperature (K)

!local variables
real(dp) :: Vcmax               !maximum rate of Rubisco activity (umol m-2 s-1)
integer  :: k                   !
real(dp) :: small               !
real(dp) :: dratio              !
real(dp) :: upper_bound         !
real(dp) :: lower_bound         !
real(dp) :: Anet                !net photosynthesis (umol CO2 m-2 s-1)
real(dp) :: Rd                  !leaf respiration (umol m-2 s-1)
real(dp) :: rtemp               !leaf temperature for respiration calculation
real(dp) :: A_gross             !net photosynthesis integrated over timestep !(umol CO2 m-2 s-1)
real(dp) :: adtmm

!point pointers
c4                => prm%c4
ratio             => veg%ratio
leaftemp          => sv(j)%temp_veg
Patm              => sv(j)%Patm

!---------

small = 0.01_dp * max_gt

do k = 1,20

  dratio = upper_bound - lower_bound
  ratio(pft) = lower_bound + dratio / 2._dp

  ! Calculate gross photosynthesis for this value of ratio
  call photosynthesis(j,pft,APAR,n,Vcmax,A_gross)

  rtemp = leaftemp(2,pft) - Tfreeze

  ! Calc rate of mitochondrial respiration (c3) CLM 3.0 formulation
  if (.not. c4(pft)) then  
    
    Rd = 0.015 * Vcmax * ft(ftbase,rtemp)   !(umol m-2 s-1) 
    
  else

    Rd = 0.02  * Vcmax * ft(ftbase,rtemp)   !(umol m-2 s-1)  
    
  end if

  ! Net leaf level photosynthesis (umol CO2 leaf m-2 s-1)
  Anet = A_gross - Rd

  ! Net photosynthesis at the leaf level (g C leaf m-2 timestep-1)
  gpp_temp = Anet * umol2g * dt 

  ! Net photosynthesis in terms of gas loss (mm s-1) at the leaf scale level
  ! this is the absolute amount of gas that passes through the stomata.
  adtmm = ABS(Anet * umol2g / c12mass * Rgas * leaftemp(2,pft) / Patm)

  ! Relate transpiration to photosynthetic uptake (mm s-1)

  gt_temp = drCO2H2O * adtmm  / (co2(1) * 1.e-6 * (1._dp - ratio(pft))) * lai  !

  ! If this value is close to the maximum transpiration continue, else iterate to a solution

  if ((ratio(pft) - lower_bound) < 1.e-3 .or. (upper_bound - ratio(pft)) < 1.e-3) exit

  if (gt_temp < max_gt - small) then

    lower_bound = ratio(pft)

  else if (gt_temp > max_gt + small) then

    upper_bound = ratio(pft)

  else !convergence, optimal ratio found

    exit 

  end if
  
end do

end subroutine optgpp

!-----------------------------

subroutine photosynthesis(j,pft,APAR,n,Vcmax,A_gross)

! Calculates the gross leaf photosynthesis. Two options are available here. The flag
! do_CLM_photosynthesis in the joboptions file switches between the CLM 4.0 photosynthesis
! scheme and the LPJ scheme (Haxeltine & Prentice 1996)
! Joe Melton 08.2010

use iovariables, only : do_CLM_photosynthesis

implicit none

! arguments
integer, intent(in)  :: j   !grid index
integer, intent(in)  :: pft !pft
real(dp), intent(in) :: APAR                    !visible only shortwave rad absorbed by vegetation 
                                                ! averaged for leaf area (W m-2) sunlight or shaded
integer, intent(in) :: n                        !1 = sunlight leaves calc, 2 = shaded leaves calc        
real(dp), intent(out) :: Vcmax                   !maximum rate of Rubisco activity (umol m-2 s-1)                              
real(dp), intent(out) :: A_gross                    !actual rate of gross photosynthesis (g C m-2 s-1)

! pointers
real(dp), pointer, dimension(:) :: ratio        !the lambda ratio of stomatal opening (fraction)
real(dp), pointer :: Patm                       !atmospheric pressure (Pa)
real(dp), pointer :: tlim_top                   !high temperature limit for CO2 unptake
real(dp), pointer :: tlim_bot                   !low temperature limit for CO2 uptake
real(dp), pointer :: topt_bot                   !lower range of temperature optimum for photosynthesis
real(dp), pointer :: topt_top                   !upper range of temperature optimum for photosynthesis
real(dp), pointer :: lambdam                    !lambda ratio of optimal stomatal opening
logical, pointer  :: c4                         !
real(dp), pointer :: temp_veg                   !vegetation temperature (K)
real(dp), pointer, dimension(:) :: CO2comp                    !

!local variables
real(dp) :: Ci                          !CO2 intercellular pressure
real(dp) :: leaftemp                    !temperature of the leaf (deg C)
real(dp) :: APAR_temp                   !temporary variable.
real(dp) :: Wc                          !rate of gross photosynthesis when the photosynthetic enzyme&
                                        ! system (RuBP) is limiting (umol CO2 m-2 s-1)
real(dp) :: Wl                          !light limited rate of gross photosynthesis (umol CO2 m-2 s-1)
real(dp) :: k1,k2,k3,low,high           !intermediate variables
real(dp) :: tstress                     !temperature inhibition function
real(dp) :: phipi                       !scalar to reduce the photosynthesis rate below its optimum
                                        ! value, captures effect of reduced intercellular pressure
real(dp) :: c1,c2                       !functions from Eqn. 14 & 15 in Haxeltine & Prentice 1996
real(dp) :: Kc                          !Michaelis-Menten co-efficients (umol mol-1)
real(dp) :: Ko                          !mmol mol-1
real(dp) :: pO2                         !O2 partial pressure (Pa)
real(dp) :: b_term                      !
real(dp) :: bot                         !used in the Michaelis-Menten co-efficients
real(dp) :: leafT                       !leaf temp - rt
real(dp) :: wj                          !
real(dp) :: we                          !
real(dp) :: baser                       !

!parameters
real(dp), parameter :: phot2J = 4.6_dp          !conversion from photons to J (umol photons / J)
real(dp), parameter :: alphac4 = 0.053_dp       !C4 intrinsic quantum efficiency
real(dp), parameter :: alphac3 = 0.08_dp        !C3 intrinsic quantum efficiency
real(dp), parameter :: bc4 = 0.02_dp            !leaf respiration as fraction of vmax for C4 plants
real(dp), parameter :: bc3 = 0.015_dp           !leaf respiration as fraction of vmax for C3 plants
real(dp), parameter :: theta = 0.7_dp           !colimitation (shape) parameter
real(dp), parameter :: slo2    = 20.9_dp        !O2 partial pressure in kPa
real(dp), parameter :: a = 42.75                !parameters used in the Ko and Kc calculations
real(dp), parameter :: b = 37830._dp
real(dp), parameter :: c = 404.9
real(dp), parameter :: g = 79430._dp
real(dp), parameter :: e = 278.4
real(dp), parameter :: f = 36380._dp
real(dp), parameter :: rt = 298._dp             !room temp in degrees kelvin
real(dp), parameter :: tmc3 = 45._dp            !maximum temperature for C3 photosynthesis
real(dp), parameter :: tmc4 = 55._dp            !maximum temperature for C4 photosynthesis
real(dp), parameter, dimension(npft) :: quanteffic = [ 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.4 ] 
real(dp), parameter :: maxO2CO2 = 0.21          !maximum rates of oxygenation to carboxylation (Farquhar and Caemmerer 1982)

!point pointers
CO2comp                   => veg%CO2comp
ratio                     => veg%ratio
temp_veg                  => sv(j)%temp_veg(2,pft)
Patm                      => sv(j)%Patm
tlim_top                  => prm(pft)%tlim_top
tlim_bot                  => prm(pft)%tlim_bot
topt_bot                  => prm(pft)%topt_bot
topt_top                  => prm(pft)%topt_top
lambdam                   => prm(pft)%lambdam
c4                        => prm(pft)%c4

!-------------------------------------------------------------------
!begin calculations

!Put veg temp into celcius
leaftemp = temp_veg - Tfreeze
leafT = temp_veg - rt !leaf temp minus room temp.

! I suspect that CLM photosynthesis works for sunlight/shaded but that LPJ formulation will
! not. To use the LPJ formulation likely best to formulate a 'big leaf' model. JM 08.03.2011

if (do_CLM_photosynthesis) then

    ! find the Michaelis-Menten constants (Pa) for
    ! CO2
    baser = 2.1
    Kc = 30. * ft(baser,leaftemp)
    
    ! and O2
    baser = 1.2
    Ko = 30000. * ft(baser,leaftemp)

    ! The leaf internal O2 partial pressure (Pa) 
    pO2 = 0.209 * Patm

    !The CO2 compensation point is (Pa)
    CO2comp(pft) = 0.5 * Kc / Ko * maxO2CO2 * pO2 

    ! intercellular CO2 partial pressure in Pa
    Ci = ratio(pft) * Patm * co2(1) * 1.e-6

    ! Use the CLM 4.0 Vcmax formulas
       call calc_Vcmax(j,pft,n,Vcmax)

    if (.not. c4) then

       ! find the RuBP carobylase (Rubisco) limited rate of carboxylation
       ! (umol CO2 m-2 s-1)
       Wc = Vcmax * max(0._dp,(ci - CO2comp(pft))) / (ci + Kc*(1. + pO2 / Ko))

       ! find the maximum rate of carboxylation allowed by the capacity to regenerate
       ! RuBP (the light limited rate)(umol CO2 m-2 s-1)
       wj = max(0._dp,(ci - CO2comp(pft)))* APAR * phot2J * quanteffic(pft) / (ci + 2._dp * CO2comp(pft))

       ! find the export limited rate of carboxylation for C3 plants and the PEP
       ! carboxylase limited rate of carboxylation for C4 plants (umol CO2 m-2 s-1)
       we = 0.5 * Vcmax

    else !C4 photosynthesis

       ! find the RuBP carobylase (Rubisco) limited rate of carboxylation
       ! (umol CO2 m-2 s-1)
       Wc = Vcmax

       ! find the maximum rate of carboxylation allowed by the capacity to regenerate
       ! RuBP (the light limited rate)(umol CO2 m-2 s-1)
       wj = APAR * phot2J * quanteffic(pft)

       ! find the export limited rate of carboxylation for C3 plants and the PEP
       ! carboxylase limited rate of carboxylation for C4 plants (umol CO2 m-2 s-1)
       we = 4000._dp * Vcmax * Ci / Patm 

    end if

    ! leaf level gross photosynthesis (umol CO2 m-2 s-1)
    A_gross = min(wc,wj,we) 
   
else !LPJ photosynthesis

    APAR_temp = APAR * phot2J * 1.e-6  !convert Watts/m2 (J/s/m2) to mol photons m-2 s-1

    !calculate temperate inhibition function
                if (leaftemp < tlim_top) then
                   k1 = 2._dp * log(0.010101_dp) / (tlim_bot - topt_bot)
                   k2 = (tlim_bot + topt_bot) / 2._dp
                   low = 1._dp / (1._dp + exp(k1 * (k2 - leaftemp)))
                   k3 = log(99._dp) / (tlim_top - topt_top)
                   high = 1._dp - 0.01_dp * exp(k3 * (leaftemp - topt_top))
                   tstress = (low * high)
                else
                   tstress = 0._dp
                end if

                if (tstress < 1.e-2) tstress = 0._dp

            !Non-water-stressed intercellular CO2 partial pressure in Pa
            !Eqn 7, Haxeltine & Prentice 1996
            !Eqn 8, Haxeltine & Prentice 1996

            Ci = co2(1) * 1.e-6 * Patm

            !pO2 = slo2 * 1.e6 / Patm * ratio
            pO2 = slo2 * 1.e3

    !---Optimal photosynthesis--
    !first calculate catalytic capacity of rubisco, Vcmax, assuming optimal
    !(non-water-stressed) value for ratio, i.e. lambdam

        if (.not.c4) then  !C3 photosynthesis

            bot = rt * (Rgas * 1.e-3) * temp_veg         !used in the Michaelis-Menten co-efficients

            !Find the in-vivo temperature dependence of the Michaelis-Menten coefficients of Rubisco, Kc and Ko
            !calculate Kc, (relation originally from Bernacchi et al. 2001)
              Kc = c * exp(g * leafT / bot) * 1.e-6 * Patm        !(II eqn 5)

            !calculate Ko, (from Bernacchi et al 2001)
              Ko = e * exp(f * leafT / bot) * 1.e-3 * Patm   !(II eqn 6)

            !calculate the temperature dependence of the CO2 compensation point (from Bernacchi et al 2001)
              CO2comp(pft) = a * exp(b * leafT / bot) * 1.e-6 * Patm !(II eqn 12) (Pa)

            !Calculation of C1C3, Eqn 4, Haxeltine & Prentice 1996
            !Notes: - there is an error in this equation in the above paper (missing
            !          2.0* in denominator) which is fixed here (see Eqn A2, Collatz
            !          et al 1991)

            c1 = tstress * alphac3 * ((Ci - CO2comp(pft)) / (Ci + 2.0 * CO2comp(pft)))

            ! High temperature inhibition modelled primarily by suppression of LUE
            ! by decreased relative affinity of rubisco for CO2 relative to O2 with
            ! increasing temperature, but we also implement a step function to
            ! prohibit any C3 photosynthesis above 45 degrees (Table 3.7, Larcher
            ! 1983)

            if (leaftemp > tmc3) c1 = 0._dp

            !Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996

            c2 = (Ci - CO2comp(pft)) / (Ci + Kc * (1._dp + pO2 / Ko))

            b_term = bc3   !Choose C3 value of b_term for Eqn 10, Haxeltine & Prentice 1996

        else !C4

            !Specify C1, C2
            !Eqns 14,15, Haxeltine & Prentice 1996

            c1 = tstress * alphac4

             !High-temperature inhibition modelled conservatively as a step function
             !prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)
            if (leaftemp > tlim_top) c1 = 0._dp

            c2 = 1._dp

            b_term = bc4

        end if

            !Use the CLM Vcmax formulas
            call calc_Vcmax(j,pft,n,Vcmax)

            !Now use this Vcmax value to calculate actual photosynthesis

          if (.not.c4) then  !C3 photosynthesis

              Ci = ratio(pft) * Patm * co2(1) * 1.e-6

            !Recalculation of C1C3, C2C3 with actual Ci

            c1 = tstress * alphac3 * ((Ci - CO2comp(pft)) / (Ci + 2._dp * CO2comp(pft)))

            if (leaftemp > tmc3) c1 = 0._dp  !high-temperature inhibition

            c2 = (Ci - CO2comp(pft)) / (Ci + kc * (1._dp + pO2 / ko))

          else  !C4 photosynthesis

            !find phipi, parameter accounting for effect of reduced intercellular CO2
            !concentration on photosynthesis.Eqn 14,16, Haxeltine & Prentice 1996
            !Fig 1b, Collatz et al 1992
            phipi = min(ratio(pft) / lambdam, 1._dp)

            !recalculate c1 with the effect of reduced intercellular CO2 concentration
            c1 = tstress * phipi * alphac4

            if (leaftemp > tmc4) c1 = 0._dp  !high-temperature inhibition

          end if

    !Calculation of PAR-limited photosynthesis rate, Wl, umol CO2 /m2/ s
    !Eqn 3, Haxeltine & Prentice 1996
          Wl = c1 * APAR_temp * c12mass / umol2g

    !Calculation of rubisco-activity-limited photosynthesis rate Wc, umol CO2/m2/ s
    !Eqn 5, Haxeltine & Prentice 1996  
          Wc = c2 * Vcmax 

    !Calculation of daily gross photosynthesis, Agd, umol CO2 /m2/timestep
    !Eqn 2, Haxeltine & Prentice 1996
    !Note: - there is an error in this equation in the above paper (missing
    !        theta in 4*theta*Wl*Wc term) which is fixed here

    if (Wl < 0._dp .or. Wc <= 0._dp) then
          A_gross = 0._dp  !do not allow negative gross photosynthesis
    else
          A_gross = (Wl + Wc - sqrt((Wl + Wc)**2 - 4._dp * theta * Wl * Wc))&
                 / (2._dp * theta)
    end if

end if !CLM vs. LPJ
 
end subroutine photosynthesis

!--------------------------------------------------

subroutine calc_Vcmax(j,pft,n,Vcmax)

!Calculates Vcmax using CLM 4.0 formulation.

implicit none

integer, intent(in) :: j        !grid index
integer, intent(in) :: pft      !pft
integer, intent(in) :: n        !1 is sunlight calc, 2 is shaded
real(dp), intent(out) :: Vcmax  !maximum rate of Rubisco activity (umol m-2 s-1)

!pointers
real(dp), pointer :: F_lnr
real(dp), pointer :: f_N
real(dp), pointer :: leafcn                     !leaf carbon to N ratio (g C g-1 N)
real(dp), pointer :: betaT                      !soil moisture function limiting transpiration (fraction)
real(dp), pointer :: maxdl                      !longest day length of the year for this cell (s)
real(dp), pointer :: temp_veg                   !vegetation temperature (K)
real(dp), pointer :: sla_sun      !specific leaf area sunlight (m2 gC-1)
real(dp), pointer :: sla_sha      !specific leaf area shaded (m2 gC-1)
!real(dp), pointer :: fsun         !sunlit fraction of canopy
!real(dp), pointer :: fsha         !shaded fraction of canopy


!parameters
real(dp), parameter :: F_nr  = 7.16          !mass ratio of total rubisco molec. mass to N in rubisco (g Rub / gN in Rub)
real(dp), parameter :: a_r25  = 60._dp       !specific activity of rubisco (umol CO2 g-1 Rub s-1)

!variables
real(dp) :: f_DYL                  !function to scale Vcmax for daylength and seasonal variation
real(dp) :: f_Tv                   !function to scale Vcamx by temp and mimic thermal breakdown of metabolic processes
real(dp) :: Na                     !area-based N conc.(gN m-2 leaf area)
real(dp) :: Vcmax25                !Vcmax at 25 deg C
real(dp) :: leaftemp               !temperature of the leaf (deg C)
real(dp) :: sla_temp

!point pointers
F_lnr                           => prm(pft)%F_lnr
f_N                             => prm(pft)%f_N
leafcn                          => prm(pft)%leafcn
betaT                           => sv(j)%betaT(pft)
maxdl                           => sv(j)%maxdl
temp_veg                        => sv(j)%temp_veg(2,pft)
sla_sun                         => veg%sla_sun(pft)
sla_sha                         => veg%sla_sha(pft)
!fsun                            => pft_state(pft)%fsun 
!fsha                            => pft_state(pft)%fsha 

!--------------------

!Put veg temp into celcius
leaftemp = temp_veg - Tfreeze

        !find the mean instantaneous SLA
        
        ! NOTE: when it enters this at night, it is presently using the sunlight sla.
        ! this might need to be changed to the shaded value. JM 08.03.2011
        
        !sla_temp = sla_sun * fsun + sla_sha * fsha  !FLAG JM 04.03.2011
        if (n == 1) then
        sla_temp = sla_sun
        else
        sla_temp = sla_sha
        end if

        !find area-based N conc.(gN m-2 leaf area)
        Na = 1._dp / (leafcn * sla_temp)

        !find Vcmax at 25 deg C
        Vcmax25 = Na * F_lnr * F_nr * a_r25

        !function to scale Vcmax for daylength and seasonal variation
        f_DYL = min(1._dp,(max(0.01_dp,(dtime*dtime)/(maxdl*maxdl))))

        !function to scale Vcamx by temp and mimic thermal breakdown of metabolic processes
        f_Tv = 1._dp /(1._dp + exp ((-220000._dp + 710._dp * temp_veg) / (0.001 * Rgas * temp_veg)))

        !find Vcmax adjusting to temp, thermal breakdown of metabolic processes, daylength, and soil water.
        Vcmax = Vcmax25 * ft(2.4_dp,leaftemp) * f_Tv * BetaT * f_DYL * f_N  

end subroutine calc_Vcmax

!-----------------------------

real(dp) function ft(base,temp) !Q10 relationship for leaf respiration

implicit none

real(dp), intent(in) :: base !
real(dp), intent(in) :: temp !temperature (C)

ft = base**(0.1 * (temp - 25._dp))

end function ft

!-----------------------------------------------------------

end module gppmod
