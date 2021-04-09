module latentsensiblefluxes

implicit none

public :: latentsens

contains

!-----------------------------------------------------------------------
subroutine latentsens(j)

! Calculate the 'equilibrium' potential evapotranspiration and sensible
! heat flux for bare soil, ground under canopy, or vegetated surfaces

!NOTE: At present, ARVE uses the CLM method for calculation of sensible
!heat fluxes but uses a variant of the Prescott equation (Prentice et al.
!2003) for latent heat flux. This approach for the latent heat fluxes was chosen 
!as we do not presently take in the relative humidity of the location. One 
!weakness of this approach is that dew drop is not possible. The code is set up
!to allow an easy change for inclusion of a more robust approach to latent heat fluxes.
!Further improvements to ARVE should include RH in the weather generator . JM 24.05.2011

use arveparams,  only : dp,minplant,rh,lvap,grav,rwv,pliq,pice,pi,Cair,npft,sigsb
use statevars,  only : sv,dt
use esatdesdT, only : esat,desdT
use soilstate_vars, only : Swg_b,Swg_v,Lwg_b,Lwg_v,Lwv,Swv,ratm,TairK,Beta_soil,dqsath_dT,dqgdTg &
                                ,qg,eatm,Tair,Eg_b,dEgdt_b,Hg_b,dHgdt_b,theta,Hg_v,dEgdt_v,dHgdt_v,dLgdt_v, &
                                lam,Eveg,qs,Eg_v,surf
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j                                !grid-cell index

!local variables
integer :: pft
real(dp) :: gamma                                       !psychrometer constant (Pa K-1)                                            
real(dp) :: ss                                          !rate of increase of saturated vapor pressure with temperature (Pa K-1)   
real(dp) :: gamma_soil                                  !psychrometer constant (Pa K-1) for soil temperature                                            
real(dp) :: ss_soil                                     !rate of increase of saturated vapor pressure with (soil) temperature (Pa K-1)   
 real(dp) :: qatm                                        !atmos. specific humidity (kg kg-1)                                             
real(dp) :: qsath                                       !saturation specific humidity (kg kg-1)                                        
real(dp) :: alpha                                       !alpha factor for soil humidity calculation (combined)                     
real(dp) :: asoil                                       !alpha factor for soil humidity calculation (soil)                         
real(dp) :: r_stomata_sun                               !sunlight stomatal resistance (s/mm)                                           
real(dp) :: r_stomata_sha                               !shaded stomatal resistance (s/mm)                                             
real(dp) :: surf_wet                                    !surface soil wetness (fraction)                                            
real(dp) :: sPsi                                        !surface soil water matric potential (mm)                                       
real(dp) :: theta_1                                     !soil surface volumetric water content                                       
real(dp) :: prop                                        !theta_1 / theta_fc                                                             
real(dp) :: esat_val                                    !value of esat from called function                                         
real(dp) :: Evpot                                       ! potential evaporation from wet foliage per unit wetted area                  
real(dp) :: cah                                         !sensible heat conductances (canopy air to atm)                          
real(dp) :: cgh                                         !        "                (ground to canopy)                             
real(dp) :: cvh                                         !        "                (leaf to canopy)                               
real(dp) :: caw                                         !latent heat conductances (canopy air to atm)                            
real(dp) :: cgw                                         !        "                (ground to canopy)                             
real(dp) :: cvw                                         !        "                (leaf to canopy)                               
real(dp) :: Edem_canopy                                 !demand for evaporation from the wetted canopy                                  
real(dp) :: canopy_avail                                !water in the canopy that is available for evaporation                          
real(dp) :: fsno_filt                                   !snow cover of litter                                                      
real(dp) :: r_litter                                    !litter resistance (s m-1)                                                  
real(dp) :: litt_eff                                    !find effective litter area index (m2 m-2) (area not covered by snow)       

!pointers
real(dp), pointer, dimension(:,:) :: temp_veg           !timestep average vegetation temperature (K)
real(dp), pointer, dimension(:) :: pet                  !potential evapotranspiration for canopy (mm s-1)
real(dp), pointer, dimension(:) :: red_pet              !potential evapotranspiration for canopy reduced by canopy evaporation (mm s-1)
real(dp), pointer :: gpet                               !potential evapotranspiration for bare ground with no veg (mm s-1)
real(dp), pointer, dimension(:) :: underpet             !potential evapotranspiration for ground beneath canopy (mm s-1)
real(dp), pointer, dimension(:) :: Tsoil                !soil temperature (K)
real(dp), pointer, dimension(:) :: betaT                !soil moisture function limiting transpiration (fraction)
real(dp), pointer, dimension(:) :: lai_sno              !leaf area index (with snow burial)
real(dp), pointer, dimension(:) :: stem_sno             !stem "                "        "
real(dp), pointer, dimension(:) :: lai_sun              !leaf area index (sunlight)
real(dp), pointer, dimension(:) :: lai_sha              !leaf are index (shaded)
integer, pointer :: snl                                 !soil layers
real(dp), pointer :: Patm                               ! atmospheric pressure (Pa)
real(dp), pointer                       :: fsnow        !fraction of the gridcell covered by snow (fraction)
real(dp), pointer :: zsno                               !snow depth (m)
real(dp), pointer :: raw_b                              !resistance to latent heat transfer (bare ground and atm)(s m-1)
real(dp), pointer, dimension(:) :: raw                  !        "        "        "        (canopy and atm)(s m-1)
real(dp), pointer, dimension(:) :: rawp                 !        "        "        "        (canopy and ground)(s m-1)
real(dp), pointer :: rah_b                              !resistance to sensible heat transfer (bare ground and atm)(s m-1)
real(dp), pointer, dimension(:) :: rah                  !        "        "        "        (canopy and atm)(s m-1)
real(dp), pointer, dimension(:) :: rahp                 !        "        "        "        (canopy and ground)(s m-1)
real(dp), pointer, dimension(:) :: rdbp                 !fraction of potential evapotranspiration
real(dp), pointer, dimension(:) :: rb                   !leaf boundary layer resistance (s m-1)
real(dp), pointer, dimension(:) :: Wcan                 !water in the canopy (mm)
real(dp), pointer, dimension(:) :: f_wet                !fraction of the canopy that is wet
real(dp), pointer, dimension(:) :: f_dry                !fraction of the canopy that is dry
real(dp), pointer, dimension(:) :: gt_sun               !canopy conductance (from last time step) (mm/s) sunlight
real(dp), pointer, dimension(:) :: gt_sha               !canopy conductance (from last time step) (mm/s) shaded
real(dp), pointer :: theta_fc                           !the volumetric water content at field capacity (fraction)
real(dp), pointer :: ustar_b                            !friction velocity estimated from scalar wind speed (m s-1)
real(dp), pointer, dimension(:) :: rdp_dry              !term for contributions by wet leaves and transpiration, limited by avail H2O and pot. evapo.
real(dp), pointer, dimension(:) :: Psat                 !
real(dp), pointer, dimension(:) :: Tsat                 !
real(dp), pointer, dimension(:) :: dz                   !
real(dp), pointer, dimension(:) :: Bexp                 !
real(dp), pointer, dimension(:) :: Wliq                 !
real(dp), pointer, dimension(:) :: Wice                 !
real(dp), pointer, dimension(:) :: can_evap             !canopy evaporation (mm s-1)

!parameter
real(dp), parameter :: asnow = 1._dp                    !alpha factor for soil humidity calculation (snow)
real(dp), parameter :: s_psi_max = -1.e8                !limit for surface soil matric potential (mm)
real(dp), parameter :: dz_litter = 0.05_dp              !assumed typical litter depth (m)
real(dp), parameter :: leaf_lit = 0.5_dp                !assumed effective litter area index (m2 m-2)

!point pointers
ustar_b    => surf%ustar_b
temp_veg   => sv(j)%temp_veg
pet        => veg%pet
red_pet    => veg%red_pet
underpet   => veg%underpet
Tsoil      => sv(j)%Tsoil
snl        => sv(j)%snl
gpet       => veg%gpet
Patm       => sv(j)%Patm
raw_b      => surf%raw_b
rah_b      => surf%rah_b
raw        => surf%raw
rawp       => surf%rawp
rah        => surf%rah
rahp       => surf%rahp
rdbp       => surf%rdbp
Wcan       => veg%Wcan
f_wet      => veg%f_wet
f_dry      => veg%f_dry
lai_sun    => veg%lai_sun
lai_sha    => veg%lai_sha
gt_sun     => veg%gt_sun
gt_sha     => veg%gt_sha
betaT      => sv(j)%betaT
Psat       => sv(j)%Psat
Tsat       => sv(j)%Tsat
dz         => sv(j)%dz
Bexp       => sv(j)%Bexp    
Wliq       => sv(j)%Wliq    
Wice       => sv(j)%Wice    
theta_fc   => sv(j)%theta_fc
rb         => surf%rb
rdp_dry    => surf%rdp_dry
can_evap   => veg%can_evap
lai_sno    => sv(j)%lai_sno
stem_sno   => sv(j)%stem_sno
zsno       => sv(j)%zsno
fsnow      => surf%fsnow

!-----------------------

! Calculate the latent heat flux for bare ground first

 ! Initial pre-calcs.
        esat_val = esat(Tsoil(snl+1))

        qsath = 0.622 * esat_val/ (Patm - 0.378 * esat_val)  !CLM 4.0 eqn 5.142

        dqsath_dT = 0.622 * Patm * desdT(Tsoil(snl+1)) / ((Patm - 0.378 * esat_val)*(Patm - 0.378 * esat_val))  !CLM 4.0 eqn 5.143

   !r_litter is the litter resistance, introduced in Sakaguchi and Zeng 2009.

        !find snow cover of litter
        fsno_filt = zsno / dz_litter

        !find effective litter area index (m2 m-2) (area not covered by snow)
        litt_eff = leaf_lit * (1._dp - min(fsno_filt, 1._dp))

        !litter resistance (s m-1) CLM 4.0 Eqn 5.106
        r_litter = 1._dp / (0.004 * ustar_b) * (1._dp - exp(-litt_eff))

   ! Find the surface wetness and surface soil water matric potential (CLM 5.65-5.66)
        surf_wet = 1._dp / (dz(1) * Tsat(1)) * ((Wliq(1) / pliq) + (Wice(1) / pice))
        surf_wet = min(surf_wet, 1._dp)
        surf_wet = max(surf_wet, 0.01_dp)
        sPsi = max(Psat(1) * surf_wet**(-Bexp(1)), s_psi_max)

        asoil = exp((sPsi * grav) / (1.e3 * Rwv * Tsoil(1)))  !CLM 5.64
        alpha = asoil * (1._dp - fsnow) + asnow * fsnow   !CLM 5.63

        qg = qsath * alpha

        ! eatm (atmospheric vapour pressure) is precalculated in calcshortstep
        qatm = (0.622 * eatm) / (Patm - 0.378 * eatm)  !CLM

        ! add limiting condition to prevent large increase (decreases) in qg for
        ! small increases (decreases) in soil moisture in very dry soils (CLM 3 tech note pg .58)
        if (qsath > qatm .and. qatm > qg) then
                qg = qatm
        end if

        !-
        ! CLM 4.0 added in a new Beta function to represent the molecular diffusion process
        ! from the soil pore to the surface within the very dry part of the soil (ref: Sakaguchi
        ! and Zeng 2009). This will act to limit evaporation.

        ! soil surface volumetric water content (CLM 4 eq. 5.69)
        theta_1 = min(1._dp,1._dp / dz(1) * ((Wliq(1) / pliq) + (Wice(1) / pice)))
        theta_1 = max(theta_1, 0.01_dp)

        ! the volumetric water content at field capacity (theta_fc) is derived by assuming a hydraulic
        ! conductivity of 0.1mm/day and inverting the hydraulic conductivity function (CLM 4 eq. 5.70)
        ! this is precalculated in soilinit

        if (theta_1 >= theta_fc .or. (qatm - qg) > 0._dp) then
                Beta_soil = 1._dp
        else if (theta_1 < theta_fc) then
                prop = min(theta_1 / theta_fc, 1._dp)
                prop = max(prop, 0.01_dp)
                Beta_soil = 0.25 * (1._dp - fsnow) * (1._dp - cos(pi * prop))**2 + fsnow
        end if

        ! Find the latent heat flux from bare soil (gpet)

        !----------Prescott        
        !psychrometer constant, weakly temperature dependent (uses degrees C, but units of gamma are Pa K-1)
        gamma_soil = 65.05_dp + Tsoil(snl+1) * 0.064_dp

        !find the rate of increase of saturated vapor pressure with temperature (Pa K-1)
        ss_soil = desdT(Tsoil(snl+1))

        gpet = max((ss_soil / (ss_soil + gamma_soil)) * (Swg_b - Lwg_b) / lvap, 0._dp)  !(mm s-1)
        !----------
        
        !Not used at present.
        ! gpet = -ratm * Beta_soil * (qatm - qg) / raw_b  !CLM 4.0 Eqn 5.62
        !---
        
        Eg_b = gpet

        dqgdTg = alpha * dqsath_dT  !CLM 4.0 Eqn 5.74

        dEgdt_b = Beta_soil * ratm * dqgdTg / raw_b   !CLM 4.0 Eqn 5.73

    ! Calculate sensible heat flux between the air and soil surface (CLM eqn 5.60, CLM 4.0 Eqn 5.61)
        Hg_b = -ratm * Cair * (theta - Tsoil(snl+1)) / rah_b

        dHgdt_b = ratm * Cair / rah_b

! Now find the latent and sensible heat fluxes for vegetation and the ground below vegetation
do pft = 1,npft
    if (lai_sno(pft) + stem_sno(pft) == 0._dp) then !no vegetation

            !no canopy hence 0
            pet(pft) = 0._dp
            underpet(pft) = 0._dp
            red_pet(pft) = 0._dp
            can_evap(pft) = 0._dp
            Hg_v(pft) = 0._dp
            dEgdt_v(pft) = 0._dp
            dHgdt_v(pft) = 0._dp
            dLgdt_v(pft) = 0._dp

            cycle

    else if (lai_sno(pft) + stem_sno(pft) > minplant) then !vegetated surfaces

        ! Calculate vegetation potential evapotranspiration using the LPJ relation

            !----------Prescott
            !psychrometer constant, weakly temperature dependent (uses degrees C, but units of gamma are Pa K-1)
            gamma = 65.05_dp + Tair * 0.064_dp

            !find the rate of increase of saturated vapor pressure with temperature (Pa K-1)
            ss = desdT(TairK)

            ! variant of the Prescott equation (Prentice et al. 1993)
            pet(pft) = max((ss / (ss + gamma)) * (Swv(pft) - Lwv(pft)) / lvap, 0._dp)  !(mm s-1)
            !----------

            ! Reduce the PET by the amount of water that can be evaporated from the wetted portion of the canopy.
                     if (pet(pft) > 0._dp .and. Wcan(pft) > 0._dp) then

                            ! The water flux from wetted leaves and stems can be found from the amount of water in the canopy
                             canopy_avail = Wcan(pft) / dt !(mm s-1)

                             ! The PET is decreased by the fraction of leaf and stem area that is wet and non-transpiring
                             red_pet(pft) = pet(pft) * f_dry(pft)

                             ! The demand for evaporation from the wetted canopy is then
                             Edem_canopy = pet(pft) * f_wet(pft)

                             if (canopy_avail > Edem_canopy) then

                                 can_evap(pft) = Edem_canopy !(mm s-1)

                             else ! the evaporative demand is more than what is available in the canopy

                                 can_evap(pft) = canopy_avail   !(mm s-1)

                             ! the excess demand goes back to red_pet
                             red_pet(pft) = red_pet(pft) + Edem_canopy - canopy_avail

                             end if

                       else if (pet(pft) > 0._dp .and. Wcan(pft) == 0._dp) then !evaporation is not reduced and canopy is dry

                             can_evap(pft) = 0._dp
                             red_pet(pft) = pet(pft)

                   else if (pet(pft) <= 0._dp) then

                  ! Condensing conditions can not exist as we have no way of knowing the amount of water to condense
              
                        can_evap(pft) = 0._dp
                        red_pet(pft) = 0._dp

                       end if !condensing/evap

        ! Calculate the latent and sensible heat flux below vegetation

             ! use eatm, qatm,qsath from the soil temp as we do not uniquely calc the veg temp

            ! potential evaporation from wet foliage per unit wetted area
            ! this uses the last timestep's qs value.qs is canopy specific humidity. Somewhat
            ! meaningless at present as we do not track canopy humidity...
            Evpot = - ratm * (qs - qsath) / rb(pft)                !CLM 4.0 Eqn 5.95

            !calc fraction of potential evapotranspiration !CLM 4.0 Eqn 5.94
            if (Evpot > 0._dp .and. betaT(pft) > 1.e-10) then

                 !NOTE: this is using the stomatal resistance from the time step before.
                 !FLAG CHECK if this is ok! JM 21.12.2010
                if (gt_sun(pft) > 0. .and. gt_sha(pft) > 0. .and. lai_sno(pft) > 0._dp) then   
                 r_stomata_sun = 1._dp / gt_sun(pft)
                 r_stomata_sha = 1._dp / gt_sha(pft)
              
                rdp_dry(pft) = (f_dry(pft) * raw_b / lai_sno(pft)) * (lai_sun(pft) / (raw_b + r_stomata_sun) &
                                                + lai_sha(pft) / (raw_b + r_stomata_sha))
               
                rdbp(pft) = min(f_wet(pft) + rdp_dry(pft), (Evpot * rdp_dry(pft) + Wcan(pft) / dt) / Evpot)

                else
                   rdbp(pft) = f_wet(pft)
                end if

            else if (Evpot > 0._dp .and. betaT(pft) <= 1.e-10) then
               rdbp(pft) = min(f_wet(pft),(Evpot * rdp_dry(pft) + Wcan(pft) / dt) / Evpot)
            else
               rdbp(pft) = 1._dp
            end if

            !calculate the potential evapotranspiration for that skin temperature of the ground below the vegetation
            
            !calculate the conductances to latent heat transfer
            caw = 1._dp / raw(pft)
            cgw = Beta_soil / (rawp(pft) + r_litter)        !CLM 4.0 changes this to add in rlitter (see eqn. 5.93)
            cvw = (lai_sno(pft) + stem_sno(pft)) / raw_b * rdbp(pft)

            !----------Prescott
            underpet(pft) = max((ss_soil / (ss_soil + gamma_soil)) * (Swg_v(pft) - Lwg_v(pft)) / lvap, 0._dp)  !(mm s-1)
            !----------
            
            !NOT presently used---
            !canopy specific humidity (CLM 4.0 Eqn 5.90)
            !qs = caw * qatm + cgw * qg + cvw * qsath / (caw + cvw + cgw)  
                                  
            !PET under the canopy                      
            !underpet(pft) = -ratm * (caw * qatm + cvw * qsath - (caw + cvw) * qg) * cgw / (caw + cvw + cgw)  !CLM 5.97
            !---

            Eg_v(pft) = underpet(pft)

           !find the sensible heat conductances
                       !from the canopy air to the atmosphere
                       cah = 1. / rah(pft)    !eq. 5.76

                       !from the ground to the canopy
                       cgh = 1. / rahp(pft)

                       !from the leaf to the canopy
                       cvh = (lai_sno(pft) + stem_sno(pft)) / rb(pft)

              !calculate the sensible heat flux from the ground
              Hg_v(pft) = -ratm * Cair * (cah * theta + cvh * temp_veg(2,pft) - &
                    (cah + cvh) * Tsoil(snl+1)) * (cgh / (cah + cvh + cgh)) !eq. 5.80

              dEgdt_v(pft) = Beta_soil * ratm / (rawp(pft) + r_litter) * (caw + cvw) / (caw + cvw + cgw) * dqgdTg  !CLM 4.0 Eqn 5.111

              dHgdt_v(pft) = ratm * Cair / rahp(pft) * (cah + cvh) / (cah + cvh + cgh)  ! CLM 4.0 Eqn 5.110

              ! this presently uses the dqsath_dT value for the soil, not veg temp.
              dLgdt_v(pft) = lam * ratm * (caw + cgw) * cvw / (caw + cvw + cgw) * dqsath_dT

    else !lai_sno + stem_sno is > 0._dp but less than minplant. Vegetation present but not affecting heat fluxes.Need a red_pet for GPP.
            !assign it the ground value.

            red_pet(pft) = max(0._dp,gpet)  !this translates into demand so should be positive.
            underpet(pft) = 0._dp
            Eg_v(pft) = 0._dp 
            can_evap(pft) = 0._dp
            Hg_v(pft) = 0._dp
            dEgdt_v(pft) = 0._dp
            dHgdt_v(pft) = 0._dp
            dLgdt_v(pft) = 0._dp

    end if
end do !pft loop

end subroutine latentsens

end module latentsensiblefluxes
