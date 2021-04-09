subroutine radiativeflux(j,i)

! calculates the heat flux from the atmosphere to the land surface
! also the sunlight and shaded plant absorbed photosynthetically active
! radiation. Based upon CLM 4.0 and other sources as listed
! coded by JM May 9 08

use arveparams,                 only : dp,Tfreeze,npft,band,sigsb,Esnow,Esoil,pi,minplant
use statevars,                  only : sv,dt,dayl
use soilstate_vars,             only : asurf_g,TairK,Swv,Swg_b,Swg_v,Lwv,Lwg_b,Lwg_v,dLgdT_b,&
                                        Latm,Lwvd,Eveg,Egrnd,surf
use metvarsmod,                 only : met,atm
use pft_state_variables,        only : veg

implicit none

! arguments
integer, intent(in) :: j             !index of the current gridcell
integer, intent(in) :: i             !intra-day iteration count

! pointers
logical, pointer                        :: polar                !true if polar night
real(dp), pointer                       :: cldf                 !cloud cover fraction (1 is fully overcast)
integer,  pointer                       :: snl                  !index value of top snow layer (negative is more layers)
real(dp), pointer                       :: rZ0                  !zenith angle at solar noon (radians)
real(dp), pointer                       :: zen                  !cosine of the solar zenith angle (rad)
real(dp), pointer                       :: Tdew                 !air dew point temperature
real(dp), pointer                       :: absrad_annual        !annual total radiation absorbed (J m-2)
real(dp), pointer, dimension(:)         :: sw                   !downwelling shortwave (W m-2)
real(dp), pointer, dimension(:,:)       :: temp_veg             !temperature of the vegetation (K)
real(dp), pointer, dimension(:,:)       :: fabi                 !diffuse flux absorbed by vegetation
real(dp), pointer, dimension(:,:)       :: fabd                 !direct flux absorbed by vegetation
real(dp), pointer, dimension(:,:)       :: fabi_vir             !diffuse flux absorbed by vegetation (virtual leaves)
real(dp), pointer, dimension(:,:)       :: fabd_vir             !direct flux absorbed by vegetation (virtual leaves)
real(dp), pointer, dimension(:,:)       :: ftid                 !downward diffuse fluxes per unit incident direct beam
real(dp), pointer, dimension(:,:)       :: ftii                 !downward diffuse fluxes per unit incident diffuse beam
real(dp), pointer, dimension(:)         :: Tsoil                !soil temperature (K)
real(dp), pointer, dimension(:)         :: Tsoiln               !soil temperature (K) previous time step
real(dp), pointer, dimension(:,:)       :: ftdd                 !down direct flux below veg per unit dir flx
real(dp), pointer, dimension(:)         :: lai_sno              !instantaneous leaf area index per pft(m2 m-2)
real(dp), pointer, dimension(:)         :: lai_vir              !virtual leaf area index per pft(m2 m-2) (for phenology)
real(dp), pointer, dimension(:)         :: lai_sun              !sunlight leaf area index per pft (m2 m-2)
real(dp), pointer, dimension(:)         :: lai_sha              !shaded leaf area index per pft (m2 m-2)
real(dp), pointer, dimension(:)         :: stem_sno             !instantaneous stem area index per pft(m2 m-2)
real(dp), pointer, dimension(:)         :: APAR_sun             !visible only shortwave rad absorbed by vegetation averaged for leaf area (W m-2) sunlight 
real(dp), pointer, dimension(:)         :: APAR_sha             !visible only shortwave rad absorbed by vegetation averaged for leaf area (W m-2) shaded
real(dp), pointer, dimension(:,:)       :: scatcf_tot           !fraction of intercepted radiation that is scattered (0 to 1)
real(dp), pointer, dimension(:,:)       :: scatcf_tot_vir       !fraction of intercepted radiation that is scattered (0 to 1)(virtual leaves)
real(dp), pointer, dimension(:)         :: fsun                 !sunlit fraction of canopy
real(dp), pointer, dimension(:)         :: fsha                 !shaded fraction of canopy
real(dp), pointer, dimension(:)         :: bigK                 !the optical depth of direct beam per unit leaf and stem area
logical, pointer, dimension(:)          :: virtual_leaf         !true if the leaf being calculated is a virtual one to determine phenology.
real(dp), pointer                       :: fsnow                !fraction of the gridcell covered by snow (fraction)

! local variables
integer :: k,pft                                !counter
real(dp), dimension(band,npft) :: Svt           !temporary variable
real(dp), dimension(band) :: Sgt_v,Swg          !SW absorbed by ground below vegetation, and bare ground
real(dp) :: plantLS                             !total LAI + stem indices.
real(dp), dimension(band) :: Swdir              !incident direct beam solar fluxes (W m-2)
real(dp), dimension(band) :: Swdif              !incident diffuse beam solar fluxes (W m-2)
real(dp) :: D                                   !dew point depression (Tdew-TairK)
real(dp) :: Sw_noondir                          !solar noon downwelling shortwave direct (W m-2)
real(dp) :: Sw_noondif                          !solar noon downwelling shortwave diffuse (W m-2)
real(dp) :: Sw_dir_i                            !downwelling shortwave direct (W m-2) for this short timestep
real(dp) :: Sw_dif_i                            !downwelling shortwave diffuse (W m-2) for this short timestep
real(dp) :: phi_dirbeam_dir                     !unscattered direct beam absorbed by the canopy (visible)
real(dp) :: phi_dirbeam_dif                     !scattered direct beam absorbed as diffuse radiation by the canopy (visible)
real(dp) :: phi_difbeam                         !incoming diffuse radiation absorbed by the canopy (visible)

!parameters
real(dp), parameter :: mu_avg = 1._dp            !average inverse optical depth for longwave rad.
real(dp), parameter :: la = 10.77_dp             !coefficients for the correction of downwelling longwave dependent on cld cover
real(dp), parameter :: lb =  2.34_dp             !
real(dp), parameter :: lc =-18.44_dp             !
real(dp), parameter :: visnir_part = 0.5_dp      !partition of visible and near infrared (0.5 and 0.5)

! point pointers
snl             => sv(j)%snl
cldf            => met(j,0)%cldf
sw              => surf%ra_dp
fabi            => veg%fabi
fabd            => veg%fabd
fabi_vir        => veg%fabi_vir
fabd_vir        => veg%fabd_vir
ftid            => veg%ftid
ftii            => veg%ftii
stem_sno        => sv(j)%stem_sno
lai_sno         => sv(j)%lai_sno
lai_vir         => sv(j)%lai_vir
Tsoil           => sv(j)%Tsoil
Tsoiln          => sv(j)%Tsoiln
temp_veg        => sv(j)%temp_veg
polar           => sv(j)%polar
rZ0             => surf%rZ0
zen             => surf%zen
ftdd            => veg%ftdd
Tdew            => atm%Tdew
APAR_sun        => veg%APAR_sun
APAR_sha        => veg%APAR_sha
lai_sun         => veg%lai_sun
lai_sha         => veg%lai_sha
scatcf_tot      => veg%scatcf_tot
scatcf_tot_vir  => veg%scatcf_tot_vir
fsun            => veg%fsun
fsha            => veg%fsha
bigK            => veg%bigK
virtual_leaf    => sv(j)%virtual_leaf
absrad_annual   => surf%absrad_annual
fsnow           => surf%fsnow

!--------

! begin calculations

! Initial default values

    plantLS = 0._dp
    APAR_sun(1:npft) = 0._dp
    APAR_sha(1:npft) = 0._dp
    Swv(1:npft) = 0._dp
    Swg_b = 0._dp
    Swg_v(1:npft) = 0._dp
    Eveg(1:npft) = 0._dp
    Lwv(1:npft) = 0._dp             
    Lwvd(1:npft) = 0._dp            
    Lwg_v(1:npft) = 0._dp           
    
! This routine has been adapted from the CLM routine in parts with other references as noted.

      if (i == 1 .and. (dayl(1) > 0._dp) .or. polar) then !only do SW during day.

       ! Partition the total day downwelling shortwave direct and diffuse components into the amount delivered in this short timestep
       ! From Wang et al. 2002 Ecol. Modelling 00 1-14. JM Oct 3 08. Regressions from Figure 2.

       ! Find the downwelling shortwave at solar noon for both diffuse and direct components
       Sw_noondir = sw(1) * (2.28_dp - 1.1_dp * rZ0 + 0.8_dp * rZ0**2 - 0.23_dp * rZ0**3)
       Sw_noondif = sw(2) * (1.73_dp - 0.81_dp * rZ0 + 0.58_dp * rZ0**2 - 0.16_dp * rZ0**3)

       Sw_dir_i = Sw_noondir * cos((acos(zen) - rZ0) / (0.5_dp * pi - rZ0) * 0.5_dp * pi) * zen / cos(rZ0)  !Eqn 5
       Sw_dif_i = Sw_noondif * cos((acos(zen) - rZ0) / (0.5_dp * pi - rZ0) * 0.5_dp * pi)  !Eqn 7

       ! Split solar flux into vis and near infrared (it is 50% in each)
       Swdir(1) = visnir_part * Sw_dir_i
       Swdif(1) = visnir_part * Sw_dif_i

       Swdir(2) = visnir_part * Sw_dir_i
       Swdif(2) = visnir_part * Sw_dif_i

       !increment the annual total radiation total
       absrad_annual = absrad_annual + (Sw_dir_i + Sw_dif_i) * dt
  
       ! Shortwave radiative flux to the soil and canopy surface

               ! Solar radiation is conserved as the sum of the diffuse and direct shortwave equaling Sv + Sg + the reflected solar
               ! radiation

               do pft = 1,npft
                 
                plantLS = lai_sno(pft) + stem_sno(pft)

                   
                   ! Find the total solar radiation absorbed by vegetation and ground. Virtual leaves will not be considered here.
                   
                   if (plantLS > minplant) then !if there is 'real' vegetation present.

                        ! Total solar radiation absorbed by the vegetation (both direct and diffuse)
                         Svt(1:band,pft) = max(0._dp, (Swdir(1:band) * fabd(1:band,pft) + Swdif(1:band) * fabi(1:band,pft)))   !eq. 4.3
                         Swv(pft) = sum(Svt(1:band,pft))

                        ! Total solar radiation absorbed by the ground (under vegetation)
                         Sgt_v(1:band) = max(0._dp,(Swdir(1:band) * ftdd(1:band,pft) * (1. - asurf_g(1,1:band)) + (Swdir(1:band)&
                                         *  ftid(1:band,pft) + Swdif(1:band) * ftii(1:band,pft)) * (1. - asurf_g(2,1:band)))) !eq. 4.4
                         Swg_v(pft) = sum(Sgt_v(1:band))

                  end if

                  ! Find the total absorbed photosynthetically active radiation for both real and virtual leaves
                  
                  ! Store the radiation absorbed in visible for photosynthesis routine for 'real' leaves
                  if (lai_sno(pft) > 0._dp) then

                            k = 1  !1 is the visible band
                            phi_dirbeam_dir = Swdir(k) * (1._dp - exp(-bigK(pft)*plantLS)) * (1._dp - scatcf_tot(k,pft))!CLM 4.0 Eqn 4.9
                            phi_dirbeam_dif = max(0._dp, (Swdir(k) * fabd(k,pft) - phi_dirbeam_dir))        !CLM 4.0 Eqn 4.10
                            phi_difbeam     = Swdif(k) * fabi(k,pft)                                        !CLM 4.0 Eqn 4.11

                            ! visible only shortwave rad absorbed by vegetation SUNLIGHT leaves averaged for the ground area(W m-2)
                            if (lai_sun(pft) > 0._dp)   APAR_sun(pft) = (phi_dirbeam_dir + phi_dirbeam_dif * fsun(pft) &
                                                           + phi_difbeam * fsun(pft)) &
                                                        * (lai_sno(pft) / plantLS) / lai_sun(pft) !CLM 4.12   
                                                        
                            ! visible only shortwave rad absorbed by vegetation SHADED leaves averaged for the ground area(W m-2)                            
                            if (lai_sha(pft) > 0._dp) APAR_sha(pft) = (phi_dirbeam_dif * fsha(pft) &
                                                            + phi_difbeam * fsha(pft)) &
                                                        * (lai_sno(pft) / plantLS)/ lai_sha(pft) !CLM 4.13
                                                                  
                  ! Find APAR_sun for virtual leaves (for phenology), all else set to 0
                  else if (virtual_leaf(pft)) then

                            k = 1  !1 is the visible band
                            phi_dirbeam_dir = Swdir(k) * (1._dp - exp(-bigK(pft) * (lai_vir(pft) + stem_sno(pft)))) * (1._dp - scatcf_tot_vir(k,pft))!CLM 4.0 Eqn 4.9
                            phi_dirbeam_dif = max(0._dp, (Swdir(k) * fabd_vir(k,pft) - phi_dirbeam_dir))        !CLM 4.0 Eqn 4.10
                            phi_difbeam     = Swdif(k) * fabi_vir(k,pft)                                        !CLM 4.0 Eqn 4.11

                            ! visible only shortwave rad absorbed by vegetation SUNLIGHT leaves averaged for the ground area(W m-2)
                            if (lai_vir(pft) > 0._dp) APAR_sun(pft) = (phi_dirbeam_dir + phi_dirbeam_dif * fsun(pft) &
                                                           + phi_difbeam * fsun(pft)) *&
                                                             (lai_vir(pft)/ (lai_vir(pft) + stem_sno(pft)))&
                                                             / lai_vir(pft)! CLM 4.0 Eqn 4.12                  
                            
                  end if  !lai if loop
                  
                end do  !pft loop

                ! total solar radiation absorbed by the ground (bare ground, no canopy above)
                Swg(1:band) = max(0._dp, (Swdir(1:band) * (1. - asurf_g(1,1:band)) + Swdif(1:band) * (1. - asurf_g(2,1:band))))  !eq. 4.5
                Swg_b = sum(Swg(1:band))

      end if

!--------
! Longwave fluxes

! Downwelling atmospheric longwave flux

              ! the following equations for downwelling longwave come from
              ! Josey et. al 2003, JGR Atmospheres 108, doi:10.1029/2002JC001418  eqn. 14

              D = Tdew - TairK
              Latm = sigsb * (TairK + la * cldf**2 + lb * cldf + lc + 0.84_dp * (D + 4.01_dp))**4

              ! calc ground emissivity
              Egrnd = (Esnow * fsnow) + (Esoil * (1. - fsnow))


! the following calculations are required each iteration

  ! Upwelling longwave fluxes
        do pft = 1,npft
        
            plantLS = lai_sno(pft) + stem_sno(pft)

            if (plantLS > minplant) then !if there is vegetation present
                  
                ! calc veg emissivity
                Eveg(pft) = 1._dp - exp(-(plantLS) / mu_avg)  !CLM eqn 4.25

                !find the downward longwave below vegetation
                Lwvd(pft) = ((1._dp - Eveg(pft)) * Latm) + Eveg(pft) * sigsb * (temp_veg(1,pft)**4) &
                                + 4._dp * Eveg(pft) * sigsb * (temp_veg(1,pft)**3) * (temp_veg(2,pft) - temp_veg(1,pft))   !4.17 (CLM 4.0 Eqn 4.21)

                !net longwave flux for vegetation
                Lwv(pft) = (2._dp - Eveg(pft) * (1._dp - Egrnd)) * Eveg(pft) * sigsb * (temp_veg(2,pft)**4) &
                                - Eveg(pft) * Egrnd * sigsb * (Tsoiln(snl+1)**4) &
                                - Eveg(pft) * (1._dp + (1._dp - Egrnd) * (1._dp - Eveg(pft))) * Latm   !eq. 4.19 (CLM 4.0 Eqn 4.23)

                !find the NET longwave flux for ground (positive towards atmosphere) (CLM 4.0 Eqn 4.22)
                Lwg_v(pft) = Egrnd * sigsb * Tsoiln(snl+1)**4 -    &                                            !upwelling longwave
                             Egrnd * Lwvd(pft) -  &                                                             !downwelling below veg
                             4._dp * Egrnd * sigsb * Tsoiln(snl+1)**3 * (Tsoil(snl+1) - Tsoiln(snl+1))         !rate of change

          end if
        end do
        
           !Longwave fluxes on bare ground
             Lwg_b = Egrnd * sigsb * Tsoiln(snl+1)**4 -    &                                    !upwelling longwave
                     Egrnd * Latm + &                                                                 !downwelling
                     4._dp * Egrnd * sigsb * Tsoiln(snl+1)**3 * (Tsoil(snl+1) - Tsoiln(snl+1))         !rate of change

             dLgdT_b = -4._dp * Egrnd * sigsb * Tsoiln(snl+1)**3

           !the net longwave flux for ground (both vegetated and not) is calculated in surfaceheatflux

end subroutine radiativeflux
