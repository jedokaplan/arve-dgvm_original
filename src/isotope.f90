!SUBROUTINE ISOTOPE
!This subroutine is for calculating the total fractionation of 13C
!as it goes from free air (as 13CO2) to fixed carbon in the leaf.
!Based upon the model used by Lloyd and Farquhar (1994).
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code
!also adapted to be timestep indepedent

subroutine isotope(pft,j,i)

use arveparams,only : dp,npft,ncvar,Tfreeze,Rgas,minplant
use statevars, only : sv,co2,dayl
use pftparametersmod, only : prm
use pft_state_variables, only : veg
use metvarsmod, only : met

implicit none

!PARAMETERS
real(dp), parameter :: a  = 4.4    !fractionation against 13CO2 during diffusion in free air (per mil) (orig. Craig 1953)
real(dp), parameter :: es = 1.1    !equilibration fractionation as CO2 enters solution (per mil)(Mook et al 1974)
real(dp), parameter :: a1 = 0.7    !fractionation against 13CO2 during diffusion in water (per mil) (O'Leary 1984)
real(dp), parameter :: b  = 27.5   !discrimination against 13CO2 during photosynthetic fixation of CO2
real(dp), parameter :: e  = 0._dp   !fractionation associated with respiration.
real(dp), parameter :: f  = 8._dp   !value used in Lloyd & Farquhar (from Rooney 1988). fractionation associated with
                                   !photorespiration
real(dp), parameter :: b3 = 30._dp  !13C discrimation during CO2 fixation by Rubisco (29per mil) with respect to dissolved
                                   !CO2 also allow for isotope effect of diffusion against 13CO2 in water (1.1 permil)
!ARGUMENTS
integer, intent(in)               :: pft        !pft
integer, intent(in)               :: i                !timestep
integer, intent(in)               :: j                !gridcell

real(dp), pointer, dimension(:,:) :: annual_gpp
real(dp), pointer, dimension(:,:)   :: dgpp
real(dp), pointer, dimension(:,:) :: dresp
real(dp), pointer                  :: temp        !temperature for that timestep (deg C)
real(dp), pointer                  :: leaftemp   !timestep average veg temp (K)
real(dp), pointer                 :: stem_sno   !leaf area index  (m2 m-2)
real(dp), pointer                 :: lai_sno    !stem area index  (m2 m-2)
real(dp), pointer, dimension(:)   :: gamma      !CO2 compensation point in the absence of mitochondrial respiration
real(dp), pointer, dimension(:)   :: CO2ratio   !ratio of the mole fraction of CO2 in the substomatal cavity (intercellular mole fraction of
                                                       !CO2) to the ambient mole fraction of CO2 in the atmosphere surrounding the leaf
logical, pointer                  :: c4                !true if pft is c4

!LOCAL
real(dp) :: deltaa      !total discrimination against 13CO2 during CO2 assimilation by plants
real(dp) :: b4          !discrimination by enzyme PEP-c during photosynthesis carboxylation
real(dp) :: catm        !ambient mole fraction of CO2 in the air surrounding the leaf
real(dp) :: k
real(dp) :: rdi
real(dp) :: q
real(dp) :: r
real(dp) :: s
real(dp) :: t
real(dp) :: phi         !ratio of the rate of CO2 leakage from the bundle sheath cells to the rate of PEP carboxylation
real(dp) :: rdtemp      !temporary rd
real(dp) :: leafT       !local leaf temperature (K)

!point pointers
CO2ratio   => veg%ratio
annual_gpp => veg%annual_gpp
c4         => prm(pft)%c4
dgpp       => veg%dgpp
dresp      => veg%dresp
gamma      => veg%CO2comp
temp       => met(j,0)%temp(i)
leaftemp   => sv(j)%Tveg_ave(2,pft)
stem_sno   => sv(j)%stem_sno(pft)
lai_sno    => sv(j)%lai_sno(pft)

!--------------------------------
!begin calculations

  if (dgpp(pft,1) > 0._dp) then

    !if the vegetation amount is enough to enact vegheatflux, use the vegetation temperature calculated therein
    !other wise take the timestep average temperature and cool it by 5% (leaves cooler due to latent heat flux also less heat
    !capacity)
            if ((lai_sno + stem_sno) > minplant) then
                leafT = leaftemp
            else
                leafT = temp !set to air temperature
            end if

     if (CO2ratio(pft) < 0.05_dp) CO2ratio(pft) = 0.05_dp

     !define fractionation parameters

      if (dresp(pft,1) <= 0._dp) then
         rdtemp = 0.01_dp
      else
         rdtemp = dresp(pft,1)
      end if

     !gamma is the CO2 compensation point (umol mol-1). Calculated in photosynth (CO2comp), we now use that version -JM
     !gamma = 1.54 * (leafT - Tfreeze) !old way.

     rdi = rdtemp / dayl(1)   !NOTE: This is not setup quite correctly but since e (fractination associated with respiration) is taken to
                                !be 0 this is okay to keep as present. IF this value is revised then check original paper and setup
                                !accordingly. JM - June 13 08
     k = rdi / 11._dp

     catm = co2(1) !umol mol-1

     b4 = (26.19_dp - (9483._dp / leafT)) * 1.e-3 !from Henderson et al 1992, Lloyd & Farquhar 1994 eqn 2b

     if (.not. c4) then

        !calculate the fractionation -> C3 pathway
        q = a * (1._dp - CO2ratio(pft) + 0.025_dp)
        r = 0.075_dp * (es + a1)
        s = b * (CO2ratio(pft) - 0.1_dp)
        t = ((e * rdi / k) + f * gamma(pft)) / catm
        deltaa = q + r + s - t
        dgpp(pft,2) = deltaa - co2(2)

     else

        !phi = 0.2   !suggested by Lloyd & Farquhar (1994)

        !NEW WAY: calculate PHI for C4 pathway to account for temperature dependence.
        !to replace old subroutine (calcphi.f90), linear regression
         !based upon Fig 4 (table 1) in Henderson et al 1992 -JM
        phi = (-0.0045 * (leafT - Tfreeze)) + 0.3319
          phi = max(0.18,real(phi))
           phi = min(0.3,real(phi))

        !calculate fractionation -> C4 pathway
        deltaa = a * (1._dp - CO2ratio(pft) + 0.0125_dp) + 0.0375_dp * (es + a1) + (b4 + (b3 - es &
                   - a1) * phi) * (CO2ratio(pft) - 0.05_dp)
        dgpp(pft,2) = deltaa - co2(2)

     end if

  else
     dgpp(pft,2) = 0._dp
  end if

return

end subroutine isotope
