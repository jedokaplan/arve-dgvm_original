module pftinitassign

implicit none

public :: pftassign

contains

subroutine pftassign

!called by arve_dgvm. Calculates lm_sapl,sla,rm_sapl,sm_sapl initial values
!coded by Joe Melton 2010.

use arveparams, only : dp,npft,pi,ncvar,OMorgC,k5,k6,frac_cr_tot,reinickerp,&
                        allom2,allom3
use statevars, only : sv,co2
use pftparametersmod, only : prm
use iovariables, only : do_allocation_daily
use pft_state_variables, only : veg

implicit none

!ARGUMENTS
logical, pointer, dimension(:)    :: tree               !true is pft is a tree
real(dp), pointer, dimension(:)   :: sapl_mass          !sapling total structural carbon (gC/ ind)
real(dp), pointer, dimension(:)   :: leaflongevity      !leaf longevity (years)
logical, pointer, dimension(:)    :: c4                 !plants with C4 (true) or C3 (false) photosynthetic pathway
real(dp), pointer, dimension(:)   :: sla                !PFT specific leaf area (m2/gC)
real(dp), pointer, dimension(:,:) :: lm_sapl            !initial (sapling) leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_sapl            !initial (sapling) sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_sapl            !initial (sapling) heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_sapl            !initial (sapling) root mass (gC)
real(dp), pointer, dimension(:)   :: latosa             !leaf area to sapwood area (m2 m-2)
real(dp), pointer, dimension(:)   :: sapl_NSC_frac      !sapling non-structural carbon as a fraction of total sapling mass 
real(dp), pointer, dimension(:,:) :: NSC_sapl           !sapling non-structural carbon store (gC / ind)
real(dp), pointer, dimension(:)   :: allom1             !Packing constrain from Reinicke's rule exponent
real(dp), pointer, dimension(:)   :: wooddens           !density of wood (sap + heart) (g m-3)

!parameters
real(dp), parameter :: theta = 0.85!0.75_dp             !allometric exponent (3/4)

real(dp), parameter, dimension(npft) :: sapl_lai = (/ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.001, 0.001 /)   !sapling initial lai (LPJ)
real(dp), parameter, dimension(npft) :: sapl_mass_lpj = (/ 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0 /)  !sapling initial mass (LPJ)
real(dp), parameter, dimension(npft) :: leafrootratio = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75 /)!non-water stressed leaf to root ratio (LPJ)

!LOCAL VARIABLES
integer  :: pft
real(dp) :: x
real(dp) :: stem_diam
real(dp) :: height_sapl
real(dp) :: lm_temp
real(dp) :: rm_temp
real(dp) :: lmtorm

!point pointers
lm_sapl                 => veg%lm_sapl
sm_sapl                 => veg%sm_sapl
hm_sapl                 => veg%hm_sapl
rm_sapl                 => veg%rm_sapl
sla                     => veg%sla
latosa                  => prm%latosa
tree                    => prm%tree
sapl_mass               => prm%sapl_mass
leaflongevity           => prm%leaflongevity
c4                      => prm%c4
sapl_NSC_frac           => prm%sapl_NSC_frac
NSC_sapl                => veg%NSC_sapl
allom1                  => prm%allom1
wooddens                => prm%wooddens

!--------------------------

!begin calculations

if (do_allocation_daily) then

      do pft = 1,npft

        !Calculate specific leaf area (SLA) for each PFT from leaf longevity
        !Include conversion (OMorgC) from m2/g(dry wt) to m2/gC
        !SLA in m2/gC, leaf_longevity in years
        !This SLA is different from the LPJ one. I took the same figure but used graphclick to get 
        !the true equation of the line from Reich et al. 1997 PNAS Figure 1f. Relation is for the global
        !dataset

        sla(pft) = 1.e-4 / OMorgC * 267.9 * (leaflongevity(pft) * 12._dp)**(-0.39)

           !---
           !Define initial mass structure

           if (tree(pft)) then  !woody PFTs

              !FLAG! SOME of these can be further parameterized to angiosperm vs. conifer

              x = sapl_mass(pft) / 1.e3 / OMorgC  !sapling total structual carbon mass in kg dry matter

              !from eqn 2 Niklas and Spatz PNAS 2004
              lm_temp = 0.137 * x**theta  !lm in kg dry matter

              lm_sapl(pft,1) = lm_temp * 1.e3 * OMorgC  !in g C

              !from Table 1 Enquist and Niklas Science 2002
              rm_temp = (lm_temp / 0.4)**(1./0.79)  !rm in kg dry matter

              rm_sapl(pft,1) = rm_temp * 1.e3 * OMorgC  !in g C

              ! Assign carbon to the sapling NSC stores !FLAG not sure which is best
!              NSC_sapl(pft,1) = lm_sapl(pft,1) * sapl_NSC_frac(pft)
              NSC_sapl(pft,1) = lm_sapl(pft,1) * prm(pft)%NSC_opt  

              !remained is sapwood as saplings/seedlings have no heartwood                   
              sm_sapl(pft,1) = sapl_mass(pft) - lm_sapl(pft,1) - rm_sapl(pft,1)

              hm_sapl(pft,1) = 0._dp  

    ! JM 21.02.2011
    ! Try new stem diam calculation from Niklas and Spatz 2004          
    !          stem_diam = (4.d0 * (sm_sapl(pft,1) + hm_sapl(pft,1)) * lm_sapl(pft,1) * sla(pft) / (pi * sm_sapl(pft,1) * latosa(pft)))**0.5  

              stem_diam = ((sm_sapl(pft,1) + hm_sapl(pft,1)) * 1.e-3 /(202.3 * k5))**(0.375)
    !--

              !from Niklas and Spatz PNAS 2004 eqn 5
              height_sapl = k5 * stem_diam**(2./3.) - k6  !height in m

           else

              !grass no sap or heartwood, remaining structural carbon evenly split between roots and leaves
              lm_sapl(pft,1) = 0.001 / sla(pft)
              rm_sapl(pft,1) = lm_sapl(pft,1)  !FLAG FLAG test JM 01.03.2011
              !lm_sapl(pft,1) = sapl_mass(pft) * 0.5_dp
              !rm_sapl(pft,1) = sapl_mass(pft) * 0.5_dp

              NSC_sapl(pft,1) = lm_sapl(pft,1) * prm(pft)%NSC_opt  

              sm_sapl(pft,1) = 0.
              hm_sapl(pft,1) = 0.

          end if  !tree/grass

    !write(*,'(i5,6f12.4)')pft,lm_sapl(pft,1),rm_sapl(pft,1),sm_sapl(pft,1),NSC_sapl(pft,1),stem_diam,height_sapl

    end do ! pft loop

else  !annual C allocation scheme. Retain original LPJ formulation

    ! Annual allocation does not track NSC    
      NSC_sapl(:,:) = 0._dp

    do pft = 1,npft
    ! Calculate specific leaf area (SLA) for each PFT from leaf longevity
    ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
    ! Equation based on Reich et al 1997, Fig 1f:

    ! SLA in m2/gC, leaf_longevity in years

    sla(pft) = 2.0e-4 * exp(6.15 - 0.46 * log(leaflongevity(pft) * 12.d0))

!    sla(pft) = 1.e-4 / OMorgC * 267.9 * (leaflongevity(pft) * 12._dp)**(-0.39)


    ! Define initial mass structure

    if (tree(pft)) then  !woody PFTs

    !  Calculate leafmass for a sapling individual
    !   (1) lai = leafmass * sla / (crown area)
    !   (2) (leaf area) = latosa * (sapwood xs area)
    !          (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
    !   (3) (crown area) = allom1 * (stem diameter) ** reinickerp
    !          (Reinickes theory)
    !  From (1),
    !   (4) leafmass = lai * (crown area) / sla
    !  From (1) & (3),
    !   (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
    !  From (2),
    !   (6) leafmass = latosa * (sapwood xs area) / sla
    !   (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
    !  From (6) and (7),
    !   (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
    !  From (8),
    !   (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
    !  (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
    !  Define x,
    !  (11) x = [ (sapwood diameter)+(heartwood diameter) ] / 
    !	    (sapwood diameter)
    !  From (10) & (11),
    !  (12) (stem diameter) = x * (sapwood diameter)
    !  From (5), (9) & (12),
    !  (13) leafmass = lai * allom1 * x**reinickerp * 
    !		 (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
    !  From (13),
    !  (14) leafmass = [ lai * allom1 * x**reinickerp *
    !         (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(1-1/reinickerp)

      x = sapl_mass_lpj(pft)

      lm_sapl(pft,1) = (sapl_lai(pft) * allom1(pft) * x**reinickerp * (4.d0 * sla(pft) / pi / latosa(pft))**(reinickerp*0.5) &
                       / sla(pft))**(1.d0 - 1.d0 / reinickerp)  !eqn 14
      lm_sapl(pft,2) = 17.8 - co2(2)   ! initial 13C value from Llyod & Farquhar, 1994
      lm_sapl(pft,3) = 0.d0

    !  Calculate sapling stem diameter
    !  From (9) & (12),
    !  (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

      stem_diam = x * (4.d0 * lm_sapl(pft,1) * sla(pft) / pi / latosa(pft))**0.5  !Eqn 15

    !  Calculate sapling height
    !  (16) height = allom2 * (stem diameter)**allom3            (source?)

      height_sapl = allom2 * stem_diam**allom3   !Eqn 16

    !  Calculate sapling sapwood mass
    !  (17) (sapwood volume) = height * (sapwood xs area)
    !  (18) (sapwood xs area) = leafmass * sla / latosa
    !  From (17) & (18),

    ! (19) (sapwood volume) = height * leafmass * sla / latosa
    !  (20) (sapwood mass) = (wood density) * (sapwood volume)
    !  From (19) & (20),
    !  (21) (sapwood mass) = (wood density) * height * leafmass * sla /
    !         latosa

      sm_sapl(pft,1) = wooddens(pft) * height_sapl * lm_sapl(pft,1) * sla(pft) / latosa(pft)   !Eqn 21
      sm_sapl(pft,2) = lm_sapl(pft,2)   ! 13C value in permille
      sm_sapl(pft,3) = lm_sapl(pft,3)

    !  Calculate sapling heartwood mass
    !  From (11),
    !  (22) (heartwood mass) = (x-1) * (sapwood mass)

      hm_sapl(pft,1) = (x - 1.d0) * sm_sapl(pft,1)  !Eqn 22
      hm_sapl(pft,2) = sm_sapl(pft,2)   ! 13C value in permille
      hm_sapl(pft,3) = sm_sapl(pft,3)

    else ! grass PFT

      lm_sapl(pft,1) = sapl_lai(pft) / sla(pft)

    !  Set initial 13C values for saplings, grass

      if (c4(pft)) then   !C4 plants
        lm_sapl(pft,2) = 3.6 - co2(2)  !from lloyd & farquhar,1994           
        lm_sapl(pft,3) = 0.d0
      else                           !C3 plpants
        lm_sapl(pft,2) = 17.8 - co2(2)  !from lloyd & farquhar,1994          
        lm_sapl(pft,3) = 0.d0
      end if

      sm_sapl(pft,2) = 0.d0             ! no sapwood and hartwood
      hm_sapl(pft,2) = 0.d0             ! for grass PFT
      sm_sapl(pft,3) = 0.d0             ! no sapwood and hartwood
      hm_sapl(pft,3) = 0.d0             ! for grass PFT

    end if

    !Calculate sapling or initial grass rootmass
    !(23) lmtorm = (leafmass) / (rootmass)

    lmtorm = leafrootratio(pft) 
    rm_sapl(pft,1) = (1.d0 / lmtorm) * lm_sapl(pft,1)  !From Eqn 23
    rm_sapl(pft,2) = lm_sapl(pft,2)       ! 13C value in permille
    
    end do !pft

end if

return

end subroutine pftassign

end module pftinitassign
