module calcshortstep

implicit none

public :: shorttimestep

contains

!-----------------------------------------------------------------------
subroutine shorttimestep(yr,j,i)

! All calculations that happen on the shorter timestep are called from this subroutine
! the length of the short timestep is set in arveparams (shortstep)

use pft_state_variables, only : veg
use statevars,  only : sv,dt,counter,counter_lim,doy,dtime,dayl
use soiltemperaturemod, only : surfaceheatflux,soiltemperature
use diurnaltempmod, only : instantemp
use soilstate_vars,  only : TairK,Tair,eatm,emb,ratm,theta,lam
use arveparams,  only : dp,Tfreeze,rh,Rda,Cair,npft,shortstep,lsub,lvap
use esatdesdT, only : esat,desdT
use gppmod, only : gpp
use leafarea, only : sla_twoleaf
use coszen, only : diszen
use increment, only : fast_increment
use canopywater, only : canopy_wet
use metvarsmod, only : atm

implicit none

!arguments
integer, intent(in)  :: yr                              !year (not strictly required)
integer, intent(in)  :: j                               !index of the current gridcell
integer, intent(in)  :: i                               !intra-day iteration count

!pointers
logical, pointer :: polar                               !true if polar night
integer,  pointer :: snl                                !index value of top snow layer (negative is more layers)
real(dp), pointer :: Patm                               !atmospheric pressure (Pa)
real(dp), pointer :: Patm30                             !Pressure at reference height (30m above ground) (Pa)
real(dp), pointer, dimension(:) :: lai_sno              !instantaneous leaf area index (m2 m-2)
logical, pointer, dimension(:) :: present               !true if pft is present in grid cell
real(dp), pointer :: Tdew                               !air dew point temperature (K)
real(dp), pointer, dimension(:,:) :: Tveg_ave           !timestep average vegetation temperature (K)
real(dp), pointer, dimension(:) :: supply_ts            !water supply to plants  timestep average (mm s-1)
real(dp), pointer, dimension(:) :: demand_ts            !plant water demand  timestep average (mm s-1)
real(dp), pointer, dimension(:) :: daet                 !timestep actual evapotranspiration (mm s-1)
real(dp), pointer, dimension(:) :: Wcan                 !water in canopy (mm)
real(dp), pointer, dimension(:) :: f_wet                !fraction of canopy that is wet
real(dp), pointer, dimension(:) :: f_dry                !fraction of canopy that is dry
real(dp), pointer, dimension(:) :: Wliq                 !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice                 !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: stem_sno             !instantaneous stem area index (m2 m-2)
real(dp), pointer, dimension(:) :: can_evap             !canopy evaporation (mm s-1)
real(dp), pointer, dimension(:,:) :: temp_veg           !vegetation temperature (K)
real(dp), pointer, dimension(:) :: Tsoil                !soil temperature (K)

!local variables
integer  :: pft                                         !counters
real(dp) :: tempVnew                                    !new veg temperature (K)

!point pointers
Patm            => sv(j)%Patm
snl             => sv(j)%snl
lai_sno         => sv(j)%lai_sno
present         => sv(j)%present
polar           => sv(j)%polar
Tdew            => atm%Tdew
Patm30          => sv(j)%Patm30
Tveg_ave        => sv(j)%Tveg_ave
supply_ts       => veg%supply_ts
demand_ts       => veg%demand_ts
daet            => veg%daet
f_wet           => veg%f_wet
f_dry           => veg%f_dry
Wcan            => veg%Wcan
Wliq            => sv(j)%Wliq
Wice            => sv(j)%Wice
stem_sno        => sv(j)%stem_sno
can_evap        => veg%can_evap
temp_veg        => sv(j)%temp_veg
Tsoil           => sv(j)%Tsoil
!----------

! sets the timestep to shortstep for the following surfacefluxloop
  counter = 1
  counter_lim = int(dtime / shortstep)
  dt = dtime / real(counter_lim)
  
! Short time step loop begins
surfacefluxloop : do        

            ! exit if it has completed the necessary number of short timesteps for the daily timestep
            if (counter == counter_lim + 1)  exit

    ! Find the instaneous temperature for this 'shortstep' timestep. (calls instantemp from diurnaltempmod)
      tair = instantemp(i,j,counter)
      TairK = tair + Tfreeze    ! change the air temp into degrees Kelvin

    ! Set vegetation temperature. NOTE: Presently we do not calculate the vegetation temperature
    ! assign all veg the same temp. -> average of soil and air temp.
    tempVnew = (TairK + Tsoil(snl+1)) / 2._dp
    temp_veg(:,:) = eoshift(temp_veg(:,:),1,tempVnew)

    ! State of the atmosphere calculations

            ! Atmospheric vapor pressure (Pa)
              eatm = rh / 100._dp * esat(TairK)  !CLM p.162

            ! Atmospheric vapor pressure (millibar)
              emb = eatm * 0.01_dp

            ! Dewpoint temperature
              Tdew = 34.07_dp + 4157._dp / log(2.1718e8 / emb)

            ! Density of moist air (kg m-3)
              ratm = (Patm - 0.378_dp * eatm) / (Rda * TairK)  !eqn 5.8

            ! Atmospheric potential temperature (K)
              theta = TairK * ((Patm/Patm30)**(Rda / Cair))        !eq. 12.1 p. 161

    ! Determine lam for the conditions (done once per short timestep and stored in soilstate_vars)
            if (Wliq(snl+1) == 0._dp .and. Wice(snl+1) > 0._dp) then
              lam = lsub
            else
              lam = lvap
            end if

          if (i == 1 .and. (dayl(1) > 0._dp) .or. polar) then !only do these routines during day.

                     ! Find the cosine of the zenith angle for this short time-step
                        call diszen(j)

                    ! Presently the albedo is calculated on a daily timestep. This is to speed up the model.
                    ! However, it is possible to do this on the short timestep. To do so, you must uncomment this
                    ! subroutine below. Change the mu in albedo to zen from avg_cZ. JM 02.08.2010.
                       call albedo(j)

           end if

                    ! Find the SLA, LAI_sun, LAI_sha for this timestep
                     do pft = 1,npft
                      if (present(pft)) then

                        call sla_twoleaf(j,pft)

                      end if
                     end do

    ! Call the hydraulic water supply to determine the amount of water that is available for this timestep    
       ! NOTE! This is NOT used! This hydraulic routine is designed for use in a cohort mode similar to LPJ-GUESS
       ! this subroutine does not work with avg individual models. Mostly a result of the strange tree shapes that
       ! avg individual models have (tall skinny trees have trouble getting water into the leaves). Water supply to 
       ! trees in ARVE is presently done like LPJ-DGVM where the wetness of the soil layers forms the constraint. If the 
       ! allometric relations are changed, we could try this again, otherwise keep commented out. JM 25.10.10
              !  do pft = 1,npft
              !    if (lai_sno(pft) > 0._dp .and. present(pft)) then
              !      call hydraulic(j,pft)
              !
              !  end do

    ! Find the resistances for sensible and latent heat transfer
            call resistance(j)

    ! Calculate the radiation flux (Swv,Swg,Lwg,Lwv,Latm) for this ground/veg temperature
            call radiativeflux(j,i)

    ! Calculate the heat flux from the atmosphere to the land surface 
            call surfaceheatflux(j)
                     !calls
                     !-> latentsens

    ! Call GPP with the red_pet (=demand) for this timestep to do photosynthesis
            call gpp(j,i)
                !contains
                !->optgpp
                !->photosynthesis
                !->Vcmax
                
                   
    ! Update the amount of water in the canopy for evaporative loss
          do pft = 1,npft
            if (can_evap(pft) /= 0._dp .and. lai_sno(pft) + stem_sno(pft) > 0._dp) then  !only update if there was a change in canopy water.

              call canopy_wet(j,pft)

            else if (lai_sno(pft) + stem_sno(pft) == 0._dp) then  !reset values if no stem/leaf area

              f_wet(pft) = 0._dp
              f_dry(pft) = 1._dp 
              Wcan(pft)  = 0._dp

            end if
          end do

    ! Call the soil temperature routine (part of short timestep loop)
              call soiltemperature(j)
                    !calls:
                    !-> tridiag.f90
                    !-> soilphasechg.f90

    ! Increment the short timestep values to daily
            call fast_increment(j)

    ! Increment the counter and loop again
    counter = counter + 1

end do surfacefluxloop

        ! Create some day/night average values
          Tveg_ave(2,:) = Tveg_ave(2,:) / real(counter_lim)   !average timestep vegetation temperature (K)

          demand_ts(:) = demand_ts(:) / dtime
          supply_ts(:) = supply_ts(:) / dtime
          daet(:) = daet(:) / dtime

end subroutine shorttimestep

end module calcshortstep
