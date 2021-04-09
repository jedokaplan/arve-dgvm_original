module increment

! This module handles all the incrementing of faster timestep values
! to values of the timestep above.

use arveparams, only : dp,npft,ncvar,climbuf,dayspy,nl,daysec,sp
use statevars, only : sv,dt,counter_lim,ov,dtime,doy
use pft_state_variables, only : veg
use soilstate_vars, only : fsnow_mo,surf
use metvarsmod,       only : met,dm,atm
use pftparametersmod, only : prm
use iovariables,        only : do_allocation_daily

implicit none

public :: fast_increment
public :: daily_increment
public :: annual_increment

contains

!----------------------------

subroutine fast_increment(j)

! Increments fast time step values to daily ones.

! arguments
integer, intent(in) :: j

! local variables
integer :: pft

! pointers
logical, pointer, dimension(:) :: present       !true if pft present
real(dp), pointer, dimension(:) :: daet                       !daily actual evapotranspiration (mm)
real(dp), pointer, dimension(:,:) :: dgpp       !daily gross primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:) :: aet                        !daily actual evapotranspiration (mm s-1)
real(dp), pointer, dimension(:,:) :: sgpp       !gross primary productivity (gC m-2 short timestep-1)
real(dp), pointer, dimension(:) :: dpet         !daily potential evapotranspiration (mm s-1)
real(dp), pointer, dimension(:) :: pet          !potential evapotranspiration (mm s-1)
real(dp), pointer, dimension(:) :: supply       !plant water supply (mm/s)
real(dp), pointer, dimension(:) :: demand       !plant water demand (mm/s)
real(dp), pointer, dimension(:) :: supply_ts    !timestep water supply (mm)
real(dp), pointer, dimension(:) :: demand_ts    !time step water demand (mm)
real(dp), pointer, dimension(:) :: Tsoil        !soil temperature (K)
real(dp), pointer, dimension(:) :: Tsoil_ave    !average day/night soil temperature (K)
real(dp), pointer, dimension(:) :: underpet     !potential evapotranspiration for ground beneath canopy (mm s-1)
real(dp), pointer :: gpet                       !bare ground evapotranspiration (mm s-1)
real(dp), pointer, dimension(:) :: can_evap     !canopy evaporation (mm s-1)
real(dp), pointer               :: lastpet      !potential evapotranspiration from last day (mm/timestep)
real(dp), pointer               :: qseva_tot    !timestep soil evaporation flux (mm)
real(dp), pointer               :: qsubl_tot    !timestep soil sublimation flux (mm)                     
real(dp), pointer               :: qsdew_tot    !timestep soil surface dew (mm)                          
real(dp), pointer               :: qfrost_tot   !timestep soil surface frost (mm)         
real(dp), pointer               :: qseva        !soil evaporation flux (mm s-1)                 
real(dp), pointer               :: qsubl        !soil sublimation flux (mm s-1)                 
real(dp), pointer               :: qsdew        !soil surface dew (mm s-1)                      
real(dp), pointer               :: qfrost       !soil surface frost (mm s-1) 
                   

! point pointers
present         => sv(j)%present
pet             => veg%pet
dpet            => veg%dpet
gpet            => veg%gpet
underpet        => veg%underpet
lastpet         => sv(j)%lastpet
supply          => veg%supply
demand          => veg%demand
supply_ts       => veg%supply_ts
demand_ts       => veg%demand_ts
can_evap        => veg%can_evap
Tsoil           => sv(j)%Tsoil
Tsoil_ave       => sv(j)%Tsoil_ave
daet            => veg%daet
aet             => veg%aet    
dgpp            => veg%dgpp
sgpp            => veg%sgpp
qseva_tot       => surf%qseva_tot
qsubl_tot       => surf%qsubl_tot
qsdew_tot       => surf%qsdew_tot
qfrost_tot      => surf%qfrost_tot
qseva           => surf%qseva
qsubl           => surf%qsubl
qsdew           => surf%qsdew
qfrost          => surf%qfrost

!----------

! begin calcs

do pft = 1,npft
  if (present(pft)) then
  
        ! increment the daily gpp (gC m-2 d-1)
        dgpp(pft,1) = dgpp(pft,1) + sgpp(pft,1)

        ! increment daily aet (this stores mm, it is converted to mm s-1 in calcshortstep)
        daet(pft) = daet(pft) + aet(pft) * dt

        ! increment daily demand and supply (in mm)
        supply_ts(pft) = supply_ts(pft) + supply(pft) * dt
        demand_ts(pft) = demand_ts(pft) + demand(pft) * dt

        ! increment the daily PET per pft (in mm)
        !dpet(pft) = dpet(pft) + pet(pft) * dt

  end if
end do

        !increment gpet (in mm)
        !last_gpet(1) = last_gpet(1) + gpet * dt

        !timestep average soil temperature.        
        Tsoil_ave(1:nl) = Tsoil_ave(1:nl) + Tsoil(1:nl)
        
           ! increment the lastpet (only for day values) (in mm/day)
       !if (i == 1) then  !FLAG is this correct?? JM 16.04.2011
        lastpet = lastpet + dt * (gpet + sum(pet(:)) + sum(underpet(:))) !JM 07.03.2011 removed fpc_grid, already done in gpp
       !end if
       
       qseva_tot = qseva_tot + qseva * dt
       qsubl_tot = qsubl_tot + qsubl * dt
       qsdew_tot = qsdew_tot + qsdew * dt
       qfrost_tot = qfrost_tot + qfrost * dt

        ! reset short time step values to 0
        pet(:) = 0._dp
        aet(:) = 0._dp
        sgpp(:,:) = 0._dp
                     
end subroutine fast_increment

!---------------------------------------------------------------------------------------------------
subroutine daily_increment(j,i)

! Increments daily values to annual ones.

! arguments
integer, intent(in) :: j
integer, intent(in) :: i

! local variables
integer :: pft
real(dp) :: smp_inst            !soil matric potential in the root zone per pft, instant value.

! pointers
logical, pointer, dimension(:)          :: present      !true if pft present
real(dp), pointer, dimension(:,:)       :: annual_npp    !annual net primary productivity (g C m-2 yr-1)  
real(dp), pointer, dimension(:,:)       :: annual_resp   !                                                
real(dp), pointer, dimension(:,:)       :: annualgpp    !cumulative net photosynthesis                   
real(dp), pointer, dimension(:)         :: aaet         !annual actual evapotranspiration (mm/yr)
integer,  pointer, dimension(:)         :: pstate       !instantaneous phenological state flag (0=no leaves,
                                                        !1=full leaf out,2=leafing out,3=senescing)
real(dp), pointer, dimension(:)         :: wscal        !water stress factor
real(dp), pointer, dimension(:)         :: aleafdays    !days with leaf cover
real(dp), pointer, dimension(:)         :: awscal       !water stress factor
real(dp), pointer, dimension(:)         :: arh          !annual heterotrophic respiration (gC m-2 yr-1)
real(dp), pointer, dimension(:)         :: drh          !timestep heterotrophic respiration (gC/m2 timestep-1)
real(dp), pointer, dimension(:)         :: daet         !daily actual evapotranspiration (mm)
real(dp), pointer, dimension(:,:)       :: dnpp         !net primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: dgpp         !gross primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: dresp        !maintenance respiration (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: gresp        !growth respiration (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: Tveg_ave     !average vegetation temperature (K)
real(dp), pointer, dimension(:)         :: supply_ts    !water supply to plants  timestep average
real(dp), pointer, dimension(:)         :: demand_ts    !plant water demand  timestep average
real(dp), pointer, dimension(:)         :: fpc_grid     !instantaneous gridcell foliar projective cover (FPC) for gridcell (m2 m-2)
real(dp), pointer, dimension(:)         :: Tsoil_ave    !
real(dp), pointer, dimension(:)         :: dturnover_ind !daily total turnover of living biomass per individual (gC)
real(dp), pointer, dimension(:)         :: turnover_ind  !
real(dp), pointer, dimension(:)         :: reprod       !
real(dp), pointer, dimension(:)         :: dreprod      !
real(dp), pointer                       :: fsnow        !fraction of the gridcell covered by snow (fraction)
real(dp), pointer, dimension(:)         :: gdd_bioclim  !for use in bioclimatic
real(dp), pointer, dimension(:)         :: betaT        !soil moisture function limiting transpiration (fraction)
real(dp), pointer, dimension(:)         :: dsw_ra       !
real(dp), pointer, dimension(:,:)       :: smp_ra       !soil matric potential in the root zone, running average over nts timesteps
real(dp), pointer, dimension(:)         :: minT_ra      !
real(dp), pointer                       :: dsw_b        !total surface downwelling shortwave (kJ m-2 day-1)
real(dp), pointer, dimension(:)         :: gddbase      !base temperature for growing degree days
real(dp), pointer                       :: tmin         !24 hour minimum temperature (C)
real(dp), pointer, dimension(:)         :: gddsum       !accumulated GDD sum for this PFT
real(dp), pointer                       :: tmean        !mean temperature of 24 hour period (C)
real(dp), pointer                       :: zsno         !snow depth (m)
real(dp), pointer                       :: zsno_mo      !monthly cumulative snow depth (m) 
real(dp), pointer                       :: mo_runoff_tot !monthly total overland runoff (mm)
real(dp), pointer, dimension(:)         :: fsnow_ra     !fractional snow cover running values array
real(dp), pointer, dimension(:)         :: Psi          !soil water potential at layer midpoint (mm)
real(dp), pointer, dimension(:,:)       :: rootfracl    !root fraction in that soil layer
real(dp), pointer                       :: qseva_tot    !timestep soil evaporation flux (mm)
real(dp), pointer                       :: qsubl_tot    !timestep soil sublimation flux (mm)             
real(dp), pointer                       :: qsdew_tot    !timestep soil surface dew (mm)                  
real(dp), pointer                       :: qfrost_tot   !timestep soil surface frost (mm)   
real(dp), pointer                       :: qover        !liquid water surface runoff (kg m-2 sec-1)   
integer, pointer  :: gnl                                !index value of lowest 'soil' layer

! point pointers
present         => sv(j)%present
annual_npp      => veg%annual_npp
annual_resp     => veg%annual_resp
annualgpp       => veg%annual_gpp
aaet            => veg%aaet
pstate          => sv(j)%pstate
wscal           => veg%wscal
aleafdays       => sv(j)%aleafdays
awscal          => veg%awscal
arh             => veg%arh
drh             => veg%drh
Tveg_ave        => sv(j)%Tveg_ave
supply_ts       => veg%supply_ts
demand_ts       => veg%demand_ts
fpc_grid        => veg%fpc_grid
daet            => veg%daet
Tsoil_ave       => sv(j)%Tsoil_ave
turnover_ind    => veg%turnover_ind
reprod          => veg%reprod
gddbase         => prm%gddbase
dresp           => veg%dresp   
dsw_b           => met(j,0)%dswb        !NOTE: this is purposefully the kJ/m2/day value
gdd_bioclim     => veg%gdd_bioclim
betaT           => sv(j)%betaT  
smp_ra          => sv(j)%smp_ra   
dsw_ra          => sv(j)%dsw_ra     
minT_ra         => sv(j)%minT_ra    
tmin            => met(j,0)%tmin    
gddsum          => sv(j)%gddsum     
tmean           => met(j,0)%temp24  
zsno            => sv(j)%zsno
dreprod         => veg%dreprod
gresp           => veg%gresp
dgpp            => veg%dgpp
dnpp            => veg%dnpp
dturnover_ind   => veg%dturnover_ind
zsno_mo         => surf%zsno_mo
mo_runoff_tot   => surf%mo_runoff_tot
fsnow           => surf%fsnow
fsnow_ra        => sv(j)%fsnow_ra
Psi             => sv(j)%Psi
rootfracl       => sv(j)%rootfracl
gnl             => sv(j)%gnl
qseva_tot       => surf%qseva_tot
qsubl_tot       => surf%qsubl_tot
qsdew_tot       => surf%qsdew_tot
qfrost_tot      => surf%qfrost_tot
qover           => surf%qover
     
!--------

! begin calcs

if (.not. do_allocation_daily .and. i == 1) then

 ! dsw (this is an averaged value to prevent one really cloudy day from prematurely causing senescence)
 ! accumulates KJ/m2/day
 dsw_ra(:) = eoshift(dsw_ra(:),1,dsw_b)

 ! minimum temperature of the day
 minT_ra(:) = eoshift(minT_ra(:),1,tmin)
 
 ! fractional snow cover at day end
 fsnow_ra(:) = eoshift(fsnow_ra(:),1,fsnow)

end if

do pft = 1,npft

   ! accumulate the gdd_bioclim for this pft (only once per day). daily C accounting scheme does this in phenology.
   if (i == 1 .and. .not. do_allocation_daily) gdd_bioclim(pft) = gdd_bioclim(pft) + max(tmean - gddbase(pft),0._dp)

  if (present(pft)) then  
     
    if (i == 1) then  !daytime increments
        
             ! increment today's values to the running averages arrays. These are to prevent a 
             ! single erratic day to cause undue influence increment the gddsum if the day is 
             ! warm enough. Since this is the daily mean, only do increment during day.
             gddsum(pft) = gddsum(pft) + max(tmean - gddbase(pft),0._dp)

          ! find the water scalar            
          if (pstate(pft) /= 0) then
            
             if (demand_ts(pft) > 0._dp) then
                wscal(pft) = min(1._dp,supply_ts(pft)/demand_ts(pft))
             else ! no demand so wscal is 1
                wscal(pft) = 1._dp
             end if
             
             smp_inst = sum(Psi(1:gnl) * rootfracl(1:gnl,pft))
                          
             ! increment water scalar on days with leaf cover (used in annual C allocation scheme)
             awscal(pft) = awscal(pft) + wscal(pft)

             ! increment count of days with leaf cover
             aleafdays(pft) = aleafdays(pft) + 1._dp

          else if (pstate(pft) == 0) then
        
               wscal(pft) = min(1._dp,betaT(pft))
 
               smp_inst = sum(Psi(1:gnl) * rootfracl(1:gnl,pft))

          end if  !pstate
          
         ! water scalar running avg increment.
           smp_ra(pft,:) = eoshift(smp_ra(pft,:),1,smp_inst)
        
        else !night time increments
        
        ! increment annual aet (mm)
        aaet(pft) = aaet(pft) + daet(pft)
        
        !increment the annual gpp (gC m-2 yr-1)
        annualgpp(pft,1) = annualgpp(pft,1) + dgpp(pft,1)
        annualgpp(pft,2:ncvar) = annualgpp(pft,2:ncvar) + dgpp(pft,2:ncvar) * dgpp(pft,1)

        !increment annual npp (gC m-2 yr-1)
        annual_npp(pft,1) = annual_npp(pft,1) + dnpp(pft,1)

        !increment annual respiration (includes both auto and growth)
        annual_resp(pft,1) = annual_resp(pft,1) + dresp(pft,1) + gresp(pft,1)

        !balance for isotopes of C (running flux-weighted total of isotope ratio)
        annual_resp(pft,2:ncvar) = annual_resp(pft,2:ncvar) + dresp(pft,2:ncvar) * (dresp(pft,1) + gresp(pft,1))
        if (dnpp(pft,1) > 0) annual_npp(pft,2:ncvar) = annual_npp(pft,2:ncvar) + dnpp(pft,2:ncvar) * dnpp(pft,1)

        if (do_allocation_daily) then

          !increment reproduction cost
          reprod(pft) = reprod(pft) + dreprod(pft)    

          ! increment the annual turnover (gC yr-1)
          turnover_ind(pft) = turnover_ind(pft) + dturnover_ind(pft)      

        end if  

        !reset the daily values to 0
        dgpp(pft,1) = 0._dp
        dnpp(pft,1) = 0._dp
        dresp(pft,1) = 0._dp
        gresp(pft,1) = 0._dp
        dturnover_ind(pft) = 0._dp
        dreprod = 0._dp
        wscal(pft) = 0._dp
           
        end if
            
    end if  !present
end do  !pft loop

       ! shift the average timestep vegetation temperature
       Tveg_ave(:,:) = eoshift(Tveg_ave(:,:),1,0._dp)
       
       ! increment the heterotrophic respiration values
       arh(1) = arh(1) + drh(1)
       arh(2:ncvar) = arh(2:ncvar) + drh(2:ncvar) * drh(1)
       
       ! increment the monthly runoff (mm)
       mo_runoff_tot = mo_runoff_tot + qover * dtime

       !increment fractional coverage by snow and snow depth
       fsnow_mo = fsnow_mo + fsnow * dtime/daysec
       zsno_mo = zsno_mo + zsno * dtime/daysec       

       !some resets
       daet(1:npft) = 0._dp
       demand_ts(1:npft) = 0._dp
       supply_ts(1:npft) = 0._dp
       qseva_tot = 0._dp
       qsubl_tot = 0._dp
       qsdew_tot = 0._dp
       qfrost_tot = 0._dp

end subroutine daily_increment

!-----------------------------------------------------------------------------------

subroutine annual_increment(j)

implicit none

! arguments
integer, intent(in) :: j   !gridcell index

! local variables
integer :: k,pft
integer :: yrs

! pointers
logical, pointer, dimension(:)          :: present              !true if pft is present in gridcell
real(dp), pointer, dimension(:)         :: aleafdays            !count of days with foliage
real(dp), pointer, dimension(:)         :: awscal               !annual water stress factor
real(dp), pointer, dimension(:)         :: nind                 !indiv/m2
real(dp), pointer, dimension(:,:)       :: annual_npp            !annual net primary productivity (g C m-2 yr-1)
real(dp), pointer, dimension(:,:)       :: annual_resp           !
real(dp), pointer, dimension(:,:)       :: bm_inc               !annual biomass increment (gC/m2) (minus reproduction costs)
real(dp), pointer, dimension(:)         :: arh                  !annual heterotrophic respiration (gC m-2 yr-1)
real(dp), pointer, dimension(:,:)       :: annual_gpp                 !annual gross primary productivity (g C m-2 yr-1)
integer, pointer, dimension(:)          :: yrpres               !number of years the pft was present
real(sp), pointer                       :: max_ald              !maximal active layer depth for the year (m)
real(dp), pointer, dimension(:)         :: ald                  !active layer depth (m)
real(dp), pointer, dimension(:)         :: reprod               !summed reproductive cost (gC/m2)
real(dp), pointer, dimension(:)         :: lm_tot_max           !
real(dp), pointer, dimension(:,:)       :: lm_max               !
real(dp), pointer                       :: MAT                  !mean annual temperature (C)
real(dp), pointer                       :: MAT_d                !stores running total of MAT (C)

! point pointers
present         => sv(j)%present
aleafdays       => sv(j)%aleafdays
awscal          => veg%awscal
bm_inc          => veg%bm_inc
annual_npp      => veg%annual_npp
annual_resp     => veg%annual_resp
arh             => veg%arh
annual_gpp      => veg%annual_gpp
yrpres          => sv(j)%yrpres
nind            => sv(j)%nind
max_ald         => ov(j)%max_ald
ald             => surf%ald
reprod          => veg%reprod
lm_tot_max      => veg%lm_tot_max
lm_max          => sv(j)%lm_max
MAT             => atm%MAT
MAT_d           => atm%MAT_d

!-----------------------------

!begin calcs

  ! Copy NPP to bm_inc, which will store running total of C production remaining for growth
  if (do_allocation_daily) then
  
   ! remove the amount allocated to repoduction (calc'd in daily_allocation)
    do k = 1,ncvar
      bm_inc(1:npft,k) = annual_npp(1:npft,k) - reprod(1:npft) 
    end do
 
   ! reset the max leaf mass to the maximal leaf mass reached this year
   lm_max(1:npft,1) = lm_tot_max(1:npft)

   lm_tot_max(1:npft) = 0._dp
 
  else !annual C allocation scheme
  
      bm_inc(1:npft,1:ncvar) = annual_npp(1:npft,1:ncvar)

  end if
   
   do pft = 1,npft

     if (present(pft)) then !if the PFT is present in the grid cell

        if (aleafdays(pft) > 0._dp) then
                awscal(pft) = min(1._dp,awscal(pft) / aleafdays(pft))
        else
                awscal(pft) = 0._dp
        end if

         ! increment the year present counter for this pft
         yrpres(pft) = yrpres(pft) + 1

         if (yrpres(pft) < climbuf) then  !averaging period is climbuf
           yrs = yrpres(pft)
         else
           yrs = climbuf
         end if

     end if

     ! make adjustments to isotope pools
      if (annual_resp(pft,1) > 1._dp) then !JM 01.03.2011
        annual_resp(pft,2:ncvar) = annual_resp(pft,2:ncvar) / annual_resp(pft,1)
      else
        annual_resp(pft,2:ncvar) = 0._dp
      end if

      if (annual_npp(pft,1)  > 1._dp) then  !JM 01.03.2011
         annual_npp(pft,2:ncvar) = annual_npp(pft,2:ncvar) / annual_npp(pft,1)
       else
         annual_npp(pft,2:ncvar) = 0._dp
       end if

      if (annual_gpp(pft,1) > 1._dp) then  !JM 01.03.2011
         annual_gpp(pft,2:ncvar) = annual_gpp(pft,2:ncvar) / annual_gpp(pft,1)
       else
         annual_gpp(pft,2:ncvar) = 0._dp
       end if

  end do

  !adjust annual heterotrophic respiration isotope pools
  if (arh(1) > 1._dp) then  !JM 01.03.2011
   arh(2:ncvar) = arh(2:ncvar) / arh(1)
  end if

  !find the mean annual temperature
  MAT = MAT_d / dayspy

  MAT_d = 0._dp
  
  !find the active layer depth (deepest depth of melt in the year)
  max_ald = real(maxval(ald(:)))
  
end subroutine annual_increment

end module increment
