module canopywater
!determines the partitioning of precipitation into canopy intercepted, canopy drip,
!and canopy throughfall. Also determines the wetted fraction of the canopy. Used in the demand calculatation
!Coded May 7 08 JM from CLM tech note ch.7.1  also CLM 3.5 Documentation

implicit none

public :: canopy_interception
public :: canopy_wet

contains

subroutine canopy_interception(j,pft)
!determines the partitioning of precipitation into canopy intercepted, canopy drip,
!and canopy throughfall.

use arveparams,only : dp
use statevars, only : sv,dtime
use soilstate_vars, only : surf
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j !gridcell index
integer, intent(in) :: pft !plant functional type

!pointers
real(dp), pointer, dimension(:) :: lai_sno
real(dp), pointer, dimension(:) :: stem_sno
real(dp), pointer, dimension(:) :: Wcan_old  !canopy water from previous timestep (mm)
real(dp), pointer, dimension(:) :: Wcanmax   !maximum quantity of water the canopy can hold (mm)
real(dp), pointer, dimension(:) :: Wcan      !canopy water (mm)
real(dp), pointer, dimension(:) :: f_wet     !fraction of the canopy that is wet
real(dp), pointer, dimension(:) :: f_dry     !fraction of the canopy that is dry
real(dp), pointer               :: qliq      !liquid precipitation (kg m-2 sec-1)
real(dp), pointer               :: qsno      !snowfall (kg m-2 sec -1)
real(dp), pointer, dimension(:) :: qgrnd_l   !total rate of liquid precip reaching the ground under canopy(mm s-1)
real(dp), pointer, dimension(:) :: qgrnd_s   !total rate of solid precip reaching the ground under canopy(mm s-1)

!parameters
real(dp), parameter :: pmax = 0.1       !max storage of water in canopy (mm)
real(dp), parameter :: alpha = 0.25     !factor scales the interception from point to grid cell.

!variables
real(dp) :: qintr   !water intercepted by the canopy (mm s-1)
real(dp) :: qthru_l !liquid water throughfall the canopy (mm s-1)
real(dp) :: qthru_s !snow throughfall the canopy (mm s-1)
real(dp) :: qdrip_l !canopy drip (liquid) (mm s-1)
real(dp) :: qdrip_s !canopy drip (snow) (mm s-1)
real(dp) :: Wcan_i  !canopy water after accounting for inteception (mm)
real(dp) :: fintr   !fraction intercepted

!point pointers
lai_sno          => sv(j)%lai_sno
stem_sno         => sv(j)%stem_sno
Wcan_old         => sv(j)%Wcan_old
Wcanmax          => veg%Wcanmax
Wcan             => veg%Wcan
f_wet            => veg%f_wet
f_dry            => veg%f_dry
qliq             => surf%qliq
qsno             => surf%qsno
qgrnd_l          => surf%qgrnd_l
qgrnd_s          => surf%qgrnd_s

!--------------------

!calculations begin

if (lai_sno(pft) + stem_sno(pft) > 0._dp .and. (qliq + qsno > 0._dp .or. Wcan_old(pft) > 0._dp)) then  !7.12

    !find maximum quantity of water the canopy can hold
    Wcanmax(pft) = pmax * (lai_sno(pft) + stem_sno(pft))  !7.8

    !precipitation intercepted by the canopy
    fintr = alpha * (1. - exp(-0.5 * (lai_sno(pft) + stem_sno(pft))))  !eqn A1 CLM 3.5 Documentation
    qintr = (qliq + qsno) * fintr

    !canopy water after accounting for interception
    Wcan_i = max((Wcan_old(pft) + qintr * dtime), 0._dp)   !7.7

    !precipitation that falls through the canopy
    qthru_l = qliq * (1._dp - fintr)   !7.3
    qthru_s = qsno  * (1._dp - fintr)   !7.4

    !canopy drip
    if (Wcan_i > Wcanmax(pft) .and. (qliq + qsno) > 0._dp) then

            qdrip_l = max(((Wcan_i - Wcanmax(pft)) / dtime * qliq / (qliq + qsno)), 0._dp)  !7.5
            qdrip_s = max(((Wcan_i - Wcanmax(pft)) / dtime * qsno / (qliq + qsno)), 0._dp)  !7.6
    else
              qdrip_l = 0._dp
              qdrip_s = 0._dp
    end if

    !find the total rate of precipitation reaching the ground
    qgrnd_l(pft) = qthru_l + qdrip_l  !7.10 add qgrnd_l to the soil water

    qgrnd_s(pft) = qthru_s + qdrip_s  !7.11 add qgrnd_s  to the snow pack.

    !update the canopy water (but not including dew or evaporation!!)
    Wcan(pft) = max((Wcan_old(pft) + dtime * (qintr - qdrip_l - qdrip_s)), 0._dp)  !7.9

    !find wetted fraction of the canopy for surface albedo calculations
    !this does not account for dew/evap over the course of the timestep.
    !that happens in canopy_wet

          f_wet(pft) = min(1._dp,(Wcan(pft) / Wcanmax(pft))**(2./3.))

          f_dry(pft) = (1. - f_wet(pft)) * lai_sno(pft) / (lai_sno(pft) + stem_sno(pft))  !7.13

    ! move Wcan_old to the new Wcan value
    Wcan_old(pft) = Wcan(pft)

else
  Wcan(pft) = 0._dp
  f_wet(pft) = 0._dp
  
    if (lai_sno(pft) > 0.d0) then !able to transpire
    f_dry(pft) = 1._dp 
    else
    f_dry(pft) = 0._dp   
    end if
    
  qgrnd_l(pft) = 0._dp
  qgrnd_s(pft) = 0._dp
end if

end subroutine canopy_interception

!-----

subroutine canopy_wet(j,pft)
!Update the water in the canopy after vegheatflux subroutine. Also update the wetted fraction of the canopy

use arveparams,only : dp
use statevars, only : sv,dt
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j                !gridcell index
integer, intent(in) :: pft              !plant functional type

!pointers
real(dp), pointer, dimension(:) :: lai_sno                  !instantaneous leaf area index (m2 m-2)
real(dp), pointer, dimension(:) :: stem_sno                 !instantaneous stem area index (m2 m-2)
real(dp), pointer, dimension(:) :: Wcan_old                 !canopy water from previous timestep (mm)          
real(dp), pointer, dimension(:) :: Wcan                     !canopy water (mm)                                 
real(dp), pointer, dimension(:) :: f_wet                    !fraction of the canopy that is wet                
real(dp), pointer, dimension(:) :: f_dry                    !fraction of the canopy that is dry                
real(dp), pointer, dimension(:) :: Wcanmax                  !maximum quantity of water the canopy can hold (mm )
real(dp), pointer, dimension(:) :: can_evap                 !canopy evaporation (mm s-1)                       
real(dp), pointer, dimension(:) :: xs_can_drip              !excess canopy water (mm)

!point pointers
lai_sno         => sv(j)%lai_sno
stem_sno        => sv(j)%stem_sno
Wcan_old        => sv(j)%Wcan_old
Wcan            => veg%Wcan
f_wet           => veg%f_wet
f_dry           => veg%f_dry
Wcanmax         => veg%Wcanmax
can_evap        => veg%can_evap
xs_can_drip     => veg%xs_can_drip

!--------------------

!calculations begin
xs_can_drip(pft) = 0._dp

!update the canopy water for evaporation/condensation over the course of the timestep.
Wcan(pft) = max((Wcan(pft)  - can_evap(pft) * dt), 0._dp)  !7.9

!Wcan can not become greater than Wcanmax. Assume that extra water passes to ground
!NOTE: this assumes it is water and not snow.
  if (Wcan(pft) > Wcanmax(pft)) then
        xs_can_drip(pft) = Wcan(pft) - Wcanmax(pft)
        Wcan(pft) = min(Wcan(pft),Wcanmax(pft))
  end if

!find wetted fraction of the canopy for plant water demand

      f_wet(pft) = min(1._dp,(Wcan(pft) / Wcanmax(pft))**(2./3.))    !7.12

        ! for stability, if the f_wet is less than 0.01 then call it dry.
        if (f_wet(pft) < 0.01_dp) f_wet(pft) = 0.d0

      f_dry(pft) = (1. - f_wet(pft)) * lai_sno(pft) / (lai_sno(pft) + stem_sno(pft))  !7.13

!set up next timesteps canopy water
Wcan_old(pft) = Wcan(pft)

end subroutine canopy_wet

end module canopywater
