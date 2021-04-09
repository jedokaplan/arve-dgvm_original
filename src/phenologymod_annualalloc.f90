module phenologymod_annualalloc

! Annual C allocation scheme phenology subroutine
! Instantaneous calculation of phenological state for summergreen and drought-deciduous plants
! routine is based upon several sources (listed in subroutines). JM Oct 20 08

use arveparams,                 only : dp,daysec,nts,Tfreeze,psimax
use statevars,                  only : sv,doy,dayl,dtime
use pftparametersmod,           only : prm
use metvarsmod,                 only : met,dm

implicit none

public  :: phenology_annualalloc
private :: summerphen
private :: rainphen
private :: grassphen
private :: everphen

!-------------Module level declarations

!input pointers
logical,  pointer               :: tree         !true if pft is a tree
integer,  pointer               :: leafout      !status of leaf out (1) able to leaf out or presently leafed out, (2) leaf shed due to cold/daylight
                                                ! (3) leaf shed due to moisture stress, (4) leaf shed due to leaf longevity
integer,  pointer               :: phentype     !phenological type: evergreen (1), summergreen (2), raingreen (3), any type (4)
real(dp), pointer               :: gddsum       !accumulated GDD sum for this PFT
real(dp), pointer               :: phenramp     !summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
real(dp), pointer               :: longevity    !leaf longevity (years)
real(dp), pointer               :: drydrop      !minimum water scalar before leaves are dropped
real(dp), pointer, dimension(:) :: dsw_ra       !arrays where the values for the running averages are stored (dsw)
real(dp), pointer, dimension(:,:) :: smp_ra     !soil matric potential in the root zone running average array
real(dp), pointer, dimension(:) :: minT_ra      !(min air temp)
real(dp), pointer, dimension(:) :: fsnow_ra     !fractional snow cover running values array
real(dp), pointer               :: chillday     !winter chilling days (# of days with daily temp <= 5C) reset in Aug (N. Hemi)
real(dp), pointer               :: lat          !latitude of vegetation
real(dp), pointer               :: tmean        !mean temperature of 24 hour period (C)
real(dp), pointer               :: betaT        !soil moisture function limiting transpiration (fraction)
real(dp), pointer               :: gddleaf      !growing degree-day value at leaf out initiation
real(dp), pointer               :: p_wet_m      !precipitation of the wettest month
real(dp), pointer               :: p_dry_m      !precipitation of the driest month
real(dp), pointer               :: p_mo_avg     !mean annual monthly precipitation
integer,  pointer               :: pstate       !instantaneous phenological state flag (0=no leaves, 1=full leaf out,2=leafing out,3=senescing)
real(dp), pointer               :: phdays       !number of days the PFT has been in the given phenological state
real(dp), pointer               :: phfrac       !instantaneous phenological state fraction (0=no leaves, 1=full leaf out)

!local variables
real(dp) :: smp_av      !average of nts for wscal
real(dp) :: dsw_av      !average of nts for dsw
real(dp) :: minT_av     !average of nts for min temp.
real(dp) :: fsnow_av    !average of nts for fsnow
real(dp) :: ttr         !thermal time requirement (degree days with temperature >5C) before leafout can occur
real(dp) :: dsw_slope   !the direction of solar insolation (increasing/decreasing cloud cover proxy) FLAG if clouds are
                        !fixed they can be used instead JM OCt 30 08
integer  :: dper        !nts / 2 = half of the number of timesteps to take a running dsw average

!parameters
real(dp), parameter :: dorm = 20._dp     !days of imposed dormancy after leaf lose due to leaf age being more than leaf longevity
real(dp), parameter :: dsen =   22._dp   !number of days for senescence for summergreen
!real(dp), parameter :: lmin = 10000._dp  !low light intensity at which senesence is triggered (KJ m-2 d-1)
real(dp), parameter :: alpha = 56.51_dp  !coefficients for the ttr calculation
real(dp), parameter :: beta = 677.24_dp
real(dp), parameter :: gamma = -0.0232_dp
real(dp), parameter :: rdsen =   20._dp   !number of days for senescence for raingreen/herbaceous
real(dp), parameter :: minairT = -2._dp   !minimum air temperature cut-off for grass (deg C)

contains

!---------------

subroutine phenology_annualalloc(j,pft)

implicit none

! arguments
integer, intent (in) :: j               !gridcell index
integer, intent (in) :: pft             !plant functional type

! point pointers
pstate                    => sv(j)%pstate(pft)
gddsum                    => sv(j)%gddsum(pft)
lat                       => sv(j)%lat
leafout                   => sv(j)%leafout(pft)
tmean                     => met(j,0)%temp24
tree                      => prm(pft)%tree
phentype                  => prm(pft)%phentype
longevity                 => prm(pft)%leaflongevity
betaT                     => sv(j)%betaT(pft)
gddleaf                   => sv(j)%gddleaf(pft)

! variables only found in this annual allocation version of phenology.
phdays                    => sv(j)%phdays(pft)
smp_ra                    => sv(j)%smp_ra
dsw_ra                    => sv(j)%dsw_ra
minT_ra                   => sv(j)%minT_ra
fsnow_ra                  => sv(j)%fsnow_ra
phenramp                  => prm(pft)%phenramp
drydrop                   => prm(pft)%drydrop
chillday                  => sv(j)%chilld
phfrac                    => sv(j)%phfrac(pft)

!--------------

! begin calculations

        ! find some running averages
        dsw_av = sum(dsw_ra(:)) / real(nts)
        smp_av = sum(smp_ra(pft,:)) / real(nts)
        minT_av = sum(minT_ra(:)) / real(nts)
        fsnow_av = sum(fsnow_ra(:)) / real(nts)
        
   !select type of phenology and calculate
   select case (phentype)

   case(1)  !evergreen phenology
     call everphen(j,pft)

   case(2) !summergreen phenology
     call summerphen()

   case(3) !raingreen phenology
     call rainphen()

   case(4) !herbaceous PFTs
     call grassphen() 
     
   end select

end subroutine phenology_annualalloc

!--------------------------------------------------------------------------------------------------
subroutine everphen(j,pft)

use metvarsmod, only : dm

implicit none

! local variables
integer, intent(in) :: j
integer, intent(in) :: pft
real(dp) ::  enter_dry_season         !precipitation amount for entry of dry season (mm)
real(dp) ::  enter_wet_season         !precipitation amount for entry of wet season (mm)

! parameters
real(dp), parameter :: dsen =   60._dp  !number of days for leaf abscission for evergreen
real(dp), parameter :: BB   =  300._dp  !number of gdd's (5 degree base, should be changed) before budburst can begin

! pointers
real(dp), pointer :: prec_ra

! point pointer
prec_ra => dm(j)%prec_ra
p_wet_m => dm(j)%p_wet_m
p_dry_m => dm(j)%p_dry_m
p_mo_avg=> dm(j)%p_mo_avg

!-------------------
! calculations begin

! The tropical evergreen routine below is not tested and should not be used.
! pft /= 0 prevents it from being used.
   
       pstate = 1
       phfrac = 1._dp
       phdays = phdays + 1._dp

    return  ! do not use the code below...
 
    if (pft /= 1) then  !tropical evergreen

    !Tropical rainforest vegetation exhibits seasonal swings in LAI (Myneni et al. PNAS 2007)
    !Evergreens (only in tropics! PFT 1) should account for the leaf flushing with dry season
    !and leaf abcission with the wet season, total swings ~25%). The logic is that with the decreased
    !light intensity due to cloudy conditions in the wet season, the plant can minimize its metabolic
    !costs by abcissing some leaves and flushing when light intensity is higher.

      dper = nint(0.5 * real(nts))

      !find the slope of the dsw for the nts period.
      dsw_slope = sum(dsw_ra(1:dper)) / sum(dsw_ra(dper+1:nts))
      enter_dry_season = (p_mo_avg + p_dry_m) / 2._dp
      enter_wet_season = (p_mo_avg + p_wet_m) / 2._dp

      select case (pstate)

      case(0)  !for this evergreen, 0 is reduced leaf cover. Leaf abscission occurs in rainy season to avoid high
              !respiration costs.
       if (prec_ra <= enter_dry_season .and. dsw_slope <= 1._dp) then
          pstate = 2 !leafing out
          phdays = 0._dp
          gddsum = 0._dp
       else
          phdays = phdays + 1._dp
       end if

      case(1)  !full leaf
        if (prec_ra > enter_wet_season.and. dsw_slope > 1._dp) then   !move into wet season, loss of ~25% of leaves to avoid the high respiration costs
          pstate = 3 !senescing
          phdays = 0._dp
        else
          phdays = phdays + 1._dp
        end if

      case(2)  !leafing out
        phfrac = min((gddsum / (phenramp + BB) + 0.7_dp), 1._dp)
        if (phfrac == 1._dp) then
          pstate = 1 !full leaf out
          phdays = 0._dp
        else
          phdays = phdays + 1._dp
        end if

      case(3)  !senescing
        phfrac = max(1._dp - (phdays / (2._dp * dsen)) * 0.3_dp, 0.7_dp)
        if (phfrac == 0.75_dp) then
          pstate = 0  !reduced leaf cover is achieved.
          phdays = 0._dp
          gddsum = 0._dp
        else
          phdays = phdays + 1._dp
        end if

      end select
     end if
     
end subroutine everphen

!-----------------------

subroutine summerphen()

! In Summerphen, leaves are lost or gained on the basis of temperature or a combination of temperature
! and decreasing daylight. Also possible is leaf longevity

implicit none

!--------------
!calculations begin

! set conditions for leaf out.

if (leafout == 2) then ! if the leaves were originally lost due to cold/daylight apply conditions below

   ! The chilling day requirement is reset on Aug 1st (N. Hemi) and the ttr is reset January 1 (N. Hemi)

   if (lat >= 0._dp) then !N. Hemi
      if (doy == 1) then    ! January 1st
           gddsum = 0._dp
           leafout = 1
      else if (doy == 213) then ! August 1st
           chillday = 0._dp
      end if
   else ! S. Hemi
      if (doy == 182) then ! July 1st
           gddsum = 0._dp
           leafout = 1
      else if (doy == 32) then ! February 1st
           chillday = 0._dp
      end if
   end if

else if (leafout == 4) then  ! if leaves were lost due to leaf longevity, impose dormancy.
        if (phdays > dorm) leafout = 1
end if

phdays = phdays + 1._dp
      
select case (pstate)

case(0)  !no leaves

     ! for a summerphen to be able to put on leaves again, it must fulfill a chilling requirement
     ! this prevents the plant from reacting to erratic warm periods in fall. Also as the chilling requirement is
     ! fulfilled this reduces the amount of gdd that is required for the plant to green up in the spring
     ! from Zhang et al. GRL 2007 doi:10.1029/2007GL031447 and Zhang et al. GCB 2004 doi:10.1111/j.1365-2486.200400784.x

     if (tmean <= 5._dp) then  ! if the day's temperature is less than 5 C then add it as a chillday
       chillday = chillday + 1._dp
     end if

     ! find the thermal time requirement (degree days with temperature >5C)
     ttr = alpha + beta * exp(gamma * chillday)

     ! If: 1) the growing degree days has satisfied ttr, 2) it is not below freezing,
     ! 3)not barren of water (soil water is available to roots; betaT > 0), and 4) it is a new season (leafout = 1):
     ! then leafout!
     
     if (gddsum >= ttr .and. leafout == 1 .and. betaT > 0._dp) then
       pstate = 2   !leafing out
       phdays = 0._dp
       gddleaf = gddsum
     end if

case(1)  !full leaf

    ! leaves are full unless they are lost due to 1) cold temps (less than 0C), 2) low insolation
    ! (~photoperiod with cool temps (<11.15C)) AND the daylength is getting shorter OR leaf longevity
    ! is too long

    if ((minT_av < -2._dp .or. (dayl(1) < 39300 .and. minT_av <= 11.15_dp)) .and. dayl(1) >= dayl(2)) then !cold
      pstate = 3 !senescing
      phdays = 0._dp
      leafout = 2
    else if (phdays >= (365.0 * longevity)) then !old leaves
      pstate = 3 !senescing
      phdays = 0._dp
      leafout = 4
    end if

case(2)  !leafing out

    ! the speed of leafing out is determined by the gddsum. Presently we have gddsum as only the number of degrees above
    ! gddbase for each day. NOTE: This does not work well for areas with a long period of cool spring temperatures. We could better
    ! use the Baskerville-Emin method. However, our present method is retained as all our pft parameter values were determined using the old
    ! technique. So... we retain the averaging method for now. JM Oct 20 08

    ! NOTE: in some warmer locations, the chilling day requirement is low so the ttr is high, and can be higher than the
    ! phenramp. To prevent instantaneous leafout, an additional 50 degree days is required above ttr for leaf out to occur in
    ! in instances like this. JM

      if (gddsum > phenramp) then
        phfrac = min(gddsum / (gddleaf + 50._dp), 1._dp)
      else
        phfrac = min(gddsum / phenramp, 1._dp)
      end if

    if (phfrac == 1._dp) then
      pstate = 1 !full leaf out
      phdays = 0._dp
    end if

case(3)  !senescing

    phfrac = max(1._dp - phdays / dsen, 0._dp)

    if (phfrac <= 0._dp) then
      pstate = 0  !no leaves (done senescence)
      phdays = 0._dp
    end if

end select

end subroutine summerphen

!--------------------------------------------------------------------------------------------------
subroutine rainphen()
!Drought phenology and net phenology for today. Drought deciduous PFTs shed their leaves when their
!water scalar falls below their PFT specific minimum value (drydrop). Leaves are replaced with a phenological ramp
!once the 10% + minimum water scalar is exceeded.

implicit none

!--------------
!calculations begin

!set conditions for leaf out.

!if the leaves were originally lost due to leaf longevity, impose
if (leafout == 4 .and. pstate == 0 .and. phdays > dorm) then
    leafout = 1
else if (leafout == 3) then !no imposed dormancy, able to leaf out again once smp_av is high enough.
    leafout = 1
end if

phdays = phdays + 1._dp
      
select case (pstate)

case(0)  !no leaves
      
    if (smp_av > psimax .and. leafout == 1) then  !takes more water to start leafing out than to lose leaves 
      pstate = 2 !leafing out
      phdays = 0._dp
    end if

case(1)  !full leaf

    if (smp_av <= drydrop) then   !water below water threshold
      phdays = 0._dp
      pstate = 3 !senescing
      leafout = 3
    else if  (phdays >= (365.0 * longevity)) then  !leaves have reached end of life cycle
      phdays = 0._dp
      pstate = 3 !senescing
      leafout = 4
    end if

case(2)  !leafing out

    phfrac = max(min(phdays / rdsen, 1._dp), 0._dp)

    if (phfrac == 1._dp) then
      pstate = 1 !full leaf out
      phdays = 0._dp
    end if

case(3)  !senescing

    phfrac = max(1._dp - phdays / rdsen, 0._dp)

    if (phfrac <= 0._dp) then
      pstate = 0  !no leaves (done senescence)
      phdays = 0._dp
      gddsum = 0._dp
    end if

end select

end subroutine rainphen

!--------------------------------------------------------------------------------------------------
subroutine grassphen()

implicit none

!--------------
!calculations begin

phdays = phdays + 1._dp

select case (pstate)

case(0)  !no leaves 

      ! if it is not covered with snow or barren of water, then leaf out
      if (smp_av > psimax .and. fsnow_av <= 0.1_dp .and. minT_av > minairT) then
              !    write(*,'(a,es12.4,f12.4,i5,2f12.4)')'start leaf out',smp_av,fsnow_av,doy,phdays,sv(1)%zsno
        pstate = 2 !leafing out
        phdays = 0._dp
        gddleaf = gddsum
      end if

case(1)  !full leaf

      !moisture stress causes leaf drop or full burial
      if (smp_av <= drydrop .or. fsnow_av >= 0.5 .or. minT_av < minairT) then  
                !  write(*,'(a,es12.4,f12.4,i5,3f12.4)')'start drop',smp_av,fsnow_av,doy,phdays,sv(1)%zsno,minT_av 
          leafout = 3
          pstate = 3 !senescing
          phdays = 0._dp
      end if

case(2)  !leafing out

       !in some warmer locations, the chilling day requirement is low so the ttr is high, and can be higher than the
       !phenramp. To prevent instantaneous leafout, an additional 50 degree days is required above ttr for leaf out to occur.

      if (gddsum > phenramp - 50.) then
        phfrac = min(gddsum / (gddleaf + 50._dp), 1._dp)
      else
        phfrac = min(gddsum / phenramp, 1._dp)
      end if

      ! to avoid becoming fully leafed for a quick event, we can return to a no leaf state
      if ((smp_av < drydrop .or. fsnow_av >= 0.9  .or. minT_av < minairT) .and. phdays < 7._dp) then
    !write(*,'(a,es12.4,f12.4,i5,2f12.4)')'return to sens',smp_av,fsnow_av,doy,phdays,minT_av
        pstate = 3 !senescing
        phdays = 0._dp
      end if

      if (phfrac == 1._dp) then  !leaf fraction is one
    !write(*,'(a,es12.4,f12.4,i5,2f12.4)')'full leaves',smp_av,fsnow_av,doy,phdays,sv(1)%zsno
        pstate = 1 !full leaf out
        phdays = 0._dp
      end if

case(3)  !senescing
        
      !if during a sensecence, conditions improve take advantage and leaf out.
      if (smp_av > psimax .and. fsnow_av < 0.5 .and. minT_av > minairT .and. phdays < 7._dp) then
     ! write(*,'(a,es12.4,f12.4,2i5,f12.4)')'return to leaf',smp_av,fsnow_av,doy,nint(phdays),minT_av
        pstate = 2 !leafing out
        phdays = 0._dp
        gddleaf = gddsum
      end if

      phfrac = max(1._dp - phdays / rdsen, 0._dp)

      if (phfrac <= 0._dp) then !no leaves (done senescence)
        !         write(*,'(a,es12.4,f12.4,i5,2f12.4)')'no leaves',smp_av,fsnow_av,doy,phdays,sv(1)%zsno
        pstate = 0  
        phdays = 0._dp
      end if

end select

end subroutine grassphen

end module phenologymod_annualalloc
