subroutine calcdaily(j,yr)

! Daily calculations in arve-dgvm
! expects the meteorological variables in the data structure met (module metvarsmod)
! to be set for this day and gridcell (performed in the driver that calls calcdaily)

use arveparams,         only : dp,npft,Tfreeze,nl,shortstep,minplant,ns,hr1,hr23
use metvarsmod,         only : met
use statevars,          only : sv,counter,counter_lim,dt,gridded,doy,dtime,dayl,ov
use soilhydrologymod,   only : newsnowflux,soilwaterflux,snowdynamics
use soiltemperaturemod, only : soilthermalprop,surfaceheatflux
use phenologymod,       only : phenology
use phenologymod_annualalloc, only: phenology_annualalloc
use soilstate_vars,     only : surf,hs  !just to print out to screen
use surfrad,            only : simplesurfrad,surface_rad
use increment,          only : daily_increment
use leafarea,           only : leafburial
use calcshortstep,      only : shorttimestep
use iovariables,        only : do_allocation_daily
use calcannual,         only : roots
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j                        !gridcell index
integer, intent(in) :: yr                       !NOTE brought in for testing. Not really needed JM 28.04.2010

!pointers
real(dp), pointer :: prec                       !precipitation (mm)  !just for printing out
real(dp), pointer :: avg_cZ                     !cosine of the zenith angle (rad) of the incident beam
logical, pointer, dimension(:) :: present       !true if pft present
integer,  pointer :: snl                        !index value of top snow layer (negative is more layers)
logical, pointer  :: polar                      !true if polar day
real(dp), pointer :: maxdl                      !longest day length of the year for this cell (s)
real(dp), pointer, dimension(:) :: Tsoil_ave    !avg. soil temp across night/day (K)

!local variables
real(dp), dimension(2) :: tstep                  !timestep temporary variable
integer :: i                                     !day/night index counter
integer :: pft                                   !plant functional type counter

!point pointers
snl             => sv(j)%snl
prec            => met(j,0)%prec
present         => sv(j)%present
polar           => sv(j)%polar
avg_cZ          => surf%avg_cZ
maxdl           => sv(j)%maxdl
Tsoil_ave       => sv(j)%Tsoil_ave

!---------------------
! Begin calculations

! Calculate the timestep for day/night loop (limited to 1 hr minimum)
        if (dayl(1) > hr1) then
          if (dayl(1) < hr23) then
            tstep(1) = dayl(1)
            tstep(2) = 86400._dp - dayl(1)  ! normal day
            polar = .false.
          else
            tstep(1) = hr23                ! polar day
            tstep(2) = hr1
            polar = .true.
          end if
        else
          tstep(1) = hr1                   ! polar night
          tstep(2) = hr23
          polar = .false.
        end if

! DAYTIME/NIGHTTIME loop.
! This loop must be performed sequentially - no parallelization!

do i = 1,2  !1=day, 2=night

                 !if (yr > 1) then
                 !write(*,'(a,i5,2f10.4,l5,2i5,f10.4)')'calcdaily:',j,sv(j)%lon,sv(j)%lat,sv(j)%peat,doy,i,sv(j)%zsno
                 !if (yr > 1 .and. doy > 110) then
                 ! read(*,*)
                 !end if
                 !end if
                
  dtime = tstep(i)

  ! Record the longest day of the year for this gridcell
    if (yr == 1) maxdl = max(dtime,maxdl)

  ! Find shortwave radiation -> Different for gridded versus point simulations!
    if (gridded) then !for ARVE-grid:

     ! Calculate the attentuation of light from dsw_t as well as the partitioning between diffuse and direct
     ! downwelling radiation at the surface (output in W m-2, alse sets met%dswb in kJ m-2 d-1)
       call surface_rad(j,i)

    else !for ARVE-Point:

     ! Calculate the partitioning of station measured light into direct and diffuse beams (output in W m-2)
       call simplesurfrad(j,i)

    end if
    
    if (.not. do_allocation_daily) then  !--------for annual allocation only--
    ! Calculate phenology for annual allocation scheme
      if (i == 2) then  ! Only once per day
        do pft = 1,npft
          if (present(pft)) then  ! Only if the PFT is present in the grid cell
            call phenology_annualalloc(j,pft)
          end if
        end do
      end if
      
    end if  !-----------------------------
    

  ! Calculate snow and rain flux to surface also calculates
  ! canopy water (interception, throughfall and canopy drip)
    call newsnowflux(j,i)
                !calls
                !->canopy_interception

  ! Find the leaf area, stem, and foliar protective cover sum for all pfts accounting for snow burial
    if (any(present)) then
      call leafburial(j)
    end if

  ! Calculate soil thermal properties
    call soilthermalprop(j)

  ! Calculate the albedo of the canopy and soil
  !  call albedo(j)  !FLAG could put in a joboptions flag to make this daily?
    
  !>>=============
  ! Calculate soil temperature, soil phase change and surface/veg heat fluxes (will run on short sub-timestep)
  ! also calculates GPP (photosynthesis), radiative fluxes
  ! surface resistances. All calculations are performed per PFT.

    call shorttimestep(yr,j,i)

  !>>=============

  Tsoil_ave(1:nl) = Tsoil_ave(1:nl) / counter_lim

  do pft = 1,npft
    if (present(pft)) then  ! Only if the PFT is present in the grid cell

      ! Calculate autotrophic respiration
      call auto_resp(j,i,pft)
   
      ! Calculate the carbon isotope fractionation in plants
      ! NOTE this is turned off, the changes in the isotope tracking scheme have not been fixed up
      ! for the changes due to the daily carbon allocation. Thus the isotope values from ARVE will
      ! not be correct until fixed so calculation is only allowed on the annual C alloc scheme. JM 05.01.2011.
      if (.not. do_allocation_daily) then
        call isotope(pft,j,i)
      end if

    end if   !end pft present
  end do   !end pft loop

  if (do_allocation_daily) then  !-------------- for daily allocation only!
  
          if (i == 2) then  !only during night

          ! Calculate phenology
            call phenology(j)

          ! Calculate daily carbon allocation, turnover, reproduction C loss, and sapwood to heartwood conversion
            call daily_allocation(j)  

          ! Calculate change in root distribution 
            call roots(j)

          end if
          
  end if  !-------------------------------------------------

  ! Calculate soil water flux and depth of the active layer
    call soilwaterflux(j)

  ! Calculate soil and litter respiration
    if (any(present)) then
      call het_resp(j,yr)
    end if

  ! Calculate methane emissions and the associated isotope ratio
    !call methane(j)

  ! If snow exists calculate snow layer compaction, division, and combination
    if (-snl > 0) then
      call snowdynamics(j)

    end if
   
               ! write(87,'(3i5,6f10.4)') yr, doy,i, met(j,0)%temp(i),met(j,0)%tmax,met(j,0)%tmin,met(j,0)%prec,met(j,0)%cldf,sv(j)%zsno
               ! call flush(87)

                !                if (i == 1) then  !daytime timestep

                                !'(i5,6f12.6)'
                 !               write(15,*)doy,met(j,0)%tmax,met(j,0)%tmin,sv(j)%Psi(1:5),sv(j)%Wliq(1:5)

                  !              end if
                   !             call flush(15)
                                !write(*,*)'=============================  done',doy,i
                                !call flush_(35)

                    !            write(25,*)met(j,0)%temp(i),met(j,0)%tmax,met(j,0)%tmin,sv(j)%zsno
                     !           call flush(25)
                                !-------
                                
       22 format(3i5,20es12.3)

       !Print out some PFT variable info to file  
     !  if (i == 2) then                       
     !  do l = 1,npft
     !    if (sv(j)%present(l)) then
     !       write(17,22)yr,doy,l,sv(j)%crownarea(l),sv(j)%rootdepth(l),sv(j)%lm_max(l,1),sv(j)%rm_ind(l,1)&
     !                             ,sv(j)%lm_ind(l,1),sv(j)%litter_ag_fast(l,1),sv(j)%sm_ind(l,1)&
     !                             ,sv(j)%hm_ind(l,1),ov%acflux_fire(1),sv(j)%lai_ind(l),sv(j)%height(l),sv(j)%NSC_storage(l,1)&
     !                             ,veg%dgpp(1,l),veg%dnpp(1,l),veg%dresp(1,l)&
     !                             ,sv(j)%lai_sno(l),veg%supply_ts(l),veg%demand_ts(l),surf%avg_cZ,veg%scatcf_tot(1,l)

     !      end if
     !   end do
     !   end if
     !   call flush(17)
   

  ! Do some incrementing of daily calculated values 
  call daily_increment(j,i)

                                !-----Print some output to screen
                                !write(55,'(3i5,3f10.4)')yr,doy,i,sv(j)%last_gpet(2),met(j,0)%temp(i),sv(j)%zsno
                               ! write(55,'(3i5,3f10.4)')yr,doy,i,sv(j)%tsoil_ave(1)-Tfreeze,sv(j)%Tveg_ave(1,6)-Tfreeze,met(j,0)%temp(i)
                               ! call flush(55)

  ! Reset daily average soil temp (used in autoresp)
  Tsoil_ave(1:nl) = 0._dp

 end do  !day/night loop

end subroutine calcdaily
