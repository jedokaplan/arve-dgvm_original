subroutine dailygrid_driver(yr,ncells)

! ARVE-Grid
! Central driver called from arve_dgvm
! from dailygrid_driver the daily and annual subroutines are called
! this is for gridded data (point simultations use dailysite_driver)

use arveparams,     only : sp,dp,mdays,nl,Tfreeze,daysec,imdays
use statevars,      only : sv,ov,doy,dayl
use iovariables,    only : ibuf,spinupyears,dosoilm_assim
use metvarsmod,     only : atm,met,dm,atm
use insolationmod,  only : orbital_parameters,dayins
use weathergenmod,  only : weathergen,rmsmooth
use diurnaltempmod, only : diurnaltemp
use pft_state_variables, only : veg,initpftvars
use airmass,        only : get_airmass
use coszen,         only : avgZen
use soilstate_vars, only : fsnow_mo,surf,init_surf_state_vars
use iovariables, only : soil            !for testing purposes
use udunits                             !not connected up fully -Feb 16 08
use calcannual, only : annual_calcs

use soilmoistureforcing, only : anom_in,nudgesoil

implicit none

!arguments
integer, intent(in)             :: yr                  !year of simulation
integer, intent(in)             :: ncells              !number of cells in simulation

!pointers
integer, pointer                :: gnl          !index value of lowest 'soil' layer                                                 
real(dp), pointer, dimension(:) :: temp         !vector for mean temperature of the day and night timesteps: 1=day, 2=night (deg C) 
real(dp), pointer               :: sdecl        !solar declination (deg)                                                            
real(dp), pointer               :: sdecl_n      !the next day's solar decl (deg)                                                    
real(dp), pointer               :: dsw_t        !integrated calculated downwelling shortwave radiation top of atmos(KJ m-2 day)     
real(dp), pointer               :: dsw_tn       !the next day's dsw_t                                                               
real(dp), pointer               :: prec         !24 hour total precipitation (mm)                                                   
real(dp), pointer               :: Pjj          !precipitation equitability index                                                   
real(dp), pointer, dimension(:) :: met_res      !residuals, used in weathergenmod
real(dp), pointer               :: p_wet_m      !precipitation of the wettest month                  
real(dp), pointer               :: p_dry_m      !precipitation of the driest month                   
real(dp), pointer               :: p_mo_avg     !mean annual monthly precipitation                   
real(dp), pointer               :: prec_ra      !precipitation running sum over 30 days (mm)         
real(dp), pointer, dimension(:) :: mprec        !quasi-monthly precipitation total (mm)
real(dp), pointer, dimension(:) :: tmp          !monthly mean temp (C) 
real(dp), pointer, dimension(:) :: dtr          !monthly mean diurnal temperature range
real(dp), pointer, dimension(:) :: pre          !monthly precip (mm)
real(dp), pointer, dimension(:) :: cld          !monthly mean cloud fraction
real(dp), pointer, dimension(:) :: dtmx         !tmax (daily smoothed vectors (365 days))           
real(dp), pointer, dimension(:) :: dtmn         !tmin (daily smoothed vectors (365 days))           
real(dp), pointer, dimension(:) :: dcld         !cloud fraction (daily smoothed vectors (365 days)) 
real(dp), pointer               :: zsno_mo      !monthly cumulative snow depth (m) 
real(dp), pointer, dimension(:) :: daytempv     !vector of daytime averaged temperatures
real(dp), pointer               :: twm          !maximum monthly temperature for the year being simulated
real(dp), pointer               :: tminm        !minimum monthly temperature for the year being simulated
real(dp), pointer, dimension(:) :: runtemp      !holds one month of temp data for mean monthly calc.
real(dp), pointer               :: aprec        !annual precipitation (mm)
real(dp), pointer               :: MAT          !mean annual temperature (C)
real(dp), pointer               :: MAT_d        !stores running total of MAT (C)
real(dp), pointer               :: mo_runoff_tot !monthly total overland runoff (mm)
        
!local variables
integer                         :: j,l            !counters
integer                         :: m            !month
integer                         :: mo           !month to pass to weathergen
integer                         :: dom          !day of month
integer, dimension(1)           :: wm           !warmest month
integer, dimension(1)           :: cm           !coldest month
integer, dimension(1)           :: wet_m        !wettest month
integer, dimension(1)           :: dry_m        !driest month
real(dp)                        :: p_wm         !precipitation of the warmest month
real(dp)                        :: p_cm         !precipitation of the coldest month
real(dp), dimension(size(ibuf(1)%tmp)) :: tmn
real(dp), dimension(size(ibuf(1)%tmp)) :: tmx
real(dp), dimension(365)        :: dmsoilanom   !smooth daily interpolation of soil moisture anomaly
real(sp)                        :: dfm          !1/number of days in the month
real(sp)                         :: Sbar        !mean soil wetness fraction over all soil layers
real(sp)                        :: Wmax         !Tpor * dzmm
real(sp)                        :: Wtot         !sum of Wliq + Wice

!local variables for output
!real(dp), dimension(365)    :: tvect
real(sp), dimension(365)        :: zsno_y
real(sp), dimension(365)        :: drh_y
real(sp), dimension(nl,365)     :: stemp
real(sp), dimension(nl,365)     :: sliq_y
real(sp), dimension(nl,365)     :: tliq_y
real(sp), dimension(nl,365)     :: wliq_y
real(sp), dimension(nl,365)     :: wice_y

!parameter
real(dp), parameter             :: missing = -999.9_dp  !missing value for climate files

!---------------------------------------------

! Calculate orbital parameters for the current year (eccentricity, precession, obliquity, perihelion)
! not dependendent on grid location. Expects year as integer

   call orbital_parameters(0)

!Parallel grid loop

do j = 1,ncells

  !Check if it is a valid cell, if not cycle

        if (sv(j)%valid_cell) then  ! the check for valid SOIL cells took place in arve_dgvm.f90

          !Check if it is a valid CLIMATE cell
          
          !temperature is missing or cloud cover is negative JM (17.10.2010)

          if (ibuf(j)%tmp(1) < missing+1. .or. ibuf(j)%cld(1) < 0._dp) then
            sv(j)%valid_cell = .false.
            cycle
          end if

        else  !Not a valid cell, cycle
          cycle
        end if

        write(*,'(a,2i8,2f8.2)')'grid:',j,yr,sv(j)%lon,sv(j)%lat

  ! Point pointers with the cell index, j, and fix any potential problems.

        ! Monthly input values (12 mos)
        tmp => ibuf(j)%tmp   !mean temperature
        dtr => ibuf(j)%dtr   !diurnal temperature range
        cld => ibuf(j)%cld   !cloud cover
        pre => ibuf(j)%pre   !precip

        ! Correct the input values for potential inconsistencies
        dtr = max(dtr,0._dp)

        ! Calculate min and max temperature
        tmn = tmp - 0.5_dp * dtr
        tmx = tmp + 0.5_dp * dtr

        ! Convert cloud cover (in percent) to a fraction
        cld = 0.01_dp * cld

        ! Point pointers for pseudo-daily smoother vectors (365 days)
        dtmx => dm(j)%tmx
        dtmn => dm(j)%tmn
        dcld => dm(j)%cld

        prec_ra  => dm(j)%prec_ra
        mprec    => dm(j)%mprec
        p_wet_m  => dm(j)%p_wet_m
        p_dry_m  => dm(j)%p_dry_m
        p_mo_avg => dm(j)%p_mo_avg

        ! Point some more pointers
        temp    => met(j,0)%temp
        sdecl   => surf%sdecl
        sdecl_n => sv(j)%sdecl_n
        dsw_t   => surf%dsw_t
        dsw_tn  => sv(j)%dsw_tn
        prec    => met(j,0)%prec
        Pjj     => atm%Pjj
        gnl     => sv(j)%gnl
        met_res => met(j,0)%met_res
        zsno_mo => surf%zsno_mo
        daytempv=> atm%daytempv
        twm     => atm%twm
        tminm   => atm%tminm
        runtemp => atm%runtemp
        aprec   => atm%aprec
        MAT     => atm%MAT
        MAT_d   => atm%MAT_d
        mo_runoff_tot => surf%mo_runoff_tot
        
        !Initialize
        runtemp =    0._dp
        twm     = -100._dp
        tminm   =  100._dp
        aprec   =    0._dp
        doy     =    0
        met_res =    0._dp

  !-----------------
  
  ! Begin calculations

  ! Calculate daylength and top of atmosphere shortwave for first day

  call dayins(j,doy)

  ! Shift values to previous days position (only happens to prepare for 1st day)
  dayl(1) = dayl(2)
  sdecl = sdecl_n
  dsw_t = dsw_tn

  ! Print out some stuff if desired.
  goto 120
          if (yr == 1) then
              write(*,*)soil(j)%maxd
              write(*,'(a4,7a12)')'','zipos','sand','silt','clay','rock','sorg','bulk'
           do l = 1,nl
            write(*,'(i4,7f12.4)')l,sv(j)%zipos(l),sv(j)%sand(l),sv(j)%silt(l),sv(j)%clay(l),sv(j)%rock(l)*100.,sv(j)%sorg(l),sv(j)%bulk(l)
           end do
           read(*,*)
          end if
  120 continue
  
         ! output some yearly climate data
         !write(0,*)
         !write(0,'(12f7.1)')tmp
         !write(0,'(12f7.1)')pre
         !write(0,'(12f7.1)')tmx
         !write(0,'(12f7.1)')tmn
         !write(0,'(12f7.1)')cld
         !write(0,'(12f7.1)')ibuf(j)%wet
         !read(*,*)


  ! Create a time series of one year of pseudo-daily means for min and max temp and cloud cover

  if (dm(j)%init) then

    ! For first year of model run, use the first month as boundary condition

    call rmsmooth(tmx,imdays,[tmx(1),tmx(12)],dtmx)
    call rmsmooth(tmn,imdays,[tmn(1),tmn(12)],dtmn)
    call rmsmooth(cld,imdays,[cld(1),cld(12)],dcld)
    
    MAT = 0._dp
    MAT_d = 0._dp

    ! Because the diurnal temperature calculation requires temperature of the previous day,
    ! for the very first day of the first year of the model run, we precalculate a day's weather

    call weathergen(j,1,1)
    dm(j)%init = .false.

  else

    call rmsmooth(tmx,imdays,[dtmx(365),tmx(12)],dtmx)
    call rmsmooth(tmn,imdays,[dtmn(365),tmn(12)],dtmn)
    call rmsmooth(cld,imdays,[dcld(365),cld(12)],dcld)

  end if

    !Min/max monthly mean temp accounting

    twm   = maxval(tmp)
    tminm = minval(tmp)    

  ! Set up the precipitation equitability index for swrad.f90. This uses the month with the maximum and minimum
  ! monthly temperature. If planning to couple ARVE into a GCM then this will need to be changed as
  ! the hottest/coldest months will not be known a priori for each year. This set up is fine in ARVE stand-alone
  ! JM Sept 30 2008

  wm   = maxloc(tmp)
  cm   = minloc(tmp)

  p_wm = pre(wm(1))
  p_cm = pre(cm(1))

  !Find wettest and driest months

  wet_m = maxloc(pre)
  dry_m = minloc(pre)
  p_wet_m = pre(wet_m(1))
  p_dry_m = pre(dry_m(1))

  ! Annual mean monthly precipitation (for annual allocation phenology)  
  
  p_mo_avg = sum(pre)/12.

  ! Calculation of the "precipitation equitability index"

  if (p_wm + p_cm > 0.) then
    Pjj = 2._dp * (p_wm - p_cm) / (p_wm + p_cm)
    Pjj = max(Pjj,0._dp)
  else
    Pjj = 0._dp
  end if

  !--------------------------------
  !calculation for soil moisture anomaly experiment (none-base ARVE code)
  if (yr > spinupyears .and. dosoilm_assim) then
    !if it is the first year of the transient, interpolate the climatological mean soil moisture to daily
!    if (yr == spinupyears + 1) call rmsmooth(sv(j)%meansoilwater,imdays,[sv(j)%meansoilwater(12),sv(j)%meansoilwater(1)],sv(j)%dmsoilm)
    !on every year, interpolate the soil anomaly to smooth daily values
    call rmsmooth(anom_in(1,1,:),imdays,[anom_in(1,1,12),anom_in(1,1,1)],dmsoilanom)
  end if
  !----------------------------------
  
  !Start daily loop

  doy = 1
  do m = 1,12

    dfm = 1. / mdays(m)
    Sbar = 0.

    do dom = 1,imdays(m)

     !--------------------------------
     ! if (yr > spinupyears .and. dosoilm_assim) then  !only calculate if not in spinup and desired
     !    call nudgesoil(j,sv(j)%dmsoilm(doy),dmsoilanom(doy))
     ! end if
     !--------------------------------
      
      !Calculate top of the atmosphere solar radiation, daylength; sets: sv%dayl, sv%dsw_t (KJ m-2 day)

      call dayins(j,doy+1)

      !Calculate weather for this day: min/max temp, precip, surface sw, and cloud cover
      !generate following day's weather. This is still entered if the doweather flag is false.

           if (dom == imdays(m)) then !since weathergen is one day ahead, move the month to the next when you are on
             if (m /= 12) then        ! the last day of the present month.
               mo = m + 1               
             else 
               mo = 12
             end if
           else 
             mo = m
           end if 

      call weathergen(j,mo,doy)

      !Calculate integrated daytime and nighttime temperatures, sets dm(j)%temp(1:2)

      call diurnaltemp(j)

      !Add day to vector of DAYTIME averaged temperatures

      daytempv(doy) = temp(1)

      !Find the diurnally averaged cosine of the solar zenith angle

      call avgZen(j)

      !Calculate elevation-corrected daily mean and instantaneous air mass

      call get_airmass(j)

      !Calculate daily biophysics and biogeochemistry
      
      call calcdaily(j,yr)

      !Increment yearly precipitation and MAT, and shift the dayl,sdecl,dsw_t

      aprec = aprec + prec
      MAT_d = MAT_d + (dayl(1) * temp(1) + (daysec - dayl(1)) * temp(2)) / daysec  !weighted value of day and night temps
      mprec(:) = eoshift(mprec(:),1,prec)
      prec_ra = sum(mprec(:)) !find the running 30 day sum

      !Move the 'next day' values to present day in preparation for another daily cycle

      dayl(1) = dayl(2)
      sdecl = sdecl_n
      dsw_t = dsw_tn

      !Accumulate daily soil stats into annual array

      drh_y(doy)  = veg%drh(1)
      zsno_y(doy) = sv(j)%zsno

      stemp(:,doy) = sv(j)%Tsoil(1:nl) - Tfreeze
      sliq_y(:,doy)  = sv(j)%Tliq(1:gnl) / sv(j)%Tsat(1:gnl)
      tliq_y(:,doy)  = sv(j)%Tliq(1:gnl)
      wliq_y(:,doy)  = sv(j)%Wliq(1:gnl)
      wice_y(:,doy)  = sv(j)%Wice(1:gnl)

      !-----------
      !mean soil moisture over all layers (soilassim code)

      Wmax = sum(sv(j)%Tsat(1:3) * sv(j)%dzmm(1:3))  !total soil porosity (mm)
      Wtot = sum(sv(j)%Wliq(1:3) + sv(j)%Wice(1:gnl))  !total soil water content (mm) (liq + ice)
      Sbar = Sbar + Wtot / Wmax * dfm  !mean 

      !-----------

      !Increment day counter
      doy = doy + 1

    end do !day in month loop
    
    ! output some monthly mean values    
    ov(j)%runoff(m) = real(mo_runoff_tot * dfm)    
    ov(j)%mw1(m) = real(Sbar)
    !ov(j)%albedo(m,1) = albedo_mo_tot(1) * dfm ! vis
    !ov(j)%albedo(m,2) = albedo_mo_tot(2) * dfm ! nir
    ov(j)%snow_frac(m) = real(fsnow_mo * dfm)
    ov(j)%snow_depth(m) = real(zsno_mo * dfm)
    
    !reset values for coming month
    mo_runoff_tot = 0._dp
    runtemp = 0._dp
    fsnow_mo = 0._dp
    zsno_mo = 0._dp
   ! albedo_mo_tot = 0._dp
    
  end do !month loop

  !If you wish to rewind the text file output, do so here.
  !rewind(87) !FLAG
  !rewind(15) !FLAG
  !rewind(25) !FLAG
  !rewind(35) !FLAG
  !rewind(12) !FLAG

        !Do annual vegetation dynamics 
        call annual_calcs(j,yr)
        

  !tally running mean for soil moisture if in spinup (none-base ARVE code)
!  if (yr <= spinupyears) then
!    sv(j)%meansoilwater = (ov(j)%mw1 + real(yr-1) * sv(j)%meansoilwater) / real(yr)
!  end if
  !print soil variables to output in netcdf
  !ov(j)%stemp(:,:)=stemp(:,:)

  !reset the pft vars
  call initpftvars()
  call init_surf_state_vars()

end do !done year for this grid cell.

end subroutine dailygrid_driver
