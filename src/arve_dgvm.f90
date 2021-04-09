program arve_dgvm

!to call the program the command line needs to be in the form:
!arvegrid [joboptions file] [lon]/[lat] [outputfile]

!requires libs
!-L/usr/local/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lSystemStubs -lm -lz

use typeSizes
use netcdf
use netcdf_error
use iovariables,         only : ibuf,inputclimlen,timebuflen,climate_path1,soil_path,slope_path,endy, &
                                latvect,lonvect,endx,cntx,cnty,srtx,srty,climate_path2,dotransient, &
                                inputlonlen,inputlatlen,niv,vname,sname,climfid,va,spinupyears,transientyears, &
                                soilfid,vsname,nsiv,soil,slope,siv,slopefid,CO2file, &
                                dosoilm_assim,dovegetation,doweather,do_allocation_daily,climateyears, &
                                do_CLM_photosynthesis,dosoilm_assim,spinupCO2_12C,dospinup,cycle_spinup,cycle_years, &
                                save_spinup_file,save_transient_file,save_spinup_state,save_transient_state,fromsavedstate, &
                                saved_model_state_file,spinupCO2_13C,spinupCO2_14C,outputfile,fmaxfid,fmax_path,climatemonths
use arveparams,          only : dp,npft,mdays,Pstd
use statevars,           only : sv,ov,initstatevars,area,gridded,CO2
use pftparametersmod,    only : initpftpars
use pft_state_variables, only : initpftvars,veg
use metvarsmod,          only : met,dm
use pftinitassign,       only : pftassign
use soilmoistureforcing, only : initsoilmfile,getsoilanom
use soilinit,            only : soillayers,soilprep
use input_dataset_drivers, only : initclimate,soildriver,slopedriver,openCO2file,readCO2file
use udunits  !not connected up fully -Feb 16 08
use soilstate_vars,     only : surf,init_surf_state_vars

implicit none

!------
        
!local variables
integer :: memcheck
integer :: ncells
character(60) :: status_msg
integer(2), dimension(:,:,:), allocatable :: var_in
real(dp), dimension(:,:,:), allocatable :: realvar_in
!integer :: dimid
integer :: varid
integer :: i,j,x,y,year
integer :: nyears
integer :: tpos
integer :: tlen
integer :: px,py
real(dp) :: yr_tot
integer :: mo

!=============================================================================

write(0,*)
write(0,*)'***--- Starting up ARVE-DGVM v.0.95 ---***'
write(0,*)'------        arve.epfl.ch          ------'
write(0,*)

! Open text files for diagnostic output

!open(05,file='diag_test.dat',status='unknown')
!open(12,file='subdaily.dat',status='unknown')
!open(15,file='diag_unknown.dat',status='unknown')
!open(16,file='pft_output.dat',status='unknown')
!open(17,file='pft_daily.dat',status='unknown')
!open(25,file='diag_soiltemp.dat',status='unknown')
!open(35,file='diag_soilmoisture.dat',status='unknown')
!open(45,file='diag_heatest.dat',status='unknown')
!open(55,file='diag_latent.dat',status='unknown')
!open(65,file='annual_meteo.dat',status='unknown')
!open(87,file='test_meteo.dat',status='unknown') 
!open(99,file='/Volumes/arve/shared/06_Joe/allometry/lpj_allom_-90_45.dat',status='old')

! Begin by setting up the job

call initjob()

! Open the binary files for save state or to load a saved state
if (fromsavedstate) open(97,file=saved_model_state_file,status='old',form='unformatted')
if (save_spinup_state) open(98,file=save_spinup_file,status='unknown',form='unformatted')
if (save_transient_state) open(99,file=save_transient_file,status='unknown',form='unformatted')

!status = utopen('/usr/local/etc/udunits.dat')

! Open the soil initial conditions files and allocate
! the soils input matrix (and lat and lon vect). Also does soil depth
! allocates lonvect, latvect and soil%

call soildriver(soil_path)

! Set up boundaries

srtx = max(1,srtx)
srty = max(1,srty)

cntx = min(cntx,inputlonlen)
cnty = min(cnty,inputlatlen)

endx = srtx + cntx - 1
endy = srty + cnty - 1

! Set up for the fractional slope and maximum saturated fraction of grid cells.

call slopedriver(slope_path,fmax_path)

! Set up the climate

if (dospinup)then
        call initclimate(climate_path1)
else  !doing only transient
        call initclimate(climate_path2)
end if

! Memory check and find number of cells to model

memcheck = (inputclimlen*timebuflen*3*4 + inputlonlen*inputlatlen*2*2*4) / 1048576
ncells = min(cntx*cnty,inputclimlen)

! Allocate variable arrays

allocate(sv(ncells))
allocate(met(ncells,-1:1))   !met has three positions, -1 = previous day, 0 = current day, 1 = next day
allocate(dm(ncells))
allocate(ov(ncells))

! Initialize the state variables

call initstatevars()
call init_surf_state_vars()

! Open the CO2 file

call openCO2file(CO2file)

! Set up year loop
if (dotransient) then

  !if the transient is starting from a pre-saved model state, read that information into sv now.
  if (fromsavedstate) then
     write(0,'(a,a)')'State vars are being read-in from file:',trim(saved_model_state_file)
     write(0,*)
     read(97)sv 
     close(97)
  end if
  
  if (.not. dospinup) spinupyears = 0
  
  nyears = spinupyears + transientyears
  
else !only do a spinup
 
  nyears = spinupyears
  transientyears = 0

end if

! Output initial conditions to screen

write(0,'(a,i4,a,i4,a,i4,a)')'Your input files have ',inputlonlen,' columns and ',inputlatlen,' rows with a buffer of ',timebuflen,' months'
write(0,'(a,i,a,i4,a)')'Number of cells in climate vector:',inputclimlen,'. Climate data set has ',climateyears,' years total'
write(0,'(a,i4,a,i4,a,i4,a,i6)')'Total memory requirement of input data: ',memcheck,' MB. Cells to calculate: ',cntx,' by ',cnty,' = ',ncells
write(0,'(a,i4,a)')'I will simulate ',nyears,' years total' 
if (dospinup) write(0,'(a,i4,a,f8.2,a,a)')'The spinup is ',spinupyears,' years with a CO2 of ',spinupCO2_12C,' and climate is: ',trim(climate_path1)
if (cycle_spinup) write(0,'(a,i4,a)')'The spinup will cycle over the first',cycle_years,' years of my spinup climate file'
if (dotransient) write(0,'(a,i4,a,a)')'The transient is ',transientyears,' years. The transient climate file is: ',trim(climate_path2)
if (dotransient) write(0,'(a,a)')'and the transient CO2 is: ',CO2file
if (save_spinup_state) write(0,'(a,a)')'Spinup state is saved to: ',save_spinup_file
if (save_transient_state) write(0,'(a,a)')'Transient state is saved to: ',save_transient_file
write(0,'(a,a)')'The NetCDF output file is:',outputfile
write(0,'(a,a)')'soil file is:',soil_path
write(0,'(a,l,a,l,a,l)')'Vegetation: ',dovegetation,'  Weather generator: ',doweather,'  Daily C allocation: ',do_allocation_daily
write(0,'(a,l,a,l)')'Soil anomaly exp: ',dosoilm_assim,'  CLM photosynthesis: ',do_CLM_photosynthesis
write(0,'(a,2f8.2,a,2i4,a)')'Starting calculation at:',lonvect(srtx),latvect(srty),' (',srtx,srty,')'
write(0,*)

! Allocate the soil data arrays to read in the soil data

allocate(var_in(cntx,cnty,timebuflen))
allocate(realvar_in(cntx,cnty,5))

!-----------------

! Get soil data

do i = 1,nsiv

  status = nf90_inq_varid(soilfid,vsname(i),varid)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(soilfid,varid,realvar_in,start=[srtx,srty,1],count=[cntx,cnty,5])
  if (status /= nf90_noerr) call handle_err(status)

  do y = 1,cnty
    do x = 1,cntx
      j = x + cntx * (y-1)

      ! Set up soil data
      select case (i)
      case(1)
        soil(j)%sdto = realvar_in(x,y,:)
      case(2)
        soil(j)%stpc = realvar_in(x,y,:)
      case(3)
        soil(j)%clpc = realvar_in(x,y,:)
      case(4)
        soil(j)%totc = realvar_in(x,y,:)
      case(5)
        soil(j)%bulk = realvar_in(x,y,:)
      case(6)
        soil(j)%cfrag = realvar_in(x,y,:)
      case(7)
        soil(j)%tawc = realvar_in(x,y,:)
      end select

    end do
  end do
end do

deallocate(realvar_in)
allocate(realvar_in(cntx,cnty,1))

! Assign dataset value of soil max depth

status = nf90_inq_varid(soilfid,'maxd',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(soilfid,varid,realvar_in,start=[srtx,srty,1],count=[cntx,cnty,1])
if (status /= nf90_noerr) call handle_err(status)

do y = 1,cnty
  do x = 1,cntx

    j = x + cntx * (y-1)
    soil(j)%maxd = (realvar_in(x,y,1))

  end do
end do

deallocate(realvar_in)
allocate(realvar_in(cntx,cnty,1))

!---------------------------------------
! Find the fraction of the gridcell below slope tolerances, gridcell size and the mean elevation

do i = 1,siv

  status = nf90_inq_varid(slopefid,sname(i),varid)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(slopefid,varid,realvar_in,start=[srtx,srty,1],count=[cntx,cnty,1])
  if (status /= nf90_noerr) call handle_err(status)

  do y = 1,cnty
    do x = 1,cntx
      j = x + cntx * (y-1)

      ! Assign topographic data from dataset
      select case (i)
      case(1)
        slope(j)%elevation = (realvar_in(x,y,1))
      case(2)
        slope(j)%sf2 = max(0._dp,(realvar_in(x,y,1)))
      case(3)
        slope(j)%sf075 = max(0._dp,(realvar_in(x,y,1)))
      case(4)
        slope(j)%sf05 = max(0._dp,(realvar_in(x,y,1)))
      end select

      ! Find the grid cell size
      sv(j)%grid_area =  area(sv(j)%lat,inputlonlen,inputlatlen)

    end do
  end do

end do

deallocate(realvar_in)
allocate(realvar_in(cntx,cnty,1))

!------------------------------------------------
! Take in the maximum fractional saturated area

status = nf90_inq_varid(fmaxfid,'FMAX',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(fmaxfid,varid,realvar_in,start=[srtx,srty,1],count=[cntx,cnty,1])
if (status /= nf90_noerr) call handle_err(status)

  do y = 1,cnty
    do x = 1,cntx
      j = x + cntx * (y-1)
         soil(j)%fmax = (realvar_in(x,y,1))
    end do
  end do

deallocate(realvar_in)

!------------------------------------------------
! Initializations of pft parameters

call initpftpars()
call initpftvars()

! This is a gridded data run so set the signpost to true

gridded = .true.

if (cycle_spinup) then
  ! when cycle_spinup is set, set the length of the climate file to 
  ! the smaller of the actual length of the climate file or cycle_years * 12 months
  tlen = min(climatemonths,cycle_years * 12)
else
  !otherwise tlen is just the length of the climate file.
  tlen = climatemonths  
end if

! Create netcdf output file

call netcdf_create()

!----------------------------
!***soil moisture assimilation -Ben Cook**
if (dosoilm_assim) then
     call initsoilmfile()
end if
!----------------------------

! Initialize soil and PFT properties

do y = 1,cnty
  do x = 1,cntx
    j = x + cntx * (y-1)

    ! Check if it is a valid cell, if not cycle

    if (soil(j)%sdto(1) < 0._dp .or. soil(j)%maxd <= 0._dp) then  ! if the upper layer has no defined soil properties then don't calc.

      sv(j)%valid_cell = .false.
      cycle

    else

      ! Assign lm_sapl, rm_sapl, sla, sm_sapl
      call pftassign

      ! Assign soil layers and allocate the soil hydrological variables
      call soillayers(j)

      ! Intialize and calculate soil properties
      call soilprep(j)

      ! Elevation is used to get the atmospheric pressure (Diehl 1925 (from Lloyd and Farquhar 1994))
      sv(j)%Patm = Pstd*(1._dp - slope(j)%elevation /44308.)**5.2568

      ! Also calculate the reference height pressure (30m above ground) (same ref as above)
      sv(j)%Patm30 = Pstd * (1._dp - (slope(j)%elevation + 30._dp) / 44308.)**5.2568
      
      !Check if this 'valid' grid cell has a valid fmax value from the dataset
      !if not, assigned it an arbitrary value of 0.35. This will not happen in many
      !instances (JM 13.05.2011)
      if (soil(j)%fmax == -999.99) soil(j)%fmax = 0.35

      ! This is the first year of the simulations so init is true.
      dm(j)%init = .true.

    end if !valid cell loop

  end do
end do

!----------
! YEARLOOP

do year = 1,nyears

  ! set tpos, which determines the year of the climate dataset to load
  if (year < spinupyears+1) then  ! for the spin up
    tpos = mod(1 + 12 * (year - 1),tlen)
    if (tpos == 0) tpos = tlen
  else ! and for the transient
    tpos = mod(1 + 12 * ((year - spinupyears) - 1),tlen)
    if (tpos == 0) tpos = tlen
  end if
  
  ! Make internal model 'write' statements to output on the screen.

  write(status_msg,*)'working on year:',year

  call overprint(status_msg)

  ! Transient - load separate climate input file for transient
  if (year == spinupyears+1 .and. dotransient .and. dospinup) then
    
    write(0,*)'loading transient climate input'
        
    tpos = 1
    
    ! Deallocate old variables
    deallocate(ibuf)
    
    ! set up the transient climate file.
    call initclimate(climate_path2)

  ! Spinup using the first cycle_years of climate repeatedly
  else if (year == cycle_years + 1 .and. year < spinupyears .and. dospinup .and. cycle_spinup) then
 
   tpos = 1
       	
  end if
  
  !-----------------
  ! Get climate data

  do i = 1,niv
  
    status = nf90_inq_varid(climfid,vname(i),varid)
    if (status /= nf90_noerr) call handle_err(status)
        
    status = nf90_get_var(climfid,varid,var_in,start=[srtx,srty,tpos],count=[cntx,cnty,12])
    if (status /= nf90_noerr) call handle_err(status)

    do y = 1,cnty
      do x = 1,cntx
        j = x + cntx * (y-1)

        if (.not. sv(j)%valid_cell) cycle  ! do not get climate for this gridcell, it is not valid

        px = srtx + x-1
        py = srty + y-1
        sv(j)%lon = lonvect(px)
        sv(j)%lat = latvect(py)

        ! Set up yearly climate data
        select case (i)
        case(1)
          ibuf(j)%cld = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
        case(2)
          ibuf(j)%dtr = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
        case(3)
          ibuf(j)%pre = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
        case(4)
          ibuf(j)%tmp = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
        case(5)
          ibuf(j)%wet = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
        end select

        !FLAG FLAG this is a temporary fix for problems with dataset pre amounts
        !JM 09.05.2011
        do mo = 1,12
        yr_tot = sum(ibuf(j)%pre)
        if (ibuf(j)%pre(mo) > 0.9 * yr_tot) ibuf(j)%pre(mo) = ibuf(j)%pre(mo) * 0.01
        end do
        !----------------------------------
        
      end do

    end do
  end do

  !--------
  ! Read in CO2 values
  if (dospinup .and. year < spinupyears+1) then
      CO2(1) = spinupCO2_12C
      CO2(2) = spinupCO2_13C
      CO2(3) = spinupCO2_14C
  else !read in transient CO2 from file    
      call readCO2file
  end if    

  if (year > spinupyears .and. dosoilm_assim) then

    ! Get the soil moisture anomaly for this year
    call getsoilanom(year-spinupyears)

  end if

  !--------------------------------------------------------
  ! CORE: Call the driver for daily (and annual) subroutines

  call dailygrid_driver(year,ncells)

  ! Output
  call netcdf_output(year)
  
  ! Save state if done spinup and you want to save the spinup state
  if (year == spinupyears .and. save_spinup_state) then
  write(*,*)'writing to spinup save-file',save_spinup_file
     write(98)sv
     call flush(98)
     close(98)
  end if
  
  ! Save state if done the transient and you want to save the transient state   
  if (year == spinupyears+transientyears .and. save_transient_state) then
     write(*,*)'writing to transient save-file',save_transient_file
     write(99,*)sv
     call flush(99)
     close(99)
  end if   
  
end do  ! YEARLOOP

call netcdf_close()

write(0,*)
write(0,*)'***--- Simulation is complete ---***'
write(0,*)'***---      That was fun.     ---***'
write(0,*)'***---         Goodbye.       ---***'
write(0,*)

end program arve_dgvm
