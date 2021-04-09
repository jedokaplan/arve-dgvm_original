program arvepoint

!THIS HAS NOT BEEN UPDATED RECENTLY, check dailygrid_driver for likely needed changes. JM 03.03.2011
!does not have fmax set up yet.

!Use to run ARVE for points using weather station data. Requires the use of sitedatamod to open
!the station weather information. Added May 27th 08 by JM and JOK


!New version 9 June 2009
!lon, lat, soils, elevation and slope data are read from joboptions namelist
!climate data expects a value for sun fraction - if you don't have it specify a flag (-1) and then it will be estimated by the program

use typeSizes
use netcdf
use netcdf_error
use iovariables,        only : soil,CO2file,slope
use arveparams,         only : sp,dp,npft,imdays,Tfreeze,nl,Pstd,daysec
use metvarsmod,         only : met,dm,atm
use statevars,          only : sv,ov,initstatevars,area,gridded,doy,dayl
use pftparametersmod,   only : initpftpars
use sitedatamod,        only : openfile,readfile,tsoil,date,year
use input_dataset_drivers, only : openCO2file,readCO2file
use daylengthmod,       only : daylength
use diurnaltempmod,     only : diurnaltemp
use airmass,            only : get_airmass
use pft_state_variables,only : initpftvars,veg
use soilstate_vars,     only : Fhti,hs,surf
use pftinitassign,      only : pftassign
use coszen,             only : avgZen
use soilinit,           only : soillayers,soilprep
use udunits  !not connected up fully -Feb 16 08
use calcannual,         only : annual_calcs

!----------------------

implicit none

!parameters

integer, parameter :: ns = 2 !number of soil layers in input data set

!local variables

character(100) :: jobfile
character(100) :: climatefile

integer :: j            !pixel index number
integer :: l
integer :: m,dom
integer :: runlength    !length of the simulation
integer :: yr         !year of simulation

real(dp) :: calcrad     !possible calc shortwave in case of missing measured value
real(dp) :: sunf        !sun fraction for the make-up shortwave calc.


!pointers
integer,  pointer :: gnl        !index of the lowest soil layer
real(dp), pointer :: lon        !longitude (degree)
real(dp), pointer :: lat        !latitude (degrees)
real(dp), pointer :: sdecl      !solar declination (degrees)
real(dp), pointer :: sdecl_n    !next days solar declination (degrees)
real(dp), pointer :: dswb       !downwelling shortwave (kJ m-2 d-1)
real(dp), pointer :: avg_cZ     !diurnally averaged cos of the zenith angle (rad) of the incident beam
real(dp), pointer :: Patm       !atmospheric pressure (Pa)
real(dp), pointer :: Patm30     !atmospheric pressure at reference height (30m) (Pa)
real(dp), pointer :: prec_ra         ! precipitation running sum over 30 days (mm)
real(dp), pointer, dimension(:) :: mprec        !quasi-monthly precipitation total (mm)
real(dp), pointer, dimension(:)   :: daytempv   !vector of daytime averaged temperatures
real(dp), pointer                 :: twm        !maximum monthly temperature for the year being simulated
real(dp), pointer                 :: tminm      !minimum monthly temperature for the year being simulated
real(dp), pointer, dimension(:)   :: runtemp    !holds one month of temp data for mean monthly calc.
real(dp), pointer                 :: aprec      !annual precipitation (mm)
real(dp), pointer                 :: MAT        !mean annual temperature (C)
real(dp), pointer                 :: MAT_d      !stores running total of MAT (C)
real(dp), pointer                 :: realsoilT  !used when weather station data is used. Can read in the real soil T measured at the site.

type inputvalues
  real(sp) :: lon
  real(sp) :: lat
  real(sp) :: maxd
  real(sp) :: elevation
  real(sp) :: sf2
  real(sp) :: sf075
  real(sp) :: sf05
  real(sp), dimension(ns) :: sdto  !sand (all values %)
  real(sp), dimension(ns) :: stpc  !silt
  real(sp), dimension(ns) :: clpc  !clay
  real(sp), dimension(ns) :: totc
  real(sp), dimension(ns) :: bulk
  real(sp), dimension(ns) :: cfrag
  real(sp), dimension(ns) :: tawc
end type inputvalues

type(inputvalues) :: ival

!real(dp), dimension(365)    :: tvect
real(sp), dimension(365)    :: zsno_y
real(sp), dimension(365)    :: drh_y
real(sp), dimension(17,365) :: stemp
real(sp), dimension(17,365) :: sliq_y
real(sp), dimension(17,365) :: tliq_y
real(sp), dimension(17,365) :: wliq_y
real(sp), dimension(17,365) :: wice_y

namelist /sitedata/ ival,climatefile,co2file,runlength

!------------------------------------------------------------------------------

!Begin run setup

!set pixel number (always 1 in ARVEpoint)
j = 1

!This is a point simulation so set the signpost to false
gridded = .false.

!allocate variable arrays
allocate(sv(1))
allocate(met(1,-1:1))
allocate(dm(1))
allocate(ov(1))

!allocate soils structure
allocate(soil(1))

allocate(soil(1)%sdto(5))
allocate(soil(1)%stpc(5))
allocate(soil(1)%clpc(5))
allocate(soil(1)%totc(5))
allocate(soil(1)%bulk(5))
allocate(soil(1)%cfrag(5))
allocate(soil(1)%tawc(5))

!allocate slope structure
allocate(slope(1))

!some intializations
runtemp =    0._dp  !metvars
twm     = -100._dp
tminm   =  100._dp
aprec   =    0._dp

!point pointers
lon     => sv(j)%lon
lat     => sv(j)%lat
gnl     => sv(j)%gnl
sdecl   => surf%sdecl
sdecl_n => sv(j)%sdecl_n
dswb    => met(j,0)%dswb
avg_cZ  => surf%avg_cZ
Patm    => sv(j)%Patm
Patm30  => sv(j)%Patm30
prec_ra => dm(j)%prec_ra
mprec   => dm(j)%mprec
daytempv=> atm%daytempv
twm     => atm%twm
tminm   => atm%tminm
runtemp => atm%runtemp
aprec   => atm%aprec
MAT     => atm%MAT
MAT_d   => atm%MAT_d
realsoilT => surf%realsoilT

!-----------------------
!the getarg intrinsic reads character strings off the command line

call getarg(1,jobfile)  !this is the joboptions, lon, lat, elevation, soil and slope file

open(99,file=jobfile)

read(99,nml=sitedata)

close(99)

lon = ival%lon
lat = ival%lat

! NOTE! specific options for these runs, may not be correct for yours!
!because we have only two layers in the soil dataset by ARVE expects 5 values, distribute the input values over the expected points
! Presently set up for 1 layer surface to 30 cm and 2nd layer down.
soil(1)%sdto(1:2) = ival%sdto(1)
soil(1)%sdto(3:5) = ival%sdto(2)

soil(1)%stpc(1:2) = ival%stpc(1)
soil(1)%stpc(3:5) = ival%stpc(2)

soil(1)%clpc(1:2) = ival%clpc(1)
soil(1)%clpc(3:5) = ival%clpc(2)

soil(1)%totc(1:2) = ival%totc(1) * 10._dp  !ARVE expects soil organic matter as organic carbon (g/kg) and converts it
soil(1)%totc(3:5) = ival%totc(2) * 10._dp  !to organic matter (mass percent). Multiply this input (org C mass %) by 10
                                          !so that it is correct for soil init.
soil(1)%bulk(1:2) = ival%bulk(1)
soil(1)%bulk(3:5) = ival%bulk(2)

soil(1)%cfrag(1:2) = ival%cfrag(1)
soil(1)%cfrag(3:5) = ival%cfrag(2)

soil(1)%tawc(1:2) = ival%tawc(1)
soil(1)%tawc(3:5) = ival%tawc(2)

soil(1)%maxd  = ival%maxd

slope(1)%elevation = ival%elevation
slope(1)%sf2   = ival%sf2
slope(1)%sf075 = ival%sf075
slope(1)%sf05  = ival%sf05

!----
!open the climate data input file
call openfile(climatefile)

!open the CO2 data input file
call openCO2file(CO2file)

!initialize the state variables
call initstatevars()

!find the grid cell size
sv(j)%grid_area = 1.0 !set to one square meter, because anything else doesn't really make sense for the point version

!-------------------------------------------------------------

!initializations of pft parameters and vars
call initpftpars()
call initpftvars()

!assign lm_sapl, rm_sapl, sla, sm_sapl
call pftassign

!assign soil layers and allocate the soil hydrological variables
call soillayers(j)

!intialize and calculate soil properties
call soilprep(j)

!allocate dimension of soil output variables (dependent on the soil depth)
!GND FLAG, maybe do this in the future, leave for now JM 14.06.2010

!start the daylength buffer
doy = 0  !day zero
call daylength(doy,sunf,calcrad)
dayl(1) = dayl(2)
sdecl = sdecl_n

!elevation is used to get the atmospheric pressure
Patm = Pstd * (1._dp - slope(j)%elevation / 44308)**5.2568  !(Diehl 1925 (from Lloyd and Farquhar 1994))

!also calculate the reference height pressure (30m above ground)
Patm30 = Pstd * (1._dp - (slope(j)%elevation + 30._dp) / 44308)**5.2568  !(Diehl 1925 (from Lloyd and Farquhar 1994))

!print some output on the soils
!do l = 1,nl
! write(*,'(i4,6f12.4)')l,sv(j)%sand(l),sv(j)%silt(l),sv(j)%clay(l),sv(j)%rock(l)*100.,sv(j)%sorg(l),sv(j)%bulk(l)
!end do
!write(*,*)'soildep',soil(j)%maxd


!OUTPUT

!open text files for diagnostic output
open(05,file='diag_test.dat',status='unknown')
open(15,file='diag_unknown.dat',status='unknown')
open(25,file='diag_soiltemp.dat',status='unknown')
open(35,file='diag_soilmoisture.dat',status='unknown')
!open(45,file='diag_heatest.dat',status='unknown')
open(55,file='diag_latent.dat',status='unknown')


!     !open the diagnostic netcdf output file
!     status = nf90_open('arve_diag_soil.nc',nf90_write,ncid)
!     if (status /= nf90_noerr) call handle_err(status)
!
!     !write the values for soil depth and layer thickness
!
!     status = nf90_inq_varid(ncid,'depth',varid)
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_put_var(ncid,varid,-sv(j)%zipos(1:nl))
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_inq_varid(ncid,'dz',varid)
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_put_var(ncid,varid,sv(j)%dz(1:nl))
!     if (status /= nf90_noerr) call handle_err(status)

MAT = 0._dp
MAT_d = 0._dp

!-------------------------------------------------------------
!start the year loop

do yr = 1,runlength

  write(0,*)'year:',year,j

  !get this years CO2 values
  call readCO2file

  !run ARVE loop for one year
  doy = 1
  do m = 1,12
    do dom = 1,imdays(m)

      !get one day of climate data (Note: most station data records downwelling shortwave as only a total value.
      !in ARVE-point this is inputed as direct downwelling shortwave. In simplesurfrad, it is partitioned into
      !diffuse and direct components.If your dataset is already partitioned. This will need to be accounted for. JM Oct 08.
      call readfile()

      !calculate daylength
      call daylength(doy+1,sunf,calcrad)

      !calculate day and night temperature
      call diurnaltemp(j)

      !find the diurnally averaged cosine of the solar zenith angle
      call avgZen(j)

      !add day to vector of daytime averaged temperatures
      daytempv(doy) = met(j,0)%temp(1)

      !if measured value is missing and the sun should rise, then use a calculated dsw value
      if (dswb < 0._dp .and. avg_cz > 0._dp) then
        dswb = calcrad * 10000. !use calculated rad, convert kj cm-2 d-1  to kj m-2 day-1
      else if (dswb < 0._dp .and. avg_cz == 0._dp) then !sun is below horizon and dswb is set to less than zero, make it 0.
        dswb = 0._dp
      end if

      !calculate elevation-corrected daily mean and instantaneous air mass
      call get_airmass(j)

      !calculate daily biophysics and biogeochemistry
      call calcdaily(j,yr)

      !increment the running precip total and mean annual temperature
      MAT_d = MAT_d + ((dayl(1) * met(j,0)%temp(1)) + ((daysec - dayl(1)) * met(j,0)%temp(2))) / daysec
      aprec = aprec + met(1,0)%prec
      mprec(:) = eoshift(mprec(:),1,met(1,0)%prec)
      prec_ra = sum(mprec(:)) !find the running 30 day sum

      !move the 'next day' values to present day in preparation for another daily cycle
      dayl(1) = dayl(2)
      sdecl = sdecl_n

      !accumulate daily soil stats into annual array

      drh_y(doy)  = veg%drh(1)
      zsno_y(doy) = sv(j)%zsno

      stemp(:,doy) = sv(j)%Tsoil(1:nl) - Tfreeze
      sliq_y(:,doy)  = sv(j)%Tliq(1:nl) / sv(j)%Tsat(1:nl) !GND -- likely change this to gnl
      tliq_y(:,doy)  = sv(j)%Tliq(1:nl)        !GND -- likely change this to gnl
      wliq_y(:,doy)  = sv(j)%Wliq(1:nl)        !GND -- likely change this to gnl
      wice_y(:,doy)  = sv(j)%Wice(1:nl)        !GND -- likely change this to gnl

   do l = 1, nl
     write(35,'(2i5,8f10.4)') yr, doy, -sv(j)%zpos(l), (sv(j)%Tliq(l)/sv(j)%Tsat(l)), sv(j)%Tsat(l), sv(j)%Tliq(l), sv(j)%Tice(l), sv(j)%Tsat(l)-sv(j)%Tice(l), sv(j)%Tsoil(l) - Tfreeze,-sv(j)%zw
   end do

        call flush(35)

   !write(*,*)'watertable',sv(j)%zw

      !increment day counter
      doy = doy + 1

    end do !day in month loop

    !min/max monthly mean temp accounting
    !FLAG THESE ARE NOT DONE PROPERLY!!! THESE MUST BE CHANGED TO USE ARVE-POINT! jm 03.05.2011
    twm = max(twm,(runtemp(1)/dom))
    tminm = min(tminm,(runtemp(2)/dom))
    runtemp = 0._dp

  end do !month loop

  !calculate annual values (allocation, mortality, etc.)
  call annual_calcs(j,yr)

   !rewind(05)
  !rewind(15)
  !rewind(25)
  !rewind(35)
  !rewind(45)

  !--------------------
  !write out the year of the model run's worth of data to the soil diagnostic output file
!
!     dsrt = 1 + 365 * (yr - 1)
!     dcnt = 365
!
!     doy = 365 * (year - 2)
!
!     do i = 1,365
!        tvect(i) = doy + i - 1
!     end do
!
!     status = nf90_inq_varid(ncid,'time',varid)
!     status = nf90_put_var(ncid,varid,tvect,start=[dsrt])
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_inq_varid(ncid,'zsno',varid)
!     status = nf90_put_var(ncid,varid,zsno_y,start=[dsrt])
!     if (status /= nf90_noerr) call handle_err(status)
!
!    ! status = nf90_inq_varid(ncid,'wliq',varid)
!    ! status = nf90_put_var(ncid,varid,wliq_y,start=[1,dsrt])
!    ! if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_inq_varid(ncid,'stemp',varid)
!     status = nf90_put_var(ncid,varid,stemp,start=[1,dsrt])
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_inq_varid(ncid,'sliq',varid)
!     status = nf90_put_var(ncid,varid,sliq_y,start=[1,dsrt])
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_inq_varid(ncid,'tliq',varid)
!     status = nf90_put_var(ncid,varid,tliq_y,start=[1,dsrt])
!     if (status /= nf90_noerr) call handle_err(status)
!

  !--------------------

end do   !end of the year loop


!close(15)
close(25)
close(35)
!close(45)
close(55)

!status = nf90_close(ncid)
!if (status /= nf90_noerr) call handle_err(status)


end program arvepoint
