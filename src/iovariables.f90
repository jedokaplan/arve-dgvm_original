module iovariables

!This module stores the input/output variables. It is a bit messy at the moment.
! Needs more comments on what is what!

use arveparams, only : sp,dp,ncvar

implicit none

! Input from the job options file

!Input/Output Files
character(180) :: climate_path1         !spin up climate path
character(180) :: climate_path2         !transient climate path
character(180) :: soil_path             !soil textures dataset path
character(180) :: slope_path            !ground slope dataset path
character(180) :: fmax_path              !maximum fractional saturated area dataset path
character(180) :: save_spinup_file      !file name where the spinup-state values are written
character(180) :: save_transient_file   !file name where the transient-state values are written
character(180) :: landfile              !fractional land file, not strictly required
character(180) :: icefile               !ice file, not strictly required
character(180) :: co2file               !transient CO2 file path
character(180) :: saved_model_state_file!model state file from previous run, used to restart at same conditions
character(80)  :: outputfile             !output NetCDF file, read in from command line

!Switches
logical :: save_spinup_state            !true if the model state should be saved after the spinup
logical :: save_transient_state         !true if the model state should be saved after the transient
logical :: dosoilm_assim                ! BIC, true if run in soil moisture assimilation mode
logical :: dospinup                     !true if you want to run a spinup, otherwise it goes right to transient
logical :: cycle_spinup                 !cycle over the first cycle_years years of the climate file
logical :: dotransient                  !true if the transient simulation is desired, other only the spinup
logical :: fromsavedstate               !true if transient simulation uses a saved equilbrium state, else it uses state from spinup or bare start
logical :: doweather                    !true if the weather generator should be used, otherwise pseudo daily smoothed values
logical :: dovegetation                 !true if vegetation should be simulated, otherwise bare ground simulation
logical :: do_allocation_daily          !false if allocation is only annually using the LPJ allocation scheme. true gives daily allocation
logical :: do_CLM_photosynthesis        !true does CLM 4.0 photosynthesis, false uses LPJ formulation

!Variables
real(dp) :: spinupco2_12C               !spin up 12C- CO2 mixing ratio (ppm)
real(dp) :: spinupco2_13C               !spin up 13C- CO2 mixing ratio (ppm)
real(dp) :: spinupco2_14C               !spin up 14C- CO2 mixing ratio (ppm)
integer :: spinupyears                  !years to run spin up   
integer :: transientyears               !years to run transient
integer :: cycle_years                  !number of years to cycle over from start of climate dataset
integer :: climateyears                 !number of years in the climate dataset
integer :: climatemonths                !number of months in the climate dataset

!---------------------------------------------------

! Variables to be read in from the different input files (climate, soil, slope)
integer, parameter :: niv = 5
character(3), dimension(niv), parameter :: vname = [ 'cld','dtr','pre','tmp','wet' ]

integer, parameter :: nsiv = 7
character(5), dimension(nsiv), parameter :: vsname = [ 'sdto ','stpc ','clpc ','totc ','bulk ','cfrag','tawc ' ]

integer, parameter :: siv = 4
character(9), dimension(siv), parameter :: sname = [ 'elevation','sf2      ','sf075    ','sf05     ' ]

real(sp) :: fill_value = -9999.         !value given for empty fields in the NetCDF output files

integer :: ncid
integer :: slopefid
integer :: soilfid
integer :: co2fid
integer :: climfid
integer :: fmaxfid

integer :: inputlatlen
integer :: inputlonlen
integer :: inputclimlen
integer :: inputsoillen

integer :: timebuflen  !the number of months of climate data to store in the input buffer

!maxibufsize sets the maximum number of gridcells the input and output buffers can hold.
!A maximum sized box of 259200 pixels (global 0.5 deg resolution), has some 67420 valid
!(non-water non-ice pixels). Global soils data and 12 months of climate data results in
!a memory requirement of ca. 13 MB

!integer, parameter :: maxibufsize = 300000

real(dp), dimension(4) :: bounds

!----------------
!CLIMATE

type inputbuffer
  !Monthly values read in for arve_grid
  real(dp), allocatable, dimension(:) :: cld  !cloud cover (%)
  real(dp), allocatable, dimension(:) :: dtr  !diurnal temp range(C)
  real(dp), allocatable, dimension(:) :: pre  !total precip. (mm)
  real(dp), allocatable, dimension(:) :: tmp  !mean temp (C)
  real(dp), allocatable, dimension(:) :: wet  !Number of days in month with precipitation
  !note min and max temp have been removed as these are just a function of the primary variables tmp and dtr
end type inputbuffer

type(inputbuffer), allocatable, save, target, dimension(:) :: ibuf

!----------
!SOIL

type soiltype

  logical :: water                      !I don't think this is used JM 22.03.2011

  real(dp), allocatable, dimension(:) :: sdto  !sand (mass %)
  real(dp), allocatable, dimension(:) :: stpc  !silt (mass %)
  real(dp), allocatable, dimension(:) :: clpc  !clay(mass %)
  real(dp), allocatable, dimension(:) :: totc  !organic carbon content (g kg-1)
  real(dp), allocatable, dimension(:) :: bulk  !bulk density (kg dm-3)
  real(dp), allocatable, dimension(:) :: cfrag !coarse fragments > 2mm (volume %)
  real(dp), allocatable, dimension(:) :: tawc  !available h2o capacity (cm m-1)
  real(dp) :: maxd  !soil max depth (cm) (up to 250cm)
  
  !Maximum saturated fraction from topographic index calculations
  !this comes from a different dataset than that of the soil variables
  real(dp) :: fmax            !gridcell maximum saturated fraction

end type soiltype

type(soiltype),    allocatable, save, target, dimension(:) :: soil

!----------------
!SLOPE

type slopetype

  real(dp)  :: elevation     !gridcell mean elevation above mean sea level (m)
  real(dp)  :: sf2           !fraction of the grid that is less 0.002 m m-1 gradient
  real(dp)  :: sf075         !fraction of the grid that is less 0.00075 m m-1 gradient
  real(dp)  :: sf05          !fraction of the grid that is less 0.0005 m m-1 gradient

end type slopetype

type(slopetype),   allocatable, save, target, dimension(:) :: slope

!----------------


integer :: index         !the index of the current gridcell in the input buffer

real(dp), allocatable, dimension(:) :: lonvect
real(dp), allocatable, dimension(:) :: latvect

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty
integer :: endx
integer :: endy

type var_atts
  real(sp) :: scale_factor
  real(sp) :: add_offset
end type var_atts

type(var_atts), dimension(niv) :: va

end module iovariables
