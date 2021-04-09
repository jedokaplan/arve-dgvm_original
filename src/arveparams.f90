module arveparams

implicit none

!type

integer, parameter :: sp = selected_real_kind(6)
integer, parameter :: dp = selected_real_kind(13)  


!constants

real(sp), dimension(12), parameter :: mdays  = [ 31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
integer,  dimension(12), parameter :: imdays = [ 31, 28, 31 ,30, 31, 30, 31, 31, 30, 31, 30, 31 ]
integer,  dimension(12), parameter :: midday = [ 16, 44, 75,105,136,166,197,228,258,289,319,350 ]
real(dp), parameter :: dayspy = 365._dp       !number of days in one year
real(dp), parameter :: daysec    = 86400._dp  !number of seconds in 24hrs
real(dp), parameter :: hr1  =  3600._dp       !one hour in seconds
real(dp), parameter :: hr23 = 82800._dp       !23 hours in seconds

!parameters

integer, parameter  :: npft     =  9             !number of pfts
integer, parameter  :: ncvar    =  3             !1:total carbon, 2:13C, 3:14C
integer, parameter  :: climbuf  = 20             !years in the climate buffer
integer, parameter  :: band     = 2              !1 is visible, 2 is near infrared
integer, parameter  :: numtimestep = 2           !number of major timesteps in 24 hours (we have smaller iterations in some areas)
integer, parameter  :: nts = 14                  !number of days to take a running average.MUST BE AN EVEN NUMBER
real(dp), parameter :: shortstep =  3600._dp    !the length of the short soil routine timestep. (seconds)
real(dp), parameter :: minplant  =     0.05_dp  !threshold plant and stem area for energy calculations
real(dp), parameter :: g_m       =     3.26_dp  !minimum stomatal conductance (mm s-1)(Gerten et al. 2004 J Hydrology)
real(dp), parameter :: aq_thick  = 25000._dp    !thickness of the aquifer in mm
real(dp), parameter :: allom2 =  40._dp         !original LPJ value (Sitch et al. 2003; 40)
real(dp), parameter :: allom3 =   0.67_dp       !value from Ben Smith (LPJ-GUESS) original LPJ value (Sitch et al. 2003; 0.5)
real(dp), parameter :: bet_com = 0.192_dp       !beta allometry relation for all datasets relating leaf mass to non-photosynthetic mass
                                                        !Niklas & Enquist, PNAS 2001 pg. 2923, fig 3A
real(dp), parameter :: alph_com = 0.76_dp       !alpha allometry relation for all data sets, Niklas & Enquist, PNAS 2001
real(dp), parameter :: k5 = 34.64               !k5 constant from Niklas & Spatz PNAS 2004 Fig 1
real(dp), parameter :: k6 = 0.475               !k6 constant from Niklas & Spatz PNAS 2004 Fig 1
real(dp), parameter :: frac_cr_tot = 0.2        !percent of the above ground biomass that goes to coarse roots (see refs in Ryan et al.
                                                ! Adv Ecol Res 1997
real(dp), parameter :: reinickerp = 1.6_dp      !Reinicke's rule exponent (LPJ original formulation)
real(dp), parameter :: reprod_cost = 0.1        !proportion of NPP lost to reproduction (Harper 1977)

!parameters used in the soil physics and surface energy balance routines

integer, parameter  :: nl     = 17         !number of soil layers
integer, parameter  :: snomax = 5          !maximum number of snow layers
integer, parameter  :: ns     = 1 - snomax !index of top snow layer
real(dp), parameter :: zsno_max = 4._dp   !maximum snow depth (m)

!fixed CLM parameters

real(dp), parameter :: pi    = 3.14159265358979323846_dp !26433 83279 50288 41971 69399 37510 (unitless)
real(dp), parameter :: twopi = 2._dp * pi
real(dp), parameter :: grav  = 9.80616_dp       !Gravitational acceleration  (m s-2)
real(dp), parameter :: Pstd  = 101325._dp       !Standard pressure           (Pa)
real(dp), parameter :: sigsb = 5.67e-8          !Stefan-Boltzmann constant   (W m-2 K-4)
real(dp), parameter :: kbz   = 1.38065e-23      !Boltzmann constant          (J K-1 molecule-1)
real(dp), parameter :: avo   = 6.0221415e26     !Avogadro's number           (molecule kmol-1)
real(dp), parameter :: Rgas  = avo * kbz        !Universal gas constant      (J K-1 kmol-1)
real(dp), parameter :: MWda  = 28.966_dp        !Molecular weight of dry air (kg kmol-1)
real(dp), parameter :: Rda   = Rgas / MWda      !Dry air gas constant        (J K-1 kg-1)
real(dp), parameter :: MWwv  = 18.016_dp        !Mol. weight of water vapor  (kg kmol-1)
real(dp), parameter :: Rwv   = Rgas / MWwv      !Water vapor gas constant    (J K-1 kg-1)
real(dp), parameter :: vkarm = 0.4_dp           !von Karman constant (unitless)
real(dp), parameter :: visco = 1.5e-5           !kinematic viscosity of air (m2 s-1)

real(dp), parameter :: Tfreeze = 273.15_dp      !freezing temperature of freshwater (K) 

real(dp), parameter :: pliq  = 1000._dp         !density of water (kg m-3)
real(dp), parameter :: pice  =  917._dp         !density of ice (kg m-3)

real(dp), parameter :: Cair   =     1.25e3!1.00464e3  !Heat capacity of dry air (J kg-1 K-1) (CLM parameterization)
real(dp), parameter :: Cliq   =     4.18800e3   !Heat capacity of liquid water  (J kg-1 K-1)
real(dp), parameter :: Cice   =     2.11727e3   !Heat capacity of ice (typical) (J kg-1 K-1)

real(dp), parameter :: lvap  = 2.501e6          !Latent heat of vaporization (J kg-1)
real(dp), parameter :: Lf    = 3.337e5          !Water latent heat of fusion (J kg-1)
real(dp), parameter :: lsub  = lvap + Lf        !Latent heat of sublimation  (J kg-1)

real(dp), parameter :: Kliq   =     0.57_dp     !Thermal conductivity of liquid water (typical) (W m-1 K-1)
real(dp), parameter :: Kice   =     2.2_dp      !Thermal conductivity of ice (typical) (W m-1 K-1)
real(dp), parameter :: Kair   =     0.0243_dp   !Thermal conductivity of dry air (W m-1 K-1)

!other fixed parameters and conversion factors

real(dp), parameter :: g   =  6.6742e-8         !the Earth's gravitational constant (cm3 g-1 s-2)
real(dp), parameter :: c12mass = 12.0107        !atomic mass of carbon (amu)
real(dp), parameter :: n14mass = 14.0067        !atomic mass of nitrogen (amu)
real(dp), parameter :: o16mass = 15.9994        !atomic mass of oxygen (amu)
real(dp), parameter :: CO2mass = 40.0208        !atomic mass of CO2 (amu)

real(dp), parameter :: radius  = 6378.137_dp     !km, WGS-84 spherical approximation
real(dp), parameter :: solarc = 1360.            !Solar constant (1360 W/m2)
real(dp), parameter :: d2r    = pi / 180.        !conversion factor from degrees to radians
real(dp), parameter :: a2s    = 1. / (300. * pi) !conversion factor from solar angular to seconds
real(dp), parameter :: umol2g  = 1.e-6 * c12mass !convert umol CO2 to grams C

!derived CLM parameters

real(dp), parameter :: Timp = 0.05_dp           !volumetric liquid water content at which the soil is impermeable (fraction)
real(dp), parameter :: Pmax =-1.e8              !maximum allowed value for soil matric potential (Pa)
real(dp), parameter :: psimax = -1.5e5          !wilting point potential of leaves (CLM 8.11)
real(dp), parameter :: Esoil  = 0.96_dp         !thermal emissivity of bare soil (fraction)
real(dp), parameter :: Ewater = 0.96_dp         !thermal emissivity of water (fraction)
real(dp), parameter :: Esnow  = 0.97_dp         !thermal emissivity of snow (fraction)
real(dp), parameter :: minsno = 0.05_dp         !minimum snow depth to consider burial of vegetation (m)

real(dp), parameter :: Tstune = 0.34_dp         !Tuning factor to turn first layer T into surface T (CLM parameterization, pg 88).
real(dp), parameter :: CNfac  = 0.5_dp          !Crank Nicholson factor between 0 and 1
real(dp), parameter :: oneminusCNfac = 1._dp - CNfac

real(dp), parameter :: pdrysnow =  50.00000_dp  !density of dry snow falling at under -15 deg C (kg m-3)
real(dp), parameter :: pwetsnow = 169.15775_dp  !density of wet snow falling at over 2 deg C (kg m-3) (50._dp + 1.7_dp * 17._dp**1.5, CLM eqn 7.18)
real(dp), parameter :: Tc = 2.5_dp              !critical threshold temperature separating rain from snow (deg C)

real(dp), parameter :: lb = 1.e-5               !baseflow parameter (mm s-1) CLM eqn. 7.117
real(dp), parameter :: kd = 0.04_dp             !saturated soil hydraulic cond. contributing to baseflow (mm s-1) CLM eqn. 7.118

real(dp), parameter :: Wpond = 10._dp           !max. quantity of water for ponding (kg m-2)
real(dp), parameter :: wpondflat = 100._dp      !max water for ponding on flat land (mm) !FLAG not used JM 07.04.2011

real(dp), parameter :: alpha = 3._dp            !adjustable scale dependent parameter (aquifer, soilwaterflux) CLM 4.0 p.155

real(dp), parameter :: fdecay = 2.5             !decay factor (m-1)(aquifer and infiltration)

real(dp), parameter :: Smin  = 0.033_dp         !irreducible water saturation of snow (fraction)

real(dp), parameter :: z0mg  = 0.01_dp          !momentum roughness length for bare soil (m, CLM eqn 3.49)
        
real(dp), parameter :: rh  = 67._dp             !arbitrary relative humidity

!other parameters (Soil Properties)
real(dp), parameter :: OMorgC = 0.5800464               !g org C =sorg/ 1.724 (conversion factor from Nelson & Sommers 1996)
real(dp), parameter :: peatlim= 25._dp                  !percent organic matter in soil where the soil is treated as peat
real(dp), parameter :: soilbulk = 2650._dp              !soil bulk density (kg m-3)
real(dp), parameter :: omd = 1.3                        !(g cm-3) particle density of org matter (from Hillel 1982)
real(dp), parameter :: sandd = 2.128_dp                 ! (g cm-3) particle density of sand(0.8 x quartz (2.66))
real(dp), parameter :: siltd = 2.35_dp                  !(average soil density...   FLAG!!!must look up a better value for silt!!!
real(dp), parameter :: clayd = 2.65_dp                  ! (g cm-3) particle density of other minerals
real(dp), parameter :: rockd = 2.70_dp                  ! (g cm-3) particle density of coarse frags
real(dp), parameter :: sandquartz = 0.65                ! fraction of quartz in sand on average (geology.uprm.edu/Morelock/terrigenous.html)
real(dp), parameter :: rockquartz = 0.25                ! fraction of quartz in rock on average (estimate)
real(dp), parameter :: bedrock_tc = 3.7                 ! bedrock thermal conductivity (W m-1 K-1)
real(dp), parameter :: Corg =  2.51e6                   !(1e6)Heat capacity of soil organic matter (typical) (J m-3 K-1)
real(dp), parameter :: Csand = 2.128e6                  !(1e6) Volumetric heat capacity of sand (J m-3 K-1)
real(dp), parameter :: Csilt = 2.439e6                  !(1e6) Volumetric heat capacity of silt (J m-3 K-1) from Ren et al SSSAJ 67 (6) 2003
real(dp), parameter :: Cclay = 2.8385e6                 !(1e6) Volumetric heat capacity of clay (J m-3 K-1)
real(dp), parameter :: Crock = 2.01e6                   !(1e6) Volumetric heat capacity of quartz,other minerals (J m-3 K-1)
real(dp), parameter :: Kom = 0.25_dp                    !thermal cond of organic matter (W m-1 K-1)
real(dp), parameter :: Kom_dry = 0.05d0                 !thermal cond of dry organic matter (W m-1 K-1) (Farouki, 1981)
real(dp), parameter :: Kquartz = 8.0_dp                 !thermal cond. of quartz (W m-1 K-1) (applicable across more temps
                                                        ! than previous value)
real(dp), parameter :: Kmineral = 2.5_dp                !thermal cond. of other minerals (clay) (W m-1 K-1)

!-----------------------------------------------------------------------------------------------

end module arveparams

!note from wikipedia: air conductivity 0.024, snow (typical) 0.11, soil (typical) 0.17-1.13
!1.9e6 J m-3 K-1 = 0.45 cal cm-3 degC-1
!3.2217e-6 W m-1 K-1 = 7.694902e-9 cal cm-1 sec-1 degC-1
