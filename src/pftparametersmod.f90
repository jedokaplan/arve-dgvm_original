module pftparametersmod

! Vegetation is represented by a combination of the following plant functional
! types (PFTs)

!       1. tropical broadleaved evergreen tree 
!       2. tropical broadleaved raingreen tree
!       3. temperate needleleaved evergreen tree 
!       4. temperate broadleaved evergreen tree 
!       5. temperate broadleaved summergreen tree
!       6. boreal needleleaved evergreen tree 
!       7. boreal needleleaved summergreen tree
!       8. C3 perennial grass
!       9. C4 perennial grass

use arveparams, only : dp,npft

implicit none

public :: initpftpars,prm

!PFT PARAMETERS

type pftparams

  ! basic plant information
  logical  :: c4                     !plants with C4 (true) or C3 (false) photosynthetic pathway
  logical  :: tree                   !true is pft is a tree
  real(dp) :: maxcrowna              !tree maximum crown area (m2)  
  integer  :: leaftype               !leaf type: broadleaved (1), needleleaved (2) or grass (3)
  integer  :: phentype               !phenology type: evergreen (1), summergreen (2), raingreen (3),any type (4)
  real(dp) :: wooddens               !density of wood (sapwood and heartwood) (g m-3)
  
  ! allometry   
  real(dp) :: latosa                 !leaf area to sapwood area (m2 m-2) (Ben Smith LPJ-GUESS)  
  real(dp) :: allom1                 !Packing constrain from Reinicke's rule exponent. Values from Ben Smith (LPJ-GUESS) 
  real(dp) :: alph_lmrm              !scaling exponent for rm_ind to lm_ind (Enquist & Niklas 2002 table 1)  !FLAG duplicate of those in daily_alloc
  real(dp) :: beta_lmrm              !allometric constant for rm_ind to lm_ind (Enquist & Niklas 2002 table 1)  !FLAG duplicate of those in daily_alloc
  real(dp) :: alph_stemlm            !scaling exponent for lm_ind to stem mass (Enquist & Niklas 2002 table 1)  !FLAG NOT USED?
  real(dp) :: beta_stemlm            !allometric constant for lm_ind to stem mass (Enquist & Niklas 2002 table 1)  !FLAG NOT USED?
  real(dp) :: sapl_mass              !sapling mass (gC)
  real(dp) :: sapl_NSC_frac          !sapling non-structural carbon as a fraction of total sapling mass 
  real(dp) :: NSC_opt                !allometric plant non-structural carbon multiple of maximal leaf mass 

  ! photosynthesis and respiration
  real(dp) :: gmin                   !canopy conductance component (gmin, mm/s) not associated with photosynthesis (Haxeltine & Prentice 1996, Table 4)
  real(dp) :: respcoef               !maintenance respiration coefficient
  real(dp) :: tlim_bot               !low temperature limit for CO2 uptake
  real(dp) :: topt_bot               !lower range of temperature optimum for photosynthesis
  real(dp) :: topt_top               !upper range of temperature optimum for photosynthesis
  real(dp) :: tlim_top               !high temperature limit for CO2 unptake
  real(dp) :: lambdam                !lambda ratio of optimal stomatal opening
  real(dp) :: psi_c                  !soil water potential (mm) whtn stomata are fully closed 
  real(dp) :: psi_o                  !soil water potential (mm) whtn stomata are fully open
  real(dp) :: Nmax                   !maximum foliar N content (mg/g) (Haxeltine & Prentice 1996a, Fig 4)
  real(dp) :: leafcn                 !leaf C:N mass ratio
  real(dp) :: sapwoodcn              !sapwood C:N mass ratio
  real(dp) :: rootcn                 !root C:N mass ratio
  real(dp) :: F_lnr                  !fraction of leaf N in rubisco
  real(dp) :: f_N                    !nitrogen availability factor (CLM 4.0 Table 8.1)

  ! turnover
  real(dp) :: leafturnover           !leaf turnover period (years to completely turnover the leaf mass)
  real(dp) :: leaflongevity          !leaf longevity (years)
  real(dp) :: sapturnover            !sapwood turnover period (sapwood converted to heartwood) (years)
  real(dp) :: frootQ10               !fine root turnover 'Q10' relation with temperature (deg C) (Gill & Jackson 2000)

  ! phenology
  real(dp) :: leaf_thres             !fraction of maximal leaf mass at which the plant switches to full leafout allocation pattern
  real(dp) :: phenramp       !summergreen phenology ramp, GDD5 requirement to grow full leaf canopy (annual C allocation scheme)
  real(dp) :: drydrop        !water scalar value at which leaves shed by drought deciduous PFT (annual C allocation scheme)
  
  ! fire
  real(dp) :: flam                   !flammability threshold
  real(dp) :: fireres                !fire resistance index
  real(dp) :: CH4ef                  !Emission factor for CH4 production from fire (g CH4/ kg C per PFT) (Andreae & Merlet GBC 2001)          
  real(dp) :: COef                   !Emission factor for CO production from fire (g / kg C per PFT) (Andreae & Merlet GBC 2001)              
  real(dp) :: VOCef                  !Emission factor for VOC production from fire (g/ kg C per PFT) (Andreae & Merlet GBC 2001)              
  real(dp) :: TPMef                  !Emission factor for total particulate matter production from fire (g/ kgC)(Andreae & Merlet GBC 2001)   
  real(dp) :: NOxef                  !Emission factor for NOx production from fire (g/ kg C per PFT) (Andreae & Merlet GBC 2001)              
  
  ! Bioclimatic limits
  real(dp) :: mintcm                 !minimum coldest monthly mean temperature
  real(dp) :: maxtcm                 !maximum coldest monthly mean temperature
  real(dp) :: gddmin                 !minimum growing degree days (at or above 5 deg C)
  real(dp) :: maxtwm                 !upper limit of temperature of the warmest month
  real(dp) :: gddbase                !base temperature for growing degree days

  ! Plant optical properties from surfaceheatflux (table 3.1 CLM Tech note)
  real(dp) :: lfang                  !leaf angle for pft                        
  real(dp) :: alpha_s1               !stem reflectances (visible)               
  real(dp) :: alpha_s2               !stem reflectances (nir)                   
  real(dp) :: alpha_l1               !leaf reflectances (visible)               
  real(dp) :: alpha_l2               !leaf reflectances (nir)                   
  real(dp) :: tau_l1                 !leaf transmittance (visible)              
  real(dp) :: tau_l2                 !leaf transmittance (nir)                  
  real(dp) :: tau_s1                 !stem transmittance (visible)              
  real(dp) :: tau_s2                 !stem transmittance (nir)                  
  real(dp) :: fi1                    !relates leaf angles to projected area !NOTE any changes to fi1, means f2 and mu_bar must be recalculated!    
  real(dp) :: fi2                    !  !NOTE any changes to fi1, means f2 and mu_bar must be recalculated!
  real(dp) :: mu_bar                 !average inverse optical depth per unit leaf and stem area !see NOTE above          

  ! TWO-LEAF parameters (from CLM 3.5 Thornton and Zimmerman 2007)
  real(dp) :: m                      !linear slope coefficient (dSLA/dLAI, projected area basis [m^2/gC])
  real(dp) :: sla_top                !SLA at the top of the canopy (m2 / gC)

  ! Resistance
  real(dp) :: Rzmom                  !Ratio of momentum roughness length to canopy top height (CLM Tech Note)
  
  ! Roots.               Values are from Arora & Boer, Earth Interactions, Vol 7 paper no. 6. 2003
  real(dp) :: abar                  !parameter representing mean root distribution profile
  real(dp) :: Bbar                  !average standing root biomass (kg m-2)
  real(dp) :: littleb               !parameter representing variable root distribution profile (=abar * Bbar exp alpha)  

end type pftparams

type(pftparams), dimension(npft), target :: prm

contains

!-------------------------------------------------------------------------------------------------------

subroutine initpftpars()

implicit none

!initialize the values

! Basic plant information

prm(1:8)%c4 =  .false.
prm(9)%c4   =  .true.

prm(1:7)%tree = .true.
prm(8:9)%tree = .false.

prm%maxcrowna = [ 15.0,15.0,15.0,15.0,15.0,15.0,15.0, 0.0, 0.0 ]

prm%leaftype = [ 1,1,2,1,1,2,2,3,3 ]

prm%phentype = [ 1,3,1,1,2,1,2,4,4 ]

prm%wooddens = [ 3.0e5, 3.0e5, 2.0e5, 3.0e5, 3.0e5, 2.0e5, 2.0e5, 1., 1. ]
!prm(1:9)%wooddens = 2.0e5

! Allometry

!prm(1:9)%latosa = 8000.  !LPJ orig value
prm%latosa = [ 4000., 4000., 3000., 4000., 4000., 3000., 3000., 0., 0.] !changed values to be consistent with lit estimates. JM 02.09.2010
!prm%latosa = [ 8000., 8000., 6000., 8000., 8000., 6000., 6000., 0., 0.]       !see Ben Smith's compendium of LPJ-GUESS parameters

!prm%allom1 = [ 250.0,250.0,150.0,250.0,250.0,150.0,250.0, 0.0, 0.0 ]  !old.
prm%allom1 = [ 250.0,250.0,100.0,250.0,250.0,100.0,100.0, 0.0, 0.0 ]

prm%alph_lmrm = [ 0.76, 0.76, 0.86, 0.76, 0.76, 0.86, 0.76, 0., 0. ]

prm%beta_lmrm = [ 0.3, 0.3, 0.76, 0.3, 0.3, 0.76, 0.3, 0., 0. ]

prm%alph_stemlm = [ 0.73, 0.73, 0.78, 0.73, 0.73, 0.78, 0.73, 0., 0. ]

prm%beta_stemlm = [ 0.13, 0.13, 0.34, 0.13, 0.13, 0.34, 0.13, 0., 0. ]

!prm%sapl_mass = [ 3.0,3.0,3.0,3.0,3.0,3.0,3.0,0.07,0.07 ]  
prm%sapl_mass = [ 6.0,6.0,6.0,6.0,6.0,6.0,6.0,0.07,0.07 ]  

prm%sapl_NSC_frac = [ 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15 ]  

!prm%NSC_opt = [ 0.7, 2.0, 0.7, 0.7, 4.8, 0.7, 4.8, 0.10, 0.10 ]
prm%NSC_opt = [ 0.7, 2.0, 0.7, 0.7, 4.8, 0.7, 4.8, 2.0, 2.0 ]

! Photosynthesis and respiration

prm%gmin = [ 0.5,0.5,0.3,0.5,0.5,0.3,0.3,0.5,0.5 ]

prm%respcoef = [ 0.1, 0.1, 1.0,1.0,1.0,1.2,1.2,0.7,0.2 ]  !Rita Wania's values
!prm%respcoef = [ 0.2, 0.2, 1.2,1.2,1.2,1.2,1.2,1.2,1.2 ]  !Old LPJ values

prm%tlim_bot = [  2.0, 2.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 6.0 ]

prm%topt_bot = [ 25.0,25.0,20.0,20.0,20.0,15.0,15.0,10.0,13.0 ]

prm%topt_top = [ 30.0,30.0,30.0,30.0,25.0,25.0,25.0,30.0,45.0 ]

prm%tlim_top = [ 55.0,55.0,42.0,42.0,38.0,38.0,38.0,45.0,55.0 ]

prm%lambdam = [ 0.9,0.9,0.9,0.8,0.8,0.8,0.9,0.65,0.4 ]

prm%psi_o = [ -66000.,-35000.,-66000.,-66000.,-35000.,-66000.,-35000.,-74000.,-74000. ]

prm%psi_c = [ -255000.,-224000.,-255000.,-255000.,-224000.,-255000.,-224000.,-275000.,-275000. ]

prm%Nmax = [ 100.0,100.0,100.0,100.0,120.0,100.0,100.0,100.0,100.0 ]

!prm%leafcn = [ 29.0,29.0,29.0,29.0,29.0,29.0,29.0,29.0,29.0 ]  !LPJ values
prm%leafcn = [ 30.0,25.0,35.0,30.0,25.0,40.0,25.0,25.0,25.0 ] !CLM 4.0 values

prm%sapwoodcn = [ 330.0,330.0,330.0,330.0,330.0,330.0,330.0,0.0,0.0 ]

prm%rootcn = [ 29.0,29.0,29.0,29.0,29.0,29.0,29.0,29.0,29.0 ]

prm%F_lnr = [ 0.06, 0.09, 0.05, 0.06, 0.09, 0.04, 0.09, 0.09, 0.09 ]  !CLM 4.0 Table 8.1

prm%f_N = [ 0.83, 0.66, 0.72, 0.71, 0.64, 0.78, 0.70, 0.61, 0.64 ] !CLM 4.0 Table 8.1

! Turnover

prm%leafturnover  = [ 2.0,1.0,3.0,1.0,1.0,3.0,1.0,1.0,1.0 ]  !changed values to be consistent with lit estimates. JM 02.09.2010
                                                            !see Ben Smith's compendium of LPJ-GUESS parameters
                                                            
prm%leaflongevity = [ 2.0,0.7,3.0,1.0,0.7,3.0,0.5,1.0,1.0 ]  !Leafturnover and leaflongevity are essentially the same but for rain/summergreens.
                                                             !Note: I gave raingreens 0.7 instead of 0.5, JM 13.05.2011

prm%sapturnover = [ 20.0,20.0,20.0,20.0,20.0,20.0,20.0, 1.0, 1.0 ]

prm%frootQ10 = [ 1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.6,1.6 ]

prm%leaf_thres = [ 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50 ]

prm%phenramp = [ 1000.0,200.0,1000.0,1000.0, 200.0,1000.0, 200.0, 100.0, 100.0 ]

!prm%drydrop = [ 0.00,0.35,0.00,0.00,0.00,0.00,0.00,0.35,0.35 ]
prm%drydrop = [ 0.00, -1.5e5, 0.00, 0.00, 0.00, 0.00, 0.00, -1.5e5, -1.5e5 ]

! Fire 

prm%flam = [ 0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15 ]

prm%fireres = [ 0.12,0.50,0.12,0.50,0.12,0.12,0.12,1.00,1.00 ]

prm%CH4ef = [ 6.8, 2.2, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 2.2 ]

prm%COef  = [ 103.0, 63.0, 106.0, 106.0, 106.0, 106.0, 106.0, 106.0, 63.0 ]

prm%VOCef = [ 8.1, 3.4, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 3.4 ]

prm%TPMef = [ 8.5, 8.5, 17.6, 17.6, 17.6, 17.6, 17.6, 17.6, 8.5 ]

prm%NOxef = [ 1.85, 2.35, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.35 ]

! Bioclimatic limits

prm%mintcm = [ 15.5,15.5,-2.0,3.0,-17.0,-32.5,-1000.0,-1000.0,15.5 ]

prm%maxtcm = [ 1000.0,1000.0,22.0,18.8,15.5,-2.0,-2.0,15.5,1000.0 ]

prm%gddmin = [ 0.0,0.0,900.0,1200.0,1200.0,600.0,350.0,0.0,0.0 ]

prm%maxtwm = [ 1000.0,1000.0,1000.0,1000.0,1000.0,23.0,23.0,1000.0,1000.0 ]

prm%gddbase = [ 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 2.0, 1.0, 5.0 ]

! plant optical properties

prm%lfang = [ 0.10, 0.01, 0.01, 0.10, 0.25, 0.01, 0.25, -0.30, -0.30 ]

prm%alpha_s1 = [ 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.36, 0.36 ]

prm%alpha_s2 = [ 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58 ]

prm%alpha_l1 = [ 0.10, 0.10, 0.07, 0.10, 0.10, 0.07, 0.08, 0.11, 0.11 ]

prm%alpha_l2 = [ 0.45, 0.45, 0.35, 0.45, 0.45, 0.35, 0.40, 0.58, 0.58 ]

prm%tau_l1 = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.07 ]

prm%tau_l2 = [ 0.25, 0.25, 0.10, 0.25, 0.25, 0.10, 0.17, 0.25, 0.25 ]

prm%tau_s1 = [ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.220, 0.220 ]

prm%tau_s2 = [ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380 ]

prm%fi1 = [ 0.4334, 0.493637, 0.493637, 0.4334, 0.321125, 0.493637, 0.321125, 0.6602, 0.6602 ] 

prm%fi2 = [ 0.1168164, 0.011160702, 0.011160702, 0.1168164, 0.31374675, 0.011160702, 0.31374675, -0.2809908, -0.2809908 ]

prm%mu_bar = [ 0.980886301184341, 0.997877291640547, 0.997877291640547, 0.980886301184341, 0.963766860492100, 0.997877291640547, 0.963766860492100, 1.07731439820931, 1.07731439820931 ]    

!TWO-LEAF
prm%m = [ 0.0015, 0.004, 0.00125, 0.0015, 0.004, 0.001, 0.004, 0.0, 0.0 ]

prm%sla_top = [ 0.012, 0.03, 0.01, 0.012, 0.03, 0.008, 0.03, 0.05, 0.05 ]

! Resistance
prm%Rzmom = [ 0.075, 0.055, 0.055, 0.075, 0.055, 0.055, 0.055, 0.120, 0.120 ]

! Roots

prm%abar = [ 3.87, 3.98, 2.43, 3.46, 3.46, 5.86, 5.86, 5.86, 2.84 ]  !10: 8.99

prm%Bbar = [ 4.9, 4.1, 4.4, 4.2, 4.2, 2.9, 2.9, 1.4, 1.4 ]  !10: 1.2

prm%littleb = [ 13.80, 12.31, 7.95, 10.91, 10.91, 13.73, 13.73, 7.67, 3.72 ]  !10: 10.40

end subroutine initpftpars

end module pftparametersmod
