module soilparameters

use arveparams, only : dp,OMorgC,peatlim,clayd,soilbulk

implicit none

public :: fbulk              ! function to calculate the bulk density for organic soils and mineral soils
public :: fRock

!mostly mineral soil parameters
public :: fKsat              ! function to calculate saturated conductivity
public :: fTsat              ! function to calc. vol. water content at saturation
public :: fBexp              ! function to calc Brooks and Corey B exponent
public :: fPsat              ! function to calc saturated soil matrix potential
!peat soil parameters
public :: fKsato              ! function to calculate saturated conductivity
public :: fTsato              ! function to calc. vol. water content at saturation
public :: fBexpo              ! function to calc Brooks and Corey B exponent
public :: fPsato              ! function to calc saturated soil matrix potential

contains

!===========================++++++++++++++++++++++++

function fbulk(depth,sorg,Vcf,clay,silt)

!function to calculate the soil bulk density for org as a function of depth (kg m-3)
!Values were manipulated from Minkkinen and Laine 1998 Can. J. For. Res for peat
!For mineral soil bulk density is adapted from Heuscher et al. SSSAJ 2005

implicit none

real(dp) :: fbulk
real(dp) :: depth
real(dp) :: sorg
real(dp) :: Vcf
real(dp) :: clay
real(dp) :: silt

real(dp), parameter :: a = 250.!-0.06932325934_dp
real(dp), parameter :: b = 25.!0.600509_dp
real(dp), parameter :: c = 75.!10.0_dp
real(dp), parameter :: d = 100.!5.0_dp
real(dp), parameter :: e = 1.685_dp
real(dp), parameter :: f = 0.198_dp
real(dp), parameter :: g = 0.0079_dp
real(dp), parameter :: h = 0.00014_dp
real(dp), parameter :: i = 0.0007_dp
real(dp), parameter :: bulkfrag = 1750.   !Bulk density of gravel
!---

  if (sorg < peatlim) then

          !mineral
      if (depth < 2._dp) then

      fbulk = (e - f * (SQRT(sorg * OMorgC)) + g * clay + h * depth * 100._dp - i * silt) * 1000._dp
      else
      fbulk = (e - f * (SQRT(sorg * OMorgC)) + g * clay + h * 2._dp * 100._dp - i * silt) * 1000._dp
      end if

      if (Vcf > 0._dp) then
        !from Mehuys et al. 1975 and Berger, 1976
        fbulk = fbulk + Vcf * (bulkfrag - fbulk)
      end if

  else

          !organic
      if (depth < 0.30) then ! this is for peat soils 20kg m-3 at surface 
                             ! to 100 kgm-3 at 30 cm from Minkkinen graph inspect.
                             ! Peat bulk density from Rydin and Jeglum 2006 The biology
                             ! of Peatlands, is 0.02 g/cm3 at surface
                             ! and depth from 0.1 - 0.26 g/cm3
        fbulk = a * depth + b
        else
        fbulk = c * (depth - 0.3) + d
      end if


      if (fbulk > 260.0) then  !Value for deeper peat from Rydin and Jeglum book
          fbulk = 260.0
      end if

  end if

end function fbulk

!===========================++++++++++++++++++++++++

function fRock(Vcf,bulk)

!correction factor for rock fragments (from Skirvin correction to WEPP model)
!this function takes the volume fraction and converts it to a mass fraction for
!coarse fragments.

implicit none

!real(dp) :: fVcf
!real(dp) :: rock
real(dp) :: Vcf
real(dp) :: frock
real(dp) :: bulk

real(dp), parameter :: a = 2685.
!---

  if (Vcf > 0._dp) then
    frock = a / (bulk / Vcf - bulk + a)
  else
    frock = 0._dp
  end if

  !depricated function., presently not required. JM 01.06.2010.
  !fVcf = (rock * 0.01 * bulk) / (rock * 0.01 * bulk + a * (1. - rock * 0.01))
  !Alternative method to estimate from Brakensiek & Rawls 1994 Catena
  !this method is not presently used.
  !frock = 2. * Vcf / (1. + Vcf)

end function fRock

!===========================++++++++++++++++++++++++
function fKsat(sand,sorg)

!function to calculate saturated conductivity (mm s-1)
!related to pore geometry and particle shape also packing and structure

implicit none

!arguments
real(dp), intent(in) :: sand   !soil sand concent (percent)
real(dp), intent(in) :: sorg   !soil sorg concent (percent)

!local variables
real(dp) :: fKsat               !saturated conductivity (mm s-1)
real(dp) :: Ksat_uncon          !Ksat for the unconnected pores
real(dp) :: Ksat_min            !Ksat for the mineral portion of the soil
real(dp) :: fperc               !fraction of the grid cell that is connected through the organic matter channels
real(dp) :: Nperc               !used in fperc calculation
real(dp) :: f_uncon             !fraction of grid cell that is not connected
real(dp) :: f_om                !organic matter weight fraction

!parameters
real(dp), parameter :: fthreshold = 0.5_dp       !threshold for allowing organic 'connections' (see CLM 4.0 p. 142)
real(dp), parameter :: Beta_perc = 0.139        !
real(dp), parameter :: Ksat_om = 0.1_dp          !mm s-1 from Lawrence and Slater 2008 Climate Dyn. 30:145-160
!---

!Adopted from CLM 4.0, see pg. 142 of tech note. JM & MP 07.05.2010
  !Vereecken et al. 1990 and CLM 3.0 were also tried but gave unrealistic values.

  f_om = sorg * 0.01_dp

  Nperc = (1._dp - fthreshold)**(-Beta_perc)

  if (f_om >= fthreshold) then
   fperc = Nperc * f_om * (f_om - fthreshold)**Beta_perc
  else
    fperc = 0._dp
  end if

  f_uncon = (1._dp - fperc)

  Ksat_min = 0.0070556 * 10**(-0.884 + 0.0153 * sand)

  Ksat_uncon = f_uncon / ((1._dp - f_om) / Ksat_min + (f_om - fperc) / Ksat_om)

  fKsat = f_uncon * Ksat_uncon + (1._dp - f_uncon) * Ksat_om


end function fKsat

!===========================++++++++++++++++++++++++

function fTsat(sand,clay,sorg,Vcf,bulk,mdepth)

!function to calculate volumetric water content at saturation (mm3 mm-3) takes into
!account the clay, org, sand and rock fractions. This is from the WEPP model description
!Also for calculating CEC comes from Horn et al. J Plant Nutr Soil Sci 2005.

implicit none

real(dp), intent(in) :: sand   !soil sand concent (percent)
real(dp), intent(in) :: clay   !soil clay concent (percent)
real(dp), intent(in) :: sorg   !soil organic concent (percent)
real(dp), intent(in) :: Vcf    !fraction of coarse fragments (rocks)
real(dp), intent(in) :: bulk   !soil bulk density
real(dp), intent(in) :: mdepth  !zpos, depth of the midpoint of the soil horizon

real(dp) :: CEC                          !cation exchange capacities
real(dp) :: CECc
real(dp) :: CECr
real(dp) :: Fair
real(dp) :: fTsat
real(dp) :: f_clay,f_sand,f_sorg        !fractional weight percents

f_clay = clay * 0.01_dp
f_sand = sand * 0.01_dp
f_sorg = sorg * 0.01_dp

CEC = 0.95_dp + (0.857936_dp * sorg * 10._dp * OMorgC) + (0.053_dp * clay * 10.) ! from Horn et al. J Plant Nutr Soil Sci 2005.(in cmol/kg = meq/100g)

CECc = CEC - f_sorg * (142 + 170 * mdepth)

  if (clay > 0) then !ensures no divide by 0's also set CECr to 0 (which is not important then due to no clay)
   CECr = CECc / clay
  else
   CECr = 0
  end if

Fair = 1. - ((3.8_dp + 1.9_dp * f_clay * f_clay - 3.365_dp * f_sand + 12.6_dp * CECr * f_clay + 100. * f_sorg) * (f_sand * 0.5) * (f_sand * 0.5)) * 0.01

fTsat = (1._dp - Vcf) * (1._dp - bulk / soilbulk) * Fair

!fTsat = 0.489 - 0.00126 * (sand)   ! old CLM 3.0 function for reference. (no organic or rock)

end function fTsat

!===========================++++++++++++++++++++++++

function fBexp(clay,silt,sand,rock,sorg)

  !Based upon Bloemen, 1980 Z. Pflanzenernaehr. Bodenkd. 143, 581-605
  !Since Bloemen assumes a lower limit of 1.4 that value is used to
  !convert the 'n' value to a Bexp value suitable for use in the rest of the model
  !CLM uses a lower limit of '2'.

implicit none

real(dp) :: fBexp
real(dp), intent(in) :: clay   !soil clay concent (percent)
real(dp), intent(in) :: silt   !soil clay concent (percent)
real(dp), intent(in) :: sand   !soil sand concent (percent)
real(dp), intent(in) :: rock   !soil rock concent (mass fraction)
real(dp), intent(in) :: sorg   !soil organic concent (percent)

real(dp), dimension(4), parameter :: sclass = [ 1., 26., 1125., 126000. ]  !size classes in um 0-2,2-50,50-2000,2000-250000 (~6000) so took means
real(dp), parameter :: b = 1.4_dp
real(dp), parameter :: c = 4.536_dp
real(dp), parameter :: d = 0.75_dp
real(dp), parameter :: e = 1.6_dp

real(dp), dimension(5) :: wps
real(dp), dimension(4) :: wps_cum
real(dp), dimension(3) :: tg
real(dp), dimension(4) :: f
real(dp) :: fsum
real(dp) :: n
integer :: i
integer :: base

!---

wps(1) = clay
wps(2) = silt
wps(3) = sand
wps(4) = rock * 100.
wps(5) = sorg

!NOTE: if there is no clay in the soil, then you can get a circumstance of divide by zero in equation 8
! to prevent that, tell it to simply skip the clay (moves the start of the loops to the second consitituent, which is silt.
! JM 23.04.2010

 if (wps(1) > 0._dp) then
   base = 1
 else
   base = 2
 end if


 do i = base,4 !make the Pi in the eqn 8

  wps_cum(i) = sum(wps(1:i))

 end do

 do i = base,3 !Makes the tgi (eqn 8) for the eqn 9

  tg(i) = log10(wps_cum(i+1) / wps_cum(i)) / log10(sclass(i+1) / sclass(i))

 end do

 do i = base,3

  f(i) = wps(i+1) * tg(i)

 end do

fsum = sum(f(1:3)) / sum(wps(2:4))

if (sorg > 0._dp) then  ! equation 12

  n = (b + c * (exp(fsum) - 1._dp))- d * (fsum**e) * log10(wps(5))
  fBexp = 3._dp/ (n - 1.4_dp)

else

  n = b + c * (exp(0.3_dp * fsum) - 1._dp)  !equation 11
  fBexp = 3._dp/ (n - 1.4_dp)

end if

end function fBexp

!===========================++++++++++++++++++++++++

function fPsat(sand,rock)

!function to calculate saturated soil matric potential
!closely and directly related to largest pore size

implicit none

real(dp) :: fPsat
real(dp), intent(in) :: sand   !soil sand concent (percent)
real(dp), intent(in) :: rock   !soil sand concent (fraction)

real(dp), parameter :: a = 1.8800_dp
real(dp), parameter :: b = 0.0131_dp

!since in ARVE-DGVM we include coarse fragments, I have added in rock to this
!gravel will impact upon the largest pore size. -JM Apr 17 08

fPsat = -10._dp * 10._dp**(a - b * (sand + rock * 100._dp))

end function fPsat

!===========================++++++++++++++++++++++++

!Peat soil parameters (adapted from Bloemen 1983 Z. Pflanzenernaehr. Bodenk.)

function fKsato(bulk)

!function to calculate saturated conductivity (mm s-1) organic soil
!Eqn 4 in Bloemen paper. Expects bulk density as g cm-3 so multiple the incoming
!bulk value (kg m-3) by 1e-3

implicit none

real(dp) :: fKsato
real(dp), intent(in) :: bulk   !soil bulk density

real(dp), parameter :: a = 3.07870370e-7  !0.00266*10/86400 convert units to mm sec-1
real(dp), parameter :: b = -3.625
!--

fKsato = (a * (bulk * 0.001) ** b)    !equation 4

end function fKsato

!-----

function fTsato(sorg,bulk,Vcf,sand,clay,mdepth)

!function to calculate volumetric water content at saturation (mm3 mm-3)
!original function calculates the volume taken up by peat, thus the Tsat is the total minus
!that taken up by the peat and coarse frags.

implicit none

real(dp) :: fTsato
real(dp), intent(in) :: sorg   !soil organic content (percent)
real(dp), intent(in) :: bulk   !soil bulk density
real(dp), intent(in) :: Vcf    !fraction of coarse fragments (rocks)
real(dp), intent(in) :: sand   !soil sand concent (percent)
real(dp), intent(in) :: clay   !soil clay concent (percent)
real(dp), intent(in) :: mdepth  !zpos, depth of the midpoint of the soil horizon

real(dp), parameter :: b = 0.23584906_dp ! = 1/(clayd * 1.60)

real(dp) :: CEC                                !cation exchange capacities
real(dp) :: CECc
real(dp) :: CECr
real(dp) :: Fair
real(dp) :: f_clay,f_sand,f_sorg        !fractional weight percents
real(dp) :: minTsat                        !Tsat for the mineral components of the soil
real(dp) :: orgTsat                          !Tsat for the peat (organic) part of the soil
real(dp) :: M                !mineral fraction of the dry soil
!---

!Peat part
M = (100._dp - sorg) * 0.01

orgTsat = 1._dp - ((clayd - M) * (bulk * 0.001) * b)  !adapted equation 6

!Mineral part (as in fTsat)
!function to calculate volumetric water content at saturation (mm3 mm-3) takes into
!account the clay, org, sand and rock (Vcf) fractions. This is from the WEPP model description
!Also for calculating CEC comes from Horn et al. J Plant Nutr Soil Sci 2005.

f_clay = clay * 0.01_dp
f_sand = sand * 0.01_dp
f_sorg = sorg * 0.01_dp

CEC = 0.95_dp + (0.857936_dp * sorg * 10._dp * OMorgC) + (0.053_dp * clay * 10.) ! from Horn et al. J Plant Nutr Soil Sci 2005.(in cmol/kg = meq/100g)

CECc = CEC - f_sorg * (142 + 170 * mdepth)

  if (clay > 0) then !ensures no divide by 0's also set CECr to 0 (which is not important then due to no clay)
   CECr = CECc / clay
  else
   CECr = 0
  end if

!removed the sorg to avoid a double count.
Fair = 1. - ((3.8_dp + 1.9_dp * f_clay * f_clay - 3.365_dp * f_sand + 12.6_dp * CECr * f_clay + 100. * f_sorg) * (f_sand * 0.5)*(f_sand * 0.5)) * 0.01

minTsat = (1._dp - Vcf) * (1._dp - bulk / soilbulk) * Fair

!combine them for a weighted value
fTsato = (1. - f_sorg) * minTsat + f_sorg * orgTsat

end function fTsato

!===========================++++++++++++++++++++++++

function fBexpo(bulk)

!function to calculate Brooks & Corey B exponent
!Since Bloemen assumes a lower limit of 1.4 that value is used to
!convert the 'n' value to a Bexp value suitable for use in the rest of the model
!CLM uses a lower limit of '2'.

implicit none

real(dp) :: fBexpo
real(dp), intent(in) :: bulk   !soil bulk density

real(dp), parameter :: a = 2.54_dp
real(dp), parameter :: b = 2.42_dp
!---

fBexpo = 3._dp / ((a - b * bulk * 0.001) - 1.4_dp)   !equation 14

end function fBexpo

!===========================++++++++++++++++++++++++

function fPsato(bulk,sorg)
!function fPsato(bulk,sand,rock,sorg,Ksat)

!function to calculate saturated soil matric potential (mm)
!original calc gives cm so convert to mm by *10

implicit none

real(dp) :: fPsato
real(dp), intent(in) :: bulk   !soil bulk density
!real(dp), intent(in) :: sand   !soil sand concent (percent)
!real(dp), intent(in) :: rock   !soil sand concent (fraction)
real(dp), intent(in) :: sorg   !soil organic content (percent)
!real(dp), intent(in) :: Ksat
!real(dp) :: Psat_min
real(dp) :: fsorg

real(dp), parameter :: a = -416._dp
real(dp), parameter :: b = 1.12_dp
!real(dp), parameter :: c = 1.8800_dp
!real(dp), parameter :: d = 0.0131_dp
!real(dp), parameter :: Psat_peat = -10.3
!---

fsorg = sorg * 0.01

fPsato = (a * (0.001 * bulk) ** b) * 10._dp   !equation 9

!fPsato = -(exp((log((Ksat / 10. * 86400.) / 64565.) / -2.565)))

!since in ARVE-DGVM we include coarse fragments, I have added in rock to this
!gravel will impact upon the largest pore size. -JM Apr 17 08

!Psat_min = -10._dp * 10._dp**(c - d * (sand + rock * 100._dp))

!fPsato = (1._dp - fsorg) * Psat_min + fsorg * Psat_peat
!write(*,*)Psat_min,Psat_peat,fsorg,sand,rock

end function fPsato

end module soilparameters

