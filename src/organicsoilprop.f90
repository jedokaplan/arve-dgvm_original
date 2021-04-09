module organicsoilprop

use arveparams, only : dp

implicit none

!organic soil parameters
public :: fKsato              ! function to calculate saturated conductivity
public :: fTsato              ! function to calc. vol. water content at saturation
public :: fBexpo              ! function to calc Brooks and Corey B exponent
public :: fPsato              ! function to calc saturated soil matrix potential

!depth boundaries separating fibric, hemic and sapric layers of organic soil (NB these values may be adjusted)
real(dp), parameter, dimension(2) :: lb = [ 0.3_dp, 1.0_dp ]

contains

!Note: the functions below are extracted from Figure 2 in Letts et al., 2000
!The peat type boundaries are arbitrary, and may be adjusted

!------

real(dp) function fKsato(dz)  !saturated hydraulic conductivity (mm s-1)

real(dp), intent(in)  :: dz   !depth to soil layer midpoint (m)

if (dz < lb(1)) then
  !fibric peat
  fKsato = 0.28_dp
else if (dz < lb(2)) then
  !hemic peat
  fKsato = 0.002_dp
else
  !sapric peat
  fKsato = 0.0001_dp
end if

end function fKsato

!------

real(dp) function fTsato(dz)  !total porosity (fraction)

real(dp), intent(in)  :: dz   !depth to soil layer midpoint (m)

if (dz < lb(1)) then
  !fibric peat
  fTsato = 0.93_dp
else if (dz < lb(2)) then
  !hemic peat
  fTsato = 0.88_dp
else
  !sapric peat
  fTsato = 0.83_dp
end if

end function fTsato

!------

real(dp) function fBexpo(dz)  !B exponent (unitless)

real(dp), intent(in)  :: dz   !depth to soil layer midpoint (m)

if (dz < lb(1)) then
  !fibric peat
  fBexpo = 2.7_dp
else if (dz < lb(2)) then
  !hemic peat
  fBexpo = 6.1_dp
else
  !sapric peat
  fBexpo = 12.0_dp
end if

end function fBexpo

!------

real(dp) function fPsato(dz)  !air entry matric potential (mm)

real(dp), intent(in)  :: dz   !depth to soil layer midpoint (m)

if (dz < lb(1)) then
  !fibric peat
  fPsato = -10.3_dp
else if (dz < lb(2)) then
  !hemic peat
  fPsato = -10.2_dp
else
  !sapric peat
  fPsato = -10.1_dp
end if

end function fPsato

!------

end module organicsoilprop
