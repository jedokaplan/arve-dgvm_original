module esatdesdT

implicit none

public :: esat
public :: desdT

contains

!------------------------------------------------------------------------
function esat(temp)

  !Function to calculate saturation vapor pressure in water and ice
  !From CLM formulation, table 5.2, after Flatau et al. 1992

  use arveparams, only : dp,tfreeze

  implicit none

  real(dp) :: esat  !saturation vapor pressure (Pa)
  real(dp), intent(in) :: temp !temperature in K

  real(dp), dimension(9) :: al !coefficients for liquid water
  real(dp), dimension(9) :: ai !coefficients for ice

  real(dp), dimension(0:8) :: a !coefficients

  real(dp) :: T

  integer :: i

  al(1) = 6.11213476
  al(2) = 4.44007856e-1
  al(3) = 1.43064234e-2
  al(4) = 2.64461437e-4
  al(5) = 3.05903558e-6
  al(6) = 1.96237241e-8
  al(7) = 8.92344772e-11
  al(8) =-3.73208410e-13
  al(9) = 2.09339997e-16

  ai(1) = 6.11123516
  ai(2) = 5.03109514e-1
  ai(3) = 1.88369801e-2
  ai(4) = 4.20547422e-4
  ai(5) = 6.14396778e-6
  ai(6) = 6.02780717e-8
  ai(7) = 3.87940929e-10
  ai(8) = 1.49436277e-12
  ai(9) = 2.62655803e-15

  if (temp <= tfreeze) then   !these coefficients are for temperature values in Celcius
    a(0:8) = ai
  else
    a(0:8) = al
  end if

  T = temp - tfreeze

  esat = a(0)

  do i = 1,8
    esat = esat + a(i) * T**i
  end do

  esat = 100._dp * esat

end function esat

!------------

function desdT(temp)

  !Function to calculate the first derivative of saturation vapor pressure in water and ice vs. temperature
  !From CLM formulation, table 5.3, after Flatau et al. 1992

  use arveparams, only : dp,tfreeze

  implicit none

  real(dp) :: desdT    !derivative of saturation vapor pressure
  real(dp), intent(in) :: temp !temperature in K

  real(dp), dimension(9) :: bl !coefficients for liquid water
  real(dp), dimension(9) :: bi !coefficients for ice

  real(dp), dimension(0:8) :: b !coefficients

  real(dp) :: T

  integer :: i

  bl(1) = 4.44017302e-1
  bl(2) = 2.86064092e-2
  bl(3) = 7.94683137e-4
  bl(4) = 1.21211669e-5
  bl(5) = 1.03354611e-7
  bl(6) = 4.04125005e-10
  bl(7) =-7.88037859e-13
  bl(8) =-1.14596802e-14
  bl(9) = 3.81294516e-17

  bi(1) = 5.03277922e-1
  bi(2) = 3.77289173e-2
  bi(3) = 1.26801703e-3
  bi(4) = 2.49468427e-5
  bi(5) = 3.13703411e-7
  bi(6) = 2.57180651e-9
  bi(7) = 1.32268878e-11
  bi(8) = 3.94116744e-14
  bi(9) = 4.98070196e-17

  if (temp <= tfreeze) then
    b(0:8) = bi
  else
    b(0:8) = bl
  end if

  T = temp - tfreeze  !these coefficients are for temperature values in Celcius

  desdT = b(0)

  do i = 1,8
    desdT = desdT + b(i) * T**i
  end do

  desdT = 100._dp * desdT

end function desdT

!-------
end module esatdesdT
