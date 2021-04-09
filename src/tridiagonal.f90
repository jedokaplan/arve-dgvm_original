module tridiagonal

use arveparams, only : dp

implicit none

public :: tridiag

contains

!-------------------------------

subroutine tridiag(a,b,c,r,u)

!subroutine to solve triadiagonal system of equations

!Solves for a vector u of size N the tridiagonal linear set using given by equation (2.4.1) using a
!a serial algorithm. Input vectors b (diagonal elements) and r (right-hand side) have size N,
!while a and c (off-diagonal elements) are not defined in the first and last elements, respectively.
!Based on Numerical Recipes in F77/F90

implicit none


real(dp), dimension(:), intent(in) :: a
real(dp), dimension(:), intent(in) :: b
real(dp), dimension(:), intent(in) :: c
real(dp), dimension(:), intent(in) :: r
real(dp), dimension(:), intent(out) :: u

integer :: n
integer :: k
real(dp) :: bet
real(dp), dimension(size(b)) :: gam

!----

n = size(b)

bet = b(1)

u(1) = r(1) / bet

!decomposition and forward substitution

do k = 2,n

  gam(k) = c(k-1) / bet

  bet = b(k) - a(k) * gam(k)

  u(k) = (r(k) - a(k) * u(k-1)) / bet

end do

!backsubstitution

do k = n-1,1,-1
  u(k) = u(k) - gam(k+1) * u(k+1)
end do

end subroutine tridiag

!-------------------------------

end module tridiagonal
