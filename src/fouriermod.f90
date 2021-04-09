module fouriermod

use arveparams, only : dp,pi

implicit none

public :: fitfourier

contains

!-------------------------------

subroutine fitfourier(xdata,xbar,c,t)

implicit none

!fouri2: This subroutine fits a single cos curve to data input in the vector xdata.
!It is a general Fourier method where the phase is converted to days assuming 12
!  monthly values input (for use with estmo)
!  J.S. Auburn     2 June 1988

!Adapted from Richardson's (in par/wgen); changed so that:
!  a) only one curve is fit, rather than triplicate
!  b) amplitude calculation: c=sqrt(a*a+b*b) in place of: c=a/cos(t)  (a little better)
!  c) phase: t=atan2(-b,a) in place of t=atan(-b/a) so phase is correct in all quadrants
!  d) npts is variable rather than fixed at 13 (period is fixed at npts)

real(dp), dimension(:), intent(in) :: xdata  !the time series of the input data (of arbitrary length)

real(dp), intent(out) :: xbar  !the mean of the data series
real(dp), intent(out) :: c     !the amplitude of the curve-fit to the time series
real(dp), intent(out) :: t     !the phase of the curve-fit to the time series, in days

!local variables

integer :: npts
integer :: i

real(dp) :: suma
real(dp) :: sumb
real(dp) :: a
real(dp) :: b

real, dimension(size(xdata)) :: k

!--------

npts = size(xdata)

xbar = sum(xdata)/npts

k = [(real(i),i=1,npts)]

suma = sum((xdata - xbar) * cos(2._dp * pi * k / npts))
sumb = sum((xdata - xbar) * sin(2._dp * pi * k / npts))

a = suma * (2._dp / npts)
b = sumb * (2._dp / npts)

c = sqrt(a**2 + b**2)

!...Convert phase to days from radians,
!...and add 15 days since the first x point represents "mid-January"

t = atan2(-b,a)

t = t * (365._dp / (2._dp * pi)) + 15._dp
t = -t
t = mod(t,365._dp)

if (t < 0._dp) then
  t = 365._dp + t
end if

end subroutine fitfourier

!-------------------------------

end module fouriermod
