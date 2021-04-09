subroutine calcpixels(bounds,srtx,srty,cntx,cnty,gridres,lboundlon,uboundlat)

!FLAG this subroutine is not presently used! Check carefully before using. JM 19.03.2011

use arveparams, only : dp
use netcdf_error

implicit none

real(dp), intent(in), dimension(4) :: bounds
real(dp), intent(in), dimension(2) :: gridres
real, intent(in) :: lboundlon
real, intent(in) :: uboundlat

integer, intent(out) :: srtx
integer, intent(out) :: srty
integer, intent(out) :: cntx
integer, intent(out) :: cnty

!local variables

!these values are hardwired at the moment, but shouldn't be
!real :: grid_resx  =   60!225   !grid resolution in minutes     !FLAG Not used.
!real :: grid_resy  =   60!150   !grid resolution in minutes     !FLAG Not used.

real, dimension(2) :: pixfact

real :: minlon
real :: maxlon
real :: minlat
real :: maxlat

!-------------------------
!calculate the number of pixels to be calculated

pixfact = 60. / gridres

srtx = 1 + nint(pixfact(1) * (bounds(1) - lboundlon))
srty = 1 + nint(pixfact(2) * (uboundlat - bounds(4)))

minlon = bounds(1)
if (bounds(2) == bounds(1)) then
  maxlon = bounds(1)  + (1./pixfact(1))
else
  maxlon = bounds(2)
end if

minlat = bounds(3)
if (bounds(4) == bounds(3)) then
  maxlat = bounds(3)  + (1./pixfact(2))
else
  maxlat = bounds(4)
end if

cntx = nint(pixfact(1) * (maxlon - minlon))
cnty = nint(pixfact(2) * (maxlat - minlat))

if (cntx <= 0 .or. cnty <= 0) then
  write(0,*)cntx,cnty
  stop 'invalid coordinates'
end if

end subroutine calcpixels
