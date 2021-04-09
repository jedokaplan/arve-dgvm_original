module soilmoistureforcing

use arveparams, only : sp,dp

implicit none

!this module contains subroutines to force ARVE with externally sourced soil moisture

public :: initsoilmfile
public :: getsoilanom
public :: nudgesoil

integer :: smfid
integer :: anomid

real(sp), dimension(365) :: dsoilanom

real(dp), allocatable, dimension(:,:,:) :: anom_in

contains

!----------------------------------------------------------

subroutine initsoilmfile()

use netcdf
use typesizes
use netcdf_error
use iovariables, only : cntx,cnty,dosoilm_assim

implicit none

!integer :: ncells

character(100) :: soilmfile

!-----------------
!specify external soil moisture forcing file here

soilmfile = '/Volumes/Rossby/soilm_globetest.nc'

status = nf90_open(soilmfile,nf90_nowrite,smfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(smfid,'soilm',anomid)
if (status /= nf90_noerr) call handle_err(status)

!Allocate monthly soil moisture arrays
allocate(anom_in(cntx,cnty,12))

write(0,*)'***WARNING*** Using soil moisture external forcing for transient run'
write(0,*)trim(soilmfile)

end subroutine initsoilmfile

!----------------------------------------------------------

subroutine getsoilanom(year)

use netcdf
use typesizes
use netcdf_error
use iovariables,  only : cntx,cnty,srtx,srty

implicit none

!arguments
integer, intent(in) :: year  !needs as input the year of the transient run

!local variables
integer :: tpos

!-----

tpos = 1 + 12 * (year - 1) !some function of the year of the run

status = nf90_get_var(smfid,anomid,anom_in,start=[srtx,srty,tpos],count=[cntx,cnty,12])
if (status /= nf90_noerr) call handle_err(status)

anom_in = anom_in + 0.5

end subroutine getsoilanom

!----------------------------------------------------------

subroutine nudgesoil(j,msoil,anom)

use statevars, only : sv
use soilstate_vars, only : surf

implicit none

integer,  intent(in) :: j       !grid index
real(dp), intent(in) :: msoil   !long term mean soil moisture fraction of total porosity
real(dp), intent(in) :: anom    !anomaly forcing for soil moisture (fraction)

real(dp), pointer, dimension(:) :: Wliq    !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice    !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tsat    !soil porosity
real(dp), pointer, dimension(:) :: dzmm    !soil layer thickness

real(dp), dimension(1:sv(j)%gnl) :: icef

integer :: gnl

integer :: l

!---------------
!update soil moisture state variables

Tsat => sv(j)%Tsat
Wice => sv(j)%Wice
Wliq => sv(j)%Wliq
dzmm => sv(j)%dzmm

gnl = sv(j)%gnl

!calculate the ice/liquid fraction in each layer

do l = 1,gnl
  if (Wice(l) > 0.) then
    icef(l) = Wice(l) / (Wice(l) + Wliq(l))
  else
    icef(l) = 0.
  end if
end do

!calculate soil water content as the climatological mean * the anomaly fraction; limit the total to 1 (fully saturated soil).

!write(0,'(a,f8.3)')'anomaly: ',anom
!do l = 1,gnl
!  write(0,'(i5,4f8.3)')l,Tsat(l)*dzmm(l),msoil*Tsat(l)*dzmm(l),Wice(l),Wliq(l)
!end do

Wice(1:gnl) = min(msoil * anom, 1._dp) * Tsat(1:gnl) * dzmm(1:gnl) * icef
Wliq(1:gnl) = min(msoil * anom, 1._dp) * Tsat(1:gnl) * dzmm(1:gnl) * (1._dp - icef)


!write(0,*)'after calc'
!do l = 1,gnl
!  write(0,'(i5,4f8.3)')l,Tsat(l)*dzmm(l),msoil*Tsat(l)*dzmm(l),Wice(l),Wliq(l)
!end do
!read(*,*)

end subroutine nudgesoil

!----------------------------------------------------------

end module soilmoistureforcing
