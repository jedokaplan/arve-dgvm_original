module input_dataset_drivers

! All the input dataset drivers are kept here

implicit none

public :: initclimate
public :: soildriver
public :: slopedriver
!public :: satfracdriver
public :: openCO2file
public :: readCO2file

contains

!---

subroutine initclimate(path)

use netcdf_error
use iovariables, only : ibuf,va,timebuflen,inputclimlen,climatemonths, &
                        climateyears,inputlatlen,inputlonlen,niv,vname,climfid, &
                        latvect,lonvect,cntx,cnty,srtx,srty,bounds!,indexmat

implicit none

character(180), intent(in) :: path

character(180) :: climatefile

integer :: i

integer :: xsize
integer :: ysize
integer :: timelen

integer :: dimid
integer :: varid

integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos

!-------------------------
!generate file names. Regardless of the path name the input files must always have these names.

climatefile = trim(path)

!open climate file

status = nf90_open(climatefile,nf90_nowrite,climfid)
if (status /= nf90_noerr) call handle_err(status)

!retrieve dimensions

status = nf90_inq_dimid(climfid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(climfid,dimid,len=xsize)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(climfid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(climfid,dimid,len=ysize)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(climfid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(climfid,dimid,len=climatemonths)
if (status /= nf90_noerr) call handle_err(status)

climateyears = climatemonths / 12

inputclimlen = xsize * ysize

!retrieve scale factors and offsets

do i = 1,niv

  status = nf90_inq_varid(climfid,vname(i),varid)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_att(climfid,varid,'scale_factor',va(i)%scale_factor)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_att(climfid,varid,'add_offset',va(i)%add_offset)
  if (status /= nf90_noerr) call handle_err(status)

end do

!-------------------------
!NEW STATEMENTS ADDED 29.10.2007 - SHOULD ULTIMATELY BE MOVED SOMEWHERE ELSE. JOK
!-------------------------
!calculate the number and indices of the pixels to be calculated

deallocate(lonvect)
deallocate(latvect)

allocate(lonvect(xsize))
allocate(latvect(ysize))

status = nf90_inq_varid(climfid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(climfid,varid,lonvect)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(climfid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(climfid,varid,latvect)
if (status /= nf90_noerr) call handle_err(status)

pos = minloc(abs(lonvect - bounds(1)))

!when the compiler options were -fast, the abs above was not being calculated properly.

xpos(1) = pos(1)

pos = minloc(abs(lonvect - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(latvect - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(latvect - bounds(4)))
ypos(2) = pos(1)


srtx = minval(xpos)
srty = minval(ypos)

if (lonvect(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
cntx = 1 + abs(maxval(xpos) - srtx)

if (latvect(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
cnty = 1 + abs(maxval(ypos) - srty)

!allocate the lookup matrix for the climate vector and retrieve it
!(should be the same for all variables and must match in size to the soils data files)

if (xsize /= inputlonlen .or. ysize /= inputlatlen) then
  write(0,*)'size of soils and climate grids do not match!'
  write(*,*)'xsize',xsize,'ysize',ysize,'inputlonlen',inputlonlen,'inputlatlen',inputlatlen
  stop
end if

!allocate the vector input buffer for climate data

timelen = 12               !this version of the driver works with only one year of climate data at a time.
timebuflen = timelen * 1. !no of yrs climate data to store between reads

allocate(ibuf(inputclimlen))

do i=1,inputclimlen
 !new climate variables
  allocate(ibuf(i)%cld(timelen))
  allocate(ibuf(i)%dtr(timelen))
  allocate(ibuf(i)%pre(timelen))
  allocate(ibuf(i)%tmp(timelen))
  allocate(ibuf(i)%wet(timelen))

end do

end subroutine initclimate

!--------------------------------------------------------------------------------------------
subroutine soildriver(soil_path)
!new soil driver to take in FAO soil data. Coded Feb 08 JM

use netcdf_error
use iovariables, only : soil,soilfid,latvect,lonvect,inputsoillen,inputlonlen,inputlatlen,bounds, &
                        cntx,cnty,srtx,srty
use netcdf
use typesizes
use arveparams, only : dp

implicit none

character(180), intent(in) :: soil_path
character(100) :: soilfile

integer :: i
integer :: dimid
integer :: varid
integer :: xsize
integer :: ysize
integer :: depthint
integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos

!-------------------------
!generate file names. Regardless of the path name the input files must always have these names.
soilfile  = trim(soil_path)

!open soil files
status = nf90_open(soilfile,nf90_nowrite,soilfid)
if (status /= nf90_noerr) call handle_err(status)

!find dimensions
status = nf90_inq_dimid(soilfid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(soilfid,dimid,len=xsize)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(soilfid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(soilfid,dimid,len=ysize)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(soilfid,'depth',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(soilfid,dimid,len=depthint)
if (status /= nf90_noerr) call handle_err(status)

inputsoillen = xsize * ysize
inputlonlen = xsize
inputlatlen = ysize

!calculate the number and indices of the pixels to be calculated
allocate(lonvect(xsize))
allocate(latvect(ysize))

status = nf90_inq_varid(soilfid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(soilfid,varid,lonvect)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(soilfid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(soilfid,varid,latvect)
if (status /= nf90_noerr) call handle_err(status)

pos = minloc(abs(lonvect - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lonvect - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(latvect - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(latvect - bounds(4)))
ypos(2) = pos(1)

srtx = minval(xpos)
srty = minval(ypos)

if (lonvect(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
cntx = 1 + abs(maxval(xpos) - srtx)

if (latvect(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
cnty = 1 + abs(maxval(ypos) - srty)

!allocate the vector input buffer for soil data

allocate(soil(inputsoillen))

do i=1,inputsoillen
  allocate(soil(i)%sdto(depthint))
  allocate(soil(i)%stpc(depthint))
  allocate(soil(i)%clpc(depthint))
  allocate(soil(i)%totc(depthint))
  allocate(soil(i)%bulk(depthint))
  allocate(soil(i)%cfrag(depthint))
  allocate(soil(i)%tawc(depthint))
end do

!create soil diagnostic file here if needed.

end subroutine soildriver

!----------------------------------------------------------

subroutine slopedriver(slope_path,fmax_path)

 !new slope driver to take in the slope data for wetland location. Coded Mar 08 JM
 ! also takes in the maximal fractional saturated area values (11.05.2011 JM)

use typesizes
use netcdf
use netcdf_error
use iovariables, only : latvect,lonvect,inputsoillen,inputlonlen,inputlatlen,bounds &
                        ,cntx,cnty,srtx,srty,slope,slopefid,fmaxfid
use arveparams, only : dp

implicit none

character(180) :: slopefile
character(180), intent(in) :: slope_path

character(180) :: fmaxfile
character(180), intent(in) :: fmax_path

integer :: dimid
integer :: varid
integer :: xsize
integer :: ysize
integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos

!-------------------------

!generate file names. Regardless of the path name the input files must always have these names.
slopefile  = trim(slope_path)

fmaxfile = trim(fmax_path)

!open slope and fmax files
status = nf90_open(slopefile,nf90_nowrite,slopefid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_open(fmaxfile,nf90_nowrite,fmaxfid)
if (status /= nf90_noerr) call handle_err(status)

!find dimensions
status = nf90_inq_dimid(slopefid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(slopefid,dimid,len=xsize)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(slopefid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(slopefid,dimid,len=ysize)
if (status /= nf90_noerr) call handle_err(status)

!status = nf90_inq_dimid(slopefid,'slope',dimid)
!if (status /= nf90_noerr) call handle_err(status)

!status = nf90_inquire_dimension(slopefid,dimid,len=slopeint)
!if (status /= nf90_noerr) call handle_err(status)

!use the same dimensions as the soil dataset.
inputsoillen = xsize * ysize
inputlonlen = xsize
inputlatlen = ysize

!calculate the number and indices of the pixels to be calculated

status = nf90_inq_varid(slopefid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(slopefid,varid,lonvect)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(slopefid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(slopefid,varid,latvect)
if (status /= nf90_noerr) call handle_err(status)

pos = minloc(abs(lonvect - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lonvect - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(latvect - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(latvect - bounds(4)))
ypos(2) = pos(1)

srtx = minval(xpos)
if (lonvect(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
cntx = 1 + abs(maxval(xpos) - srtx)

srty = minval(ypos)
if (latvect(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
cnty = 1 + abs(maxval(ypos) - srty)

allocate(slope(inputsoillen))

end subroutine slopedriver

!-------------------------------------------
subroutine openCO2file(infile)

implicit none

character(80), intent(in) :: infile
!---------

open(20,file=infile,status='old') ! statue old ensures that it uses an exisiting file

!read header row
read(20,*) !header

end subroutine openCO2file

!-------------------------------
subroutine readCO2file

use statevars, only : CO2

implicit none

integer :: year

!------

read(20,*,end=7)year,CO2(1),CO2(2),CO2(3)

return

!this part handles the end-of-file condition

7 continue

rewind(20)

read(20,*) !header

end subroutine readCO2file
!--------------------------------------------

end module input_dataset_drivers
