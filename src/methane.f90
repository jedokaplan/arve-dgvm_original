subroutine methane(j)

use arveparams, only : dp
use iovariables, only : slope
use statevars, only : sv
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j !index value of the current grid cell

!local variables
integer :: i

!parameters
real(dp), parameter :: k_soil_fast10 =  0.03
real(dp), parameter :: anaerobic_f   =  0.03
real(dp), parameter :: ch4_13cfrac   = -40.

real(dp), dimension(3), parameter :: wthreshold = [0.65, 0.78, 0.92]

!pointers
real(dp), pointer, dimension(:) :: Tliq    !timestep soil liquid water content (fraction)
real(dp), pointer, dimension(:) :: drh     !daily heterotrophic respiration (gC/m2)
real(dp), pointer :: gra_dp75                   !slope of the gridcell below the threshold (fractional area) (%)
real(dp), pointer :: gra_dp50
real(dp), pointer :: grad200
wetl_area

!point pointers
Tliq    => sv(j)%Tliq
drh     => veg%drh
gra_dp75 => slope(j)%sf(1)
gra_dp50 => slope(j)%sf(2)
grad200 => slope(j)%sf(3)

write(*,*)slope(j)%elevation,slope(j)%sf2,slope(j)%sf075,slope(j)%sf05,sv(j)%grid_area
write(*,*)sv(j)%Patm,Pstd

!-------
!begin calculations
exit
!determine if gridcell is wet enough and has some heterotrophic respiration occuring
  if (Tliq(1) > 0.65 .and. drh /= 0._dp) then

  !The soil hydrology is the same for the entire gridcell (i.e. there is no spotty-ness to the soil wetness)
  !yet the slope is a fractional amount, so different parts of the gridcell will have different slope. To account
  !for this difference, we will calculate the methane and wetland area for the whole gridcell and then scale it to the
  !areas that really are wetland (i.e. are flat enough).

  !The slope fractions are cumulative, so the highest gradient (most steep) will also contain the flatter sections in
  !its value. The next gradient will not include the highest gradient but will contain any points flatter than itself.

  !determine the wetland area for this timestep.
  if (Tliq(1) >= 0.65 .and. Tliq(1) < 0.78)    slopefr = gra_dp50
  if (Tliq(1) >= 0.78 .and. Tliq(1) < 0.92)    slopefr = gra_dp75
  if (Tliq(1) >= 0.92)    slopefr = grad200

  !determine the amount of methane produced in the wetland area
  CH4 =

    wetl_area = drh *
   !   do k = 1,7
  !        if (soilwatr(i,j,m) > wthreshold(k)) then
  !          ch4flux(i,j,m) = ch4flux(i,j,m) + (anaerobic_f * hetresp(i,j,m) * slope(i,j,k))
  !        end if

  end if



end subroutine methane








!++++++++++++++++++
!this is what Jed had used for the Quest meeting. Must change siginificantly

use typeSizes
use netcdf

implicit none

integer :: status
integer :: ncid
integer :: dimid
integer :: varid
integer :: ofid

integer :: xlen
integer :: ylen
integer :: mlen
integer :: tlen

integer :: i,j,k,m,t

real(dp), allocatable, dimension(:,:,:) :: slope

real(dp), allocatable, dimension(:) :: lat
real(dp), allocatable, dimension(:,:) :: areamask
real(dp), allocatable, dimension(:,:) :: soilc
real(dp), allocatable, dimension(:,:) :: soili

real(dp), allocatable, dimension(:,:,:) :: ch4iso
real(dp), allocatable, dimension(:,:,:) :: soiltemp
real(dp), allocatable, dimension(:,:,:) :: soilwatr
real(dp), allocatable, dimension(:,:,:) :: tempresp
real(dp), allocatable, dimension(:,:,:) :: hetresp
real(dp), allocatable, dimension(:,:,:) :: ch4flux

character(20), dimension(7), parameter :: slname = &
                               [ 'gra_dp005 ','gra_dp0075','gra_dp01  ','gra_dp015 ','gra_dp02  ','gra_dp03  ','gra_dp1   ' ]

real(dp), dimension(7), parameter :: wthreshold = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]

character(30) :: fname
character(80) :: lpjfile

real(dp), parameter :: k_soil_fast10 =  0.03
real(dp), parameter :: anaerobic_f   =  0.03
real(dp), parameter :: ch4_13cfrac   =-40.

real(dp) :: arhsum
real(dp) :: ch4sum
real(dp) :: tempmn
real(dp) :: isosum

!---------------

xlen = 96
ylen = 72
mlen = 12

!---------------
!read slopes

allocate(slope(xlen,ylen,7))

do i=1,7
  fname = trim(slname(i))//'.grd'
  status = nf90_open(fname,nf90_nowrite,ncid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid,3,slope(:,:,i))
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_close(ncid)
  if (status /= nf90_noerr) call handle_err(status)
end do

!---------------
!read data for one time step and calc methane

call getarg(1,lpjfile)

status = nf90_open(lpjfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'time',dimid)
status = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

allocate(lat(ylen))

status = nf90_get_var(ncid,2,lat)
if (status /= nf90_noerr) call handle_err(status)

allocate(areamask(xlen,ylen))

do j = 1,ylen
  areamask(:,j) = area(lat(j),3.75,2.5)
end do

allocate(soilc(xlen,ylen))
allocate(soili(xlen,ylen))

allocate(soiltemp(xlen,ylen,mlen))
allocate(soilwatr(xlen,ylen,mlen))

allocate(tempresp(xlen,ylen,mlen))
allocate(hetresp(xlen,ylen,mlen))
allocate(ch4flux(xlen,ylen,mlen))
allocate(ch4iso(xlen,ylen,mlen))

do t = 1,tlen

  write(0,*)'working on: ',t

  status = nf90_get_var(ncid,7,soilwatr,start=[1,1,1,t],count=[xlen,ylen,mlen,1]) !mwatr
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(ncid,8,soiltemp,start=[1,1,1,t],count=[xlen,ylen,mlen,1]) !mtemp
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(ncid,9, soilc,start=[1,1,t],count=[xlen,ylen,1]) !mtemp
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(ncid,10,soili,start=[1,1,t],count=[xlen,ylen,1]) !mtemp
  if (status /= nf90_noerr) call handle_err(status)


  !where (soiltemp < -40.)
  !  tempresp = 0.0
  !elsewhere
  !  tempresp = exp(308.56 * ((1.0 / 56.02) - (1.0 / (soiltemp + 273.0 - 227.13)))) !Lloyd & Taylor 1994
  !end where

  !where (soilc < 0.) soilc = 0.

  !do m = 1,12
  !  hetresp(:,:,m) = soilc * (1.0 - exp(-k_soil_fast10 * tempresp(:,:,m) / 12.))
  !end do

  ch4flux = 0.

  !do j = 1,ylen
  !  do i = 1,xlen
  !    do m = 1,12
  !      do k = 1,7
  !        if (soilwatr(i,j,m) > wthreshold(k)) then
  !          ch4flux(i,j,m) = ch4flux(i,j,m) + (anaerobic_f * hetresp(i,j,m) * slope(i,j,k))
  !        end if
  !      end do
  !    end do
  !  end do
  !end do

  !do m = 1,12
  !  ch4iso(:,:,m) = ch4_13cfrac - soili
  !end do

  isosum = sum(ch4iso * ch4flux) / sum(ch4flux)

  arhsum = sum(sum(hetresp,dim=3) * areamask) * 1e-15  !Pg C
  ch4sum = sum(sum(ch4flux,dim=3) * areamask) * 1e-12  !Tg C
  tempmn = sum(soiltemp,mask=soiltemp/=-9999)/count(soiltemp/=-9999)

  write(*,*)t,arhsum,ch4sum,tempmn,isosum

end do

!-------------------------------------
contains

  subroutine handle_err(status)
  !Internal subroutine - checks error status after each netcdf call,
  !prints out text message each time an error code is returned.

  integer, intent(in) :: status

  if(status /= nf90_noerr) then
    write(0,*)trim(nf90_strerror(status))
    stop
  end if

  end subroutine handle_err


!---------

end subroutine methane
