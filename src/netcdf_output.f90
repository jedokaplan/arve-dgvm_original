subroutine netcdf_output(yr)

use typesizes
use netcdf
use netcdf_error

use iovariables
use statevars
use arveparams

implicit none

integer :: yr
integer, dimension(1) :: lyear
integer :: i,x,y
real, allocatable, dimension(:,:,:) :: ovar


!----
lyear = yr

!write the time variable

status = nf90_put_var(ncid,6,lyear,start=[lyear],count=[1])
if (status/=nf90_noerr) call handle_err(status)

!----
!write other variables

!----

allocate(ovar(cntx,cnty,12))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
    if (sv(i)%valid_cell) ovar(x,y,:) = ov(i)%mw1 !CHECK!
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,7,ovar,start=[1,1,1,lyear],count=[cntx,cnty,12,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

!FLAG mts is presently not assigned, dummy value. JM 02.03.2011
allocate(ovar(cntx,cnty,12))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = 0._dp!ov(i)%mts
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,8,ovar,start=[1,1,1,lyear],count=[cntx,cnty,12,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,npft))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = real(sv(i)%lai_max)
    i = i + 1
 end do
end do

status = nf90_put_var(ncid,9,ovar,start=[1,1,1,lyear],count=[cntx,cnty,npft,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,npft))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = ov(i)%leafondays
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,10,ovar,start=[1,1,1,lyear],count=[cntx,cnty,npft,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,npft))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = real(sv(i)%fpc_gmax)
        i = i + 1
  end do
end do

status = nf90_put_var(ncid,11,ovar,start=[1,1,1,lyear],count=[cntx,cnty,npft,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%grid_npp
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,12,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

!----

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%annual_hetresp
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,13,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

!----

allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = ov(i)%afire_frac
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,14,ovar,start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = real(sv(i)%grid_area)
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,15,ovar,start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----
allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = ov(i)%plant_carbon
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,16,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = real(sv(i)%cpool_fast)
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,17,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
        if (sv(i)%valid_cell) ovar(x,y,:) = real(sv(i)%cpool_slow)
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,18,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!----

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%acflux_fire
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,19,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,npft))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = real(sv(i)%nind)
        i = i + 1
  end do
end do

status = nf90_put_var(ncid,20,ovar,start=[1,1,1,lyear],count=[cntx,cnty,npft,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%grid_gpp
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,21,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,ncvar))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%grid_autoresp
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,22,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---
allocate(ovar(cntx,cnty,12))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%runoff
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,23,ovar,start=[1,1,1,lyear],count=[cntx,cnty,12,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,12))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%snow_frac
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,24,ovar,start=[1,1,1,lyear],count=[cntx,cnty,12,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)
!---
allocate(ovar(cntx,cnty,12))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%snow_depth
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,25,ovar,start=[1,1,1,lyear],count=[cntx,cnty,12,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%abs_rad_annual
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,26,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell)  ovar(x,y,:) = ov(i)%grndwat_depth
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,27,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---

allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell) then
         if  (ov(i)%max_ald /= 999.) ovar(x,y,:) = ov(i)%max_ald
       end if
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,28,ovar(:,:,1),start=[1,1,lyear],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)

!---
allocate(ovar(cntx,cnty,1))

ovar = fill_value

i = 1
do y = 1,cnty
  do x = 1,cntx
       if (sv(i)%valid_cell) then
         ovar(x,y,1) = 0.
         if (sv(i)%peat) ovar(x,y,1) = 1.
       end if
    i = i + 1
  end do
end do

status = nf90_put_var(ncid,29,ovar(:,:,1),start=[1,1,1],count=[cntx,cnty,1])
if (status/=nf90_noerr) call handle_err(status)

deallocate(ovar)
!--

end subroutine netcdf_output

!----------------

subroutine netcdf_close()

use typesizes
use netcdf
use netcdf_error

use iovariables
use statevars

implicit none

status = nf90_close(ncid)
if (status/=nf90_noerr) call handle_err(status)

end subroutine netcdf_close
