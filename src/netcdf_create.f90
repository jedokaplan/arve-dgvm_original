subroutine netcdf_create()

use typesizes
use netcdf
use netcdf_error

use iovariables
use statevars
use arveparams, only : dp,nl,npft

implicit none

integer :: lon,lat,layer,pft,month,time
integer :: varid
integer :: i

character(8)  :: today
character(10) :: now

real(dp), dimension(2) :: xrange
real(dp), dimension(2) :: yrange

integer, allocatable, dimension(:) :: pftnum

!----------

xrange(1) = minval(lonvect(srtx:endx))
xrange(2) = maxval(lonvect(srtx:endx)) + 0.5

yrange(1) = minval(latvect(srty:endy)) - 0.5
yrange(2) = maxval(latvect(srty:endy))

!----1

status = nf90_create(outputfile,cmode=nf90_clobber,ncid=ncid)
if (status/=nf90_noerr) call handle_err(status)

!write(*,*)'created output file ncid:',ncid

status = nf90_put_att(ncid,nf90_global,'title','ARVE-DGVM netCDF output file')
if (status/=nf90_noerr) call handle_err(status)

call date_and_time(today,now)

status = nf90_put_att(ncid,nf90_global,'timestamp',today//' '//now(1:4))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'Conventions','COARDS')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'node_offset',1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'_Format',"netCDF-4")
if (status/=nf90_noerr) call handle_err(status)

!----2

status = nf90_def_dim(ncid,'lon',cntx,lon)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'lon',nf90_float,lon,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','longitude')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degrees_east')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',xrange)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)


!----3

status = nf90_def_dim(ncid,'lat',cnty,lat)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'lat',nf90_float,lat,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','latitude')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degrees_north')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',yrange)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

!----4

status = nf90_def_dim(ncid,'layer',nl,layer)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'layer',nf90_short,layer,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','soil layer')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','layer')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----5

status = nf90_def_dim(ncid,'pft',npft,pft)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'pft',nf90_short,pft,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','Plant Functional Type')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','PFT')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----6

status = nf90_def_dim(ncid,'month',12,month)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'month',nf90_short,month,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','month')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','month')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----7

status = nf90_def_dim(ncid,'time',nf90_unlimited,time)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'time',nf90_int,time,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','time')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','years since 1900')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------
!7
status = nf90_def_var(ncid,'mw1',nf90_float,[lon,lat,month,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','mean total soil column moisture (liq+ice) fraction of total porosity')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,12,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)


!----
!8
status = nf90_def_var(ncid,'mts',nf90_float,[lon,lat,month,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','top layer soil temperature')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degC')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,12,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!----
!9
status = nf90_def_var(ncid,'lai_max',nf90_float,[lon,lat,pft,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name',',max leaf area index')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','m2 m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,npft,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!----
!10
status = nf90_def_var(ncid,'leafon_days',nf90_float,[lon,lat,pft,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','days with foliage')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','days')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,npft,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)


!----
!11
status = nf90_def_var(ncid,'cover',nf90_float,[lon,lat,pft,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','PFT cover fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,npft,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)
!----
!12
status = nf90_def_var(ncid,'grid_npp',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total gridcell npp')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!----
!13
status = nf90_def_var(ncid,'arh',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','heterotrophic respiration')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!----
!14
status = nf90_def_var(ncid,'afire_frac',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','fraction of gridcell burnt this year')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','unitless')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!----
!5
status = nf90_def_var(ncid,'grid_area',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','gridcell area')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','m2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!----
!16
status = nf90_def_var(ncid,'plant_carbon',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total plant biomass carbon')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g C/m2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!------
!17
status = nf90_def_var(ncid,'soilc_fast',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','soil carbon, fast pool')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g C m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!------
!18
status = nf90_def_var(ncid,'soilc_slow',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','soil carbon, slow pool')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g C m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!------
!19
status = nf90_def_var(ncid,'acflux_fire',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','C flux to atmosphere due to fire')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g C m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!------
!20
status = nf90_def_var(ncid,'nind',nf90_float,[lon,lat,pft,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','individual density')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','indiv/ m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,npft,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---------
!21
status = nf90_def_var(ncid,'grid_gpp',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total gridcell gpp')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!22
status = nf90_def_var(ncid,'grid_autoresp',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total gridcell autotrophic respiration')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','g m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!23
status = nf90_def_var(ncid,'runoff',nf90_float,[lon,lat,month,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total gridcell mean monthly runoff')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','mm d-1')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,12,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!24
status = nf90_def_var(ncid,'snow_fraction',nf90_float,[lon,lat,month,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total gridcell mean monthly snow cover fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,12,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---
!25
status = nf90_def_var(ncid,'snow_depth',nf90_float,[lon,lat,month,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','total gridcell mean monthly snow depth')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','m')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,12,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!26
status = nf90_def_var(ncid,'abs_rad_annual',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','annual absorbed total radiation')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','GJ m-2')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!27
status = nf90_def_var(ncid,'grndwat_depth',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','depth to groundwater')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','m')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!28
status = nf90_def_var(ncid,'max_ald',nf90_float,[lon,lat,time],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','maximum active layer depth')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','m')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,1])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

!29
status = nf90_def_var(ncid,'peatland',nf90_float,[lon,lat],varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','gridcell treated as peatland')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','none')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',fill_value)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
if (status/=nf90_noerr) call handle_err(status)

!---

status = nf90_enddef(ncid)
if (status/=nf90_noerr) call handle_err(status)

!----

!write the dimension variables (except for time)

!lonvect = lonvect + 0.25
!latvect = latvect - 0.25

allocate(pftnum(npft))

forall (i=1:npft)
  pftnum(i) = i
end forall

status = nf90_put_var(ncid,1,lonvect(srtx:endx))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,2,latvect(srty:endy))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,3,[1,2])
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,4,pftnum)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,5,[1,2,3,4,5,6,7,8,9,10,11,12])
if (status/=nf90_noerr) call handle_err(status)

!----------------

deallocate(pftnum)

end subroutine netcdf_create
