module netcdf_error

  use typeSizes
  use netcdf

  implicit none

  integer :: status

  contains

  !-------------------------------------

  subroutine handle_err(status)

  !Internal subroutine - checks error status after each netcdf call,
  !prints out text message each time an error code is returned.

  integer, intent(in) :: status

  if(status /= nf90_noerr) then
    write(0,*)'netCDF error: ',trim(nf90_strerror(status))
    stop
  end if

  end subroutine handle_err

end module netcdf_error
