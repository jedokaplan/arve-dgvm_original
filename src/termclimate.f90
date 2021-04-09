subroutine termclimate()
!going to be connencted up ? JM Feb 26 08

use netcdf_error
use iovariables

implicit none

!-------------------------
!deallocate memory

deallocate(indexmat)
deallocate(ibuf)

!-------------------------
!close save file

close(20)

!-------------------------
!close climate files

status = nf90_close(tempfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(precfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(sunpfid)
if (status /= nf90_noerr) call handle_err(status)

end subroutine termclimate
