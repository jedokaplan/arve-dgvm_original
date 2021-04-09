module sitedatamod
!Used by ARVEpoint to take in the climate data from a text file.

use arveparams, only  : mdays,dp

public :: openfile
public :: readfile

public :: tsoil,date,year
private :: simstart,month,d,tmean,tmax,tmin,prec,ra_dp,tfirst,tmeanfirst

real(dp) :: tsoil
character(10) :: date

logical :: simstart

integer :: year
integer :: month
integer :: d
real(dp) :: tmean
real(dp) :: tmax
real(dp) :: tmin
real(dp) :: prec
real(dp) :: sunp
real(dp) :: ra_dp
real(dp) :: tfirst
real(dp) :: tmeanfirst

contains

!---------------------------------
subroutine openfile(infile)

implicit none

character(90), intent(in) :: infile
!---------

open(10,file=infile,status='old')  ! statue old ensures that it uses an exisiting file

read(10,*)year,month,d,tmean,tmax,tmin,prec,ra_dp,tsoil

simstart = .true.
tfirst = tmin
tmeanfirst = tmean

backspace(10)

end subroutine openfile

!---------------------------------
subroutine readfile()

use metvarsmod, only : met,dm
use soilstate_vars, only : surf

implicit none

!-------

5 continue

read(10,*,end=10)year,month,d,tmean,tmax,tmin,prec,sunp,ra_dp,tsoil
!write(0,*)year

if (month /= 1 .and. simstart) then
  goto 5 !call readfile(day)
else
  simstart = .false.
end if


met(1,0)%temp24 = tmean
met(1,0)%prec = prec
met(1,0)%dswb = ra_dp
met(1,0)%tmax = tmax
met(1,0)%tmin = tmin
met(1,0)%cldf = 1. - sunp
surf%realsoilT = tsoil

!----
!get next day's weather
read(10,*,end=8)year,month,d,tmean,tmax,tmin,prec,sunp,ra_dp,tsoil

!write(0,*)year

dm(1)%tnext = tmin
dm(1)%tmean_next = tmean
backspace(10)

return

!----
8 continue

dm(1)%tnext = tfirst
dm(1)%tmean_next = tmeanfirst

rewind(10)

!read(10,*) !header

return

!------------
!this part handles the end-of-file condition

10 continue

rewind(10)

!read(10,*) !header

goto 5

end subroutine readfile

!--------------------------------

end module sitedatamod

