module insolationmod
!This module contains basic constants, a subroutine containing the 1978 Berger
!solution for calculating the Earth's orbital parameters, and a subroutine for
!calculating daylength and top of the atmosphere insolation
!given year, day of year, and latitude
!The year 1950 is considered ka = 0 with positive numbers going back in time

!call orbital_parameters(yr) once each year to initialize orbital parameters
!then call

!NEEDS MORE COMMENTS IN CODE!

use arveparams, only : dp

implicit none

public :: orbital_parameters
public :: dayins

!orbital parameters
type orbitpars
  real(dp) :: ecc
  real(dp) :: pre
  real(dp) :: perh
  real(dp) :: xob
end type orbitpars

type(orbitpars), save :: op

contains

!----------------------------------------------
subroutine orbital_parameters(ka)

!This routine uses the orbital solutions of Berger (1978) and is valid only
!for calculations within  +  -  1.000.000 yr centered on 1950 AD.
!For longer periods the Berger (1990) solution should be used.
!(Contact Berger for this 1990 solution).
!
!Recoded by J.O. Kaplan in 2002.
!
!Please refer to :
!  Berger A., 1978. A simple algorithm to compute long term
!                   variations of daily or monthly insolation.
!                   Contr. 18  Inst. of Astronomy and Geophysics,
!                   Universite Catholique de Louvain,
!                   Louvain - la - Neuve, Belgium
!
!  Berger A., 1978. Long term variations of daily insolation and
!                   Quaternary climatic changes.
!                   J. of Atmospheric Sciences 35, 2362 - 2367
!
!The function value returned by atan is assumed to be a real(dp)
!ranging from  - pi / 2 to pi / 2
!
!Input parameters for the orbital solution are provided in a separate module, insolaparms.f90
!----------------------------------------------

use arveparams,   only : pi,d2r
use insolaparms, only : ae,eccy,eccz,aob,obly,oblz,aop,prey,prez

implicit none

!argument
integer, intent(in)  :: ka

!parameters
integer, parameter :: nef  =  19
integer, parameter :: nob  =  47
integer, parameter :: nop  =  78
real(dp), parameter :: d2rr  =  d2r / 3600.0_dp
real(dp), parameter :: step  =  360.0_dp / 365.25_dp

!local variables
integer :: i
integer :: neff
integer :: nobb
integer :: nopp
real(dp) :: xod
real(dp) :: xop
real(dp) :: prm
real(dp) :: t
real(dp) :: xes
real(dp) :: xec
real(dp) :: arg
real(dp) :: tra
real(dp) :: rp
real(dp) :: prg
real(dp), dimension(19) :: be
real(dp), dimension(19) :: ce
real(dp), dimension(47) :: bob
real(dp), dimension(47) :: cob
real(dp), dimension(78) :: bop
real(dp), dimension(78) :: cop
real(dp) :: ecc
real(dp) :: pre
real(dp) :: perh
real(dp) :: xob

!------------------------------------------------
!   daily insolation  -  long term variation  *
!------------------------------------------------
!
!
!This program computes the total daily irradiation received at the top
!of the atmosphere for a given latitude and time in the ka (in kj m - 2).

!
!   1.earth orbital elements : eccentricity           ecc   table 1

!                            precessional parameter pre
!                              obliquity              xob   table 2
!                              general precession     prg
!                              longitude perihelion   perh  table 3


!The parameters for the orbital solution are contained in a separate module
!called insolaparms.mod. It is USEd in the header of this program.

!read amplitude a  mean rate b  phase c
!they are immediately converted in radians
!
!nef  nob  nop  may be reduced to  19  18  9
!but the input data must be changed accordingly

!---------------------------------------------
!eccentricity

be  =  eccy * d2rr
ce  =  eccz * d2r

!--------------------------------------------
!obliquity

xod  =  23.320556_dp

bob  =  obly  *  d2rr
cob  =  oblz  *  d2r

!------------------------------------------
!general precession in longitude

xop  =   3.392506_dp
prm  =  50.439273_dp

bop  =  prey  *  d2rr
cop  =  prez  *  d2r

!---------------------------------------
neff = nef
nobb = nob
nopp = nop

!   3.numerical value for ecc pre xob
!--------------------------------------
!   t is negative for the past
t = real(ka) * 1000.0_dp
xes = 0.0_dp
xec = 0.0_dp

do i = 1,neff
  arg = be(i) * t + ce(i)
  xes = xes + ae(i) * sin(arg)
  xec = xec + ae(i) * cos(arg)
end do

ecc = sqrt(xes * xes + xec * xec)
tra = abs(xec)

if(tra > 1.0e-8) then

  rp = atan(xes / xec)

  if(xec > 0.0_dp) then !line 12

    if (xes > 0.0_dp) then !line 13

      perh = rp / d2r

    else if (xes < 0.0_dp) then !line 14

      rp = rp + 2.0_dp * pi
      perh = rp / d2r

    else !line 13

      perh = rp / d2r

    end if

  else if (xec < 0.0_dp) then !line 11

    rp = rp + pi
    perh = rp / d2r

  else !line 10

    if (xes > 0.0_dp) then !line 17

      rp = pi / 2.0_dp
      perh = rp / d2r

    else if (xes < 0.0_dp) then !line 15

      rp = 1.5_dp * pi
      perh = rp / d2r

    else !line 16

      rp = 0.0_dp
      perh = rp / d2r

    end if
  end if
else
  if (xes>0.0_dp) then !line 17

    rp = pi / 2.0_dp
    perh = rp / d2r

  else if (xes<0.0_dp) then !line 15

    rp = 1.5_dp * pi
    perh = rp / d2r

  else !line 16

    rp = 0.0_dp
    perh = rp / d2r

  end if
end if

prg = prm * t

do i = 1,nop

  arg = bop(i) * t + cop(i)
  prg = prg + aop(i) * sin(arg)

end do

prg = prg / 3600.0_dp + xop
perh = perh + prg

if (perh > 0.0_dp) then !line 53

  if(perh > 360.0_dp) then

    perh = perh - 360.0_dp

  end if

else if (perh < 0.0_dp) then

  perh = perh + 360.0_dp

end if

pre = ecc * sin(perh * d2r)

xob = xod

do i = 1,nobb

  arg = bob(i) * t + cob(i)
  xob = xob + aob(i) / 3600.0_dp * cos(arg)

end do

!write out to shared values
op%ecc  = ecc
op%pre  = pre
op%perh = perh
op%xob  = xob

end subroutine orbital_parameters

!----------------------------------------------------------------------------------------------------
subroutine dayins(j,nd)

!this subroutine calculates top-of-the-atmosphere insolation given orbital
!parameters, day of the year, and latitude
!output : ww = KJ m-2 day-1  dayl = length of day (hours)
!NOTE: presently this is setup to calculate the NEXT days values (i.e. d+1)

use arveparams, only : pi,d2r
use statevars, only : sv,dayl

implicit none

!arguments
integer, intent(in) :: j   !gridcell index
integer, intent(in) :: nd  !day of the year

!pointers
real(dp), pointer :: lat                !latitude (degrees)
real(dp), pointer :: sdecl_n            !the next day's solar declination (deg)
real(dp), pointer :: dsw_tn             !the next day's dsw_t

!parameters
real(dp), parameter :: ss   = 1353.0_dp
real(dp), parameter :: step =  360.0_dp / 365.25_dp
real(dp), parameter :: tau  =   86.4_dp
real(dp), parameter :: test =    1.0e-8

!local variables
real(dp) :: adelta  !absolute value of the solar zenith angle
real(dp) :: anm
real(dp) :: anv
real(dp) :: aphi    !absolute value of the latitude
real(dp) :: at
real(dp) :: cd
real(dp) :: cp
real(dp) :: dlam
real(dp) :: dlamm
real(dp) :: ecc
real(dp) :: perh
real(dp) :: ranm
real(dp) :: ranv
real(dp) :: rau
real(dp) :: rdayl
real(dp) :: rdelta   !declination in radians
real(dp) :: rlam
real(dp) :: rphi
real(dp) :: s
real(dp) :: sd
real(dp) :: sf
real(dp) :: so
real(dp) :: spa
real(dp) :: spd
real(dp) :: stp
real(dp) :: tls
real(dp) :: tp
real(dp) :: tt
real(dp) :: xec
real(dp) :: xee
real(dp) :: xl
real(dp) :: xlam
real(dp) :: xllp
real(dp) :: xob
real(dp) :: xse
real(dp) :: delta       !solar declination (deg)
real(dp) :: ww          !top of the atmosphere radiation (KJ m-2 day-1)
integer :: d            !day to calculate

!point pointers
ecc  = op%ecc
perh = op%perh
xob  = op%xob
lat   => sv(j)%lat
dsw_tn=> sv(j)%dsw_tn
sdecl_n => sv(j)%sdecl_n

!-------------------------------------------
!correct for the next timestep adjustment
if (nd == 366) then
 d = 1
else if (nd == 0) then
 d = 1
else
 d = nd
end if
!-------------------------------------------

sf   =  tau * ss / pi
so   =  sin(xob * d2r)
xl   =  perh + 180.0_dp

xllp  =  xl * d2r
xee   =  ecc * ecc
xse   =  sqrt(1.0_dp - xee)
xlam  =  (ecc / 2.0_dp + ecc * xee / 8.0_dp) * (1.0_dp + xse) * sin(xllp) - xee / 4.0_dp * (0.5_dp + xse) &
         * sin(2.0_dp * xllp) + ecc * xee / 8.0_dp * (1.0_dp / 3.0_dp + xse) * sin(3.0_dp * xllp)

xlam  =  2.0_dp * xlam / d2r
dlamm =  xlam + (d - 80) * step
anm   =  dlamm - xl

ranm  =  anm * d2r
xec   =  xee * ecc
ranv  =  ranm + (2.0_dp * ecc - xec / 4.0_dp) * sin(ranm) + 5.0_dp / 4.0_dp * ecc * ecc *    &
             sin(2.0_dp * ranm) + 13.0_dp / 12.0_dp * xec * sin(3.0_dp * ranm)

anv   =  ranv / d2r
tls   =  anv + xl

dlam  =  tls

!----------------------------------------------------------------------------------------------------

rphi    =  lat * d2r
ranv    =  (dlam - xl) * d2r
rau     =  (1.0_dp - ecc * ecc) / (1.0_dp + ecc * cos(ranv))

s       =  sf / rau / rau
rlam    =  dlam * d2r
sd      =  so * sin(rlam)
cd      =  sqrt(1.0_dp - sd * sd)

rdelta  =  atan(sd / cd)
delta   =  rdelta / d2r
spa     =  sd * sin(rphi)

cp      =  cd * cos(rphi)
aphi    =  abs(lat)
adelta  =  abs(delta)

!singularity for aphi = 90 and delta = 0
!particular cases for lat = 0 or delta = 0

tt = abs(aphi - 90.0_dp)

if (tt <= test .and. adelta <= test) then

  dayl(2) = 0.00_dp
  ww = 0.00_dp

else if (adelta <= test) then

  dayl(2) = 12.0_dp
  ww = s * cos(rphi)

else if (aphi <= test) then

  dayl(2) = 12.0_dp
  ww = s * cos(rdelta)

else

  at = 90.0_dp - adelta
  spd = lat * delta

  if (aphi <= at) then

    tp = -spa / cp
    stp = sqrt(1.0_dp - tp * tp)
    rdayl = acos(tp)
    dayl(2) = 24.0_dp * rdayl / pi
    ww = s * (rdayl * spa + cp * stp)

  else if (spd > 0.0_dp) then

    dayl(2) = 24.00_dp
    ww = s * spa * pi

  else if (spd < 0.0_dp) then

    dayl(2) = 0.00_dp
    ww = 0.00_dp

  else

    tp =  - spa / cp
    stp = sqrt(1.0_dp - tp * tp)
    rdayl = acos(tp)
    dayl(2) = 24.0_dp * rdayl / pi
    ww = s * (rdayl * spa + cp * stp)

  end if
end if

dayl(2) = dayl(2) * 3600._dp
dsw_tn = ww
sdecl_n = delta

!write(*,*)'insol',dayl(2),lat,delta,ww,spd

end subroutine dayins

!----------------------------------------------------------------------------------------------------

end module insolationmod
