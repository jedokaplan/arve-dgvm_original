module randomdistmod

use arveparams, only : dp

implicit none

!type

integer, parameter :: lg = selected_int_kind(9) !define a long integer type (32 bits)

!functions

public :: ran_seed
public :: ran_uni
public :: ran_normal
public :: ran_gamma

private :: random_gamma1
private :: random_gamma2
private :: random_exponential

!parameters

real(dp), parameter :: one    = 1._dp
real(dp), parameter :: half   = 0.5_dp
real(dp), parameter :: vsmall = tiny(1._dp)
real(dp), parameter :: zero   = 0._dp

!global variable

integer(lg) :: seed

contains

!-------------------------------

subroutine ran_seed(put,get)

implicit none

!arguments

integer(lg), optional, intent(in)  :: put
integer(lg), optional, intent(out) :: get

!module variable
!integer(lg) :: seed !the random seed (negative integer, use once to initialize)

!----------

if (present(put)) then
  if (put > 0) then
    seed = put
  else
    seed = -put
  end if
else if (present(get)) then
  get = seed
else
  seed = -5
end if

end subroutine ran_seed

!-------------------------------

subroutine ran_uni(harvest)

!Uniform random number generator, creates values between 0 and 1, exclusive of endpoints.
!Subroutine adapted from Numerical Recipes in Fortran 90, Press et al. 1996, pg. 1142

implicit none

!argument

real(dp), intent(out) :: harvest

!module variable
!integer(lg) :: seed   !the random seed (negative integer, should be initialized once at start)

!local variables

integer(lg), parameter :: ia =      16807
integer(lg), parameter :: im = 2147483647
integer(lg), parameter :: iq =     127773
integer(lg), parameter :: ir =       2836

real(dp), save :: am

integer(lg), save :: ix = -1
integer(lg), save :: iy = -1
integer(lg), save :: k

!----------

if (seed <= 0 .or. iy < 0) then
  am = nearest(1.0, -1.0) / im
  iy = ior(ieor(888889999,abs(seed)),1)
  ix = ieor(777755555,abs(seed))
  seed = abs(seed) + 1
end if

ix = ieor(ix,ishft(ix, 13))
ix = ieor(ix,ishft(ix,-17))
ix = ieor(ix,ishft(ix,  5))

k = iy / iq

iy = ia * (iy - k * iq) - ir * k

if (iy < 0) iy = iy + im

harvest = am * ior(iand(im,ieor(ix,iy)),1)

end subroutine ran_uni

!-------------------------------

subroutine ran_normal(harvest)

!Randomly samples the normal distribution with zero mean and unit variance.
!Subroutine adapted from Numerical Recipes in Fortran 90, Press et al. 1996, pg. 1152

implicit none

!argument

real(dp), intent(out) :: harvest

!module variable
!integer(lg) :: seed   !the random seed (negative integer, should be initialized once at start)

!local variables

real(dp) :: v1
real(dp) :: v2
real(dp) :: rsq

real(dp), save :: g
logical,  save :: gaus_stored = .false.

!----------

if (gaus_stored) then
  harvest = g
  gaus_stored = .false.
else
  do
    call ran_uni(v1)
    call ran_uni(v2)

    v1 = 2. * v1 - 1.
    v2 = 2. * v2 - 1.

    rsq = v1**2 + v2**2
    if (rsq > 0. .and. rsq < 1.) exit
  end do

  rsq = sqrt(-2. * log(rsq) / rsq)

  harvest = v1 * rsq
  g = v2 * rsq
  gaus_stored = .true.

end if

end subroutine ran_normal

!----------------------------------------------------------------------------------------------------------

!The functions below are for randomly sampling the gamma distribution

!adapted from code by Alan J. Miller, URL,
!http://users.bigpond.net.au/amiller/random.html

!----------------------------------------------------------------------------------------------------------

function ran_gamma(s,b,first) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     function generates a random gamma variate.
!     calls either random_gamma1 (s > 1.0)
!     or random_exponential (s = 1.0)
!     or random_gamma2 (s < 1.0).

!     s = shape parameter of distribution (0 < real)
!     b = scale parameter

implicit none

!arguments

real(dp), intent(in) :: s       !shape parameter of the Gamma distribution (alpha, unitless)
real(dp), intent(in) :: b       !scale parameter of the Gamma distribution (Beta)
logical,  intent(in) :: first   !flag if this is the first call to the distribution
real(dp)             :: fn_val

!--------

if (s <= zero) then
  write(*, *) 'shape parameter value must be positive'
  stop
end if

if (s > one) then
  fn_val = random_gamma1(s, first)
else if (s < one) then
  fn_val = random_gamma2(s, first)
else
  fn_val = random_exponential()
end if

!scale the random variable with Beta

fn_val = b * fn_val

end function ran_gamma

!--------------------------------------------------------------------

function random_gamma1(s,first) result(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s > 1.

implicit none

!arguments

real(dp), intent(in) :: s
logical,  intent(in) :: first
real(dp)             :: fn_val

!local variables

real(dp), save  :: c
real(dp), save  :: d
real(dp)        :: u
real(dp)        :: v
real(dp)        :: x

!--------

if (first) then
  d = s - one/3.
  c = one/sqrt(9.0*d)
end if

! start of main loop
do

! generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  do
    call ran_normal(x)
    v = (one + c*x)**3
    if (v > zero) exit
  end do

! generate uniform variable u

  call ran_uni(u)
  if (u < one - 0.0331*x**4) then
    fn_val = d*v
    exit
  else if (log(u) < half*x**2 + d*(one - v + log(v))) then
    fn_val = d*v
    exit
  end if
end do

end function random_gamma1

!--------------------------------------------------------------------

function random_gamma2(s,first) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function generates a random variate in [0,infinity) from
! a gamma distribution with density proportional to
! gamma2**(s-1) * exp(-gamma2),
! using a switching method.

!    s = shape parameter of distribution
!          (real < 1.0)

implicit none

!arguments

real(dp), intent(in) :: s
logical,  intent(in) :: first
real(dp)             :: fn_val

!local variables

real(dp)       :: r
real(dp)       :: x
real(dp)       :: w
real(dp), save :: a
real(dp), save :: p
real(dp), save :: c
real(dp), save :: uf
real(dp), save :: vr
real(dp), save :: d

!--------

if (s <= zero .or. s >= one) then
  write(*, *) 'shape parameter value outside permitted range'
  stop
end if

if (first) then                        ! initialization, if necessary
  a = one - s
  p = a/(a + s*exp(-a))
  if (s < vsmall) then
    write(*, *) 'shape parameter value too small'
    stop
  end if
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*log(a)
end if

do
  call ran_uni(r)
  if (r >= vr) then
    cycle
  else if (r > p) then
    x = a - log((one - r)/(one - p))
    w = a*log(x)-d
  else if (r > uf) then
    x = a*(r/p)**c
    w = x
  else
    fn_val = zero
    return
  end if

  call ran_uni(r)
  if (one-r <= w .and. r > zero) then
    if (r*(w + one) >= one) cycle
    if (-log(r) <= w) cycle
  end if
  exit
end do

fn_val = x

end function random_gamma2

!--------------------------------------------------------------------

function random_exponential() result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function generates a random variate in [0,infinity) from
! a negative exponential distribution wlth density proportional
! to exp(-random_exponential), using inversion.

implicit none

!argument

real(dp) :: fn_val

!local variable

real(dp) :: r

!--------

do
  call ran_uni(r)
  if (r > zero) exit
end do

fn_val = -log(r)

end function random_exponential

!-------------------------------

end module randomdistmod
