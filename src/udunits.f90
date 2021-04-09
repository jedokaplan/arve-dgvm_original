module udtypes

!small module for making use of the time conversion functions available in the Unidata udunits library
!link to libudunits.a
!Jed Kaplan, December 2007

implicit none

integer, parameter :: long = selected_int_kind(13)
integer, parameter :: sp   = kind(1.0)
integer, parameter :: dp   = selected_real_kind(2*precision(1.0_sp))

type ptr
  integer(long) :: address
end type ptr

end module udtypes

!----------

module udunits

interface

integer function utopen(path)
  character(*), intent(in) :: path
end function utopen

type(ptr) function utmake()
  use udtypes
end function utmake

integer function uttime(unit)
  use udtypes
  type(ptr), intent(in) :: unit
end function uttime

integer function utorigin(unit)
  use udtypes
  type(ptr), intent(in) :: unit
end function utorigin

subroutine utclr(unit)
  use udtypes
  type(ptr), intent(out) :: unit
end subroutine utclr

subroutine utcpy(source,dest)
  use udtypes
  type(ptr), intent(in)  :: source
  type(ptr), intent(out) :: dest
end subroutine utcpy

subroutine utorig(source,amount,udresult)
  use udtypes
  type(ptr), intent(in)  :: source
  real(dp),  intent(in)  :: amount
  type(ptr), intent(out) :: udresult
end subroutine utorig

subroutine utscal(source,factor,udresult)
  use udtypes
  type(ptr), intent(in)  :: source
  real(dp),  intent(in)  :: factor
  type(ptr), intent(out) :: udresult
end subroutine utscal

subroutine utmult(term1,term2,udresult)
  use udtypes
  type(ptr), intent(in)  :: term1
  type(ptr), intent(in)  :: term2
  type(ptr), intent(out) :: udresult
end subroutine utmult

subroutine utinv(source,udresult)
  use udtypes
  type(ptr), intent(in)  :: source
  type(ptr), intent(out) :: udresult
end subroutine utinv

subroutine utdiv(numer,denom,udresult)
  use udtypes
  type(ptr), intent(in)  :: numer
  type(ptr), intent(in)  :: denom
  type(ptr), intent(out) :: udresult
end subroutine utdiv

subroutine utexp(source,power,udresult)
  use udtypes
  type(ptr), intent(in)  :: source
  integer,   intent(in)  :: power
  type(ptr), intent(out) :: udresult
end subroutine utexp

integer function utdec(spec,unit)
  use udtypes
  character, intent(in)  :: spec
  type(ptr), intent(out) :: unit
end function utdec

integer function utcaltime(value,unit,year,month,day,hour,minute,second)
  use udtypes
  real(dp),  intent(in)  :: value
  type(ptr), intent(in)  :: unit
  integer,   intent(out) :: year
  integer,   intent(out) :: month
  integer,   intent(out) :: day
  integer,   intent(out) :: hour
  integer,   intent(out) :: minute
  real,      intent(out) :: second
end function utcaltime

integer function uticaltime(year,month,day,hour,minute,second,unit,value)
  use udtypes
  integer,   intent(in)  :: year
  integer,   intent(in)  :: month
  integer,   intent(in)  :: day
  integer,   intent(in)  :: hour
  integer,   intent(in)  :: minute
  real,      intent(in)  :: second
  type(ptr), intent(out) :: unit
  real(dp),  intent(out) :: value
end function uticaltime

integer function utenc(unit,spec)
  use udtypes
  type(ptr),    intent(in)  :: unit
  character(*), intent(out) :: spec
end function utenc

integer function utcvt(from,udto,slope,intercept)
  use udtypes
  type(ptr), intent(in)  :: from
  type(ptr), intent(in)  :: udto
  real(dp),  intent(out) :: slope
  real(dp),  intent(out) :: intercept
end function utcvt

subroutine utfree(unit)
  use udtypes
  type(ptr), intent(in) :: unit
end subroutine utfree

subroutine utcls()
end subroutine utcls

end interface

end module udunits
