!> @file ascrv.f90
!> file containing the subroutine ascrv
!> @author
!> modified by Javier Burguete
!>
!> this subroutine computes shape parameters x5 and x6 for the S curve equation
!> \f\[x=\frac{y}{y+\exp(x5+x6\,y)}\f\]
!> given 2 (x,y) points along the curve. x5 is determined by solving the
!> equation with \f$x\f$ and \f$y\f$ values measured around the midpoint of the
!> curve (approx. 50% of the maximum value for \f$x\f$) and x6 is determined by
!> solving the equation with \f$x\f$ and \f$y\f$ values measured close to one of
!> the endpoints of the curve (100% of the maximum value for \f$x\f$). This
!> subroutine is called from readbsn.f90 and readplant.f90
!> @param[in] x1
!> value for x in the above equation for first datapoint, x1 should be close to
!> 0.5 (the midpoint of the curve)
!> @param[in] x2
!> value for x in the above equation for second datapoint, x2 should be close to
!> 0.0 or 1.0
!> @param[in] x3
!> value for y in the above equation corresponding to x1
!> @param[in] x4
!> value for y in the above equation corresponding to x2
!> @param[out] x5
!> 1st shape parameter for S curve equation characterizing the midpoint of the
!> curve
!> @param[out] x6
!> 2nd shape parameter for S curve equation characterizing the regions close to
!> the endpoints of the curve
subroutine ascrv(x1,x2,x3,x4,x5,x6)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes shape parameters x5 and x6 for the S curve
!!    equation x = y/(y + exp(x5 + x6*y)) given 2 (x,y) points along the curve.
!!    x5 is determined by solving the equation with x and y values measured
!!    around the midpoint of the curve (approx. 50% of the maximum value for x)
!!    and x6 is determined by solving the equation with x and y values measured
!!    close to one of the endpoints of the curve (100% of the maximum value for
!!    x) This subroutine is called from readbsn.f90 and readplant.f90

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    x1          |none          |value for x in the above equation for first
!!                               |datapoint, x1 should be close to 0.5 (the
!!                               |midpoint of the curve)
!!    x2          |none          |value for x in the above equation for second
!!                               |datapoint, x2 should be close to 0.0 or 1.0
!!    x3          |none          |value for y in the above equation correspond-
!!                               |ing to x1
!!    x4          |none          |value for y in the above equation correspond-
!!                               |ing to x2
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    x5          |none          |1st shape parameter for S curve equation
!!                               |characterizing the midpoint of the curve
!!    x6          |none          |2nd shape parameter for S curve equation
!!                               |characterizing the regions close to the
!!                               |endpoints of the curve
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    xx          |none          |temp variable, used to hold calculated
!!                               |value needed in later equations
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
   implicit none

   real*8 :: xx

   real*8, intent (in) :: x1, x2, x3, x4
   real*8, intent (out) :: x5, x6

   xx = Log(x3/x1 - x3)
   x6 = (xx - Log(x4/x2 - x4)) / (x4 - x3)
   x5 = xx + (x3 * x6)

   return
end
