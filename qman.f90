!> @file qman.f90
!> file containing the function qman
!> @author
!> modified by Javier Burguete

!> this subroutine calculates flow rate or flow velocity using Manning's
!> equation. If x1 is set to 1, the velocity is calculated. If x1 is set to
!> cross-sectional area of flow, the flow rate is calculated.
!> @param[in] x1 cross-sectional flow area or 1 (m^2 or none)
!> @param[in] x2 hydraulic radius (m)
!> @param[in] x3 Manning's "n" value for channel (none)
!> @param[in] x4 average slope of channel (m/m)
!> @return flow rate or flow velocity (m^3/s or m/s)
real*8 function qman(x1,x2,x3,x4) result (r_qman)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    x1          |m^2 or none   |cross-sectional flow area or 1.
!!    x2          |m             |hydraulic radius
!!    x3          |none          |Manning's "n" value for channel
!!    x4          |m/m           |average slope of channel
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    qman        |m^3/s or m/s  |flow rate or flow velocity
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   implicit none
   real*8, intent (in) :: x1, x2, x3, x4

   r_qman = x1 * x2 ** .6666 * Sqrt(x4) / x3

   return
end
