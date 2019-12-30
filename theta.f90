!> @file theta.f90
!> file containing the function theta
!> @author
!> modified by Javier Burguete

!> this function corrects rate constants for temperature.
!> Equation is III-52 from QUAL2E
!> @param[in] r20
!> value of the reaction rate coefficient at the standard temperature (20
!> degrees C) (1/day)
!> @param[in] thk
!> temperature adjustment factor (empirical constant for each reaction
!> coefficient) (none)
!> @param[in] tmp temperature on current day (deg C)
!> @return value of the reaction rate coefficient at the local temperature
!> (1/day)
real*8 function theta(r20, thk, tmp)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    r20         |1/day         |value of the reaction rate coefficient at
!!                               |the standard temperature (20 degrees C)
!!    thk         |none          |temperature adjustment factor (empirical
!!                               |constant for each reaction coefficient)
!!    tmp         |deg C         |temperature on current day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    theta       |1/day         |value of the reaction rate coefficient at
!!                               |the local temperature
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
   implicit none

   real*8, intent (in) :: r20, thk, tmp

   theta = r20 * thk ** (tmp - 20.)

   return
end
