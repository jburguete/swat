!> @file ee.f90
!> file containing the subroutine ee
!> @author
!> modified by Javier Burguete

!> this function calculates saturation vapor pressure at a given
!> air temperature
!> @param[in] tk mean air temperature (deg C)
!> @return saturation vapor pressure (kPa)
real*8 function ee(tk) result (r_ee)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    tk          |deg C         |mean air temperature
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ee          |kPa           |saturation vapor pressure
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   implicit none
   real*8, intent (in) :: tk

   r_ee = 0.
   if (tk + 237.3 /= 0.) then
      r_ee = (16.78 * tk - 116.9) / (tk + 237.3)
      r_ee = Exp(r_ee)
   end if

   return
end
