!> @file nuts.f90
!> file containing the subroutine nuts
!> @author
!> modified by Javier Burguete

!> this function calculates the plant stress factor caused by limited
!> supply of nitrogen or phosphorus
!> @param[in] u1 actual amount of element in plant (kg/ha)
!> @param[in] u2 optimal amount of element in plant (kg/ha)
!> @param[out] uu
!> fraction of optimal plant growth achieved where reduction is caused by plant
!> element deficiency (none)
subroutine nuts(u1, u2, uu)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    u1          |kg/ha         |actual amount of element in plant
!!    u2          |kg/ha         |optimal amount of element in plant
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    uu          |none          |fraction of optimal plant growth achieved
!!                               |where reduction is caused by plant element
!!                               |deficiency
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   implicit none

   real*8, intent (in) :: u1, u2
   real*8, intent (out) :: uu

   uu = 200. * (u1 / (u2 + .0001) - .5)

   if (uu <= 0.) then
      uu = 0.
   else
      if (uu < 99.) then
         uu = uu / (uu + Exp(3.535 - .02597 * uu))
      else
         uu = 1.
      endif
   end if

   if (u2 <= 1.e-6) uu = 1.

   return
end
