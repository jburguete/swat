!> @file albedo.f90
!> file containing the subroutine albedo
!> @author
!> modified by Javier Burguete

!> this subroutine calculates albedo in the HRU for the day
!> @param[in] j HRU number
subroutine albedo(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    laiday(:)   |m**2/m**2     |leaf area index
!!    sno_hru(:)  |mm H2O        |amount of water in snow in HRU on current day
!!    sol_alb(:)  |none          |albedo when soil is moist
!!    sol_cov(:)  |kg/ha         |amount of residue on soil surface
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    albday      |none          |albedo of ground for day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cej         |none          |constant
!!    eaj         |none          |soil cover index
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: j
   real*8, parameter :: cej = -5.e-5
   real*8 :: eaj

!! calculate albedo
   eaj = Exp(cej * (sol_cov(j) + .1))   !! equation 2.2.16 in SWAT manual

   if (sno_hru(j) <= .5) then
      !! equation 2.2.14 in SWAT manual
      albday = sol_alb(j)

      !! equation 2.2.15 in SWAT manual
      if (laiday(j) > 0.) albday = .23 * (1. - eaj) + sol_alb(j) * eaj
   else
      !! equation 2.2.13 in SWAT manual
      albday = 0.8
   end if

   return
end
