!> @file dailycn.f90
!> file containing the subroutine dailycn
!> @author
!> modified by Javier Burguete

!> calculates curve number for the day in the HRU
!> @param[in] j HRU number (none)
subroutine dailycn(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    icn         |none          |CN method flag:
!!                               |(for testing alternative methods)
!!                               |0 use traditional SWAT method which bases
!!                               |  CN on soil moisture
!!                               |1 use alternative method which bases CN on
!!                               |  plant ET
!!                               |2 use tradtional SWAT method which bases
!!                               |  CN on soil moisture but rention is adjusted
!!                               |  for mildly-sloped tiled-drained watersheds
!!    sci(:)      |none          |retention coefficient for cn method based on
!!                               |plant ET
!!    smx(:)      |none          |retention coefficient for cn method based on
!!                               |soil moisture
!!    sol_sw(:)   |mm H2O        |amount of water stored in soil profile on
!!                               |any given day
!!    sol_tmp(2,:)|deg C         |daily average temperature of second soil layer
!!    wrt(1,:)    |none          |1st shape parameter for calculation of
!!                               |water retention
!!    wrt(2,:)    |none          |2nd shape parameter for calculation of
!!                               |water retention
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cnday(:)    |none          |curve number for current day, HRU and at
!!                               |current soil moisture
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    r2          |none          |retention parameter in CN equation
!!    xx          |none          |variable used to store intermediate
!!                               |calculation result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Dmax1

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: r2, xx

   xx = wrt(1,j) - wrt(2,j) * sol_sw(j)
   if (xx < -20.) xx = -20.
   if (xx > 20.) xx = 20.

   r2 = 0.
   if (icn <= 0) then
      !! traditional CN method (function of soil water)
      if ((sol_sw(j) + Exp(xx)) > 0.001) then
         r2 = smx(j) * (1. - sol_sw(j)/(sol_sw(j)+Exp(xx)))
      end if
   else
      !! alternative CN method (function of plant ET)
      r2 = dmax1(3., sci(j))
   end if

   if (sol_tmp(2,j) <= 0.) r2 = smx(j) * (1. - Exp(- cn_froz * r2))
   r2 = dmax1(3.,r2)

   cnday(j) = 25400. / (r2 + 254.)
   sol_cnsw(j) = sol_sw(j)

   return
end
