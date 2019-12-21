!> @file water_hru.f90
!> file containing the subroutine water_hru
!> @author
!> modified by Javier Burguete

!> this subroutine compute pet and et using Priestly-Taylor and a coefficient
subroutine water_hru

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Exp
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8, parameter :: pet_alpha = 1.28
   real*8 :: d, gma, ho, tmpk
   integer :: j

   j = ihru
!! if the HRU is water compute only pet and et
!! using Priestly-Taylor and a coefficient
   albday = .08
   tmpk = tmpav(j) + 273.15
   d = Exp(21.255 - 5304. / tmpk) * 5304. / tmpk ** 2
   gma = d / (d + .68)
   ho = hru_ra(j) * (1. - albday) / 2.44
   pet_day = pet_alpha * ho * gma
   etday = .7 * pet_day

   return
end
