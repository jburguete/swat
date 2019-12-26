!> @file water_hru.f90
!> file containing the subroutine water_hru
!> @author
!> modified by Javier Burguete

!> this subroutine compute pet and et using Priestly-Taylor and a coefficient
!> @param[in] j HRU number
subroutine water_hru(j)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    d
!!    gma         |kPa/deg C     |psychrometric constant
!!    ho          |              |net radiation
!!    pet_alpha   |none          |alpha factor in Priestley-Taylor ET
!!                               |equation
!!    tmpk        |deg K         |average temperature for the day in the HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8, parameter :: pet_alpha = 1.28
   real*8 :: d, gma, ho, tmpk

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
