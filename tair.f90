!> @file tair.f90
!> file containing the function tair
!> @author
!> modified by Javier Burguete

!> this function approximates hourly air temperature from daily max and
!> min temperatures as documented by Campbell (1985)
!> @param[in] hr hour of the day (none)
!> @param[in] jj HRU number (none)
!> @return air temperature for hour in HRU (deg C)
real*8 function tair(hr, jj)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hr          |none          |hour of the day
!!    jj          |none          |HRU number
!!    tmn(:)      |deg C         |minimum temperature for the day in HRU
!!    tmp_hi(:)   |deg C         |last maximum temperature in HRU
!!    tmp_lo(:)   |deg C         |last minimum temperature in HRU
!!    tmx(:)      |deg C         |maximum temperature for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    tmp_hi(:)   |deg C         |last maximum temperature in HRU
!!    tmp_lo(:)   |deg C         |last minimum temperature in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Cos

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!!    subroutine developed by A. Van Griensven
!!    Hydrology-Vrije Universiteit Brussel, Belgium
!!    subroutine modified by SLN

   use parm
   implicit none

   integer, intent (in) ::  hr, jj

!! update hi or lo temperature depending on hour of day
   select case (hr)
    case (3)
      tmp_lo(jj) = tmn(jj)
    case (15)
      tmp_hi(jj) = tmx(jj)
   end select

!! SWAT manual equation 2.3.1
   tair = 0.5 * (tmp_hi(jj) + tmp_lo(jj) + (tmp_hi(jj) - tmp_lo(jj)&
      &* Cos(0.2618 * (hr - 15))))

   return
end
