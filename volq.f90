!> @file volq.f90
!> file containing the subroutine volq
!> @author
!> modified by Javier Burguete

!> call subroutines to calculate the current day's CN for the HRU and
!> to calculate surface runoff
!> @param[in] j HRU number (none)
subroutine volq(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ievent      |none          |rainfall/runoff code
!!                               |0 daily rainfall/curve number technique
!!                               |1 sub-daily rainfall/Green&Ampt/hourly
!!                               |  routing
!!                               |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: surq_daycn, surq_greenampt

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: j

!! Compute surface runoff for day
   select case (ievent)
    case (0)
      call surq_daycn(j)
    case (1)
      call surq_greenampt(j)
      !  call dir_rnff
      !case (3)
      !  call surq_hourly
   end select

   return
end
