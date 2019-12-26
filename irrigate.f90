!> @file irrigate.f90
!> file containing the subroutine irrigate
!> @author
!> modified by Javier Burguete

!> this subroutine applies irrigation water to HRU
!> @param[in] j HRU number (none)
!> @param[in] volmm depth irrigation water applied to HRU (mm H2O)
subroutine irrigate(j, volmm)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    volmm       |mm H2O        |depth irrigation water applied to HRU
!!    aairr(:)    |mm H2O        |average annual amount of irrigation water
!!                               |applied to HRU
!!    curyr       |none          |current year of simulation
!!    irn(:)      |none          |average annual number of irrigation
!!                               |applications in HRU
!!    nyskip      |none          |number of years to skip output summarization
!!                               |and printing
!!    sol_fc(:,:) |mm H2O        |amount of water available to plants in soil
!!                               |layer at field capacity (fc - wp)
!!    sol_nly(:)  |none          |number of soil layers
!!    sol_st(:,:) |mm H2O        |amount of water stored in the soil layer
!!                               |on any given day (less wp water)
!!    hrumono(22,:)|mm H2O       |amount of irrigation water applied to HRU
!!                               |during month
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    aairr(:)    |mm H2O        |average annual amount of irrigation water
!!                               |applied to HRU
!!    aird(:)     |mm H2O        |amount of water applied to HRU on current
!!                               |day
!!    irn(:)      |none          |average annual number of irrigation
!!                               |applications in HRU
!!    sol_st(:,:) |mm H2O        |amount of water stored in the soil layer
!!                               |on any given day (less wp water)
!!    sol_sw(:)   |mm H2O        |amount of water stored in the soil profile
!!                               |on any given day
!!    hrumono(22,:)|mm H2O       |amount of irrigation water applied to HRU
!!                               |during month
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: j
   real*8, intent (in) :: volmm

!! initialize variable for HRU
!! (because irrigation can be applied in different command loops
!! the variable is initialized here)

   aird(j) = volmm * (1. - sq_rto)
   qird(j) = volmm * sq_rto

!! summary calculations
   if (curyr > nyskip) then
      irn(j) = irn(j) + 1
      aairr(j) = aairr(j) + aird(j)
      hrumono(22,j) = hrumono(22,j) + aird(j)
   end if


   return
end
