!> @file crackflow.f90
!> file containing the subroutine crackflow
!> @author
!> modified by Javier Burguete

!> this surboutine modifies surface runoff to account for crack flow
subroutine crackflow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff for the hour in HRUS
!!    ievent      |none          |rainfall/runoff code
!!                               |0 daily rainfall/curve number technique
!!                               |1 sub-daily rainfall/Green&Ampt/hourly
!!                               |  routing
!!                               |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ihru        |none          |HRU number
!!    surfq(:)    |mm H2O        |surface runoff in the HRU for the day
!!    voltot      |mm            |total volume of cracks expressed as depth
!!                               |per unit area
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff for the hour in HRU
!!    surfq(:)    |mm H2O        |surface runoff in the HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    j           |none          |HRU number
!!    voli        |none          |volume available for crack flow
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 :: voli
   integer :: ii, j

   j = ihru

   !! subtract crack flow from surface runoff
   if (surfq(j) > voltot) then
      surfq(j) = surfq(j) - voltot
   else
      surfq(j) = 0.
   endif

   if (ievent > 0) then
      voli = voltot
      do ii = 1, nstep  !j.jeong 4/24/2009
         if (hhqday(ii) > voli) then
            hhqday(ii) = hhqday(ii) - voli
            voli = 0.
         else
            voli = voli - hhqday(ii)
            hhqday(ii) = 0.
         endif
      end do
   end if

   return
end
