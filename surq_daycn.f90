!> @file surq_daycn.f90
!> file containing the subroutine surq_daycn
!> @author
!> modified by Javier Burguete

!> predicts daily runoff given daily precipitation and snow melt
!> using a modified SCS curve number approach
!> @param[in] j HRU number (none)
subroutine surq_daycn(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    cnday(:)    |none          |curve number for current day, HRU and at
!!                               |current soil moisture
!!    fcimp(:)    |fraction      |fraction of HRU area that is classified
!!                               |as directly connected impervious
!!    iurban(:)   |none          |urban simulation code:
!!                               |0  no urban sections in HRU
!!                               |1  urban sections in HRU, simulate using USGS
!!                               |   regression equations
!!                               |2  urban sections in HRU, simulate using build
!!                               |   up/wash off algorithm
!!    precipday   |mm H2O        |precipitation for the day in HRU
!!    urblu(:)    |none          |urban land type identification number from
!!                               |urban.dat
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    surfq(:)    |mm H2O        |surface runoff for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bb          |none          |variable used to store intermediate
!!                               |calculation result
!!    cnimp       |none          |curve number for impervious areas
!!    pb          |none          |variable used to store intermediate
!!                               |calculation result
!!    r2          |none          |retention parameter in CN equation
!!    surfqimp    |mm H2O        |surface runoff from impervious area
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: j
   real*8, parameter :: cnimp = 98.
   real*8 :: bb, pb, r2, surfqimp

   r2 = 25400. / cnday(j) - 254.
   bb = .2 * r2
   pb = precipday - bb

   if (pb > 0.) then
      surfq(j) = pb * pb / (precipday + .8 * r2)
   end if


   if (iurban(j) > 0) then
      surfqimp = 0.
      r2 = 25400. / cnimp - 254.
      bb = .2 * r2
      pb = precipday - bb
      if (pb > 0.) then
         surfqimp = pb * pb / (precipday + .8 * r2)
      end if
      surfq(j) = surfq(j) * (1. - fcimp(urblu(j))) + surfqimp * fcimp(urblu(j))
   end if

   return
end
