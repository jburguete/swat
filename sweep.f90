!> @file sweep.f90
!> file containing the subroutine sweep
!> @author
!> modified by Javier Burguete

!> the subroutine performs the street sweeping operation
!> @param[in] j HRU number (none)
subroutine sweep(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j              |none          |HRU number
!!    fr_curb        |none          |availability factor, the fraction of the
!!                                  |curb length that is sweepable
!!    dirtmx(:)      |kg/curb km    |maximum amount of solids allowed to
!!                                  |build up on impervious surfaces
!!    nro(:)         |none          |sequence number of year in rotation
!!    nsweep(:)      |none          |sequence number of street sweeping
!!                                  |operation within the year
!!    sweepeff       |none          |removal efficiency of sweeping
!!                                  |operation
!!    thalf(:)       |days          |time for the amount of solids on
!!                                  |impervious areas to build up to 1/2
!!                                  |the maximum level
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    nsweep(:)    |none          |sequence number of street sweeping operation
!!                                |within the year
!!    twash(:)     |days          |time that solids have built-up on streets
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dirt         |kg/curb km    |amount of solids built up on impervious
!!                                |surfaces
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: dirt

!! calculate amount of dirt on streets prior to sweeping
   dirt = dirtmx(urblu(j)) * twash(j) / (thalf(urblu(j)) + twash(j))

!! calculate reduced amount of solid built up on impervious areas
   dirt = dirt * (1. - fr_curb * sweepeff)
   if (dirt < 1.e-6) dirt = 0.

!! set time to correspond to lower amount of dirt
   twash(j) = thalf(urblu(j)) * dirt / (dirtmx(urblu(j)) - dirt)

   nsweep(j) = nsweep(j) + 1

   return
end
