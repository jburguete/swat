!> @file vbl.f90
!> file containing the subroutine vbl
!> @author
!> modified by Javier Burguete

!> this subroutine checks the water and sediment balance for ponds
!> and reservoirs at the end of a simulation
!> @param[in] evx evaporation from water body
!> @param[in] spx seepage from water body
!> @param[in] pp precipitation on water body
!> @param[in] qin water entering water body
!> @param[in] ox water leaving water body
!> @param[inout] vx1 
!> (in) volume of water in water body at beginning of simulation\n
!> (out) dfw expressed as depth over drainage area
!> @param[inout] vy
!> (in) sediment in water body at beginning of simulation\n
!> (out) dfy expressed as loading per unit area for drainage area
!> @param[in] yi sediment entering water body
!> @param[in] yo sediment leaving water body
!> @param[in] ysx change in sediment level in water body
!> @param[in] vf volume of water in water body at end of simulation
!> @param[in] vyf sediment in water body at end of simulation
!> @param[in] aha area draining into water body
subroutine vbl(evx, spx, pp, qin, ox, vx1, vy, yi, yo, ysx, vf, vyf, aha)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    evx         |m^3 H2O       |evaporation from water body
!!    spx         |m^3 H2O       |seepage from water body
!!    pp          |m^3 H2O       |precipitation on water body
!!    qin         |m^3 H2O       |water entering water body
!!    ox          |m^3 H2O       |water leaving water body
!!    vx1         |m^3 H2O       |volume of water in water body at
!!                               |beginning of simulation
!!    vy          |metric tons   |sediment in water body at beginning
!!                               |of simulation
!!    yi          |metric tons   |sediment entering water body
!!    yo          |metric tons   |sediment leaving water body
!!    ysx         |metric tons   |change in sediment level in water body
!!    vf          |m^3 H2O       |volume of water in water body at
!!                               |end of simulation
!!    vyf         |metric tons   |sediment in water body at end of
!!                               |simulation
!!    aha         |ha            |area draining into water body
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    vx1         |mm H2O        |dfw expressed as depth over drainage
!!                               |area
!!    vy          |metric tons/ha|dfy expressed as loading per unit area
!!                               |for drainage area
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dfw         |m^3 H2O       |difference between mass balance
!!                               |calculated from watershed averages
!!                               |and actual volume of water in water
!!                               |body at end of simulation
!!    dfy         |metric tons   |difference between mass balance
!!                               |calculated from watershed averages
!!                               |and actual sediment in water body
!!                               |at end of simulation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   implicit none

   real*8, intent (in) :: evx, spx, pp, qin, ox, yi, yo, ysx
   real*8, intent (in) :: vf, vyf, aha
   real*8, intent (inout) :: vx1, vy
   real*8 :: dfw, dfy

   !! ysx undefined for reservoirs

   dfw = vx1 - evx - spx + qin + pp - ox - vf
   dfy = vy + yi - yo - ysx - vyf
   vx1 = .1 * dfw / aha
   vy = dfy / aha

   return
end
