!> @file plantmod.f90
!> file containing the subroutine plantmod
!> @author
!> modified by Javier Burguete

!> this subroutine predicts daily potential growth of total plant
!> biomass and roots and calculates leaf area index. Incorporates
!> residue for tillage functions and decays residue on ground
!> surface. Adjusts daily dry matter based on water stress.
!> @param[in] j HRU number
subroutine plantmod(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    bio_ms(:)   |kg/ha         |land cover/crop biomass (dry weight)
!!    igro(:)     |none          |land cover status code:
!!                               |0 no land cover currently growing
!!                               |1 land cover growing
!!    phubase(:)  |heat units    |base zero total heat units (used when no
!!                               |land cover is growing)
!!    phutot(:)   |heat units    |total potential heat units for year (used
!!                               |when no crop is growing)
!!    sol_rsd(:,:)|kg/ha         |amount of organic matter in the soil layer
!!                               |classified as residue
!!    tmpav(:)    |deg C         |average temperature for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    phubase(:)  |heat units    |base zero total heat units (used when no
!!                               |land cover is growing)
!!    sol_cov(:)  |kg/ha         |amount of residue on soil surface
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max
!!    SWAT: swu, grow

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j

   !! calculate residue on soil surface for current day
   sol_cov(j) = .8 * bio_ms(j) + sol_rsd(1,j)
   sol_cov(j) = Max(sol_cov(j),0.)

   !! compute plant water use and water stress
   !! compute actual plant transpiration
   if (igro(j) == 1) then
      call swu(j)
      call grow(j)
   endif
!     write (99,9994) i, hru_ra(j), bio_ms(j), laiday(j),
!    * strsw(j), strstmp(j), strsn(j), strsp(j)
!9994  format (i4,7f10.2)
   return
end
