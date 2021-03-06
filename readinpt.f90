!> @file readinpt.f90
!> file containing the subroutine readinpt
!> @author
!> modified by Javier Burguete

!> this subroutine calls subroutines which read input data for the
!> databases and the HRUs
subroutine readinpt

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    nhru          |none          |total number of HRUs in the watershed
!!    npmx          |none          |total number of pesticides modeled in
!!                                 |in watershed plus 1
!!    nope(:)       |none          |sequence number of pesticide in NPNO(:)
!!    pstflg(:)     |none          |array defining different pesticides
!!                                 |used in the watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i             |none          |number of specific reservoir or HRU
!!    irtpest       |none          |redefined to the sequence number of
!!                                 |pestcide in NPNO(:) which is to be routed
!!                                 |through the watershed
!!    ndays(:)      |julian date   |julian date for last day of preceding
!!                                 |month (where the array location is the
!!                                 |number of the month) The dates are for
!!                                 |leap years
!!    npmx          |none          |total number of pesticides modeled in
!!                                 |watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: soil_chem, soil_phys, rteinit, h2omgt_init, hydroinit,
!!          impnd_init

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer :: ii

   sol_WOC = 0.
   sol_WON = 0.

   if (irtpest > 0) irtpest = nope(irtpest)
   npmx = Sum(pstflg)         !! set equal to # pesticides modeled in
   !! watershed

   do ii = 1, nhru
      call soil_chem(ii)       !! initialize soil chemical parameters
      call soil_phys(ii)       !! initialize soil physical parameters
   end do

   call rteinit               !! initialize routing parameters
   call h2omgt_init           !! initialize water management parameters
   call hydroinit             !! initialize hydrology parameters
   call impnd_init            !! initialize impoundment parameters

   return
end
