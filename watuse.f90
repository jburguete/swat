!> @file watuse.f90
!> file containing the subroutine watuse
!> @author
!> modified by Javier Burguete

!> this subroutine removes water from appropriate source (pond,
!> shallow aquifer, and/or deep aquifer) for consumptive water use
!> @param[in] j HRU number (none)
subroutine watuse(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    deepst(:)   |mm H2O        |depth of water in deep aquifer
!!    pnd_vol(:)  |m^3 H2O       |volume of water in pond
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    wudeep(:,:) |10^4 m^3/day  |average daily water removal from the deep
!!                               |aquifer for the month for the HRU within the
!!                               |subbasin
!!    wupnd(:,:)  |10^4 m^3/day  |average daily water removal from the pond
!!                               |for the month for the HRU within the subbasin
!!    wushal(:,:) |10^4 m^3/day  |average daily water removal from the shallow
!!                               |aquifer for the month for the HRU within the
!!                               |subbasin
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    deepst(:)   |mm H2O        |depth of water in deep aquifer
!!    pnd_vol(:)  |m^3 H2O       |volume of water in pond
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cnv         |none          |conversion factor (mm/ha => m^3)
!!    sub_ha
!!    xx          |none          |auxiliar variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: cnv, sub_ha, xx

   sub_ha = da_ha * sub_fr(hru_sub(j))
   cnv = sub_ha * 10.

   pnd_vol(j) = pnd_vol(j) - wupnd(i_mo,j) * 10000.
   if (pnd_vol(j) < 0.) pnd_vol(j) = 0.

   xx = 10000. / cnv

   if (wushal(i_mo,j) < 0.) then
      rchrg_src(j) = -1. * wushal(i_mo,j) * xx
   else
      rchrg_src(j) = 0.
      shallst(j) = shallst(j) - wushal(i_mo,j) * xx
      if (shallst(j) < 0.) shallst(j) = 0.
   end if

   deepst(j) = deepst(j) - wudeep(i_mo,j) * xx
   if (deepst(j) < 0.) deepst(j) = 0.

   return
end
