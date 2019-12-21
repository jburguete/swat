!> @file resetlu.f90
!> file containing the subroutine resetlu
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU/subbasin management input file
!> (.mgt). This file contains data related to management practices used in
!> the HRU/subbasin.
subroutine resetlu

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!                                 |urban.dat
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!                                    |daily
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru
!!    j
!!    mon
!!    titldum
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   character(len=80) :: titldum
   real*8 xx
   integer :: j, mon

   open (9123,file=fname(no_lup))
   read (9123, 5101) titldum
   do j = 1, mhru
      read (9123,*,end=99) xx, hru_fr(j)
   end do

!!    reset all hru_fr variables
   do j = 1, mhru
      if (hru_fr(j) <= 0.) hru_fr(j) = .0000001
      xx = hru_fr(j)
      hru_km(j) = sub_km(hru_sub(j)) * xx
      hru_ha(j) = hru_km(j) * 100.  !MJW
      hru_dafr(j) = hru_km(j) / da_km  !MJW
      do mon = 1, 12
         wupnd(mon,j) = wupnd(mon,j) * xx
         wushal(mon,j) = wushal(mon,j) * xx
         wudeep(mon,j) = wudeep(mon,j) * xx
      end do
      pnd_psa(j) = pnd_psa(j) * xx
      pnd_esa(j) = pnd_esa(j) * xx
      pnd_pvol(j) = pnd_pvol(j) * xx
      pnd_evol(j) = pnd_evol(j) * xx
      pnd_vol(j) = pnd_vol(j) * xx
      wet_nsa(j) = wet_nsa(j) * xx
      wet_mxsa(j) = wet_mxsa(j) * xx
      wet_nvol(j) = wet_nvol(j) * xx
      wet_mxvol(j) = wet_mxvol(j) * xx
      wet_vol(j) = wet_vol(j) * xx
      hru_ha(j) = hru_km(j) * 100.
!   pot_vol(j) = 10. * pot_volmm(j) * hru_ha(j)   !! mm => m^3     NUBZ
      pot_volx(j) = pot_volxmm(j)
      pot_tile(j) = pot_tilemm(j)
   end do

5101 format (a80)
99 return
end
