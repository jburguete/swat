!> @file latsed.f90
!> file containing the subroutine latsed
!> @author
!> modified by Javier Burguete

!> this subroutine calculates the sediment load contributed in lateral flow
!> @param[in] j HRU number (none)
subroutine latsed(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    hru_km(:)   |km^2          |area of HRU in square kilometers
!!    lat_sed(:)  |g/L           |sediment concentration in lateral flow
!!    latq(:)     |mm H2O        |total lateral flow in soil profile for the
!!                               |day in HRU
!!    sedyld(:)   |metric tons   |daily soil loss caused by water erosion in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sedyld(:)   |metric tons   |daily soil loss caused by water erosion in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    qq
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: qq, xx

   !! update sediment yield for sediment in lateral flow
   qq = latq(j) + gw_q(j)
   xx = qq * hru_km(j) * lat_sed(j)
   sedyld(j) = sedyld(j) + xx
   sanyld(j) = sanyld(j) + xx * det_san(j)
   silyld(j) = silyld(j) + xx * det_sil(j)
   clayld(j) = clayld(j) + xx * det_cla(j)
   sagyld(j) = sagyld(j) + xx * det_sag(j)
   lagyld(j) = lagyld(j) + xx * det_lag(j)

   !! organic n and p in the lateral flow     - by J.Jeong BREC 2011 revised 2014
   !1mm*mg/L*1000L/m3*kg/1000000mg*10m3/(ha-mm)=0.01kg/ha
   sedorgn(j) = sedorgn(j) + qq * lat_orgn(j) / 100.
   sedorgp(j) = sedorgp(j) + qq * lat_orgp(j) / 100.

   !! bmp adjustments
   sedyld(j) = sedyld(j) * bmp_seds(j)
   sedorgp(j) = sedorgp(j) * bmp_pps(j)
   sedorgn(j) = sedorgn(j) * bmp_pns(j)

   if (sedyld(j) < 0.) sedyld(j) = 0.
   if (sanyld(j) < 0.) sanyld(j) = 0.0
   if (silyld(j) < 0.) silyld(j) = 0.0
   if (clayld(j) < 0.) clayld(j) = 0.0
   if (sagyld(j) < 0.) sagyld(j) = 0.0
   if (lagyld(j) < 0.) lagyld(j) = 0.0
   if (sedorgn(j) < 0.) sedorgn(j) = 0.0
   if (sedorgp(j) < 0.) sedorgp(j) = 0.0

   return
end
