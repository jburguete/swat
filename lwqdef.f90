!> @file lwqdef.f90
!> file containing the subroutine lwqdef
!> @author
!> modified by Javier Burguete

!> this subroutine assigns default values for the lake water quality
!> (.lwq) when the lake water quality file does not exists
!> @param[in] ii reservoir number (none)
subroutine lwqdef(ii)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii           |none          |reservoir number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chlar(:)      |none          |chlorophyll-a production coefficient for
!!                                 |reservoir
!!    lkpst_conc(:) |mg/m**3       |pesticide concentration in lake water
!!    lkpst_koc(:)  |m**3/g        |pesticide partition coefficient between
!!                                 |water and sediment in lake water
!!    lkpst_mix(:)  |m/day         |mixing velocity (diffusion/dispersion) in
!!                                 |lake water for pesticide
!!    lkpst_rea(:)  |1/day         |pesticide reaction coefficient in lake water
!!    lkpst_rsp(:)  |m/day         |resuspension velocity in lake water for
!!                                 |pesticide sorbed to sediment
!!    lkpst_stl(:)  |m/day         |settling velocity in lake water for
!!                                 |pesticide sorbed to sediment
!!    lkpst_vol(:)  |m/day         |pesticide volatilization coefficient in lake
!!                                 |water
!!    lkspst_act(:) |m             |depth of active sediment layer in lake for
!!                                 |for pesticide
!!    lkspst_bry(:) |m/day         |pesticide burial velocity in lake bed
!!                                 |sediment
!!    lkspst_conc(:)|mg/m**3       |pesticide concentration in lake bed sediment
!!    lkspst_rea(:) |1/day         |pesticide reaction coefficient in lake bed
!!                                 |sediment
!!    seccir(:)     |none          |water clarity coefficient for reservoir
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: ii

!!    set default values for parameters
   if (chlar(ii) <= 1.e-6) chlar(ii) = 1.
   if (seccir(ii) <= 1.e-6) seccir(ii) = 1.
   if (lkpst_conc(ii) <= 1.e-6) lkpst_conc(ii) = 0.
   if (lkpst_rea(ii) <= 1.e-6) lkpst_rea(ii) = 0.007
   if (lkpst_vol(ii) <= 1.e-6) lkpst_vol(ii) = 0.01
   if (lkpst_koc(ii) <= 1.e-6) lkpst_koc(ii) = 0.
   if (lkpst_stl(ii) <= 1.e-6) lkpst_stl(ii) = 1.
   if (lkpst_rsp(ii) <= 1.e-6) lkpst_rsp(ii) = 0.002
   if (lkpst_mix(ii) <= 1.e-6) lkpst_mix(ii) = 0.001
   if (lkspst_conc(ii) <= 1.e-6) lkspst_conc(ii) = 0.
   if (lkspst_rea(ii) <= 1.e-6) lkspst_rea(ii) = 0.05
   if (lkspst_bry(ii) <= 1.e-6) lkspst_bry(ii) = 0.002
   if (lkspst_act(ii) <= 1.e-6) lkspst_act(ii) = 0.030
   if (theta_n(ii) <= 0.) theta_n(ii) = 1.08
   if (theta_p(ii) <= 0.) theta_p(ii) = 1.08
   if (con_nirr(ii) <= 0.) con_nirr(ii) = 0.0
   if (con_pirr(ii) <= 0.) con_pirr(ii) = 0.0

   return
end
