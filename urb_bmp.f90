!> @file urb_bmp.f90
!> file containing the subroutine urb_bmp
!> @author
!> modified by Javier Burguete

!> this subroutine
!> @param[in] j HRU number (none)
subroutine urb_bmp(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j              |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name       |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sednppm
!!    sedppm
!!    sedpppm
!!    solnppm
!!    solpppm
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: sednppm, sedppm, sedpppm, solnppm, solpppm, xx

!! convert to ppm -> (kg/ha)*100./mm = ppm
   if (qdr(j) > 0.1) then
      xx = 100. / qdr(j)
      sedppm = 1000. * xx * sedyld(j) / hru_ha(j)
      solnppm = xx * (surqno3(j) + latno3(j) + no3gw(j))
      solpppm = xx * (surqsolp(j) + minpgw(j))
      sednppm = xx * sedorgn(j)
      sedpppm = xx * (sedorgp(j) + sedminpa(j) + sedminps(j))

      if (sedppm > sed_con (j)) then
         sedyld(j) = sed_con(j) * hru_ha(j) / xx / 1000.
      endif

      if (solnppm > soln_con(j)) then
         surqno3(j) = soln_con(j) / xx
         latno3(j) = soln_con(j) / xx
         no3gw(j) = soln_con(j) / xx
      endif

      if (solpppm > solp_con(j)) then
         surqsolp(j) = solp_con(j) / xx
         minpgw(j) = solp_con(j) / xx
      endif

      if (sednppm > orgn_con(j)) then
         sedorgn(j) = orgn_con(j) / xx
      endif

      if (sedpppm > orgp_con(j)) then
         sedorgn(j)= orgp_con(j) / xx
         sedminpa(j)= orgp_con(j) / xx
         sedminps(j)= orgp_con(j) / xx
      endif

   endif

   return
end
