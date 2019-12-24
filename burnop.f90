!> @file burnop.f90
!> file containing the subroutine burnop
!> @author
!> modified by Javier Burguete

!> this subroutine performs burning
!> @param[in] j HRU number
subroutine burnop(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    burn_frlb   |none          |fraction of biomass and residue that burn(input in
!!                               |management file) range (0 - 1.0)
!!    nro         |none          |sequence number of year in rotation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    pburn       |              |amount of phosphorus that burns - removed from plant
!!    xx          |              |auxiliar variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: pburn, xx

   xx = burn_frlb
   bio_ms(j) = bio_ms(j) * xx
   plantn(j) = plantn(j) * xx
   pburn = plantp(j) * xx
   sol_orgp(1,j) = sol_orgp(1,j) + pburn
   plantp(j) = plantp(j) - pburn
   sol_rsd(1,j) = sol_rsd(1,j) * xx
   sol_fon(1,j) = sol_fon(1,j) * xx
   sol_aorgn(1,j) = sol_aorgn(1,j) * xx
   sol_orgn(1,j) = sol_orgn(1,j) * xx

   if (cswat == 2) then
      sol_LM(1,j) = sol_LM(1,j) * xx
      sol_LS(1,j) = sol_LS(1,j) * xx
      sol_LSC(1,j) = sol_LSC(1,j) * xx
      sol_LSN(1,j) =sol_LSN(1,j) * xx
      sol_LMC(1,j) = sol_LMC(1,j) * xx
      sol_LMN(1,j) = sol_LMN(1,j) * xx
      sol_LSL(1,j) = sol_LSL(1,j) * xx

      xx = 1. - xx
      emitc_d(j) = emitc_d(j) + bio_ms(j) * xx
      emitc_d(j) = emitc_d(j) + sol_rsd(1,j) * xx
   end if

   return
end
