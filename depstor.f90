!> @file depstor.f90
!> file containing the subroutine depstor
!> @author
!> modified by Javier Burguete

!> this subroutine computes maximum surface depressional storage depth based on
!> random and oriented roughness and slope steepness
!> @param[in] j HRU number
subroutine depstor(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    hru_slp(:)  |m/m           |average slope steepness in HRU
!!    iida        |julian date   |day being simulated (current julian day)
!!    precipday   |mm H2O        |precipitation for the day in HRU
!!    ranrns_hru(:)|mm           |random roughness for a given HRU
!!    sol_bd(1,:) |Mg/m^3        |bulk density of top soil layer in HRU
!!    sol_cbn(1,:)|%             |percent organic carbon in top soil layer in HRU
!!    sol_clay(1,:)|%            |percent clay content in top soil layer in HRU
!!    sol_rsd(1,:)|kg/ha         |amount of organic matter in the top soil layer
!!                               |classified as residue in HRU
!!    sol_ori(:)  |mm            |oriented roughness (ridges) at time of a given tillage operation
!!    sol_z(1,:)  |mm            |depth to bottom of top soil layer
!!    usle_ei     |100(ft-tn in)/(acre-hr)|USLE rainfall erosion index
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    stmaxd(:)   |mm            |maximum surface depressional storage for day in a given HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cumei(:)    |Mj*mm/ha*hr   |cumulative USLE rainfall erosion index since last
!!                               |tillage operation
!!    cumrt(:)    |mm H2O        |cumulative rainfall since last tillage operation
!!    df          |none          |oriented and random roughness decay factor - based
!!                               |on cumulative EI and cumulative precipday
!!    ei          |Mj*mm/ha*hr   |USLE rainfall erosion index
!!    hru_slpp    |%             |average percent slope steepness
!!    sol_orgm    |%             |percent organic matter content in soil material
!!    sol_rrr     |cm            |random roughness after a rain event
!!    sol_orr     |cm            |oriented roughness (ridges) after a rain event
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: df, hru_slpp, sol_orgm, sol_rrr, ei, xx


!! Calculate current cummulative erosivity and rainfall
   ei = usle_ei * 18.7633
   ! itill is not initialized in any part of code. Assumed itill=0
   !if(itill(j) ==1)then
   !   cumeira(j) = cumeira(j) + ei
   !   cumei(j) = cumeira(j) - ei
   !   cumrai(j) = cumrai(j) + precipday
   !   cumrt(j) = cumrai(j) - precipday
   !end if
!! Calculate the decay factor df based on %clay and %organic matter or %organic carbon
   sol_orgm = sol_cbn(1,j) / 0.58
   xx = (0.943 - 0.07 * sol_clay(1,j) + 0.0011 * sol_clay(1,j)**2&
      &- 0.67 * sol_orgm + 0.12 * sol_orgm**2)
   if (xx > 1.) then
      df = 1.
   else
      df = Exp (xx)
   end if


!! Determine the current random and oriented roughness using cumei and cumrt and initial
!! random and oriented roughness values

   ! cumei and cumrt seem to be always 0
   !sol_rrr = 0.1 * ranrns_hru(j)&
   !   &* Exp(df * (-0.0009 * cumei(j) - 0.0007 * cumrt(j)))
   sol_rrr = 0.1 * ranrns_hru(j)

!! Compute the current maximum depressional storage using percent slope steepness
!! and current random and oriented roughness values determined above
   hru_slpp = hru_slp(j)*100
   stmaxd(j)= 0.112*sol_rrr+0.031*sol_rrr**2-0.012*sol_rrr*hru_slpp

   return
end
