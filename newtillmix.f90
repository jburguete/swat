!> @file newtillmix.f90
!> file containing the subroutine newtillmix
!> @author
!> Armen R. Kemanian,\n
!> Stefan Julich,\n
!> Cole Rossi\n
!> modified by Javier Burguete

!> this subroutine mixes residue and nutrients during tillage and
!> biological mixing.
!> Mixing was extended to all layers.
!> A subroutine to simulate stimulation of organic matter decomposition was
!> added.
!> March 2009: testing has been minimal and further adjustments are expected.
!> Use with caution!
!> @param[in] j HRU number (none)
!> @param[in] bmix biological mixing efficiency: this number is zero for tillage
!> operations (none)
subroutine newtillmix(j,bmix)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j             |none          |HRU number
!!    bmix          |none          |biological mixing efficiency: this
!!                                 |number is zero for tillage operations
!!    bactlpq(:)    |# colonies/ha |less persistent bacteria in soil solution
!!    bactlps(:)    |# colonies/ha |less persistent bacteria attached to soil
!!                                 |particles
!!    bactpq(:)     |# colonies/ha |persistent bacteria in soil solution
!!    bactps(:)     |# colonies/ha |persistent bacteria attached to soil
!!                                 |particles
!!    cnop          |none          |SCS runoff curve number for moisture
!!                                 |condition II
!!    curyr         |none          |current year of simulation
!!    deptil(:)     |mm            |depth of mixing caused by tillage
!!                                 |operation
!!    effmix(:)     |none          |mixing efficiency of tillage operation
!!    npmx          |none          |number of different pesticides used in
!!                                 |the simulation
!!    nro(:)        |none          |sequence number of year in rotation
!!    ntil(:)       |none          |sequence number of tillage operation within
!!                                 |current year
!!    nyskip        |none          |number of years to skip output printing/
!!                                 |summarization
!!    sol_actp(:,:) |kg P/ha       |amount of phosphorus stored in the
!!                                 |active mineral phosphorus pool
!!    sol_aorgn(:,:)|kg N/ha       |amount of nitrogen stored in the active
!!                                 |organic (humic) nitrogen pool
!!    sol_fon(:,:)  |kg N/ha       |amount of nitrogen stored in the fresh
!!                                 |organic (residue) pool
!!    sol_fop(:,:)  |kg P/ha       |amount of phosphorus stored in the fresh
!!                                 |organic (residue) pool
!!    sol_nh3(:,:)  |kg N/ha       |amount of nitrogen stored in the ammonium
!!                                 |pool in soil layer
!!    sol_nly(:)    |none          |number of soil layers
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool.
!!    sol_orgn(:,:) |kg N/ha       |amount of nitrogen stored in the stable
!!                                 |organic N pool
!!    sol_orgp(:,:) |kg P/ha       |amount of phosphorus stored in the organic
!!                                 |P pool
!!    sol_pst(:,:,:)|kg/ha         |amount of pesticide in layer
!!    sol_rsd(:,:)  |kg/ha         |amount of organic matter in the soil
!!                                 |classified as residue
!!    sol_solp(:,:) |kg P/ha       |amount of phosohorus stored in solution
!!    sol_stap(:,:) |kg P/ha       |amount of phosphorus in the soil layer
!!                                 |stored in the stable mineral phosphorus pool
!!    sol_z(:,:)    |mm            |depth to bottom of soil layer
!!    sumix(:)      |none          |sum of mixing efficiencies in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bactlpq(:)    |# colonies/ha |less persistent bacteria in soil solution
!!    bactlps(:)    |# colonies/ha |less persistent bacteria attached to soil
!!                                 |particles
!!    bactpq(:)     |# colonies/ha |persistent bacteria in soil solution
!!    bactps(:)     |# colonies/ha |persistent bacteria attached to soil
!!                                 |particles
!!    ntil(:)       |none          |sequence number of tillage operation within
!!                                 |current year
!!    sol_actp(:,:) |kg P/ha       |amount of phosphorus stored in the
!!                                 |active mineral phosphorus pool
!!    sol_aorgn(:,:)|kg N/ha       |amount of nitrogen stored in the active
!!                                 |organic (humic) nitrogen pool
!!    sol_fon(:,:)  |kg N/ha       |amount of nitrogen stored in the fresh
!!                                 |organic (residue) pool
!!    sol_fop(:,:)  |kg P/ha       |amount of phosphorus stored in the fresh
!!                                 |organic (residue) pool
!!    sol_nh3(:,:)  |kg N/ha       |amount of nitrogen stored in the ammonium
!!                                 |pool in soil layer
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool.
!!    sol_orgn(:,:) |kg N/ha       |amount of nitrogen stored in the stable
!!                                 |organic N pool
!!    sol_orgp(:,:) |kg P/ha       |amount of phosphorus stored in the organic
!!                                 |P pool
!!    sol_rsd(:,:)  |kg/ha         |amount of organic matter in the soil
!!                                 |classified as residue
!!    sol_solp(:,:) |kg P/ha       |amount of phosohorus stored in solution
!!    sol_stap(:,:) |kg P/ha       |amount of phosphorus in the soil layer
!!                                 |stored in the stable mineral phosphorus pool
!!    sumix(:)      |none          |sum of mixing efficiencies in HRU
!!    min_res(:)    |kg/ha         |Min residue allowed due to implementation of
!!                                 |residue managment in the OPS file.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dtil        |mm            |depth of mixing
!!    emix        |none          |mixing efficiency
!!    k           |none          |counter
!!    l           |none          |counter
!!    maxmix      |none          | maximum mixing eff to preserve specified minimum residue cover
!!    smix(:)     |varies        |amount of substance in soil profile
!!                               |that is being redistributed between
!!                               |mixed layers
!!    sol_mass
!!    sol_msm                    | sol_mass mixed
!!    sol_msn                    | sol_mass not mixed
!!    sol_thick
!!    XX
!!    WW1
!!    WW2
!!    WW3
!!    WW4
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min
!!    SWAT: tillfactor, curno

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: j
   real*8, intent (in) :: bmix
   real*8, dimension(22+npmx+12) :: smix
   real*8, dimension(sol_nly(j)) :: sol_mass, sol_msm, sol_msn, sol_thick
   real*8 :: dtil, emix, maxmix, XX, WW1, WW2, WW3, WW4
   integer :: k, l


   XX = 0.
   if (bmix > 1.e-6) then
      !! biological mixing
      emix = bmix !bmix MJW (rev 412)
      dtil = Min(sol_z(sol_nly(j),j), 50.) ! it was 300.  MJW (rev 412)
   else
      !! tillage operation
      emix = effmix(idtill)
      dtil = deptil(idtill)
   end if
   ! -------------------------------- Original D. Moriasi code replaced by code below
! Drainmod  07/2006
!      if(itill(j) == 1) then
!   cumei(j) = 0.
!   cumeira(j) = 0.
!   cumrt(j) = 0.
!        cumrai(j) = 0.
!   ranrns_hru(j) = ranrns(idtill)
!      end if
!!    Drainmod 7/2006
! --------------------------------------------------------------------
   if (idtill.GE. 1) then ! Updated dynamic depressional storage D.Moriasi 4/8/2014
      cumei(j)   = 0.
      cumeira(j) = 0.
      cumrt(j)   = 0.
      cumrai(j)  = 0.
      ranrns_hru(j) = ranrns(idtill)
   end if
! --------------------------------------------------------------------

   !!by zhang DSSAT tillage
   !!=======================
   !!    deptil(:)   |mm  |depth of mixing caused by tillage operation
   !j is hru number
   if (cswat == 2) then
      tillage_days(j) = 0
      tillage_depth(j) = dtil
      tillage_switch(j) = 1
   end if
   !!by zhang DSSAT tillage
   !!=======================


   smix = 0.
   sol_mass = 0.
   sol_thick = 0.
   sol_msm = 0.
   sol_msn = 0.

   !! incorporate bacteria - no mixing - lost from transport
   if (dtil > 10.) then
      bactpq(j) = bactpq(j) * (1. - emix)
      bactps(j) = bactps(j) * (1. - emix)
      bactlpq(j) = bactlpq(j) * (1. - emix)
      bactlps(j) = bactlps(j) * (1. - emix)
   end if

   !! calculate max mixing to preserve target surface residue MJW rev 490
   !! Assume residue in all other layers is negligible to simplify calculation and remove depth dependency
   if (min_res(j) > 1. .and. bmix < 0.001) then
      maxmix = 1 - min_res(j)/sol_rsd(1,j)
      if (maxmix <0.05)  maxmix = 0.05
      if (emix > maxmix)  emix = maxmix
   end if


   do l=1, sol_nly(j)
      if ( l == 1) then
         sol_thick(l) = sol_z(l,j)
      else
         sol_thick(l) = sol_z(l,j) - sol_z(l-1,j)
      end if

      sol_mass(l) = (sol_thick(l) / 1000.) * 10000. *&
         &sol_bd(l,j) * 1000. * (1.- sol_rock(l,j) / 100.)

   end do

   smix = 0.

   if (dtil > 0.) then
!!!  added by Armen 09/10/2010 next line only
      if (dtil < 10.0) dtil = 11.0
      do l=1, sol_nly(j)

         if (sol_z(l,j) <= dtil) then
            !! msm = mass of soil mixed for the layer
            !! msn = mass of soil not mixed for the layer
            sol_msm(l) = emix * sol_mass(l)
            sol_msn(l) = sol_mass(l) - sol_msm(l)
         else if (sol_z(l,j) > dtil.AND.sol_z(l-1,j) < dtil) then
            sol_msm(l) = emix * sol_mass(l) *&
               &(dtil - sol_z(l-1,j)) / sol_thick(l)
            sol_msn(l) =  sol_mass(l) -  sol_msm(l)
         else
            sol_msm(l) = 0.
            sol_msn(l) = sol_mass(l)
         end if

         !! calculate the mass or concentration of each mixed element
         !! mass based mixing
         WW1 = sol_msm(l)/(sol_msm(l) + sol_msn(l))
         smix(1) = smix(1) + sol_no3(l,j) * WW1
         smix(2) = smix(2) + sol_orgn(l,j) * WW1
         smix(3) = smix(3) + sol_nh3(l,j) * WW1
         smix(4) = smix(4) + sol_solp(l,j) * WW1
         smix(5) = smix(5) + sol_orgp(l,j) * WW1
         smix(6) = smix(6) + sol_aorgn(l,j) * WW1
         smix(7) = smix(7) + sol_actp(l,j) * WW1
         smix(8) = smix(8) + sol_fon(l,j) * WW1
         smix(9) = smix(9) + sol_fop(l,j) * WW1
         smix(10) = smix(10) + sol_stap(l,j) * WW1
         smix(11) = smix(11) + sol_rsd(l,j) * WW1
         smix(12) = smix(12) + sol_mc(l,j) * WW1
         smix(13) = smix(13) + sol_mn(l,j) * WW1
         smix(14) = smix(14) + sol_mp(l,j) * WW1

         !! concentration based mixing
         WW2 = XX + sol_msm(l)
         smix(15) = (XX * smix(15) + sol_cbn(l,j) * sol_msm(l)) /WW2
         smix(16) = (XX * smix(16) + sol_n(l,j) * sol_msm(l)) /WW2
         smix(17) = (XX * smix(17) + sol_clay(l,j) * sol_msm(l)) /WW2
         smix(18) = (XX * smix(18) + sol_silt(l,j) * sol_msm(l)) /WW2
         smix(19) = (XX * smix(19) + sol_sand(l,j) * sol_msm(l)) /WW2
!          smix(20) = (XX * smix(20) + sol_rock(l,j) * sol_msm(l)) / WW2
!          smix(21) = (XX * smix(21) + sol_ph(l,j) * sol_msm(l)) /WW2 !! mjw rev490
!          smix(22) = (XX * smix(22) + sol_cal(l,j) * sol_msm(l)) /WW2 !! mjw rev490
         !! mass based distribution
         do k = 1, npmx
            smix(20+k) = smix(20+k) + sol_pst(k,j,l) * WW1
         end do

         !!by zhang
         !!==============
         if (cswat == 2) then
            smix(20+npmx+1) = smix(20+npmx+1) +sol_LSC(l,j)* WW1
            smix(20+npmx+2) = smix(20+npmx+2) +sol_LSLC(l,j)* WW1
            smix(20+npmx+3) = smix(20+npmx+3) +sol_LSLNC(l,j)* WW1
            smix(20+npmx+4) = smix(20+npmx+4) +sol_LMC(l,j)* WW1
            smix(20+npmx+5) = smix(20+npmx+5) +sol_LM(l,j)* WW1
            smix(20+npmx+6) = smix(20+npmx+6) +sol_LSL(l,j)* WW1
            smix(20+npmx+7) = smix(20+npmx+7) +sol_LS(l,j)* WW1

            smix(20+npmx+8) = smix(20+npmx+8) +sol_LSN(l,j)* WW1
            smix(20+npmx+9) = smix(20+npmx+9) +sol_LMN(l,j)* WW1
            smix(20+npmx+10) = smix(20+npmx+10) +sol_BMN(l,j)* WW1
            smix(20+npmx+11) = smix(20+npmx+11) +sol_HSN(l,j)* WW1
            smix(20+npmx+12) = smix(20+npmx+12) +sol_HPN(l,j)* WW1
         end if
         !!by zhang
         !!=============

         XX = XX + sol_msm(l)
      end do

      do l=1, sol_nly(j)

         ! reconstitute each soil layer
         WW3 = sol_msn(l) / sol_mass(l)
         WW4 = sol_msm(l) / XX

         sol_no3(l,j) = sol_no3(l,j) * WW3 + smix(1) * WW4
         sol_orgn(l,j) = sol_orgn(l,j) * WW3 + smix(2) * WW4
         sol_nh3(l,j) = sol_nh3(l,j) * WW3 + smix(3) * WW4
         sol_solp(l,j) = sol_solp(l,j) * WW3 + smix(4) * WW4
         sol_orgp(l,j) = sol_orgp(l,j) * WW3 + smix(5) * WW4
         sol_aorgn(l,j) = sol_aorgn(l,j) * WW3 + smix(6) * WW4
         sol_actp(l,j) = sol_actp(l,j) * WW3 + smix(7) * WW4
         sol_fon(l,j) = sol_fon(l,j) * WW3 + smix(8) * WW4
         sol_fop(l,j) = sol_fop(l,j) * WW3 + smix(9) * WW4
         sol_stap(l,j) = sol_stap(l,j) * WW3 + smix(10) * WW4
         sol_rsd(l,j) = sol_rsd(l,j) * WW3 + smix(11) * WW4
         if (sol_rsd(l,j) < 1.e-10) sol_rsd(l,j) = 1.e-10
         sol_mc(l,j) = sol_mc(l,j) * WW3 + smix(12) * WW4
         sol_mn(l,j) = sol_mn(l,j) * WW3 + smix(13) * WW4
         sol_mp(l,j) = sol_mp(l,j) * WW3 + smix(14) * WW4

         sol_cbn(l,j) = (sol_cbn(l,j) * sol_msn(l) + smix(15)&
            &* sol_msm(l)) / sol_mass(l)
         sol_n(l,j) = (sol_n(l,j) * sol_msn(l) + smix(16)&
            &* sol_msm(l)) / sol_mass(l)
         sol_clay(l,j) = (sol_clay(l,j) * sol_msn(l) + smix(17)&
            &* sol_msm(l)) / sol_mass(l)
         sol_silt(l,j) = (sol_silt(l,j) * sol_msn(l) + smix(18)&
            &* sol_msm(l)) / sol_mass(l)
         sol_sand(l,j) = (sol_sand(l,j) * sol_msn(l) + smix(19)&
            &* sol_msm(l)) / sol_mass(l)
!  sol_rock(l,j) = (sol_rock(l,j) * sol_msn(l) + smix(20) * sol_msm(l)) / sol_mass(l)
!            sol_ph(l,j) = (sol_ph(l,j) * sol_msn(l) + smix(21)        &
!     &           * sol_msm(l)) / sol_mass(l) !! mjw rev 490 simplified, PH not linear
!            sol_cal(l,j) = (sol_cal(l,j) * sol_msn(l) + smix(22)      &
!     &           * sol_msm(l)) / sol_mass(l) !! mjw rev 490



         do k = 1, npmx
            sol_pst(k,j,l) = sol_pst(k,j,l) * WW3 + smix(20+k) * WW4
         end do

         !!by zhang
         !!=============
         if (cswat == 2) then
            sol_LSC(l,j) = sol_LSC(l,j)*WW3+smix(20+npmx+1)* WW4
            sol_LSLC(l,j) = sol_LSLC(l,j)*WW3+smix(20+npmx+2)* WW4
            sol_LSLNC(l,j) = sol_LSLNC(l,j)*WW3+smix(20+npmx+3)* WW4
            sol_LMC(l,j) = sol_LMC(l,j)*WW3 + smix(20+npmx+4)* WW4
            sol_LM(l,j) = sol_LM(l,j)*WW3 + smix(20+npmx+5)* WW4
            sol_LSL(l,j) = sol_LSL(l,j)*WW3 + smix(20+npmx+6)* WW4
            sol_LS(l,j) = sol_LS(l,j)*WW3 + smix(20+npmx+7)* WW4
            sol_LSN(l,j) = sol_LSN(l,j)*WW3 + smix(20+npmx+8)* WW4
            sol_LMN(l,j) = sol_LMN(l,j)*WW3 + smix(20+npmx+9)* WW4
            sol_BMN(l,j) = sol_BMN(l,j)*WW3 + smix(20+npmx+10)* WW4
            sol_HSN(l,j) = sol_HSN(l,j)*WW3 + smix(20+npmx+11)* WW4
            sol_HPN(l,j) = sol_HPN(l,j)*WW3 + smix(20+npmx+12)* WW4
         end if
         !!by zhang
         !!==============

      end do

      if (cswat == 1) then
         call tillfactor(j,bmix,emix,dtil,sol_thick)
      end if

      !! summary calculations
      if (curyr > nyskip) then
         sumix(j) = sumix(j) + emix
      end if

   end if

   !! perform final calculations for tillage operation

   !! count the tillage only if it is a scheduled operation biomix does not count MJW Rev 490
   if (bmix <= 1.e-6) then
      ntil(j) = ntil(j) + 1
   end if
   if (cnop > 1.e-4) call curno(cnop,j)

   !ntil(j) = ntil(j) + 1 ' orig code

   return
end
