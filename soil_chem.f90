!> @file soil_chem.f90
!> file containing the subroutine soil_chem
!> @author
!> modified by Javier Burguete

!> this subroutine initializes soil chemical properties
!> @param[in] ii HRU number
subroutine soil_chem(ii)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hrupest(:)    |none          |pesticide use flag:
!!                                 | 0: no pesticides used in HRU
!!                                 | 1: pesticides used in HRU
!!    nactfr        |none          |nitrogen active pool fraction. The fraction
!!                                 |of organic nitrogen in the active pool.
!!    npmx          |none          |number of different pesticides used in
!!                                 |the simulation
!!    npno(:)       |none          |array of unique pesticides used in watershed
!!    psp           |none          |Phosphorus availability index. The fraction
!!                                 |of fertilizer P remaining in labile pool
!!                                 |after initial rapid phase of P sorption.
!!    skoc(:)       |(mg/kg)/(mg/L)|soil adsorption coefficient normalized
!!                                 |for soil organic carbon content
!!    sol_bd(:,:)   |Mg/m**3       |bulk density of the soil
!!    sol_cbn(:,:)  |%             |percent organic carbon in soil layer
!!    sol_nly(:)    |none          |number of soil layers
!!    sol_no3(:,:)  |mg N/kg soil  |nitrate concentration in soil layer
!!    sol_orgn(:,:) |mg/kg         |organic N concentration in soil layer
!!    sol_orgp(:,:) |mg/kg         |organic P concentration in soil layer
!!    sol_pst(:,:,1)|kg/ha         |initial amount of pesticide in first layer
!!                                 |read in from .chm file
!!    sol_rsd(:,:)  |kg/ha         |amount of organic matter in the soil layer
!!                                 |classified as residue
!!    sol_solp(:,:) |mg/kg         |solution P concentration in soil layer
!!    sol_z(:,:)    |mm            |depth to bottom of soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    basminpi      |kg P/ha       |average amount of phosphorus initially in
!!                                 |the mineral P pool in watershed soil
!!    basno3i       |kg N/ha       |average amount of nitrogen initially in the
!!                                 |nitrate pool in watershed soil
!!    basorgni      |kg N/ha       |average amount of nitrogen initially in
!!                                 |the organic N pool in watershed soil
!!    basorgpi      |kg P/ha       |average amount of phosphorus initially in
!!                                 |the organic P pool in watershed soil
!!    conv_wt(:,:)  |none          |factor which converts kg/kg soil to kg/ha
!!    sol_actp(:,:) |kg P/ha       |amount of phosphorus stored in the
!!                                 |active mineral phosphorus pool
!!    sol_aorgn(:,:)|kg N/ha       |amount of nitrogen stored in the active
!!                                 |organic (humic) nitrogen pool
!!    sol_cov(:)    |kg/ha         |amount of residue on soil surface
!!    sol_fon(:,:)  |kg N/ha       |amount of nitrogen stored in the fresh
!!                                 |organic (residue) pool
!!    sol_fop(:,:)  |kg P/ha       |amount of phosphorus stored in the fresh
!!                                 |organic (residue) pool
!!    sol_hum(:,:)  |kg humus/ha   |amount of organic matter in the soil layer
!!                                 |classified as humic substances
!!    sol_kp(:,:,:) |(mg/kg)/(mg/L)|pesticide sorption coefficient, Kp; the
!!                                 |ratio of the concentration in the solid
!!                                 |phase to the concentration in solution
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool. This variable is read in as
!!                                 |a concentration and converted to kg/ha.
!!                                 |(this value is read from the .sol file in
!!                                 |units of mg/kg)
!!    sol_orgn(:,:) |kg N/ha       |amount of nitrogen stored in the stable
!!                                 |organic N pool NOTE UNIT CHANGE!
!!    sol_orgp(:,:) |kg P/ha       |amount of phosphorus stored in the organic
!!                                 |P pool NOTE UNIT CHANGE!
!!    sol_pst(:,:,:)|kg/ha         |amount of pesticide in layer NOTE UNIT
!!                                 |CHANGE!
!!    sol_solp(:,:) |kg P/ha       |amount of phosohorus stored in solution
!!                                 |NOTE UNIT CHANGE!
!!    sol_stap(:,:) |kg P/ha       |amount of phosphorus in the soil layer
!!                                 |stored in the stable mineral phosphorus
!!                                 |pool
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    actp
!!    dg          |mm            |depth of layer
!!    j           |none          |counter
!!    jj          |none          |dummy variable to hold value
!!    n           |none          |counter
!!    nly         |none          |number of soil layers
!!    RTO
!!    sol_cmass
!!    sol_mass
!!    sol_min_n
!!    sol_thick
!!    soldepth    |mm            |depth from bottom of 1st soil layer to
!!                               |the bottom of the layer of interest
!!    solp
!!    solpst      |mg/kg         |concentration of pesticide in soil
!!    SSP
!!    summinp     |kg P/ha       |amount of phosphorus stored in the mineral P
!!                               |pool in the profile
!!    sumno3      |kg N/ha       |amount of nitrogen stored in the nitrate pool
!!                               |in the soil profile
!!    sumorgn     |kg N/ha       |amount of nitrogen stored in the organic N
!!                               |pools in the profile
!!    sumorgp     |kg P/ha       |amount of phosphorus stored in the organic P
!!                               |pools in the profile
!!    wt1         |none          |converts mg/kg (ppm) to kg/ha
!!    X1
!!    xx          |none          |variable to hold value
!!    zdst        |none          |variable to hold value
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   !Frction of Mirobial Biomass, Humus Passive C pools
         !FBM = 0.0
         !IF(FBM<1.E-10)FBM=.04
         !FHP = 0.0
         !IF(FHP<1.E-10)FHP=.7-.4*EXP(-.0277*100)
         !From DSSAT
         !FBM = 0.02
         !FHP = 0.44
   integer, intent(in) :: ii
   real*8, parameter :: FBM =.04, FHP = .7 - .4 * EXP(-.0277 * 100)
   real*8 :: actp, dg, RTO, sol_cmass, sol_mass, sol_min_n,&
      &sol_thick, soldepth, solp, solpst, SSP, summinp, sumno3, sumorgn,&
      &sumorgp, wt1, xx, X1, zdst
   integer :: j, jj, n, nly

   sumno3 = 0.
   sumorgn = 0.
   summinp = 0.
   sumorgp = 0.
   nly = sol_nly(ii)

!!    calculate sol_cbn for lower layers if only have upper layer
   if (nly >= 3 .and. sol_cbn(3,ii) <= 0) then
      do j = 3, nly
         if (sol_cbn(j,ii) == 0.) then
            soldepth = sol_z(j,ii) - sol_z(2,ii)
            sol_cbn(j,ii) = sol_cbn(j-1,ii) * Exp(-.001 * soldepth)
         end if
      end do
   end if

   cmup_kgh = 0.
   cmtot_kgh = 0.
   do j = 1, nly
      sol_thick = sol_z(j,ii)
      if (j /= 1) sol_thick = sol_thick - sol_z(j-1,ii)

!! soil carbon and nitrogen
      !sol_mass = (sol_thick / 1000.) * 10000. * sol_bd(j,ii)&
      !   &* 1000. * (1 - sol_rock(j,ii) / 100.) ! removed redundant operations
      sol_mass = sol_thick * 10000. * sol_bd(j,ii) * (1 - sol_rock(j,ii) / 100.)
      sol_cmass = sol_mass * (sol_cbn(j,ii) / 100.)

      if (j == 1) cmup_kgh(ii) = sol_cmass
      cmtot_kgh(ii) = cmtot_kgh(ii) + sol_cmass
   end do


!!    calculate sol_kp as function of koc and sol_cbn
!!    and set initial pesticide in all layers equal to value given for
!!    upper layer
   if (hrupest(ii) == 1) then
      do j = 1, npmx
         jj = npno(j)
         if (jj > 0) then
            solpst = sol_pst(j,ii,1)  !!concentration of pesticide in soil
            xx = 0.
            do n = 1, nly
               dg = (sol_z(n,ii) - xx)
               xx = sol_z(n,ii)
               wt1 = sol_bd(n,ii) * dg / 100.              !! mg/kg => kg/ha
               sol_kp(j,ii,n) = skoc(jj) * sol_cbn(n,ii) / 100.
               sol_pst(j,ii,n) = solpst * wt1
            end do
         end if
      end do
   end if


!!    calculate initial nutrient contents of layers, profile and
!!    average in soil for the entire watershed
!!    convert mg/kg (ppm) to kg/ha
   xx = 0.
   sol_fop(1,ii) = sol_rsd(1,ii) * .0010 !! was 0.0003 Armen January 2009
   sol_fon(1,ii) = sol_rsd(1,ii) * .0055 !! was 0.0015 Armen January 2009
   sol_cov(ii) = sol_rsd(1,ii)
   do j = 1, nly
      dg = (sol_z(j,ii) - xx)
      wt1 = sol_bd(j,ii) * dg / 100.              !! mg/kg => kg/ha
      conv_wt(j,ii) = 1.e6 * wt1                  !! kg/kg => kg/ha

      if (sol_no3(j,ii) <= 0.) then
         zdst = Exp(-sol_z(j,ii) / 1000.)
         sol_no3(j,ii) = 10. * zdst * .7
      end if
      sol_no3(j,ii) = sol_no3(j,ii) * wt1          !! mg/kg => kg/ha
      sumno3 = sumno3 + sol_no3(j,ii)

      if (sol_orgn(j,ii) > 0.0001) then
         sol_orgn(j,ii) = sol_orgn(j,ii) * wt1      !! mg/kg => kg/ha
      else
         !! assume C:N ratio of 10:1
         sol_orgn(j,ii) = 10000. * (sol_cbn(j,ii) / 14.) * wt1  !! CN ratio changed back to 14 cibin 03022012
      end if
      sol_aorgn(j,ii) = sol_orgn(j,ii) * nactfr
      sol_orgn(j,ii) = sol_orgn(j,ii) * (1. - nactfr)
      sumorgn = sumorgn + sol_aorgn(j,ii) + sol_orgn(j,ii) + sol_fon(j,ii)

      if (sol_orgp(j,ii) > 0.0001) then
         sol_orgp(j,ii) = sol_orgp(j,ii) * wt1      !! mg/kg => kg/ha
      else
         !! assume N:P ratio of 8:1
         sol_orgp(j,ii) = .125 * sol_orgn(j,ii)
      end if

      if (sol_solp(j,ii) > 0.0001) then
         sol_solp(j,ii) = sol_solp(j,ii) * wt1      !! mg/kg => kg/ha
      else
         !! assume initial concentration of 5 mg/kg
         sol_solp(j,ii) = 5. * wt1
      end if

      !! Set active pool based on dynamic PSP MJW

      if (sol_P_model == 0) then
         !! Allow Dynamic PSP Ratio
         !! convert to concentration
         solp = sol_solp(j,ii) / conv_wt(j,ii) * 1000000.
         !! PSP = -0.045*log (% clay) + 0.001*(Solution P, mg kg-1) - 0.035*(% Organic C) + 0.43
         if (sol_clay(j,ii) > 0.) then
            psp(ii) = -0.045 * log(sol_clay(j,ii))+ (0.001 * solp)
            psp(ii) = psp(ii) - (0.035  * sol_cbn(j,ii)) + 0.43
         else
            psp(ii) = 0.4
         endif
         !! Limit PSP range
         if (psp(ii) <.05) then
            psp(ii) = 0.05
         else if (psp(ii) > 0.9) then
            psp(ii) = 0.9
         end if
      end if

      sol_actp(j,ii) = sol_solp(j,ii) * (1. - psp(ii)) / psp(ii)

      !! Set Stable pool based on dynamic coefficant
      if (sol_P_model == 0) then  !! From White et al 2009
         !! convert to concentration for ssp calculation
         actp = sol_actp(j,ii) / conv_wt(j,ii) * 1000000.
         solp = sol_solp(j,ii) / conv_wt(j,ii) * 1000000.
         !! estimate Total Mineral P in this soil based on data from sharpley 2004
         SSP = 25.044 * (actp + solp)** (-0.3833)
         !!limit SSP Range
         if (SSP > 7.) SSP = 7.
         if (SSP < 1.) SSP = 1.
         sol_stap(j,ii) = SSP * (sol_actp(j,ii) + sol_solp(j,ii))!define stableP
      else
         !! The original code
         sol_stap(j,ii) = 4. * sol_actp(j,ii)
      end if

      sol_hum(j,ii) = sol_cbn(j,ii) * wt1 * 17200.
      xx = sol_z(j,ii)
      summinp = summinp + sol_solp(j,ii) + sol_actp(j,ii) + sol_stap(j,ii)
      sumorgp = sumorgp + sol_orgp(j,ii) + sol_fop(j,ii)
   end do

   xx = hru_km(ii) / da_km
   basno3i = basno3i + sumno3 * xx
   basorgni = basorgni + sumorgn * xx
   basminpi = basminpi + summinp * xx
   basorgpi = basorgpi + sumorgp * xx

   if (cswat == 2) then
      if (rsdin(ii) > 0.) sol_rsd(1,ii) = rsdin(ii)
      do j = 1, nly
         !!kg/ha sol mass in each layer
         xx = 10000. * sol_bd(j,ii) * (1. - sol_rock(j,ii) / 100.)
         if (j == 1) then
            !sol_mass = (sol_z(j,ii)) / 1000.
            !sol_mass = sol_mass * 10000. * sol_bd(j,ii)* 1000.
            !sol_mass = sol_mass * (1- sol_rock(j,ii) / 100.)
            !removed redundant operations
            sol_mass = sol_z(j,ii) * xx
         else
            !sol_mass = (sol_z(j,ii) - sol_z(j-1,ii)) / 1000.
            !sol_mass = sol_mass * 10000. * sol_bd(j,ii)* 1000.
            !sol_mass = sol_mass * (1- sol_rock(j,ii) / 100.)
            !removed redundant operations
            sol_mass = (sol_z(j,ii) - sol_z(j-1,ii)) * xx
         end if
         !!kg/ha mineral nitrogen
         sol_min_n = sol_no3(j,ii)+sol_nh3(j,ii)

         dg = sol_z(j,ii)
         if (j > 1) dg = dg - sol_z(j-1,ii)

         sol_WOC(j,ii) = sol_mass * sol_cbn(j,ii)/100
         sol_WON(j,ii) = sol_aorgn(j,ii)+  sol_orgn(j,ii)!0.1 * sol_WOC(j,ii)


         sol_BM(j,ii)=FBM*sol_WOC(j,ii)
         sol_BMC(j,ii)=sol_BM(j,ii)
         RTO=sol_WON(j,ii)/sol_WOC(j,ii)
         sol_BMN(j,ii)=RTO*sol_BMC(j,ii)
         sol_HP(j,ii)=FHP*(sol_WOC(j,ii)-sol_BM(j,ii))
         sol_HS(j,ii)=sol_WOC(j,ii)-sol_BM(j,ii)-sol_HP(j,ii)
         sol_HSC(j,ii)=sol_HS(j,ii)
         sol_HSN(j,ii)= RTO*sol_HSC(j,ii)  !sol_aorgn(j,ii)
         sol_HPC(j,ii)=sol_HP(j,ii)
         sol_HPN(j,ii)= RTO*sol_HPC(j,ii)  !sol_orgn(j,ii)


         X1=sol_rsd(j,ii) /1000.

         sol_LM(j,ii)=500.*X1
         sol_LS(j,ii)=sol_LM(j,ii)
         sol_LSL(j,ii)=.8*sol_LS(j,ii)
         sol_LMC(j,ii)=.42*sol_LM(j,ii)

         sol_LMN(j,ii)=.1*sol_LMC(j,ii)
         sol_LSC(j,ii)=.42*sol_LS(j,ii)
         sol_LSLC(j,ii)=.8*sol_LSC(j,ii)
         sol_LSLNC(j,ii)=.2*sol_LSC(j,ii)
         sol_LSN(j,ii)=sol_LSC(j,ii)/150.
         sol_WOC(j,ii)=sol_WOC(j,ii)+sol_LSC(j,ii)+sol_LMC(j,ii)
         sol_WON(j,ii)=sol_WON(j,ii)+sol_LSN(j,ii)+sol_LMN(j,ii)

         sol_orgn(j,ii) = sol_HPN(j,ii)
         sol_aorgn(j,ii) = sol_HSN(j,ii)
         sol_fon(1,ii) = sol_LMN(j,ii) + sol_LSN(j,ii)
         sumorgn = sumorgn + sol_aorgn(j,ii) + sol_orgn(j,ii) +&
            &sol_fon(j,ii) + sol_BMN(j,ii)


      end do

   end if

   return
end
