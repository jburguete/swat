!> @file noqual.f90
!> file containing the subroutine noqual
!> @author
!> modified by Javier Burguete

!> this subroutine performs in-stream nutrient calculations. No transformations
!> are calculated. New concentrations of the nutrients are calculated based
!> on the loading to the reach from upstream.
!> @param[in] jrch reach number (none)
!> @param[in] k inflow hydrograph storage location number (none)
subroutine noqual(jrch, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    jrch         |none          |reach number
!!    k            |none          |inflow hydrograph storage location number
!!    ai0          |ug chla/mg alg|ratio of chlorophyll-a to algal biomass
!!    algae(:)     |mg alg/L      |algal biomass concentration in reach
!!    ammonian(:)  |mg N/L        |ammonia concentration in reach
!!    disolvp(:)   |mg P/L        |dissolved phosphorus concentration in reach
!!    nitraten(:)  |mg N/L        |nitrate concentration in reach
!!    nitriten(:)  |mg N/L        |nitrite concentration in reach
!!    organicn(:)  |mg N/L        |organic nitrogen concentration in reach
!!    organicp(:)  |mg P/L        |organic phosphorus concentration in reach
!!    rch_cbod(:)  |mg O2/L       |carbonaceous biochemical oxygen demand in
!!                                |reach
!!    rch_dox(:)   |mg O2/L       |dissolved oxygen concentration in reach
!!    rchwtr       |m^3 H2O       |water stored in reach at beginning of day
!!    rnum1        |none          |fraction of overland flow
!!    rtwtr        |m^3 H2O       |flow out of reach
!!    varoute(2,:) |m^3 H2O       |water
!!    varoute(4,:) |kg N          |organic nitrogen
!!    varoute(5,:) |kg P          |organic posphorus
!!    varoute(6,:) |kg N          |nitrate
!!    varoute(7,:) |kg P          |soluble phosphorus
!!    varoute(13,:)|kg            |chlorophyll-a
!!    varoute(14,:)|kg N          |ammonium
!!    varoute(15,:)|kg N          |nitrite
!!    varoute(16,:)|kg            |carbonaceous biological oxygen demand
!!    varoute(17,:)|kg O2         |dissolved oxygen
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    algae(:)    |mg alg/L      |algal biomass concentration in reach
!!    ammonian(:) |mg N/L        |ammonia concentration in reach
!!    chlora(:)   |mg chl-a/L    |chlorophyll-a concentration in reach
!!    disolvp(:)  |mg P/L        |dissolved phosphorus concentration in reach
!!    nitraten(:) |mg N/L        |nitrate concentration in reach
!!    nitriten(:) |mg N/L        |nitrite concentration in reach
!!    organicn(:) |mg N/L        |organic nitrogen concentration in reach
!!    organicp(:) |mg P/L        |organic phosphorus concentration in reach
!!    rch_cbod(:) |mg O2/L       |carbonaceous biochemical oxygen demand in
!!                               |reach
!!    rch_dox(:)  |mg O2/L       |dissolved oxygen concentration in reach
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    algcon      |mg alg/L      |initial algal biomass concentration in reach
!!    algin       |mg alg/L      |algal biomass concentration in inflow
!!    ammoin      |mg N/L        |ammonium N concentration in inflow
!!    cbodcon     |mg/L          |initial carbonaceous biological oxygen demand
!!                               |concentration in reach
!!    cbodin      |mg/L          |carbonaceous biological oxygen demand
!!                               |concentration in inflow
!!    chlin       |mg chl-a/L    |chlorophyll-a concentration in inflow
!!    disoxin     |mg O2/L       |dissolved oxygen concentration in inflow
!!    dispin      |mg P/L        |soluble P concentration in inflow
!!    nh3con      |mg N/L        |initial ammonia concentration in reach
!!    nitratin    |mg N/L        |nitrate concentration in inflow
!!    nitritin    |mg N/L        |nitrite concentration in inflow
!!    no2con      |mg N/L        |initial nitrite concentration in reach
!!    no3con      |mg N/L        |initial nitrate concentration in reach
!!    o2con       |mg O2/L       |initial dissolved oxygen concentration in
!!                               |reach
!!    orgncon     |mg N/L        |initial organic N concentration in reach
!!    orgnin      |mg N/L        |organic N concentration in inflow
!!    orgpcon     |mg P/L        |initial organic P concentration in reach
!!    orgpin      |mg P/L        |organic P concentration in inflow
!!    solpcon     |mg P/L        |initial soluble P concentration in reach
!!    wtrin       |m^3 H2O       |water flowing into reach on day
!!    wtrtot      |m^3 H2O       |inflow + storage water
!!    xx          |              |auxiliar variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: jrch, k
   real*8 :: algcon, algin, ammoin, cbodcon, cbodin, chlin, disoxin, dispin,&
      &nh3con, nitratin, nitritin, no2con, no3con, o2con, orgncon, orgnin,&
      &orgpcon, orgpin, solpcon, wtrin, wtrtot, xx

   !! initialize water flowing into reach
   wtrin = varoute(2,k) * (1. - rnum1)

   if (rtwtr / 86400. > 0.01 .and. wtrin > 0.01) then
!! concentrations
      !! initialize inflow concentrations
      if (varoute(13,k) < 1.e-6) varoute(13,k) = 0.0
      xx = 1000. * (1. - rnum1) / wtrin
      chlin = varoute(13,k) * xx
      algin = 1000. * chlin / ai0        !! QUAL2E equation III-1
      orgnin = varoute(4,k) * xx
      ammoin = varoute(14,k) * xx
      nitritin = varoute(15,k) * xx
      nitratin = varoute(6,k) * xx
      orgpin = varoute(5,k) * xx
      dispin = varoute(7,k) * xx
      if (varoute(16,k) < 1.e-6) varoute(16,k) = 0.0
      cbodin = varoute(16,k) * xx
      if (varoute(17,k) < 1.e-6) varoute(17,k) = 0.0
      disoxin = varoute(17,k) * xx

      !! initialize concentration of nutrient in reach
      if (algae(jrch) < 1.e-6) algae(jrch) = 0.0
      if (organicn(jrch) < 1.e-6) organicn(jrch) = 0.0
      if (ammonian(jrch) < 1.e-6) ammonian(jrch) = 0.0
      if (nitriten(jrch) < 1.e-6) nitriten(jrch) = 0.0
      if (nitraten(jrch) < 1.e-6) nitraten(jrch) = 0.0
      if (organicp(jrch) < 1.e-6) organicp(jrch) = 0.0
      if (disolvp(jrch) < 1.e-6) disolvp(jrch) = 0.0
      if (rch_cbod(jrch) < 1.e-6) rch_cbod(jrch) = 0.0
      if (rch_dox(jrch) < 1.e-6) rch_dox(jrch) = 0.0
      wtrtot = wtrin + rchwtr
      algcon = (algin * wtrin + algae(jrch) * rchwtr) / wtrtot
      orgncon = (orgnin * wtrin + organicn(jrch) * rchwtr) / wtrtot
      nh3con = (ammoin * wtrin + ammonian(jrch) * rchwtr) / wtrtot
      no2con = (nitritin * wtrin + nitriten(jrch) * rchwtr) / wtrtot
      no3con = (nitratin * wtrin + nitraten(jrch) * rchwtr) / wtrtot
      orgpcon = (orgpin * wtrin + organicp(jrch) * rchwtr) / wtrtot
      solpcon = (dispin * wtrin + disolvp(jrch) * rchwtr) / wtrtot
      cbodcon = (cbodin * wtrin + rch_cbod(jrch) * rchwtr) / wtrtot
      o2con = (disoxin * wtrin + rch_dox(jrch) * rchwtr) / wtrtot

      !! calculate algal biomass concentration at end of day
      algae(jrch) = algcon
      if (algae(jrch) < 1.e-6) algae(jrch) = 0.

      !! calculate chlorophyll-a concentration at end of day
      chlora(jrch) = algae(jrch) * ai0 / 1000.
      if (chlora(jrch) < 1.e-6) chlora(jrch) = 0.

!! oxygen calculations
      !! calculate carbonaceous biological oxygen demand at end
      !! of day
      rch_cbod(jrch) = cbodcon
      if (rch_cbod(jrch) < 1.e-6) rch_cbod(jrch) = 0.

      !! calculate dissolved oxygen concentration if reach at
      !! end of day
      rch_dox(jrch) = o2con
      if (rch_dox(jrch) < 1.e-6) rch_dox(jrch) = 0.
!! end oxygen calculations

!! nitrogen calculations
      !! calculate organic N concentration at end of day
      organicn(jrch) = orgncon
      if (organicn(jrch) < 1.e-6) organicn(jrch) = 0.

      !! calculate ammonia nitrogen concentration at end of day
      ammonian(jrch) = nh3con
      if (ammonian(jrch) < 1.e-6) ammonian(jrch) = 0.

      !! calculate concentration of nitrite at end of day
      nitriten(jrch) = no2con
      if (nitriten(jrch) < 1.e-6) nitriten(jrch) = 0.

      !! calculate nitrate concentration at end of day
      nitraten(jrch) = no3con
      if (nitraten(jrch) < 1.e-6) nitraten(jrch) = 0.
!! end nitrogen calculations

!! phosphorus calculations
      !! calculate organic phosphorus concentration at end of
      !! day
      organicp(jrch) = orgpcon
      if (organicp(jrch) < 1.e-6) organicp(jrch) = 0.

      !! calculate dissolved phosphorus concentration at end
      !! of day (mineral P)
      disolvp(jrch) = solpcon
      if (disolvp(jrch) < 1.e-6) disolvp(jrch) = 0.
!! end phosphorus calculations

   else
      !! all water quality variables set to zero when no flow
      algae(jrch) = 0.0
      chlora(jrch) = 0.0
      organicn(jrch) = 0.0
      ammonian(jrch) = 0.0
      nitriten(jrch) = 0.0
      nitraten(jrch) = 0.0
      organicp(jrch) = 0.0
      disolvp(jrch) = 0.0
      rch_cbod(jrch) = 0.0
      rch_dox(jrch) = 0.0
   endif

   return
end
