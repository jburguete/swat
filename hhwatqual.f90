!> @file hhwatqual.f90
!> file containing the subroutine hhwatqual
!> @author
!> modified by Javier Burguete

!> this subroutine performs in-stream nutrient transformations and water
!> quality calculations for hourly timestep
!> @param[in] jrch reach number (none)
!> @param[in] k inflow hydrograph storage location number (none)
subroutine hhwatqual(jrch, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    jrch             |none          |reach number
!!    k                |none          |inflow hydrograph storage location number
!!    ai0              |ug chla/mg alg|ratio of chlorophyll-a to algal biomass
!!    ai1              |mg N/mg alg   |fraction of algal biomass that is N
!!    ai2              |mg P/mg alg   |fraction of algal biomass that is P
!!    ai3              |mg O2/mg alg  |the rate of oxygen production per unit of
!!                                    |algal photosynthesis
!!    ai4              |mg O2/mg alg  |the rate of oxygen uptake per unit of
!!                                    |algae respiration
!!    ai5              |mg O2/mg N    |the rate of oxygen uptake per unit of NH3
!!                                    |nitrogen oxidation
!!    ai6              |mg O2/mg N    |the rate of oxygen uptake per unit of NO2
!!                                    |nitrogen oxidation
!!    algae(:)         |mg alg/L      |algal biomass concentration in reach
!!    ammonian(:)      |mg N/L        |ammonia concentration in reach
!!    bc(1,:)           |1/hr          |rate constant for biological oxidation of
!!                                    |NH3 to NO2 in reach at 20 deg C
!!    bc(2,:)           |1/hr          |rate constant for biological oxidation of
!!                                    |NO2 to NO3 in reach at 20 deg C
!!    bc(3,:)           |1/hr          |rate constant for hydrolysis of organic N
!!                                    |to ammonia in reach at 20 deg C
!!    bc(4,:)           |1/hr          |rate constant for the decay of organic P
!!                                    |to dissolved P in reach at 20 deg C
!!    chlora(:)        |mg chl-a/L    |chlorophyll-a concentration in reach
!!    dayl(:)          |hours         |day length for current day
!!    disolvp(:)       |mg P/L        |dissolved P concentration in reach
!!    frad(:,:)        |none          |fraction of solar radiation occuring
!!                                    |during hour in day in HRU
!!    hdepth(:)        |m             |depth of flow on day
!!    hhtime(:)        |hr            |flow travel time for hour
!!    hhvaroute(2,:,:) |m^3 H2O       |water
!!    hhvaroute(4,:,:) |kg N          |organic nitrogen
!!    hhvaroute(5,:,:) |kg P          |organic posphorus
!!    hhvaroute(6,:,:) |kg N          |nitrate
!!    hhvaroute(7,:,:) |kg P          |soluble phosphorus
!!    hhvaroute(13,:,:)|kg            |chlorophyll-a
!!    hhvaroute(14,:,:)|kg N          |ammonium
!!    hhvaroute(15,:,:)|kg N          |nitrite
!!    hhvaroute(16,:,:)|kg            |carbonaceous biological oxygen demand
!!    hhvaroute(17,:,:)|kg O2         |dissolved oxygen
!!    hrchwtr(ii)      |m^3 H2O       |water stored in reach at beginning of day
!!    hrtwtr(:)        |m^3 H2O       |flow out of reach
!!    hru_ra(:)        |MJ/m^2        |solar radiation for the day in HRU
!!    igropt           |none          |Qual2E option for calculating the local
!!                                    |specific growth rate of algae
!!                                    |1: multiplicative:
!!                                    | u = mumax * fll * fnn * fpp
!!                                    |2: limiting nutrient
!!                                    | u = mumax * fll * Min(fnn, fpp)
!!                                    |3: harmonic mean
!!                                    | u = mumax * fll * 2. / ((1/fnn)+(1/fpp))
!!    k_l              |MJ/(m2*hr)    |half saturation coefficient for light
!!    k_n              |mg N/L        |michaelis-menton half-saturation constant
!!                                    |for nitrogen
!!    k_p              |mg P/L        |michaelis-menton half saturation constant
!!                                    |for phosphorus
!!    lambda0          |1/m           |non-algal portion of the light extinction
!!                                    |coefficient
!!    lambda1          |1/(m*ug chla/L)|linear algal self-shading coefficient
!!    lambda2          |(1/m)(ug chla/L)**(-2/3)
!!                                    |nonlinear algal self-shading coefficient
!!    mumax            |1/hr          |maximum specific algal growth rate at
!!                                    |20 deg C
!!    nitraten(:)      |mg N/L        |nitrate concentration in reach
!!    nitriten(:)      |mg N/L        |nitrite concentration in reach
!!    organicn(:)      |mg N/L        |organic nitrogen concentration in reach
!!    organicp(:)      |mg P/L        |organic phosphorus concentration in reach
!!    p_n              |none          |algal preference factor for ammonia
!!    rch_cbod(:)      |mg O2/L       |carbonaceous biochemical oxygen demand in
!!                                    |reach
!!    rch_dox(:)       |mg O2/L       |dissolved oxygen concentration in reach
!!    rhoq             |1/hr          |algal respiration rate at 20 deg C
!!    rk(1,:)           |1/hr          |CBOD deoxygenation rate coefficient in
!!                                    |reach at 20 deg C
!!    rk(2,:)           |1/hr          |reaeration rate in accordance with Fickian
!!                                    |diffusion in reach at 20 deg C
!!    rk(3,:)           |1/hr          |rate of loss of CBOD due to settling in
!!                                    |reach at 20 deg C
!!    rk(4,:)           |mg O2/        |sediment oxygen demand rate in reach
!!                     |  ((m**2)*hr) |at 20 deg C
!!    rnum1            |none          |fraction of overland flow
!!    rs(1,:)           |m/hr          |local algal settling rate in reach at
!!                                    |20 deg C
!!    rs(2,:)           |(mg disP-P)/  |benthos source rate for dissolved P
!!                     |  ((m**2)*hr) |in reach at 20 deg C
!!    rs(3,:)           |(mg NH4-N)/   |benthos source rate for ammonia nitrogen
!!                     |  ((m**2)*hr) |in reach at 20 deg C
!!    rs(4,:)           |1/hr          |rate coefficient for organic nitrogen
!!                                    |settling in reach at 20 deg C
!!    rs(5,:)           |1/hr          |organic phosphorus settling rate in reach
!!                                    |at 20 deg C
!!    rttime           |hr            |reach travel time
!!    tfact            |none          |fraction of solar radiation that is
!!                                    |photosynthetically active
!!    tmpav(:)         |deg C         |average air temperature on current day
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
!!    soxy        |mg O2/L       |saturation concetration of dissolved oxygen
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    algcon      |mg alg/L      |initial algal biomass concentration in reach
!!    algi        |MJ/(m2*hr)    |photosynthetically active light intensity
!!                               |for hour
!!    algin       |mg alg/L      |algal biomass concentration in inflow
!!    ammoin      |mg N/L        |ammonium N concentration in inflow
!!    bc1mod      |1/day         |rate constant for biological oxidation of NH3
!!                               |to NO2 modified to reflect impact of low
!!                               |oxygen concentration
!!    bc2mod      |1/day         |rate constant for biological oxidation of NO2
!!                               |to NO3 modified to reflect impact of low
!!                               |oxygen concentration
!!    cbodcon     |mg/L          |initial carbonaceous biological oxygen demand
!!                               |concentration in reach
!!    cbodin      |mg/L          |carbonaceous biological oxygen demand
!!                               |concentration in inflow
!!    chlin       |mg chl-a/L    |chlorophyll-a concentration in inflow
!!    cinn        |mg N/L        |effective available nitrogen concentration
!!    cordo       |none          |nitrification rate correction factor
!!    disoxin     |mg O2/L       |dissolved oxygen concentration in inflow
!!    dispin      |mg P/L        |soluble P concentration in inflow
!!    f1          |none          |fraction of algal nitrogen uptake from
!!                               |ammonia pool
!!    fll         |none          |growth attenuation factor for light
!!    fnn         |none          |algal growth limitation factor for nitrogen
!!    fpp         |none          |algal growth limitation factor for phosphorus
!!    gra         |1/hr          |local algal growth rate at 20 deg C
!!    ii          |none          |counter
!!    lambda      |1/m           |light extinction coefficient
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
!!    thbc1       |none          |temperature adjustment factor for local
!!                               |biological oxidation of NH3 to NO2
!!    thbc2       |none          |temperature adjustment factor for local
!!                               |biological oxidation of NO2 to NO3
!!    thbc3       |none          |temperature adjustment factor for local
!!                               |hydrolysis of organic N to ammonia N
!!    thbc4       |none          |temperature adjustment factor for local
!!                               |decay of organic P to dissolved P
!!    thgra       |none          |temperature adjustment factor for local algal
!!                               |growth rate
!!    thour       |none          |flow duration (fraction of hr)
!!    thrho       |none          |temperature adjustment factor for local algal
!!                               |respiration rate
!!    thrk1       |none          |temperature adjustment factor for local CBOD
!!                               |deoxygenation
!!    thrk2       |none          |temperature adjustment factor for local oxygen
!!                               |reaeration rate
!!    thrk3       |none          |temperature adjustment factor for loss of
!!                               |CBOD due to settling
!!    thrk4       |none          |temperature adjustment factor for local
!!                               |sediment oxygen demand
!!    thrs1       |none          |temperature adjustment factor for local algal
!!                               |settling rate
!!    thrs2       |none          |temperature adjustment factor for local
!!                               |benthos source rate for dissolved phosphorus
!!    thrs3       |none          |temperature adjustment factor for local
!!                               |benthos source rate for ammonia nitrogen
!!    thrs4       |none          |temperature adjustment factor for local
!!                               |organic N settling rate
!!    thrs5       |none          |temperature adjustment factor for local
!!                               |organic P settling rate
!!    uu          |varies        |variable to hold intermediate calculation
!!                               |result
!!    vv          |varies        |variable to hold intermediate calculation
!!                               |result
!!    wtmp        |deg C         |temperature of water in reach
!!    wtrin       |m^3 H2O       |water flowing into reach on day
!!    wtrtot      |m^3 H2O       |inflow + storage water
!!    ww          |varies        |variable to hold intermediate calculation
!!                               |result
!!    xx          |varies        |variable to hold intermediate calculation
!!                               |result
!!    yy          |varies        |variable to hold intermediate calculation
!!                               |result
!!    zz          |varies        |variable to hold intermediate calculation
!!                               |result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Log, Min
!!    SWAT: Theta, Oxygen_saturation

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 Theta, Oxygen_saturation
   integer, intent(in) :: jrch, k
   real*8, parameter :: thbc1 = 1.083, thbc2 = 1.047, thbc3 = 1.047,&
      &thbc4 = 1.047, thgra = 1.047, thrho = 1.047, thrk1 = 1.047,&
      &thrk2 = 1.024, thrk3 = 1.024, thrk4 = 1.060, thrs1 = 1.024,&
      &thrs2 = 1.074, thrs3 = 1.074, thrs4 = 1.024, thrs5 = 1.024
   real*8 :: algcon, algi, algin, ammoin, bc1mod, bc2mod, cbodcon, cbodin,&
      &cinn, chlin, cordo, disoxin, dispin, f1, fll, fnn, fpp, gra, lambda,&
      &nh3con, nitratin, nitritin, no2con, no3con, o2con, orgncon, orgnin,&
      &orgpcon, orgpin, solpcon, thour, uu, vv, wtmp, wtrin, wtrtot, ww, xx,&
      &yy, zz
   integer :: ii

!! hourly loop
   do ii = 1, nstep
      !! initialize water flowing into reach
      xx = 1. - rnum1
      wtrin = hhvaroute(2,k,ii) * xx

      if (hrtwtr(ii) / (idt * 60.) > 0.01) then
!! concentrations
         !! initialize inflow concentrations
         chlin = 0.
         algin = 0.
         orgnin = 0.
         ammoin = 0.
         nitritin = 0.
         nitratin = 0.
         orgpin = 0.
         dispin = 0.
         cbodin = 0.
         disoxin = 0.
         cinn = 0.
         if (wtrin > 0.001) then
            xx = xx / wtrin
            chlin = 1000. * hhvaroute(13,k,ii) * xx
            algin = 1000. * chlin / ai0        !! QUAL2E equation III-1
            orgnin = 1000. * hhvaroute(4,k,ii) * xx
            ammoin = 1000. * hhvaroute(14,k,ii) * xx
            nitritin = 1000. * hhvaroute(15,k,ii) * xx
            nitratin = 1000. * hhvaroute(6,k,ii) * xx
            orgpin = 1000. * hhvaroute(5,k,ii) * xx
            dispin = 1000. * hhvaroute(7,k,ii) * xx
            cbodin = 1000. * hhvaroute(16,k,ii) * xx
            disoxin= 1000. * hhvaroute(17,k,ii) * xx
         end if

         if (chlin < 1.e-6) chlin = 0.0
         if (algin < 1.e-6) algin = 0.0
         if (orgnin < 1.e-6) orgnin = 0.0
         if (ammoin < 1.e-6) ammoin = 0.0
         if (nitritin < 1.e-6) nitritin = 0.0
         if (nitratin < 1.e-6) nitratin = 0.0
         if (orgpin < 1.e-6) orgpin = 0.0
         if (dispin < 1.e-6) dispin = 0.0
         if (cbodin < 1.e-6) cbodin = 0.0
         if (disoxin < 1.e-6) disoxin = 0.0

         !! initialize concentration of nutrient in reach
         wtrtot = wtrin + hrchwtr(ii)
         if (ii == 1) then
            algcon = (algin * wtrin + algae(jrch) * hrchwtr(ii)) / wtrtot
            orgncon = (orgnin * wtrin + organicn(jrch) * hrchwtr(ii)) / wtrtot
            nh3con = (ammoin * wtrin + ammonian(jrch) * hrchwtr(ii)) / wtrtot
            no2con = (nitritin * wtrin + nitriten(jrch) * hrchwtr(ii)) / wtrtot
            no3con = (nitratin * wtrin + nitraten(jrch) * hrchwtr(ii)) / wtrtot
            orgpcon = (orgpin * wtrin + organicp(jrch) * hrchwtr(ii)) / wtrtot
            solpcon = (dispin * wtrin + disolvp(jrch) * hrchwtr(ii)) / wtrtot
            cbodcon = (cbodin * wtrin + rch_cbod(jrch) * hrchwtr(ii)) / wtrtot
            o2con = (disoxin * wtrin + rch_dox(jrch) * hrchwtr(ii)) / wtrtot
         else
            algcon = (algin * wtrin + halgae(ii-1) * hrchwtr(ii)) / wtrtot
            orgncon = (orgnin * wtrin + horgn(ii-1) * hrchwtr(ii)) / wtrtot
            nh3con = (ammoin * wtrin + hnh4(ii-1) * hrchwtr(ii)) / wtrtot
            no2con = (nitritin * wtrin + hno2(ii-1) * hrchwtr(ii)) / wtrtot
            no3con = (nitratin * wtrin + hno3(ii-1) * hrchwtr(ii)) / wtrtot
            orgpcon = (orgpin * wtrin + horgp(ii-1) * hrchwtr(ii)) / wtrtot
            solpcon = (dispin * wtrin + hsolp(ii-1) * hrchwtr(ii)) / wtrtot
            cbodcon = (cbodin * wtrin + hbod(ii-1) * hrchwtr(ii)) / wtrtot
            o2con = (disoxin * wtrin + hdisox(ii-1) * hrchwtr(ii)) / wtrtot
         end if

         if (algcon < 1.e-6) algcon = 0.0
         if (orgncon < 1.e-6) orgncon = 0.0
         if (nh3con < 1.e-6) nh3con = 0.0
         if (no2con < 1.e-6) no2con = 0.0
         if (no3con < 1.e-6) no3con = 0.0
         if (orgpcon < 1.e-6) orgpcon = 0.0
         if (solpcon < 1.e-6) solpcon = 0.0
         if (cbodcon < 1.e-6) cbodcon = 0.0
         if (o2con < 1.e-6) o2con = 0.0
         !! calculate temperature in stream
         !! Stefan and Preudhomme. 1993.  Stream temperature estimation
         !! from air temperature.  Water Res. Bull. p. 27-45
         !! SWAT manual equation 2.3.13
         wtmp = 5.0 + 0.75 * tmpav(jrch)
         if (wtmp <= 0.) wtmp = 0.1

         !! calculate effective concentration of available nitrogen
         !! QUAL2E equation III-15
         cinn = nh3con + no3con

         !! calculate saturation concentration for dissolved oxygen
         !! QUAL2E section 3.6.1 equation III-29
         soxy = Oxygen_saturation(wtmp)
!! end initialize concentrations

!! O2 impact calculations
         !! calculate nitrification rate correction factor for low
         !! oxygen QUAL2E equation III-21
         cordo = 1.0 - Exp(-0.6 * o2con)
         !! modify ammonia and nitrite oxidation rates to account for
         !! low oxygen
         bc1mod = bc(1,jrch) * cordo
         bc2mod = bc(2,jrch) * cordo
!! end O2 impact calculations

         !! calculate flow duration
         !thour = hhtime(ii)
         !if (thour > 1.0) thour = 1.0
         thour = 1.0 ! this overwrites previous lines

!! algal growth
         !! calculate light extinction coefficient
         !! (algal self shading) QUAL2E equation III-12
         if (ai0 * algcon > 1.e-6) then
            lambda = lambda0 + (lambda1 * ai0 * algcon) + lambda2 *&
               &(ai0 * algcon) ** (.66667)
         else
            lambda = lambda0
         endif

         !! calculate algal growth limitation factors for nitrogen
         !! and phosphorus QUAL2E equations III-13 & III-14
         fnn = cinn / (cinn + k_n)
         fpp = solpcon / (solpcon + k_p)

         !! calculate hourly, photosynthetically active,
         !! light intensity QUAL2E equation III-9c
         !! Light Averaging Option # 3
         algi = frad(hru1(jrch),ii) * hru_ra(hru1(jrch)) * tfact

         !! calculate growth attenuation factor for light, based on
         !! hourly light intensity QUAL2E equation III-6a
         fll = (1. / (lambda * hdepth(ii))) *&
            &Log((k_l + algi) / (k_l + algi * (Exp(-lambda * hdepth(ii)))))

         !! calculcate local algal growth rate
         gra = 0.
         select case (igropt)
          case (1)
            !! multiplicative QUAL2E equation III-3a
            gra = mumax * fll * fnn * fpp
          case (2)
            !! limiting nutrient QUAL2E equation III-3b
            gra = mumax * fll * Min(fnn, fpp)
          case (3)
            !! harmonic mean QUAL2E equation III-3c
            if (fnn > 1.e-6 .and. fpp > 1.e-6) then
               gra = mumax * fll * 2. / ((1. / fnn) + (1. / fpp))
            else
               gra = 0.
            endif
         end select

         !! calculate algal biomass concentration at end of day
         !! (phytoplanktonic algae)
         !! QUAL2E equation III-2
         halgae(ii) = algcon + (Theta(gra,thgra,wtmp) * algcon -&
            &Theta(rhoq,thrho,wtmp) * algcon - Theta(rs(1,jrch),thrs1,wtmp)&
            &/ hdepth(ii) * algcon) * thour
         if (halgae(ii) < 1.e-6) halgae(ii) = 0.0

         !! calculate chlorophyll-a concentration at end of day
         !! QUAL2E equation III-1
         hchla(ii) = halgae(ii) * ai0 / 1000.
         if (hchla(ii) < 1.e-6) hchla(ii) = 0.0
!! end algal growth

!! oxygen calculations
         !! calculate carbonaceous biological oxygen demand at end
         !! of day QUAL2E section 3.5 equation III-26
         yy = Theta(rk(1,jrch),thrk1,wtmp) * cbodcon
         zz = Theta(rk(3,jrch),thrk3,wtmp) * cbodcon
         hbod(ii) = cbodcon - (yy + zz) * thour
         if (hbod(ii) < 1.e-6) hbod(ii) = 0.0

         !! calculate dissolved oxygen concentration if reach at
         !! end of day QUAL2E section 3.6 equation III-28
         uu = Theta(rk(2,jrch),thrk2,wtmp) * (soxy - o2con)
         vv = (ai3 * Theta(gra,thgra,wtmp) - ai4 *&
            &Theta(rhoq,thrho,wtmp)) * algcon
         ww = Theta(rk(1,jrch),thrk1,wtmp) * cbodcon
         xx = Theta(rk(4,jrch),thrk4,wtmp) / (hdepth(ii) * 1000.)
         yy = ai5 * Theta(bc1mod,thbc1,wtmp) * nh3con
         zz = ai6 * Theta(bc2mod,thbc2,wtmp) * no2con
         hdisox(ii) = o2con + (uu + vv - ww - xx - yy - zz) * thour
         if (hdisox(ii) < 1.e-6) hdisox(ii) = 0.0
!! end oxygen calculations

!! nitrogen calculations
         !! calculate organic N concentration at end of day
         !! QUAL2E section 3.3.1 equation III-16
         xx = ai1 * Theta(rhoq,thrho,wtmp) * algcon
         yy = Theta(bc(3,jrch),thbc3,wtmp) * orgncon
         zz = Theta(rs(4,jrch),thrs4,wtmp) * orgncon
         horgn(ii) = orgncon + (xx - yy - zz) * thour
         if (horgn(ii) < 1.e-6) horgn(ii) = 0.0

         !! calculate fraction of algal nitrogen uptake from ammonia
         !! pool QUAL2E equation III-18
         f1 = p_n * nh3con / (p_n * nh3con + (1. - p_n) * no3con + 1.e-6)

         !! calculate ammonia nitrogen concentration at end of day
         !! QUAL2E section 3.3.2 equation III-17
         ww = Theta(bc(3,jrch),thbc3,wtmp) * orgncon
         xx = Theta(bc1mod,thbc1,wtmp) * nh3con
         yy = Theta(rs(3,jrch),thrs3,wtmp) / (hdepth(ii) * 1000.)
         zz = f1 * ai1 * algcon * Theta(gra,thgra,wtmp)
         hnh4(ii) = nh3con + (ww - xx + yy - zz) * thour
         if (hnh4(ii) < 1.e-6) hnh4(ii) = 0.0

         !! calculate concentration of nitrite at end of day
         !! QUAL2E section 3.3.3 equation III-19
         yy = Theta(bc1mod,thbc1,wtmp) * nh3con
         zz = Theta(bc2mod,thbc2,wtmp) * no2con
         hno2(ii) = no2con + (yy - zz) * thour
         if (hno2(ii) < 1.e-6) hno2(ii) = 0.0

         !! calculate nitrate concentration at end of day
         !! QUAL2E section 3.3.4 equation III-20
         yy = Theta(bc2mod,thbc2,wtmp) * no2con
         zz = (1. - f1) * ai1 * algcon * Theta(gra,thgra,wtmp)
         hno3(ii) = no3con + (yy - zz) * thour
         if (hno3(ii) < 1.e-6) hno3(ii) = 0.0
!! end nitrogen calculations

!! phosphorus calculations
         !! calculate organic phosphorus concentration at end of
         !! day QUAL2E section 3.3.6 equation III-24
         xx = ai2 * Theta(rhoq,thrho,wtmp) * algcon
         yy = Theta(bc(4,jrch),thbc4,wtmp) * orgpcon
         zz = Theta(rs(5,jrch),thrs5,wtmp) * orgpcon
         horgp(ii) = orgpcon + (xx - yy - zz) * thour
         if (horgp(ii) < 1.e-6) horgp(ii) = 0.0

         !! calculate dissolved phosphorus concentration at end
         !! of day QUAL2E section 3.4.2 equation III-25
         xx = Theta(bc(4,jrch),thbc4,wtmp) * orgpcon
         yy = Theta(rs(2,jrch),thrs2,wtmp) / (hdepth(ii) * 1000.)
         zz = ai2 * Theta(gra,thgra,wtmp) * algcon
         hsolp(ii) = solpcon + (xx + yy - zz) * thour
         if (hsolp(ii) < 1.e-6) hsolp(ii) = 0.0
!! end phosphorus calculations

      else
         !! all water quality variables set to zero when no flow
         halgae(ii) = 0.0
         hchla(ii) = 0.0
         horgn(ii) = 0.0
         hnh4(ii) = 0.0
         hno2(ii) = 0.0
         hno3(ii) = 0.0
         horgp(ii) = 0.0
         hsolp(ii) = 0.0
         hbod(ii) = 0.0
         hdisox(ii) = 0.0
      endif

   end do
!! end hourly loop

!! set end of day concentrations
   algae(jrch) = halgae(nstep)
   chlora(jrch) = hchla(nstep)
   organicn(jrch) = horgn(nstep)
   ammonian(jrch) = hnh4(nstep)
   nitriten(jrch) = hno2(nstep)
   nitraten(jrch) = hno3(nstep)
   organicp(jrch) = horgp(nstep)
   disolvp(jrch) = hsolp(nstep)
   rch_cbod(jrch) = hbod(nstep)
   rch_dox(jrch) = hdisox(nstep)

   return
end
