!> @file pond.f90
!> file containing the subroutine pond
!> @author
!> modified by Javier Burguete

!> this subroutine routes water and sediment through ponds
!> and computes evaporation and seepage from the ponds
!> @param[in] k HRU or reach number (none)
subroutine pond(k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k           |none          |HRU or reach number
!!    bp(1,:)     |none          |1st shape parameter for pond surface area
!!                               |equation
!!    bp(2,:)     |none          |2nd shape parameter for the pond surface area
!!                               |equation
!!    chlap(:)    |none          |chlorophyll-a production coefficient for pond
!!    hru_sub(:)  |none          |subbasin in which HRU/reach is located
!!    iflod(1,:)  |none          |beginning month of non-flood season
!!    iflod(2,:)  |none          |ending month of non-flood season
!!    ipnd(1,:)   |none          |beginning month of 2nd "season" of nutrient
!!                               |settling
!!    ipnd(2,:)   |none          |ending month of 2nd "season" of nutrient
!!                               |settling
!!    i_mo        |none          |current month of simulation
!!    ndtarg(:)   |none          |number of days required to reach target
!!                               |storage from current pond storage
!!    nsetlp(1,:) |m/day         |nitrogen settling rate for 1st season
!!    nsetlp(2,:) |m/day         |nitrogen settling rate for 2nd season
!!    pet_day     |mm H2O        |potential evapotranspiration on day
!!    pnd_evol(:) |m^3 H2O       |volume of water required to fill pond
!!                               |to the emergency spillway
!!    pnd_fr(:)   |none          |fraction of HRU/subbasin area that drains
!!                               |into ponds
!!    pnd_k(:)    |mm/hr         |hydraulic conductivity through bottom of
!!                               |ponds
!!    pnd_no3(:)  |kg N          |amount of nitrate originating from surface
!!                               |runoff in pond at beginning of day
!!    pnd_no3g(:) |kg N          |amount of nitrate originating from
!!                               |groundwater in pond at beginning of day
!!    pnd_no3s(:) |kg N          |amount of nitrate originating from lateral
!!                               |flow in pond at beginning of day
!!    pnd_nsed(:) |kg/L          |normal ratio of sediment to water in pond
!!    pnd_orgn(:) |kg N          |amount of organic N originating from
!!                               |surface runoff in pond at beginning of day
!!    pnd_orgp(:) |kg P          |amount of organic P originating from
!!                               |surface runoff in pond at beginning of day
!!    pnd_psed(:) |kg P          |amount of mineral P attached to sediment
!!                               |originating from surface runoff in pond at
!!                               |beginning of day
!!    pnd_pvol(:) |m^3 H2O       |volume of water required to fill pond
!!                               |to the principal spillway
!!    pnd_sed(:)  |kg/L          |ratio of sediment to water in pond
!!    pnd_solp(:) |kg P          |amount of soluble P originating from surface
!!                               |runoff in pond at beginning of day
!!    pnd_solpg(:)|kg P          |amount of soluble P originating from
!!                               |groundwater in pond at beginning of day
!!    pnd_vol(:)  |m^3 H2O       |volume of water in pond
!!    pndflwi     |m^3 H2O       |volume of water flowing into pond on day
!!    pndsedin    |metric tons   |sediment entering pond during day
!!    psetlp(1,:) |m/day         |phosphorus settling rate for 1st season
!!    psetlp(2,:) |m/day         |phosphorus settling rate for 2nd season
!!    seccip(:)   |none          |water clarity coefficient for pond
!!    sed_stl(:)  |kg/kg         |fraction of sediment remaining suspended in
!!                               |impoundment after settling for one day
!!    sol_sumfc(:)|mm H2O        |amount of water held in the soil profile
!!                               |at field capacity
!!    sol_sw(:)   |mm H2O        |amount of water stored in the soil profile
!!                               |on any given day
!!    subp(:)     |mm H2O        |precipitation for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    pnd_chla(:) |kg chl_a      |amount of chlorophyll-a in pond at end of day
!!    pnd_no3(:)  |kg N          |amount of nitrate originating from surface
!!                               |runoff in pond at end of day
!!    pnd_no3g(:) |kg N          |amount of nitrate originating from
!!                               |groundwater in pond at end of day
!!    pnd_no3s(:) |kg N          |amount of nitrate originating from lateral
!!                               |flow in pond at end of day
!!    pnd_orgn(:) |kg N          |amount of organic N originating from
!!                               |surface runoff in pond at end of day
!!    pnd_orgp(:) |kg P          |amount of organic P originating from
!!                               |surface runoff in pond at end of day
!!    pnd_psed(:) |kg P          |amount of mineral P attached to sediment
!!                               |originating from surface runoff in pond at
!!                               |end of day
!!    pnd_seci(:) |m             |secchi-disk depth of pond
!!    pnd_sed(:)  |kg/L          |ratio of sediment to water in pond
!!    pnd_solp(:) |kg P          |amount of soluble P originating from surface
!!                               |runoff in pond at end of day
!!    pnd_solpg(:)|kg P          |amount of soluble P originating from
!!                               |groundwater in pond at end of day
!!    pnd_vol(:)  |m^3 H2O       |volume of water in pond
!!    pndev       |m^3 H2O       |evaporation from pond on day
!!    pndflwo     |m^3 H2O       |volume of water flowing out of pond on day
!!    pndpcp      |m^3 H2O       |precipitation on pond during day
!!    pndsedc     |metric tons   |net change in sediment in pond during day
!!    pndsedo     |metric tons   |sediment leaving pond during day
!!    pndsep      |m^3 H2O       |seepage from pond on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chlaco      |ppb (ug/L)    |concentration of chlorophyll-a in pond
!!    cla
!!    finsed
!!    inised
!!    lag
!!    nitrok      |none          |fraction of nitrogen in pond removed by
!!                               |settling
!!    phosk       |none          |fraction of phosphorus in pond removed by
!!                               |settling
!!    pndsa       |ha            |surface area of pond on current day
!!    remsetsed
!!    sag
!!    san
!!    sed         |kg/L          |sediment concentration in pond at beginning of
!!                               |day
!!    setsed
!!    sil
!!    susp
!!    targ        |m^3 H2O       |target storage level in pond
!!    tpco        |ppb (ug/L)    |concentration of phosphorus in pond water
!!                               |on day
!!    trappnd
!!    velofl
!!    vol         |m^3 H2O       |volume of water in pond at beginning of day
!!    xx          |none          |variable to hold intermediate calc result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: k
   real*8 :: chlaco, cla, finsed, inised, lag, nitrok, phosk, pndsa, remsetsed,&
      &sag, san, sed, setsed, sil, susp, targ, tpco, trappnd, velofl, vol, xx


   !! store initial values
   vol = pnd_vol(k)
   sed = pnd_sed(k)
   san = pnd_san(k)
   sil = pnd_sil(k)
   cla = pnd_cla(k)
   sag = pnd_sag(k)
   lag = pnd_lag(k)

   !! calculate water balance for day
   pndsa = hru_fr(k) * bp(1,k) * pnd_vol(k) ** bp(2,k)
   pndev = 10. * evpnd(k) * pet_day * pndsa
   pndsep = pnd_k(k) * pndsa * 240.
   pndpcp = subp(k) * pndsa * 10.

   !! new water volume for day
   pnd_vol(k) = pnd_vol(k) - pndsep - pndev + pndpcp + pndflwi

   if (pnd_vol(k) < 0.001) then
      !! if volume deficit in pond reduce seepage
      !! so that the pond volume is zero
      pndsep = pndsep + pnd_vol(k)
      pnd_vol(k) = 0.
      !! if seepage is less than the volume deficit, take the remainder
      !! from evaporation
      if (pndsep < 0.) then
         pndev = pndev + pndsep
         pndsep = 0.
      end if
      pnd_sed(k) = 0.
      pnd_san(k) = 0.
      pnd_sil(k) = 0.
      pnd_cla(k) = 0.
      pnd_sag(k) = 0.
      pnd_lag(k) = 0.
      pnd_solp(k) = 0.
      pnd_psed(k) = 0.
      pnd_orgp(k) = 0.
      pnd_solpg(k) = 0.
      pnd_orgn(k) = 0.
      pnd_no3(k) = 0.
      pnd_no3s(k) = 0.
      pnd_no3g(k) = 0.
      pnd_chla(k) = 0.
      pnd_seci(k) = 0.

   else

      !! compute outflow
      if (pnd_evol(k) <= 0.) then
         !! all storage over principle is outflow
         if (pnd_vol(k) <= pnd_pvol(k)) then
            pndflwo = 0.
         else
            pndflwo = pnd_vol(k) - pnd_pvol(k)
         end if
      elseif (pnd_vol(k) > pnd_evol(k)) then
         !! if emergency level is defined, anytime pond volume
         !! exceeds this level, all excess is released
         pndflwo = pnd_vol(k) - pnd_evol(k)
      else
         !! target storage based on flood season and soil water
         targ = 0.
         if (iflod(2,k) > iflod(1,k)) then
            if (i_mo > iflod(1,k) .and. i_mo < iflod(2,k)) then
               targ = pnd_evol(k)
            else
               xx = Min(sol_sw(k) / sol_sumfc(k),1.)
               targ = pnd_pvol(k) + .5 * (1. - xx) * (pnd_evol(k) -&
                  &pnd_pvol(k))
            end if
         else
            if (i_mo > iflod(1,k) .or. i_mo < iflod(2,k)) then
               targ = pnd_evol(k)
            else
               xx = Min(sol_sw(k) / sol_sumfc(k),1.)
               targ = pnd_pvol(k) + .5 * (1. - xx) * (pnd_evol(k) -&
                  &pnd_pvol(k))
            end if
         end if
         if (pnd_vol(k) > targ) then
            pndflwo = (pnd_vol(k) - targ) / ndtarg(k)
         else
            pndflwo = 0.
         end if
      end if

      !! compute new sediment concentration
      if (pndsedin < 1.e-6) pndsedin = 0.
      if (pndsa == 0.) pndsa = 0.001    !!MJW added line of code 040811
      velofl = (pndflwo / pndsa) / 10000.
      if (velofl > 1.e-6) then
         trappnd = velsetlp(k) / velofl
         if (trappnd > 1.) trappnd = 1.
         susp = 1. - trappnd
      else
         susp = 1.
      endif

      xx = pnd_vol(k)
      pnd_sed(k) = (sed * vol + susp * pndsedin) / xx
      pnd_san(k) = (san * vol + pndsanin) / xx
      pnd_sil(k) = (sil * vol + pndsilin) / xx
      pnd_cla(k) = (cla * vol + pndclain) / xx
      pnd_sag(k) = (sag * vol + pndsagin) / xx
      pnd_lag(k) = (lag * vol + pndlagin) / xx

      !! compute final pond volume
      pnd_vol(k) = pnd_vol(k) - pndflwo
      if (pnd_vol(k) < 0.) then
         pndflwo = pndflwo + pnd_vol(k)
         pnd_vol(k) = 0.
      endif

      !! compute change in sediment concentration due to settling
      if (sed_stl(k) < 1.e-6) sed_stl(k) = 0.0
      if (pnd_sed(k) > pnd_nsed(k)) then
         inised = pnd_sed(k)
         pnd_sed(k) = (pnd_sed(k) - pnd_nsed(k)) * sed_stl(k) + pnd_nsed(k)
         finsed = pnd_sed(k)
         setsed = inised - finsed

         if (pnd_lag(k) >= setsed) then
            pnd_lag(k) = pnd_lag(k) - setsed
         else
            remsetsed = setsed - pnd_lag(k)
            pnd_lag(k) = 0.
            if (pnd_san(k) >= remsetsed) then
               pnd_san(k) = pnd_san(k) - remsetsed
            else
               remsetsed = remsetsed - pnd_san(k)
               pnd_san(k) = 0.
               if (pnd_sag(k) >= remsetsed) then
                  pnd_sag(k) = pnd_sag(k) - remsetsed
               else
                  remsetsed = remsetsed - pnd_sag(k)
                  pnd_sag(k) = 0.
                  if (pnd_sil(k) >= remsetsed) then
                     pnd_sil(k) = pnd_sil(k) - remsetsed
                  else
                     remsetsed = remsetsed - pnd_sil(k)
                     pnd_sil(k) = 0.
                     if (pnd_cla(k) >= remsetsed) then
                        pnd_cla(k) = pnd_cla(k) - remsetsed
                     else
                        pnd_cla(k) = 0.
                     end if
                  end if
               end if
            end if
         end if

      end if
      !! compute sediment leaving pond
      pndsedo = pnd_sed(k) * pndflwo
      pndsano = pnd_san(k) * pndflwo
      pndsilo = pnd_sil(k) * pndflwo
      pndclao = pnd_cla(k) * pndflwo
      pndsago = pnd_sag(k) * pndflwo
      pndlago = pnd_lag(k) * pndflwo

      !! net change in amount of sediment in pond for day
      pndsedc = vol * sed + pndsedin - pndsedo - pnd_sed(k) * pnd_vol(k)

      !! determine settling rate
      !! part of equation 29.1.3 in SWAT manual
      if (i_mo >= ipnd(1,k) .and. i_mo <= ipnd(2,k)) then
         phosk = psetlp(1,k)
         nitrok = nsetlp(1,k)
      else
         phosk = psetlp(2,k)
         nitrok = nsetlp(2,k)
      endif
      xx = pndsa * 10000. / pnd_vol(k)
      phosk = phosk * xx  !setl/mean depth
      phosk = Min(phosk, 1.)
      nitrok = nitrok * xx !setl/mean depth
      nitrok = Min(nitrok, 1.)

      !! remove nutrients by settling
      !! other part of equation 29.1.3 in SWAT manual
      xx = 1. - phosk
      pnd_solp(k) = pnd_solp(k) * xx
      pnd_psed(k) = pnd_psed(k) * xx
      pnd_orgp(k) = pnd_orgp(k) * xx
      pnd_solpg(k) = pnd_solpg(k) * xx
      xx = 1. - nitrok
      pnd_orgn(k) = pnd_orgn(k) * xx
      pnd_no3(k) = pnd_no3(k) * xx
      pnd_no3s(k) = pnd_no3s(k) * xx
      pnd_no3g(k) = pnd_no3g(k) * xx

      if (pnd_vol(k) + pndflwo > 0.1) then
         tpco = 1.e+6 * (pnd_solp(k) + pnd_orgp(k) + pnd_psed(k) +&
            &pnd_solpg(k)) / (pnd_vol(k) + pndflwo)
      else
         tpco = 0.
      endif
      chlaco = 0.
      pnd_chla(k) = 0.
      if (tpco > 1.e-4) then
         !! equation 29.1.6 in SWAT manual
         chlaco = chlap(k) * 0.551 * (tpco**0.76)
         pnd_chla(k) = chlaco * (pnd_vol(k) + pndflwo) * 1.e-6
      endif
      pnd_seci(k) = 0.
      if (chlaco > 1.e-4) then
         !! equation 29.1.8 in SWAT manual
         pnd_seci(k) = seccip(k) * 6.35 * (chlaco**(-0.473))
      endif
   end if

   return
end
