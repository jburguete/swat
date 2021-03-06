!> @file pondhr.f90
!> file containing the subroutine pondhr
!> @author
!> modified by Javier Burguete

!> @param[in] j HRU or reach number (none)
!> @param[in] k current time step of the day (none)
subroutine pondhr(j,k)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes water and sediment through ponds
!!    and computes evaporation and seepage from the ponds

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU or reach number
!!    k           |none          |current time step of the day
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
!!    rainsub(:,:)|mm H2O        |precipitation for the time step during the
!!                               |day in HRU
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
!!    nitrok      |none          |fraction of nitrogen in pond removed by
!!                               |settling
!!    phosk       |none          |fraction of phosphorus in pond removed by
!!                               |settling
!!    pndsa       |ha            |surface area of pond on current day
!!    sed         |kg/L          |sediment concentration in pond at beginning of
!!                               |day
!!    targ        |m^3 H2O       |target storage level in pond
!!    tpco        |ppb (ug/L)    |concentration of phosphorus in pond water
!!                               |on day
!!    vol         |m^3 H2O       |volume of water in pond at beginning of day
!!    xx          |none          |variable to hold intermediate calc result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: j, k
   real*8 :: chlaco, nitrok, phosk, pndsa, sed, targ, tpco, vol, xx


   !! store initial values
   vol = pnd_vol(j)
   sed = pnd_sed(j)

   !! calculate water balance for day
   pndsa = bp(1,j) * pnd_vol(j) ** bp(2,j)
   pndev = 6. * pet_day * pndsa / nstep       !! urban modeling by J.Jeong
   pndsep = pnd_k(j) * pndsa * 240./ nstep    !! urban modeling by J.Jeong
   pndpcp = rainsub(j,k) * pndsa * 10.        !! urban modeling by J.Jeong

   !! new water volume for time step
   pnd_vol(j) = pnd_vol(j) - pndsep - pndev + pndpcp + pndflwi

   if (pnd_vol(j) < 0.001) then
      !! if volume deficit in pond reduce seepage so that the pond volume is zero
      pndsep = pndsep + pnd_vol(j)
      pnd_vol(j) = 0.
      !! if seepage is less than the volume deficit, take the remainder from evaporation
      if (pndsep < 0.) then
         pndev = pndev + pndsep
         pndsep = 0.
      end if
      pnd_sed(j) = 0.
      pnd_solp(j) = 0.
      pnd_psed(j) = 0.
      pnd_orgp(j) = 0.
      pnd_solpg(j) = 0.
      pnd_orgn(j) = 0.
      pnd_no3(j) = 0.
      pnd_no3s(j) = 0.
      pnd_no3g(j) = 0.
      pnd_chla(j) = 0.
      pnd_seci(j) = 0.

   else

      !! compute new sediment concentration
      pnd_sed(j) = (sed * vol + pndsedin) / pnd_vol(j)

      !! compute outflow
      if (pnd_evol(j) <= 0.) then
         !! all storage over principle is outflow
         if (pnd_vol(j) <= pnd_pvol(j)) then
            pndflwo = 0.
         else
            pndflwo = pnd_vol(j) - pnd_pvol(j)
         end if
      elseif (pnd_vol(j) > pnd_evol(j)) then
         !! if emergency level is defined, anytime pond volume
         !! exceeds this level, all excess is released
         pndflwo = pnd_vol(j) - pnd_evol(j)
      else
         !! target storage based on flood season and soil water
         xx = 0.
         if (iflod(2,j) > iflod(1,j)) then
            if (i_mo > iflod(1,j) .and. i_mo < iflod(2,j)) then
               targ = pnd_evol(j)
            else
               xx = Min(sol_sw(j) / sol_sumfc(j),1.)
               targ = pnd_pvol(j) + .5 * (1. - xx) *&
                  &(pnd_evol(j) - pnd_pvol(j))
            end if
         else
            if (i_mo > iflod(1,j) .or. i_mo < iflod(2,j)) then
               targ = pnd_evol(j)
            else
               xx = Min(sol_sw(j) / sol_sumfc(j),1.)
               targ = pnd_pvol(j) + .5 * (1. - xx) *&
                  &(pnd_evol(j) - pnd_pvol(j))
            end if
         end if
         if (pnd_vol(j) > targ) then
            pndflwo = (pnd_vol(j) - targ) / ndtarg(j)
         else
            pndflwo = 0.
         end if
      end if

      !! compute final pond volume
      pnd_vol(j) = pnd_vol(j) - pndflwo
      if (pnd_vol(j) < 0.) then
         pndflwo = pndflwo + pnd_vol(j)
         pnd_vol(j) = 0.
      endif

      !! compute change in sediment concentration due to settling
      if (pnd_sed(j) > pnd_nsed(j)) then
         pnd_sed(j) = (pnd_sed(j) - pnd_nsed(j)) * sed_stl(j) + pnd_nsed(j)
      end if

      !! compute sediment leaving pond
      pndsedo = pnd_sed(j) * pndflwo

      !! net change in amount of sediment in pond for day
      pndsedc = vol * sed + pndsedin - pndsedo - pnd_sed(j) * pnd_vol(j)

      !! determine settling rate
      !! part of equation 29.1.3 in SWAT manual
      if (i_mo >= ipnd(1,j) .and. i_mo <= ipnd(2,j)) then
         phosk = psetlp(1,j)
         nitrok = nsetlp(1,j)
      else
         phosk = psetlp(2,j)
         nitrok = nsetlp(2,j)
      endif
      xx = pndsa * 10000. / pnd_vol(k)
      phosk = phosk * xx  !setl/mean depth
      phosk = Min(phosk, 1.)
      nitrok = nitrok * xx !setl/mean depth
      nitrok = Min(nitrok, 1.)

      !! remove nutrients by settling
      !! other part of equation 29.1.3 in SWAT manual
      xx = 1. - phosk
      pnd_solp(j) = pnd_solp(j) * xx
      pnd_psed(j) = pnd_psed(j) * xx
      pnd_orgp(j) = pnd_orgp(j) * xx
      pnd_solpg(j) = pnd_solpg(j) * xx
      xx = 1. - nitrok
      pnd_orgn(j) = pnd_orgn(j) * xx
      pnd_no3(j) = pnd_no3(j) * xx
      pnd_no3s(j) = pnd_no3s(j) * xx
      pnd_no3g(j) = pnd_no3g(j) * xx

      if (pnd_vol(j) + pndflwo > 0.1) then
         tpco = 1.e+6 * (pnd_solp(j) + pnd_orgp(j) + pnd_psed(j) +&
            &pnd_solpg(j)) / (pnd_vol(j) + pndflwo)
      else
         tpco = 0.
      endif
      chlaco = 0.
      pnd_chla(j) = 0.
      pnd_seci(j) = 0.
      if (tpco > 1.e-4) then
         !! equation 29.1.6 in SWAT manual
         chlaco = chlap(j) * 0.551 * (tpco**0.76)
         pnd_chla(j) = chlaco * (pnd_vol(j) + pndflwo) * 1.e-6
      endif
      if (chlaco > 1.e-4) then
         !! equation 29.1.8 in SWAT manual
         pnd_seci(j) = seccip(j) * 6.35 * (chlaco**(-0.473))
      endif
   end if

   return
end
