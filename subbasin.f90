!> @file subbasin.f90
!> file containing the subroutine subbasin
!> @author
!> modified by Javier Burguete

!> this subroutine controls the simulation of the land phase of the
!> hydrologic cycle
!> @param[in] i current day in simulation--loop counter (julian date)
subroutine subbasin(i)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i              |julian date   |current day in simulation--loop counter
!!    auto_wstr(:)   |none          |water stress factor which triggers auto
!!                                  |irrigation
!!    bio_e(:)       |(kg/ha)/      |biomass-energy ratio
!!                   |     (MJ/m**2)|The potential (unstressed) growth rate per
!!                                  |unit of intercepted photosynthetically
!!                                  |active radiation.
!!    canev          |mm H2O        |amount of water evaporated from canopy
!!                                  |storage
!!    ep_day         |mm H2O        |actual amount of transpiration that occurs
!!                                  |on day in HRU
!!    es_day         |mm H2O        |actual amount of evaporation (soil et) that
!!                                  |occurs on day in HRU
!!    gw_q(:)        |mm H2O        |groundwater contribution to streamflow from
!!                                  |HRU on current day
!!    hru_ra(:)      |MJ/m^2        |solar radiation for the day in HRU
!!    iida           |julian date   |day being simulated (current julian date)
!!    idplt(:)       |none          |land cover code from crop.dat
!!    igro(:)        |none          |land cover status code
!!                                  |0 no land cover currently growing
!!                                  |1 land cover growing
!!    inum1          |none          |subbasin number
!!    imp_trig(:)    |none          |release/impound action code:
!!                                  |0 begin impounding water
!!                                  |1 release impounded water
!!    irrsc(:)       |none          |irrigation source code:
!!                                  |1 divert water from reach
!!                                  |2 divert water from reservoir
!!                                  |3 divert water from shallow aquifer
!!                                  |4 divert water from deep aquifer
!!                                  |5 divert water from source outside
!!                                  |  watershed
!!    iurban(:)      |none          |urban simulation code:
!!                                  |0  no urban sections in HRU
!!                                  |1  urban sections in HRU, simulate using
!!                                  |   USGS regression equations
!!                                  |2  urban sections in HRU, simulate using
!!                                  |   build up/wash off algorithm
!!    latq(:)        |mm H2O        |total lateral flow in soil profile for the
!!                                  |day in HRU
!!    nafert(:)      |none          |sequence number of auto-fert application
!!                                  |within the year
!!    nair(:)        |none          |sequence number of auto-irrigation
!!                                  |application within the year
!!    nfert(:)       |none          |sequence number of fertilizer application
!!                                  |within the year
!!    nirr(:)        |none          |sequence number of irrigation application
!!                                  |within the year
!!    nrelease(:)    |none          |sequence number of impound/release
!!                                  |operation within the year
!!    nro(:)         |none          |sequence number of year in rotation
!!    peakr          |m^3/s         |peak runoff rate
!!    pet_day        |mm H2O        |potential evapotranspiration on current
!!                                  |day in HRU
!!    phuacc(:)      |none          |fraction of plant heat units accumulated
!!    phubase(:)     |heat units    |base zero total heat units (used when no
!!                                  |land cover is growing)
!!                                  |pesticide application occurs
!!    pot_fr(:)      |km2/km2       |fraction of HRU area that drains into
!!                                  |pothole
!!    pot_vol(:)     |m**3 H2O      |current volume of water stored in the
!!                                  |depression/impounded area
!!    precipday      |mm H2O        |precipitation for the day in HRU
!!    qday           |mm H2O        |surface runoff loading to main channel from
!!                                  |HRU for day
!!    qtile          |mm H2O        |drainage tile flow in soil layer for the
!!                                  |day
!!    sci(:)         |none          |retention coefficient for CN method based
!!                                  |on plant ET
!!    sedyld(:)      |metric tons   |soil loss for day in HRU
!!    smx(:)         |none          |retention coefficient for CN method based
!!                                  |on soil moisture
!!    surfq(:)       |mm H2O        |surface runoff generated on day in HRU
!!    tmn(:)         |deg C         |minimum temperature for the day in HRU
!!    tmpav(:)       |deg C         |average temperature for the day in HRU
!!    tmx(:)         |deg C         |maximum temperature for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    albday      |none          |albedo, the fraction of the solar radiation
!!                               |reflected at the soil surface back into
!!                               |space
!!    etday       |mm H2O        |actual evapotranspiration occuring on day
!!                               |in HRU
!!    ihru        |none          |HRU number
!!    inflpcp     |mm H2O        |amount of precipitation that infiltrates
!!                               |into soil (enters soil)
!!    nafert(:)   |none          |sequence number of auto-fert application
!!                               |within the year
!!    nair(:)     |none          |sequence number of auto-irrigation
!!                               |application within the year
!!    qdfr        |none          |fraction of water yield that is surface
!!                               |runoff
!!    qdr(:)      |mm H2O        |total amount of water entering main channel
!!                               |for day from HRU
!!    sci(:)      |none          |retention coefficient for CN method based
!!                               |on plant ET
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i_wtrhru
!!    ihout1
!!    iihru
!!    iru_sub
!!    j           |none          |HRU number
!!    ovs
!!    ovsl
!!    sumdaru
!!    sumk
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max, Exp, Dmin1
!!    SWAT: sub_subbasin, varinit, water_hru, schedule_ops, albedo, solt,
!!          surface, operatn, autoirr, percmain, etpot, etact, wattable,
!!          confert, conapply, graze, plantmod, dormant, nminrl, carbon,
!!          carbon_zhang2, nitvol, pminrl, pminrl2, biozone, gwmod, gwmod_deep,
!!          washp, decay, pestlch, enrsb, pesty, orgn, orgncswat, orgncswat2,
!!          psed, nrain, nlch, solp, bacteria, urban, urbanhr, latsed, gwnutr,
!!          gw_no3, surfstor, substor, filter, buffer, filtw, grass_wway,
!!          wetland, hrupond, hrupondhr, pothole, urb_bmp, watuse, watbal,
!!          subwq, sumv, virtual, routeunit, sumhyd, routels, addh

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i
   integer, parameter :: iru_sub = 1 ! route across landscape unit
   real*8 :: ovs, ovsl, sumdaru, sumk, xx
   integer :: i_wtrhru, ihout1, j, k

   ihru = hru1(inum1)

   call sub_subbasin(ihru)

   do k = 1, hrutot(inum1)

      j = ihru


      !!    deptil(:)   |mm  |depth of mixing caused by tillage operation
      !jj is hru number
      if (cswat == 2) then
         if (tillage_switch(j) .eq. 1) then
            if (tillage_days(j) .ge. 30) then
               tillage_switch(j) = 0
               tillage_days(j) = 0
            else
               tillage_days(j) = tillage_days(j) + 1
            end if
            !tillage_depth(j) = dtil
            !tillage_switch(j) = .TRUE.
         end if
      end if



      call varinit(j)
      if (icr(j) <= 0) icr(j) = 1

      i_wtrhru = 0
      idplrot(icr(j),j) = idplt(j)
      if (idplt(j) /= 0) then
         if (cpnm(idplt(j)) == "WATR") then
            i_wtrhru = 1
         end if
      endif

      if (i_wtrhru == 1) then
         call water_hru(j)
      else

         !! Simulate land covers other than water

         !! update base zero total heat units
         if (tmpav(j) > 0. .and. phutot(hru_sub(j)) > 0.01) then
            phubase(j) = phubase(j) + tmpav(j) / phutot(hru_sub(j))
         end if

         call schedule_ops(j)

         !! calculate albedo for day
         call albedo(j)

         !! calculate soil temperature for soil layers
         call solt(j)

!       if (ipot(j) /= j .and. imp_trig(nro(j),nrelease(j),j)==1)       &  Srini pothole
!
!     &        then
         !! calculate surface runoff if HRU is not impounded or an
         !! undrained depression--
         call surface(i,j)

         !! add surface flow that was routed across the landscape on the previous day
         !!   qday = qday + surfq_ru(j)
         !!   surfq_ru(j) = 0.

         !! compute effective rainfall (amount that percs into soil)
         inflpcp = Max(0.,precipday - surfq(j))
!        end if

         !! perform management operations
         if (yr_skip(j) == 0) call operatn(j)

         if (auto_wstr(j) > 1.e-6 .and. irrsc(j) > 2) call autoirr(j)

         !! perform soil water routing
         call percmain(j)

         !! compute evapotranspiration
         call etpot(j)
         call etact(j)

         !! compute water table depth using climate drivers
         call wattable(j)

         !! new CN method
         if (icn == 1) then
            sci(j) = sci(j) + pet_day*Exp(-cncoef_sub(hru_sub(j))*sci(j)/&
               &smx(j)) - precipday + qday + qtile + latq(j) + sepbtm(j)
         else if (icn == 2) then
            sci(j) = sci(j) + pet_day*Exp(-cncoef_sub(hru_sub(j))*sci(j)/&
               &smx(j)) - precipday + qday + latq(j) + sepbtm(j) + qtile
            sci(j) = Dmin1(sci(j),smxco * smx(j))
         end if

         !! apply fertilizer/manure in continuous fert operation
         if (icfrt(j) == 1) then
            ndcfrt(j) = ndcfrt(j) + 1
            call confert(j)
         end if

         !! apply pesticide in continuous pest operation
         if (icpst(j) == 1) then
            ndcpst(j) = ndcpst(j) + 1
            call conapply(j)
         end if

         !! remove biomass from grazing and apply manure
         if (igrz(j) == 1) then
            ndeat(j) = ndeat(j) + 1
            call graze(j)
         end if

         !! compute crop growth
         call plantmod(j)

         !! check for dormancy
         if (igro(j) == 1) call dormant(j)
         !! compute actual ET for day in HRU
         etday = ep_day + es_day + canev

         !! write daily air and soil temperature file
         !! can be uncommmented if needed by user and also in readfile.f

!      write (120,12112) i,j,tmx(j),tmn(j),(sol_tmp(k,j),k=1,sol_nly(j))
!12112  format (2i4,12f8.2)

         !! compute nitrogen and phosphorus mineralization
         select case (cswat)
          case (0)
            call nminrl(j)
          case (1)
            call carbon(i, j)
          case (2)
            call carbon_zhang2(j)
         end select

         call nitvol(j)
         if (sol_P_model == 1) then
            call pminrl(j)
         else
            call pminrl2(j)
         end if

!!    compute biozone processes in septic HRUs
!!    if 1)current is septic hru and 2)  soil temperature is above zero
         if (isep_opt(j) /= 0 .and. iyr >= isep_iyr(j)) then
            if (sol_tmp(i_sep(j),j) > 0.) call biozone(j)
         endif

         !! compute ground water contribution
         call gwmod(j)
         call gwmod_deep(j)

         if (hrupest(j) /= 0) then
            !! compute pesticide washoff
            if (precipday >= 2.54) call washp(j)
            !! compute pesticide degradation
            call decay(j)
            !! compute pesticide movement in soil
            call pestlch(j)
         end if

         if (surfq(j) > 0. .and. peakr > 1.e-6) then
            if (precipday > 0.) then
               call enrsb(0, j)
               if (hrupest(j) /= 0 .and. sedyld(j) > 0.) call pesty(0, j)

               select case (cswat)
                case (0) 
                  call orgn(0, j)
                case (1)
                  call orgncswat(0, j)
                case (2)
                  call orgncswat2(0, j)
               end select

               call psed(0, j)
            end if
         end if

         !! add nitrate in rainfall to soil profile
         call nrain(j)

         !! compute nitrate movement leaching
         call nlch(j)

         !! compute phosphorus movement
         call solp(j)

         !! compute bacteria transport
         call bacteria(j)

         !! compute loadings from urban areas
         if (urblu(j) > 0) then
            if(ievent == 0) then
               call urban(j) ! daily simulation
            else
               call urbanhr(j) ! subdaily simulation J.Jeong 4/20/2009
            endif
         endif

         !! compute sediment loading in lateral flow and add to sedyld
         call latsed

         !! compute nutrient loading in groundwater flow
         call gwnutr
         call gw_no3

         !! lag nutrients and sediment in surface runoff
         call surfstor(j)

         !! lag subsurface flow and nitrate in subsurface flow
         call substor(j)

         !! compute reduction in pollutants due to edge-of-field filter strip
         if (vfsi(j) >0.)then
            call filter(i)
            if (filterw(j) > 0.) call buffer
         end if
         if (vfsi(j) == 0. .and. filterw(j) > 0.) then
            call filtw
            call buffer
         end if

         !! compute reduction in pollutants due to in field grass waterway
         if (grwat_i(j) == 1) then
            call grass_wway(j)
         end if

         !! compute water yield for HRU
         qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j)
         if (qdr(j) < 0.) qdr(j) = 0.
         if (qdr(j) > 0.) then
            qdfr = qday / qdr(j)
         else
            qdfr = 0.
         end if

         !! compute wetland processes
         call wetlan(j)

         !! compute pond processes
         if (ievent == 0) then
            call hrupond
         else
            call hrupondhr
         endif

!       Srini pothole
         if (pot_fr(j) > 0.) call pothole(i, j)

         xx = sed_con(j)+soln_con(j)+solp_con(j)+orgn_con(j)+orgp_con(j)
         if (xx > 1.e-6) then
            call urb_bmp(j)
         end if

         !! consumptive water use (ponds, shallow aquifer, deep aquifer)
         call watuse(j)

         !! perform water balance
         call watbal(j)

         !! compute chl-a, CBOD and dissolved oxygen loadings
         call subwq(j)

         !! qdayout is surface runoff leaving the hru - after wetlands, ponds, and potholes
         qdayout(j) = qday

      endif

      !! perform output summarization
      call sumv(j)

      !! summarize output for multiple HRUs per subbasin
      !! store reach loadings for new fig method
      call virtual(i,j,k)
      aird(j) = 0.

      ihru = ihru + 1
   end do

   !! route 2 landscape units
   if (ils2flag(inum1) > 0) then
      isub = inum1                        ! save the subbasin number

      !! calculate outputs from hillslope
      ihout1 = mhyd_bsn + (inum1 - 1) * 4 ! first outflow hyd number
      ihout = ihout1                      ! outflow hyd number
      inum1 = 1                           ! landscape unit number
      inum2 = isub                        ! subbasin number
      call routeunit                      ! hillslope unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout

      !! calculate outputs from valley bottom
      inum1 = 2                           ! landscape unit number
      ihout = ihout + 1                   ! outflow hyd number
      sumdaru = 0.
      do j = 1, hrutot(isub)
         sumdaru = sumdaru + hru_km(j)
      end do
      daru_km(inum2,inum1) = sumdaru
      call routeunit                      ! valley bottom unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout

      !! route output from hillslope across valley bottom
      ihout = ihout + 1                   ! outflow hyd number
      inum1 = 2                           ! valley bottom landscape unit
      inum2 = ihout1                      ! inflow hyd=outlfow from hillslope
      inum3 = isub                        ! subbasin number
      rnum1 = 1.                          ! fraction overland flow
      !! compute weighted K factor for sediment transport capacity
      sumk = 0.
      ovsl = 0.
      ovs = 0.
      do j = 1, hrutot(isub)
         sumk = sumk + usle_k(j) * hru_rufr(inum1,j)
         ovsl = ovsl + slsubbsn(j)
         ovs = ovs + hru_slp(j)
      end do
      ovsl = ovsl / hrutot(isub)
      ovs = ovs / hrutot(isub)
      ru_k(isub,inum1) = sumk
      ru_ovsl(isub,inum1) = ovsl
      ru_ovs(isub,inum1) = ovs
      ru_ktc(isub,inum1) = 50.
      ru_a(isub,inum1) = daru_km(isub,1) / ru_ovsl(isub,inum1)
      call routels(iru_sub)               ! route across valley bottom
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      inum3s(ihout) = inum3
      ihouts(ihout) = ihout

      !! add routed with valley bottom loading
      inum1 = ihout                       ! hyd from routed
      inum2 = ihout - 1                   ! hyd from loading
      ihout = ihout + 1                   ! outflow hyd number
      call addh                           ! add hyd's
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout

      !! save landscape routed output in place of subbasin output for routing
      varoute(isub,:) = varoute(ihout,:)
   end if

   return
end
