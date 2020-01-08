!> @file clicon.f90
!> file containing the subroutine clicon
!> @author
!> modified by Javier Burguete

!> this subroutine controls weather inputs to SWAT. Precipitation and
!> temperature data is read in and the weather generator is called to
!> fill in radiation, wind speed and relative humidity as well as
!> missing precipitation and temperatures. Adjustments for climate
!> changes studies are also made in this subroutine.
!> @param[in] i current day of simulation (julian date)
subroutine clicon(i)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |julian date   |current day of simulation
!!    elevb(:,:)  |m             |elevation at center of band
!!    elevb_fr(:,:)|none         |fraction of subbasin area within elevation
!!                               |band
!!    elevp(:)    |m             |elevation of precipitation gage station
!!    elevt(:)    |m             |elevation of temperature gage station
!!    hru_sub(:)  |none          |subbasin in which HRU is located
!!    huminc(:,:) |none          |monthly humidity adjustment. Daily values
!!                               |for relative humidity within the month are
!!                               |raised or lowered by the specified amount.
!!                               |(used in climate change studies)
!!    id1         |julian date   |first day of simulation in year
!!    ifirstpet   |none          |potential ET data search code
!!                               |0 first day of potential ET data located in
!!                               |  file
!!                               |1 first day of potential ET data not located
!!                               |  in file
!!    ipet        |none          |code for potential ET method
!!                               |0 Priestley-Taylor method
!!                               |1 Penman/Monteith method
!!                               |2 Hargreaves method
!!                               |3 read in daily potential ET values
!!    irgage(:)   |none          |HRU rain gage data code (gage # for rainfall
!!                               |data used in HRU)
!!    itgage(:)   |none          |HRU temperature gage data code (gage # for
!!                               |temperature data used in HRU)
!!    iyr         |year          |year currently being simulated (eg 1980)
!!    i_mo        |none          |current month of simulation
!!    nhru        |none          |number of HRUs in watershed
!!    nstep       |none          |number of lines of rainfall data for each
!!                               |day
!!    pcpsim      |none          |rainfall input code
!!                               |1 gage read for each subbasin
!!                               |2 gage simulated for each subbasin
!!    plaps(:)    |mm H2O/km     |precipitation lapse rate: precipitation
!!                               |increase due to increase in elevation
!!    radinc(:,:) |MJ/m^2        |monthly solar radiation adjustment. Daily
!!                               |radiation within the month is raised or
!!                               |lowered by the specified amount. (used in
!!                               |climate change studies)
!!    rfinc(:,:)  |%             |monthly rainfall adjustment. Daily rainfall
!!                               |within the month is adjusted to the specified
!!                               |percentage of the original value (used in
!!                               |climate change studies)
!!    rhsim       |none          |relative humidity input code
!!                               |1 measured data read for each subbasin
!!                               |2 data simulated for each subbasin
!!    slrsim      |none          |solar radiation input code
!!                               |1 measured data read for each subbasin
!!                               |2 data simulated for each subbasin
!!    tlaps(:)    |deg C/km      |temperature lapse rate: temperature increase
!!                               |due to increase in elevation
!!    tmpinc(:,:) |deg C         |monthly temperature adjustment. Daily maximum
!!                               |and minimum temperatures within the month are
!!                               |raised or lowered by the specified amount
!!                               |(used in climate change studies)
!!    tmpsim      |none          |temperature input code
!!                               |1 daily max/min read for each subbasin
!!                               |2 daily max/min simulated for each subbasin
!!    welev(:)    |m             |elevation of weather station used to compile
!!                               |weather generator data
!!    wndsim      |none          |wind speed input code
!!                               |1 measured data read for each subbasin
!!                               |2 data simulated for each subbasin
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    frad(:,:)   |none          |fraction of solar radiation occuring during
!!                               |hour in day in HRU
!!    hru_ra(:)   |MJ/m^2        |solar radiation for the day in HRU
!!    hru_rmx(:)  |MJ/m^2        |maximum solar radiation for the day in HRU
!!    ifirstpet   |none          |potential ET data search code
!!                               |0 first day of potential ET data located in
!!                               |  file
!!                               |1 first day of potential ET data not located
!!                               |  in file
!!    pcpband(:,:)|mm H2O        |precipitation for the day in band in HRU
!!    petmeas     |mm H2O        |potential ET value read in for day
!!    rainsub(:,:)|mm H2O        |precipitation for the time step during the
!!                               |day in HRU
!!    rhd(:)      |none          |relative humidity for the day in HRU
!!    subp(:)     |mm H2O        |precipitation for the day in HRU
!!    tavband(:,:)|deg C         |average temperature for the day in band in HRU
!!    tmn(:)      |deg C         |minimum temperature for the day in HRU
!!    tmpav(:)    |deg C         |average temperature for the day in HRU
!!    tmx(:)      |deg C         |maximum temperature for the day in HRU
!!    tmxband(:,:)|deg C         |maximum temperature for the day in band in HRU
!!    u10(:)      |m/s           |wind speed for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    daylbsb
!!    fradbsb(:)  |none          |hourly solar radiation fractions for subbasin
!!    ib          |none          |counter
!!    idap        |julain date   |day currently being simulated
!!    ii          |none          |counter
!!    inum3sprev  |none          |subbasin number of previous HRU
!!    iyp         |none          |year currently being simulated
!!    k           |none          |counter
!!    l           |none          |counter
!!    npcpbsb
!!    pdif        |mm H2O        |difference in precipitation for station and
!!                               |precipitation for elevation band
!!    rabsb       |MJ/m^2        |generated solar radiation for subbasin
!!    ratio       |none          |fraction change in precipitation due to
!!                               |elevation changes
!!    rbsb        |mm H2O        |generated precipitation for subbasin
!!    rhdbsb      |none          |generated relative humidity for subbasin
!!    rmxbsb      |MJ/m^2        |generated maximum solar radiation for subbasin
!!    sb          |none          |subbasin number
!!    tdif        |deg C         |difference in temperature for station and
!!                               |temperature for elevation band
!!    tmnband(:)  |deg C         |minimum temperature for the day in band in HRU
!!    tmnbsb      |deg C         |generated minimum temperature for subbasin
!!    tmxbsb      |deg C         |generated maximum temperature for subbasin
!!    u10bsb      |m/s           |generated wind speed for subbasin
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max, Min
!!    SWAT: pmeas, tmeas, smeas, hmeas, wmeas
!!    SWAT: pgen, tgen, weatgn, clgen, slrgen, rhgen, wndgen

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i
   real*8, dimension(nstep) :: fradbsb
   real*8, dimension(10) :: tmnband
   real*8 :: daylbsb, pdif, rabsb, ratio, rbsb, rhdbsb, rmxbsb, tdif, tmnbsb,&
      &tmxbsb, u10bsb
   integer :: ib, idap, ii, inum3sprev, iyp, k, l, npcpbsb, sb

   sb = hru_sub(k)

!! Precipitation: Measured !!
   if (pcpsim == 1) call pmeas(i)

!! Temperature: Measured !!
   if (tmpsim == 1) call tmeas

!! Solar Radiation: Measured !!
   if (slrsim == 1) call smeas

!! Relative Humidity: Measured !!
   if (rhsim == 1) call hmeas

!! Wind Speed: Measured !!
   if (wndsim == 1 .and. ipet == 1) call wmeas

!! Potential ET: Read in data !!
   if (ipet == 3) then
      if (ifirstpet == 0) then
         read (140,5100) petmeas
      else
         ifirstpet = 0
         do
            iyp = 0
            idap = 0
            read (140,5000) iyp, idap, petmeas
            if (iyp == iyr .and. idap == id1) exit
         end do
      end if
   end if

!! Generate Relative Humidity, Wind Speed, Radiation
!! Precipitation and Temperature if = 2
   inum3sprev = 0
   do k = 1, nhru
      tmxbsb = 0.
      tmnbsb = 0.
      rbsb = 0.
      rhdbsb = 0.
      rabsb = 0.
      rmxbsb = 0.
      daylbsb = 0.
      npcpbsb = 0
      u10bsb = 0.
      fradbsb = 0.
      !! use same generated data for all HRUs in a subbasin
      if (sb == inum3sprev .and. sb /= 0) then
         if (tmpsim == 2) then
            tmx(k) = tmxbsb
            tmn(k) = tmnbsb
         end if
         if (pcpsim == 2) then
            subp(k) = rbsb
            if (ievent > 0) then
               do l = 1, nstep
                  rainsub(k,l) = rstpbsb(l)
               end do
            end if
         end if
         if (rhsim == 2) rhd(k) = rhdbsb
         if (slrsim == 2) then
            hru_ra(k) = rabsb
            hru_rmx(k) = rmxbsb
            dayl(k) = daylbsb
            npcp(k) = npcpbsb
            do ii = 1, nstep
               frad(k,ii) = fradbsb(ii)
            end do
         end if
         if (wndsim == 2) u10(k) = u10bsb
      else
         if (pcpsim == 2) call pgen(k)
         if (tmpsim == 2) then
            call weatgn(k)
            call tgen(k)
         end if
         if (slrsim == 2) then
            call clgen(k)
            call slrgen(k)
         end if
         if (rhsim == 2) call rhgen(k)
         if (ipet == 1) then
            if (wndsim == 2) call wndgen(k)
         end if
         !! set subbasin generated values
         inum3sprev = sb
         tmxbsb = tmx(k)
         tmnbsb = tmn(k)
         rbsb = subp(k)
         if (ievent > 0) then
            rstpbsb = 0.
            do l = 1, nstep
               rstpbsb(l) = rainsub(k,l)
            end do
         end if
         rhdbsb = rhd(k)
         rabsb = hru_ra(k)
         rmxbsb = hru_rmx(k)
         daylbsb = dayl(k)
         npcpbsb = npcp(k)
         u10bsb = u10(k)
         do ii = 1, nstep
            fradbsb(ii) = frad(k,ii)
         end do
      end if
      tmpav(k) = (tmx(k) + tmn(k)) / 2.
   end do


!! Climate Change Adjustments !!
   do k = 1, nhru
      subp(k) = subp(k) * (1. + rfinc(i_mo,sb) / 100.)
      if (subp(k) < 0.) subp(k) = 0.
      if (nstep > 0) then
         do ii = 1, nstep
            rainsub(k,ii) = rainsub(k,ii) *&
               &(1. + rfinc(i_mo,sb) / 100.)
            if (rainsub(k,ii) < 0.) rainsub(k,ii) = 0.
         end do
      end if
      tmx(k) = tmx(k) + tmpinc(i_mo,sb)
      tmn(k) = tmn(k) + tmpinc(i_mo,sb)
      tmpav(k) = tmpav(k) + tmpinc(i_mo,sb)
      hru_ra(k) = hru_ra(k) + radinc(i_mo,sb)
      hru_ra(k) = Max(0.,hru_ra(k))
      rhd(k) = rhd(k) + huminc(i_mo,sb)
      rhd(k) = Max(0.01,rhd(k))
      rhd(k) = Min(0.99,rhd(k))
   end do

!! Elevation Adjustments !!
   do k = 1, nhru
      if (elevb(1,sb) > 0. .and.&
         &elevb_fr(1,sb) > 0.) then
      !! determine temperature and precipitation for individual bands
         ratio = 0.
         do ib = 1, 10
            if (elevb_fr(ib,sb) < 0.) exit
            if (tmpsim == 1) then
               tdif = (elevb(ib,sb) -&
                  &dfloat(elevt(itgage(sb)))) * tlaps(sb)&
                  &/ 1000.
            else
               tdif = (elevb(ib,sb) - welev(sb))&
                  &* tlaps(sb) / 1000.
            end if
            if (pcpsim == 1) then
               pdif = (elevb(ib,sb) -&
                  &dfloat(elevp(irgage(sb)))) * plaps(sb)&
                  &/ 1000.
            else
               pdif = (elevb(ib,sb) - welev(sb))&
                  &* plaps(sb) / 1000.
            end if
            tavband(ib,k) = tmpav(k) + tdif
            tmxband(ib,k) = tmx(k) + tdif
            tmnband(ib) = tmn(k) + tdif
            if (subp(k) > 0.01) then
               pcpband(ib,k) = subp(k) + pdif
               if (pcpband(ib,k) < 0.) pcpband(ib,k) = 0.
            end if
            ratio = ratio + pdif * elevb_fr(ib,sb)
         end do
         !! determine fraction change in precipitation for HRU
         if (subp(k) >= 0.01) then
            ratio = ratio / subp(k)
         else
            ratio = 0.
         end if
         !! determine new overall temperature and precipitation values
         !! for HRU
         tmpav(k) = 0.
         tmx(k) = 0.
         tmn(k) = 0.
         subp(k) = 0.
         do ib = 1, 10
            if (elevb_fr(ib,sb) < 0.) exit
            tmpav(k) = tmpav(k) + tavband(ib,k) * elevb_fr(ib,sb)
            tmx(k) = tmx(k) + tmxband(ib,k) * elevb_fr(ib,sb)
            tmn(k) = tmn(k) + tmnband(ib) * elevb_fr(ib,sb)
            subp(k) = subp(k) + pcpband(ib,k) * elevb_fr(ib,sb)
         end do
         if (nstep > 0) then
            do ii = 1, nstep
               if (rainsub(k,ii) > 0.01) then
                  rainsub(k,ii) = rainsub(k,ii) + ratio * rainsub(k,ii)
                  if (rainsub(k,ii) < 0.) rainsub(k,ii) = 0.
               end if
            end do
         end if
      end if
   end do

   return
5000 format (i4,i3,f5.1)
5100 format (7x,f5.1)
end
