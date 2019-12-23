!> @file readwgn.f90
!> file containing the subroutine readwgn
!> @author
!> modified by Javier Burguete

!> this subroutine reads the HRU weather generator parameters from the
!> .wgn file
!> @param[in] ii HRU number (none)
subroutine readwgn(ii)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    adj_pkr     |none          |peak rate adjustment factor in the subbasin
!!                               |Used in the MUSLE equation to account for
!!                               |impact of peak flow on erosion
!!    ffcb        |none          |initial soil water content expressed as a
!!                               |fraction of field capacity
!!    ii          |none          |HRU number
!!    idaf        |julian date   |beginning day of simulation
!!    idist       |none          |rainfall distribution code
!!                               |  0 for skewed normal dist
!!                               |  1 for mixed exponential distribution
!!    ndays(:)    |julian date   |julian date for last day of preceding
!!                               |month (where the array location is the
!!                               |number of the month) The dates are for
!!                               |leap years
!!    rndseed(:,:)|none          |random number generator seeds
!!    sub_lat(:)  |degrees       |latitude of HRU/subbasin
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    amp_r(:,:)  |none          |average fraction of total daily rainfall
!!                               |occuring in maximum half-hour period
!!                               |for month
!!    daylmn(:)   |hr            |shortest daylength occurring during the year
!!    dewpt(:,:)  |deg C         |average dew point temperature for the month
!!    dormhr(:)   |hour          |time threshold used to define dormant
!!                               |period for plant (when daylength is within
!!                               |the time specified by dl from the minimum
!!                               |daylength for the area, the plant will go
!!                               |dormant)
!!    ffc(:)      |none          |initial HRU soil water content
!!                               |expressed as fraction of field capacity
!!    irelh       |none          |irelh = 0 (dewpoint)
!!                               |      = 1 (relative humidity)
!!                               |note:  inputs > 1.0 (dewpoint)
!!                               |       inputs < 1.0 (relative hum)
!!    ireg(:)     |none          |precipitation category:
!!                               |  1 precipitation <= 508 mm/yr
!!                               |  2 precipitation > 508 and <= 1016 mm/yr
!!                               |  3 precipitation > 1016 mm/yr
!!    latcos(:)   |none          |Cos(Latitude)
!!    latsin(:)   |none          |Sin(Latitude)
!!    pcf(:,:)    |none          |normalization coefficient for precipitation
!!                               |generator
!!    pcp_stat(:,1,:)|mm/day     |average amount of precipitation falling in
!!                               |one day for the month
!!    pcp_stat(:,2,:)|mm/day     |standard deviation for the average daily
!!                               |precipitation
!!    pcp_stat(:,3,:)|none       |skew coefficient for the average daily
!!                               |precipitation
!!    phutot(:)   |heat unit     |total potential heat units for year (used
!!                               |when no crop is growing)
!!    pr_w1(:,:)  |none          |probability of wet day after dry day in month
!!    pr_w2(:,:)  |none          |probability of wet day after wet day in month
!!    pr_w3(:,:)  |none          |proportion of wet days in the month
!!    solarav(:,:)|MJ/m^2/day    |average daily solar radiation for the month
!!    tmp_an(:)   |deg C         |average annual air temperature
!!    tmpmn(:,:)  |deg C         |avg monthly minimum air temperature
!!    tmpmx(:,:)  |deg C         |avg monthly maximum air temperature
!!    tmpstdmn(:,:)|deg C        |standard deviation for avg monthly minimum air
!!                               |temperature
!!    tmpstdmx(:,:)|deg C        |standard deviation for avg monthly maximum air
!!                               |temperature
!!    welev(:)    |m             |elevation of weather station used to compile
!!                               |data
!!    wlat(:)     |degrees       |latitude of weather station used to compile
!!                               |data
!!    wndav(:,:) |m/s            |average wind speed for the month
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dl          |hour          |time threshold used to define dormant
!!                               |period for plant (when daylength is within
!!                               |the time specified by dl from the minimum
!!                               |daylength for the area, the plant will go
!!                               |dormant)
!!    j           |none          |counter
!!    k           |none          |counter
!!    lattan      |none          |Tan(Latitude)
!!    m1          |none          |array location (see definition of ndays)
!!    mdays       |none          |number of days in the month
!!    mon         |none          |monthly counter
!!    nda         |julian date   |julian date of last day in the month
!!    pcp         |mm H2O        |generated precipitation
!!    pcpd(:)     |days          |average number of days of precipitation
!!                               |in the month
!!    pcpmm(:)    |mm            |amount of precipitation in month
!!    r6          |none          |variable to hold calculation result
!!    rain_hhsm(:)|mm            |smoothed values for maximum 0.5 hour rainfall
!!    rain_yrs    |none          |number of years of recorded maximum 0.5h
!!                               |rainfall used to calculate values for
!!                               |rainhhmx(:)
!!    rainhhmx(:) |mm            |maximum 0.5 hour rainfall in month
!!                               |for entire period of record
!!    rndm1       |none          |random number between 0.0 and 1.0
!!    rnm2        |none          |random number between 0.0 and 1.0
!!    sffc
!!    sum         |none          |variable to hold summation results
!!    summm_p     |mm            |sum of precipitation over year
!!    summn_t     |deg C         |sum of mimimum temp values over year
!!    summx_t     |deg C         |sum of maximum temp values over year
!!    tav         |deg C         |average monthly temperature
!!    titldum     |NA            |title line of .wgn file (not used elsewhere)
!!    tmax        |deg C         |maximum average monthly temperature
!!    tmin        |deg C         |minimum average monthly temperature
!!    tmpsoil     |deg C         |initial temperature of soil layers
!!    x1          |none          |variable to hold calculation results
!!    x2          |none          |variable to hold calculation results
!!    x3          |none          |variable to hold calculation results
!!    xlv         |none          |variable to hold calculation results
!!    xrnd
!!    xx          |varies        |variable to hold calculation results
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sin, Cos, Tan, abs, Acos, Log, Exp, MaxVal
!!    SWAT: Aunif, Dstn1

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: ii
   real*8, dimension (12) :: pcpd, pcpmm, rain_hhsm, rainhhmx
   character (len=80) :: titldum
   real*8 :: dl, lattan, pcp, r6, rain_yrs, rndm1, rnm2, sffc, sum, summm_p,&
      &summn_t, summx_t, tav, tmpsoil, x1, x2, x3, xlv, xx
   integer :: j, k, m1, mdays, mon, nda, xrnd


   pcpd = 0.
   rainhhmx = 0.
   pcpmm = 0.

   read (114,5000) titldum
   read (114,5100) wlat(ii)
   read (114,5100) welev(ii)
   read (114,5100) rain_yrs
   read (114,5200) (tmpmx(mon,ii),mon = 1,12)
   read (114,5200) (tmpmn(mon,ii),mon = 1,12)
   read (114,5200) (tmpstdmx(mon,ii),mon = 1,12)
   read (114,5200) (tmpstdmn(mon,ii),mon = 1,12)
   read (114,5201) (pcpmm(mon),mon = 1,12)
   read (114,5200) (pcp_stat(mon,2,ii),mon = 1,12)  !pcpstd
   read (114,5200) (pcp_stat(mon,3,ii),mon = 1,12)  !pcpskw
   read (114,5200) (pr_w1(mon,ii),mon = 1,12)
   read (114,5200) (pr_w2(mon,ii),mon = 1,12)
   read (114,5200) (pcpd(mon),mon = 1,12)
   read (114,5200) (rainhhmx(mon),mon = 1,12)
   read (114,5200) (solarav(mon,ii),mon = 1,12)
   read (114,5200) (dewpt(mon,ii),mon = 1,12)
   read (114,5200) (wndav(mon,ii),mon = 1,12)

!! determine if input for dewpt is relative humidity
   do mon = 1,12
      if (dewpt(mon,ii) > 1.0  .or. dewpt(mon,ii) < 0.0) then
         irelh(ii) = 0
      end if
   end do

!! variables needed for radiation calcs.
   x2 = 0.0
   if (sub_lat(ii) < 1.e-4) sub_lat(ii) = wlat(ii)
   xx = sub_lat(ii) / 57.296
   !!convert degrees to radians (2pi/360=1/57.296)
   latsin(ii) = Sin(xx)
   latcos(ii) = Cos(xx)
   lattan = Tan(xx)
!! calculate minimum daylength
!! daylength=2*acos(-tan(sd)*tan(lat))/omega
!! where solar declination, sd, = -23.5 degrees for minimum daylength in
!!                      northern hemisphere and -tan(sd) = .4348
!!       dabsolute value is taken of tan(lat) to convert southern hemisphere
!!                      values to northern hemisphere
!!       the angular velocity of the earth's rotation, omega, = 15 deg/hr or
!!                      0.2618 rad/hr and 2/0.2618 = 7.6394
   x1 = .4348 * abs(lattan)
   if (x1 < 1.) x2 = Acos(x1)
   !!x1 will be >= 1. if sub_lat > 66.5 or < -66.5
   daylmn(ii) = 7.6394 * x2

!! calculate day length threshold for dormancy
   if (dorm_hr < -1.e-6) then
      if (abs(sub_lat(ii)) > 40.) then
         dl = 1.
      else if (abs(sub_lat(ii)) < 20.) then
         dl = -1.
      else
         dl = (abs(sub_lat(ii)) - 20.) / 20.
      end if
   else
      dl = dorm_hr
   end if



!! calculate smoothed maximum 0.5hr rainfall amounts
   rain_hhsm = 0.
   rain_hhsm(1) = (rainhhmx(12) + rainhhmx(1) + rainhhmx(2)) / 3.
   do mon = 2, 11
      rain_hhsm(mon) = (rainhhmx(mon-1) + rainhhmx(mon) + rainhhmx(mon+1)) / 3.
   end do
   rain_hhsm(12) = (rainhhmx(11) + rainhhmx(12) + rainhhmx(1)) / 3.


!! calculate missing values and additional parameters
   summx_t = 0.
   summn_t = 0.
   summm_p = 0.
   !tmin = 100. ! not used
   !tmax = 0. ! not used
   pcpdays(ii) = 0.
   do mon = 1, 12
      mdays = ndays(mon+1) - ndays(mon)
      tav = (tmpmx(mon,ii) + tmpmn(mon,ii)) / 2.
      !if (tav > tmax) tmax = tav ! not used
      !if (tav < tmin) tmin = tav ! not used
      summx_t = summx_t + tmpmx(mon,ii)
      summn_t = summn_t + tmpmn(mon,ii)

      !! calculate total potential heat units
      if (tav > 0.) phutot(ii) = phutot(ii) + tav * mdays

      !! calculate values for pr_w if missing or bad
      if (pr_w2(mon,ii) <= pr_w1(mon,ii).or.pr_w1(mon,ii) <= 0.) then
         if (pcpd(mon) < .1) pcpd(mon) = 0.1
         pr_w1(mon,ii) = .75 * pcpd(mon) / mdays
         pr_w2(mon,ii) = .25 + pr_w1(mon,ii)
      else
         !! if pr_w values good, use calculated pcpd based on these values
         !! using first order Markov chain
         pcpd(mon) = mdays * pr_w1(mon,ii) /&
         &(1. - pr_w2(mon,ii) + pr_w1(mon,ii))

      end if

      !! calculate precipitation-related values
      if (pcpd(mon) <= 0.) pcpd(mon) = .001
      pr_w3(mon,ii) = pcpd(mon) / mdays
      pcp_stat(mon,1,ii) = pcpmm(mon) / pcpd(mon)
      if (pcp_stat(mon,3,ii) < 0.2) pcp_stat(mon,3,ii) = 0.2
      summm_p = summm_p + pcpmm(mon)
      pcpdays(ii) = pcpdays(ii) + pcpd(mon)
   end do

   tmp_an(ii) = (summx_t + summn_t) / 24.

   !! calculate initial temperature of soil layers
   if (idaf > ndays(2)) then
      do mon = 2, 12
         m1 = mon + 1
         nda = ndays(m1) - 1
         if (idaf <= nda) exit
      end do
   else
      mon = 1
   end if
   tmpsoil = (tmpmx(mon,ii) + tmpmn(mon,ii)) / 2.

   xrnd = rndseed(idg(3),ii)
   rndm1 = Aunif(xrnd)
   do mon = 1, 12
      !! calculate precipitation correction factor for pcp generator
      if (idist == 0) then
         rnm2 = 0.
         xlv = 0.
         pcp = 0.
         sum = 0.
         r6 = pcp_stat(mon,3,ii) / 6.
         do j = 1, 1000
            rnm2 = Aunif(xrnd)
            xlv = (Dstn1(rndm1,rnm2) -r6) * r6 + 1
            rndm1 = rnm2
            xlv = (xlv**3 - 1.) * 2 / pcp_stat(mon,3,ii)
            pcp = xlv * pcp_stat(mon,2,ii) + pcp_stat(mon,1,ii)
            if (pcp < 0.01) pcp = 0.01
            sum = sum + pcp
         end do
         if (sum > 0.) then
            pcf(mon,ii) = 1000. * pcp_stat(mon,1,ii) / sum
         else
            pcf(mon,ii) = 1.
         end if
      end if

      !! calculate or estimate amp_r values
      if (rain_yrs < 1.0) rain_yrs = 10.
      x1 = .5 / rain_yrs
      x2 = x1 / pcpd(mon)
      x3 = rain_hhsm(mon) / Log(x2)
      if (pcp_stat(mon,1,ii) > 1.e-4) then
         amp_r(mon,ii) = adj_pkr * (1. - Exp(x3 / pcp_stat(mon,1,ii)))
      else
         amp_r(mon,ii) = 0.95
      end if
      if (amp_r(mon,ii) < .1) amp_r(mon,ii) = .1
      if (amp_r(mon,ii) > .95) amp_r(mon,ii) = .95
   end do

!!    determine precipitation category (ireg initialized to category 1)
   xx = summm_p
   if (summm_p > 508.) ireg(ii) = 2
   if (summm_p > 1016.) ireg(ii) = 3

!!    set fraction of field capacity in soil
   if (ffcb <= 0.) then
      sffc = summm_p / (summm_p + Exp(9.043 - 0.002135*summm_p))
      !!S-curve equation Jeff made up.
   else
      sffc = ffcb
   end if

!! assign HRU values
   do j = 1, hrutot(ii)
      ihru = nhru + j
      do k = 1, sol_nly(ihru)
         sol_tmp(k,ihru) = tmpsoil
      end do
      ffc(ihru) = sffc
      dormhr(ihru) = dl
   end do

   close (114)
   return
5000 format (a)
5100 format (12x,f7.2)
5200 format (12f6.2)
5201 format (12f6.1)
end
