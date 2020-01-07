!> @file readsub.f90
!> file containing the subroutine readsub
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU/subbasin general input file
!> (.sub). This file contains data related to general processes modeled
!> at the HRU/subbasin level.
!> @param[in] i subbasin number (none)
subroutine readsub(i)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |subbasin number
!!    da_km       |km2           |area of the watershed in square kilometers
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_k(1,:)    |mm/hr         |effective hydraulic conductivity of tributary
!!                               |channel alluvium
!!    ch_l(1,:)    |km            |longest tributary channel length in subbasin
!!    ch_n(1,:)    |none          |Manning's "n" value for the tributary channels
!!    ch_s(1,:)    |m/m           |average slope of tributary channels
!!    ch_w(1,:)    |m             |average width of tributary channels
!!    cncoef_sub  |              |soil water depletion coefficient used
!!                               | in the new (modified curve number method)
!!                               | same as soil index coeff used in APEX
!!                               | range: 0.5 - 2.0
!!    co2(:)      |ppmv          |CO2 concentration
!!    driftco(:)  |none          |coefficient for pesticide drift directly
!!                               |onto stream
!!    elevb(:,:)  |m             |elevation at the center of the band
!!    elevb_fr(:,:)|none         |fraction of subbasin area within elevation
!!                               |band (the same fractions should be listed for
!!                               |all HRUs within the subbasin)
!!    epco(:)     |none          |plant water uptake compensation factor (0-1)
!!    esco(:)     |none          |soil evaporation compensation factor (0-1)
!!    harg_petco(:)              |coefficient related to radiation used
!!                               | in hargreaves eq (range: 0.0019 - 0.0032)
!!    hru_km(:)   |km**2         |area of HRU in square kilometers
!!    huminc(:,:) |none          |monthly humidity adjustment. Daily values
!!                               |for relative humidity within the month are
!!                               |rasied or lowered by the specified amount.
!!                               |(used in climate change studies)
!!    ifld(:)     |none          |number of HRU (in subbasin) that is a
!!                               |floodplain
!!    ihgage(:)   |none          |subbasin relative humidity data code
!!    ipot(:)     |none          |number of HRU (in subbasin) that is ponding
!!                               |water--the HRU that the surface runoff from
!!                               |current HRU drains into. This variable is
!!                               |used only for rice paddys or closed
!!                               |depressional areas
!!    irgage(:)   |none          |subbasin rain gage data code
!!    irip(:)     |none          |number of HRU (in subbasin) that is a
!!                               |riparian zone
!!    isgage(:)   |none          |subbasin radiation gage data code
!!    itgage(:)   |none          |subbasin temp gage data code
!!    iwgage(:)   |none          |subbasin wind speed gage data code
!!    plaps(:)    |mm H2O/km     |precipitation lapse rate: precipitation
!!                               |change due to change in elevation
!!    radinc(:,:) |MJ/m^2        |monthly solar radiation adjustment. Daily
!!                               |radiation within the month is raised or
!!                               |lowered by the specified amount. (used in
!!                               |climate change studies)
!!    rfinc(:,:)  |%             |monthly rainfall adjustment. Daily rainfall
!!                               |within the month is adjusted to the specified
!!                               |percentage of the original value (used in
!!                               |climate change studies)
!!    sub_smfmn(:)|mm/deg C/day  |min melt rate for snow during year (Dec 21)
!!                               | for subbasin(:)
!!                               | (range: -5.0/5.0)
!!    sub_smfmx(:)|mm/deg C/day  |max melt rate for snow during year (June 21)
!!                               | for subbasin(:)
!!                               | where deg C refers to the air temperature
!!                               | SUB_SMFMX and SMFMN allow the rate of snow
!!                               | melt to vary through the year.  These
!!                               | parameters are accounting for the impact
!!                               | of soil temperature on snow melt.
!!                               | (range: -5.0/5.0)
!!                               | where deg C refers to the air temperature
!!    sub_sftmp(:)|deg C         |Snowfall temperature for subbasin(:)
!!                               | Mean air temperature at which precip is
!!                               | equally likely to be rain as snow/freezing
!!                               | rain (range: -5.0/5.0)
!!    sub_smtmp(:)|deg C         |Snow melt base temperature for subbasin(:)
!!                               | mean air temperature at which snow melt will
!!                               | occur (range: -5.0/5.0)
!!    sno_hru(:)  |mm H2O        |amount of water stored as snow
!!    snoeb(:,:)  |mm H2O        |initial snow water content in elevation band
!!    sub_elev(:) |m             |average elevation of subbasin
!!    sub_km(:)   |km**2         |area of subbasin in square kilometers
!!    sub_lat(:)  |degrees       |latitude of subbasin
!!    tlaps(:)    |deg C/km      |temperature lapse rate: temperature change
!!                               |due to change in elevation
!!    tmpinc(:,:) |deg C         |monthly temperature adjustment. Daily maximum
!!                               |and minimum temperatures within the month are
!!                               |raised or lowered by the specified amount
!!                               |(used in climate change studies)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_ls
!!    chmfile     |NA            |HRU soil chemical data
!!    eof         |none          |end of file flag (=-1 if eof, else =0)
!!    gwfile      |NA            |HRU groundwater data
!!    hrufile     |NA            |name of HRU general data file
!!    if          |none          |number of HRU in subbasin that is floodplain
!!    ils
!!    ip          |none          |number of HRU in subbasin that is pothole
!!    ir          |none          |number of HRU in subbasin that is riparian zone
!!    j           |none          |counter
!!    jj          |none          |variable to take special HRUs into account
!!    k
!!    kk
!!    l           |none          |HRU number
!!    mgtfile     |NA            |HRU management data
!!    mon         |none          |monthly counter
!!    opsfile     |NA            |name of operation schedule file for Phil G.
!!    pndfile     |NA            |subbasin impoundment file (.pnd)
!!    sdrfile     |NA            |name of subbasin drainage file
!!    septfile    |NA            |septic input data file (.sep)
!!    snofile
!!    solfile     |NA            |HRU soil data
!!    sumebfr     |none          |total area of subbasin modeled in
!!                               |elevation bands
!!    titldum     |NA            |title line of .sub file (not used)
!!    wgnfile     |NA            |name of weather generator data file
!!    wusfile     |NA            |subbasin water use file (.wus)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: caps, readsno, readsepticbz, readsdr, readhru, readchm, readmgt,
!!          readsol, readgw, readops, readwgn, readpnd, readwus

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: i
   character (len=80) :: snofile, titldum
   real*8, dimension(10) :: ssnoeb
   character (len=13) :: chmfile, gwfile, hrufile, mgtfile, opsfile, pndfile,&
      &sdrfile, septfile, solfile, wgnfile, wusfile
   real*8 :: ch_ls, sno_sub, sumebfr
   integer :: eof, ils, j, jj, k, kk, l, mon

   lu_nodrain = "    "
   wgnfile = ""
   pndfile = ""
   wusfile = ""
   eof = 0
   ssnoeb = 0.
   sno_sub = 0.
   ch_ls = 0.
   jj = 1
   !ip = 0
   !if = 0
   !ir = 0

   read (101,5100) titldum
   read (101,*) sub_km(i)
   if (isproj == 3) then
      read (101,5101) harg_petco(i), cncoef_sub(i), sub_smfmx(1,i),&
         &sub_smfmn(1,i), sub_sftmp(1,i), sub_smtmp(1,i), sub_timp(1,i)
      do j = 2, 10
         sub_smfmx(j,i) = sub_smfmx(1,i)
         sub_smfmn(j,i) = sub_smfmn(1,i)
         sub_sftmp(j,i) = sub_sftmp(1,i)
         sub_smtmp(j,i) = sub_smtmp(1,i)
         sub_timp(j,i) = sub_timp(1,i)
      end do
   else
      read (101,5100) titldum
   end if
   read (101,5100) titldum
   read (101,*) sub_lat(i)
   read (101,*) sub_elev(i)
   read (101,*) irgage(i)
   read (101,*) itgage(i)
   read (101,*) isgage(i)
   read (101,*) ihgage(i)
   read (101,*) iwgage(i)
   read (101,5300) wgnfile
   call caps(wgnfile)
   open (114,file=wgnfile)
   read (101,*) fcst_reg(i)
   read (101,5100) titldum
   read (101,5100) titldum
   read (101,5200) (elevb(j,i), j = 1, 10)
   read (101,5100) titldum
   read (101,5200) (elevb_fr(j,i), j = 1, 10)
   read (101,5100) titldum
   read (101,5200) (ssnoeb(j), j = 1, 10)
   read (101,*) plaps(i)
   read (101,*) tlaps(i)
   read (101,*) sno_sub
   read (101,5100) titldum
   read (101,*) ch_ls
   read (101,*) ch_s(1,i)
   read (101,*) ch_w(1,i)
   read (101,*) ch_k(1,i)
   read (101,*) ch_n(1,i)
   read (101,5100) titldum
   read (101,5300) pndfile
   call caps(pndfile)
   open (104,file=pndfile)
   read (101,5100) titldum
   read (101,5300) wusfile
   call caps(wusfile)
   open (105,file=wusfile)
   read (101,5100) snofile
   if (snofile /='             ' .or. snofile /= 'Climate Change')then
      if (snofile /='             ') then
         if (snofile /= 'Climate Change') then
            open (113,file=snofile)
            call caps (snofile)
            call readsno(i)
         endif
      endif
   endif
   read (101,*) co2(i)
   read (101,5100) titldum
   read (101,5200) (rfinc(i,mon),mon = 1,6)
   read (101,5100) titldum
   read (101,5200) (rfinc(i,mon),mon = 7,12)
   read (101,5100) titldum
   read (101,5200) (tmpinc(i,mon),mon = 1,6)
   read (101,5100) titldum
   read (101,5200) (tmpinc(i,mon),mon = 7,12)
   read (101,5100) titldum
   read (101,5200) (radinc(i,mon),mon = 1,6)
   read (101,5100) titldum
   read (101,5200) (radinc(i,mon),mon = 7,12)
   read (101,5100) titldum
   read (101,5200) (huminc(i,mon),mon = 1,6)
   read (101,5100) titldum
   read (101,5200) (huminc(i,mon),mon = 7,12)
   read (101,5100) titldum
!! read HRU input data
   read (101,*) hrutot(i)
   read (101,5100) titldum
   read (101,5100) titldum
   read (101,5100) titldum
   read (101,5100) titldum
   read (101,5100) titldum
   read (101,5100) titldum
   read (101,5100) titldum
   !!General HRUs
   read (101,5100) titldum
   do j = jj, hrutot(i)
      l = nhru + j
      if (j == 1) hru1(i) = l
!         ipot(l) = ip
!         ifld(l) = if
!         irip(l) = ir
      chmfile = ""
      hrufile = ""
      mgtfile = ""
      solfile = ""
      gwfile = ""
      opsfile = ""
      septfile = ""
      sdrfile = ""
      read (101,5300) hrufile, mgtfile, solfile, chmfile, gwfile,&
         &opsfile, septfile, sdrfile, ils2(l)
      call caps(hrufile)
      call caps(mgtfile)
      call caps(solfile)
      call caps(chmfile)
      call caps(gwfile)
      if (septfile /='             ') then
         call caps (septfile)
         open (172,file=septfile, status='old')
         isep_hru(l) = 1
         call readsepticbz(l)
      end if
      if (sdrfile /= '             ') then
         call caps(sdrfile)
         open (112,file=sdrfile)
         call readsdr(l)
         if (sdrain(l) <= 0.) sdrain(l) = sdrain_bsn   !! todd campbell 06/07/18
      else
         if (re(l) <= 0.) re(l) = re_bsn
         if (sdrain(l) <= 0.) sdrain(l) = sdrain_bsn
         if (drain_co(l) <= 0.) drain_co(l) = drain_co_bsn
         if (pc(l) <= 0.) pc(l) = pc_bsn
         if (latksatf(l) <= 0.) latksatf(l) = latksatf_bsn
         if (sstmaxd(l) <= 0.) sstmaxd(l) = sstmaxd_bsn
         !           sdrain(l) = 0.      !!!! nbs 11/25/15
      end if

      open (106,file=chmfile)
      open (107,file=solfile)
      open (108,file=hrufile)
      open (109,file=mgtfile)
      open (110,file=gwfile)
      call readhru(i, l)
      call readchm(l)
      call readmgt(l)
      call readsol(l)
      call readgw(i, l)
      if (opsfile /= '             ') then
         call caps(opsfile)
         open (111,file=opsfile)
         call readops(l)
      end if

      ! set up variables for landscape routing
!          if (ils_nofig == 1) then
      if (ils2(l) == 0) then
         ils = 1
      else
         ils = 2
         ils2flag(i) = 1
      end if
      daru_km(i,ils) = daru_km(i,ils) + hru_fr(j) * sub_km(i)
!         end if

      ! estimate drainage area for urban distributed bmps in hectares - jaehak
      if (urblu(l)>0) then
         kk=1
         bmpdrain(l) = 1
         do while(lu_nodrain(kk) .ne. "    ")
            if (urbname(urblu(l)) .eq. lu_nodrain(kk)) then
               bmpdrain(l) = 0
               exit
            end if
            kk = kk + 1
            if (kk>30) exit
         end do
         if(bmpdrain(l)==1) then
            sub_ha_imp(i) = sub_ha_imp(i) + hru_ha(l) * fimp(urblu(l))
            sub_ha_urb(i) = sub_ha_urb(i) + hru_ha(l)
         end if
      end if

      ! HRU selection criteria for Irrigation by retention-irrigation basins
      if (num_ri(i)>0) then
         if (sol_z(sol_nly(l),l)>300& !!    - soil thickness > 12 inches
            &.AND.sol_k(1,l)>0.76&       !    - permeability > 0.03 inches/hr (=0.76mm/hr)
            &.AND.hru_slp(l)<0.1&        !!    - hru slope < 10%
            &.AND.urblu(l)>0) then !urban LU

            ri_luflg(l) = 1 !irrigate HRU
         end if

         do kk=1,num_noirr(i)
            if (urbname(urblu(l)).eq.ri_nirr(i,kk)) then
               ri_luflg(l) = 0 !exclude these land uses from irrigation
            end if
         end do

         if (ri_luflg(l) == 1) then
            ri_subkm(i) = ri_subkm(i) + hru_km(l) * (1.-fimp(urblu(l))) !km2
         end if
      end if

      ! estimate impervious cover in the upstream drainage area for on-line bmps
      if (iurban(l) > 0) then
         subdr_ickm(i) = subdr_ickm(i) + hru_km(l) * fimp(urblu(l))
      end if
      ! estimate average Curve Number for the subbasin
      sub_cn2(i) = sub_cn2(i) + cn(2,l) * hru_fr(l)
   end do      ! hru loop

   !! set up routing unit fractions for landscape routing
   do j = jj, hrutot(i)
      l = nhru + j
      if (ils2(l) == 0) then
         ils = 1
      else
         ils = 2
      end if

   end do
   if (ils == 2) then
      do j = jj, hrutot(i)
         hru_rufr(ils,l) = hru_fr(l) * sub_km(i) / daru_km(i,ils)
      end do
   end if

!! commented the following statements and moved above in the 'else'
!! where it reads the sdrfile.  Jeff should check.
!!    set default values
!      do l = jj, hrutot(i)
!        if (re(l) <= 0.) re(l) = re_bsn
!   if (sdrain(l) <= 0.) sdrain(l) = sdrain_bsn
!   if (drain_co(l) <= 0.) drain_co(l) = drain_co_bsn
!   if (pc(l) <= 0.) pc(l) = pc_bsn
!        if (latksatf(l) <= 0.) latksatf(l) = latksatf_bsn
!   if (sstmaxd(l) <= 0.) sstmaxd(l) = sstmaxd_bsn
!     end do

   !     estimate drainage area for urban on-line bmps in square km
   !subdr_km(i) = subdr_km(i) + sub_km(i)


!!    set default values
   if (sub_km(i) <= 0.) sub_km(i) = 1.
   if (harg_petco(i) < 1.e-6) harg_petco(i) = .0023
   if (cncoef_sub(i) <= 1.e-6) then
      if (cncoef > 1.e-6) then
         cncoef_sub(i) = cncoef
      else
         cncoef_sub(i) = 1.
      endif
   endif

   if (fcst_reg(i) <= 0.) fcst_reg(i) = 1
   if (co2(i) <= 0.) co2(i) = 330.
   if (ch_s(1,i) <= 0.) ch_s(1,i) = .0001
   if (ch_n(1,i) <= 0.005) ch_n(1,i) = 0.005
   if (ch_n(1,i) >= 0.70) ch_n(1,i) = 0.70
   do j = 1, 10
      if (sub_smtmp(j,i) < 1.e-6) sub_smtmp(j,i) = smtmp
      if (sub_sftmp(j,i) < 1.e-6) sub_sftmp(j,i) = sftmp
      if (sub_smfmx(j,i) < 1.e-6) sub_smfmx(j,i) = smfmx
      if (sub_smfmn(j,i) < 1.e-6) sub_smfmn(j,i) = smfmn
      if (sub_timp(j,i) < 1.e-6) sub_timp(j,i) = timp
   end do

!!    check elevation band fractions
   sumebfr = 0.
   do j = 1, 10
      sumebfr = sumebfr + elevb_fr(j,i)
   end do
   if (sumebfr > 1.e-5) then
      if (sumebfr < .99) write (24,1000) i
   end if

!!    This equation given to us by EPA, in the process of getting reference
!   sdrift = .01 * (10.**(-.00738 * (7.62 * ch_w(1,i)) - 2.5889) + .2267) / 2.
!   not used

!! assign subbasin values to HRUs where needed
   do j = 1, hrutot(i)
      l = nhru + j
      hru_sub(l) = i
!!   hru_seq = sequential hru number within the subbasin
      hru_seq(l) = j
      hrugis(l) = subgis(i)
      do k = 1, 10
         snoeb(k,l) = ssnoeb(k)
      end do
      sno_hru(l) = sno_sub
      ch_l(1,l) = ch_ls
!      driftco(l) = sdrift ! not used
   end do

!! calculate watershed land area
   da_km = da_km + sub_km(i)

!!read in weather generator parameter values
   call readwgn(i)
   plaps(i) = plaps(i) / pcpdays(i)
!!read in subbasin impoundment parameter values
   call readpnd(i)
!!read in subbasin water use parameter values
   call readwus(i)

!! sediment delivery ration for the subbasin..... urban modeling by J.Jeong
   dratio(i) = 0.42 * sub_km(i) ** (-0.125)
   if (dratio(i) > 0.9) dratio(i) = 0.9

   close (101)
   return
1000 format ('ERROR: Elevation Band Fractions in Subbasin ',i4,&
   &' do not add up to 100% of subbasin area!')
5100 format (a)
5101 format (f8.4,f4.2,5f8.3)
5200 format (10f8.1)
5300 format (8a13,i6)
end
