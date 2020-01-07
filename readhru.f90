!> @file readhru.f90
!> file containing the subroutine readhru
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU general input file (.hru).
!> This file contains data related to general processes modeled
!> at the HRU level.
!> @param[in] i subbasin number (none)
!> @param[in] j HRU number (none)
subroutine readhru(i, j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |subbasin number
!!    j           |none          |HRU number
!!    ifld(:)     |none          |number of HRU (in subbasin) that is a
!!                               |floodplain
!!    ipot(:)     |none          |number of HRU (in subbasin) that is ponding
!!                               |water--the HRU that the surface runoff from
!!                               |current HRU drains into. This variable is
!!                               |used only for rice paddys or closed
!!                               |depressional areas
!!    irip(:)     |none          |number of HRU (in subbasin) that is a
!!                               |riparian zone
!!    da_km       |km2           |area of the watershed in square kilometers
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    canmx(:)    |mm H2O        |maximum canopy storage
!!    cf          |              |this parameter controls the response of
!!                               |decomposition to the combined effect of soil
!!                               |temperature and moisture.
!!    cfh         |              |Maximum humification rate``
!!    cfdec       |              |the undisturbed soil turnover rate under
!!                               |optimum soil water and temperature. Increasing
!!                               |it will increase carbon and organic N decomp.
!!    dis_stream(:) | m          |average distance to stream
!!    epco(:)     |none          |plant water uptake compensation factor (0-1)
!!    erorgn(:)   |none          |organic N enrichment ratio, if left blank
!!                               |the model will calculate for every event
!!    erorgp(:)   |none          |organic P enrichment ratio, if left blank
!!                               |the model will calculate for every event
!!    esco(:)     |none          |soil evaporation compensation factor (0-1)
!!    hru_fr(:)   |km2/km2       |fraction of subbasin area contained in HRU
!!    hru_ha(:)   |ha            |area of HRU in hectares
!!    hru_km(:)   |km**2         |area of HRU in square kilometers
!!    hru_slp(:)  |m/m           |average slope steepness
!!    lat_sed(:)  |g/L           |sediment concentration in lateral flow
!!    lat_ttime(:)|days          |lateral flow travel time
!!    n_reduc     |              |nitrogen uptake reduction factor (not currently used; defaulted 300.)
!!    n_ln        |dimensionless |power function exponent for calculating nitrate concentration in
!!                                 subsurface drains (1.0 - 3.0)
!!    n_lnco      |dimensionless |coefficient for power function for calculating nitrate concentration
!!                                 in subsurface drains (0.5 - 4.0)
!!    ov_n(:)     |none          |Manning's "n" value for overland flow
!!    pot_fr(:)   |km2/km2       |fraction of HRU area that drains into pothole
!!    pot_k       |(mm/hr)       |hydraulic conductivity of soil surface of pothole
!!                   [defaults to conductivity of upper soil (0.01--10.) layer]
!!    pot_solp(:) |1/d           | soluble P loss rate in the pothole (.01 - 0.5)
!!    pot_tile(:) |m3/s          |average daily outflow to main channel from
!!                               |tile flow if drainage tiles are installed in
!!                               |pothole (needed only if current HRU is IPOT)
!!    pot_vol(:)  |mm            |initial volume of water stored in the
!!                               |depression/impounded area (read in as mm
!!                               |and converted to m^3) (needed only if current
!!                               |HRU is IPOT)
!!    pot_volx(:) |mm            |maximum volume of water stored in the
!!                               |depression/impounded area (read in as mm
!!                               |and converted to m^3) (needed only if current
!!                               |HRU is IPOT)
!!    r2adj       |dimensionless |curve number retention parameter adjustment factor to
!!                                 adjust surface runoff for flat slopes (0.5 - 3.0)
!!    rsdin(:)    |kg/ha         |initial residue cover
!!    slsoil(:)   |m             |slope length for lateral subsurface flow
!!    slsubbsn(:) |m             |average slope length for subbasin
!!    surlag      |days          |Surface runoff lag time.
!!                               |This parameter is needed in subbasins where
!!                               |the time of concentration is greater than 1
!!                               |day. SURLAG is used to create a "storage" for
!!                               |surface runoff to allow the runoff to take
!!                               |longer than 1 day to reach the subbasin outlet
!!    usle_ls(:)  |none          |USLE equation length slope (LS) factor
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag (=-1 if eof, else =0)
!!    epcohru     |none          |plant water uptake compensation factor (0-1)
!!    escohru     |none          |soil evaporation compensation factor (0-1)
!!    evpot       |none          |pothole evaporation coefficient
!!    fld_fr      |km2/km2       |fraction of HRU area that drains into floodplain (not used)
!!    n_lag       |dimensionless |lag coefficient for calculating nitrate concentration in subsurface
!!                                 drains (0.001 - 1.0)
!!    pot_no3l    |1/day         |nitrate decay rate in impounded area
!!    pot_nsed    |mg/L          |normal sediment concentration in impounded
!!                               |water (needed only if current HRU is IPOT)
!!    rip_fr      |km2/km2       |fraction of HRU area that drains into riparian
!!                               |zone (not used)
!!    sin_sl      |none          |Sin(slope angle)
!!    titldum     |NA            |title line of .sub file (not used)
!!    xm          |none          |exponential in equation to calculate
!!                               |USLE LS
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Sin, Atan

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: i, j
   character (len=80) :: titldum
   real*8 :: epcohru, escohru, evpot, fld_fr, n_lag, pot_no3l, pot_nsed,&
      &pot_solpl, pot_volmm, r2adj, rip_fr, sin_sl, xm, xx
   integer :: eof



   eof = 0
   escohru = 0.
   epcohru = 0.
   pot_nsed = 0.

   do
      read (108,5100) titldum
      read (108,*) hru_fr(j)
      read (108,*) slsubbsn(j)
      read (108,*) hru_slp(j)
      read (108,*) ov_n(j)
      read (108,*) lat_ttime(j)
      read (108,*) lat_sed(j)   !read in in mg/L
      read (108,*) slsoil(j)
      read (108,*,iostat=eof) canmx(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) escohru
      if (eof < 0) exit
      read (108,*,iostat=eof) epcohru
      if (eof < 0) exit
      read (108,*,iostat=eof) rsdin(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) erorgn(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) erorgp(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_fr(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) fld_fr ! not used
      if (eof < 0) exit
      read (108,*,iostat=eof) rip_fr ! not used
      if (eof < 0) exit
      read (108,5100,iostat=eof) titldum
      if (eof < 0) exit
!      if (ipot(j) == j) then   Srini pothole
      read (108,*,iostat=eof) pot_tilemm(j)    !!NUBZ
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_volxmm(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_volmm
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_nsed
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_no3l
      if (eof < 0) exit
!        read (108,5100,iostat=eof) titldum
!        if (eof < 0) exit
!        read (108,5100,iostat=eof) titldum
!        if (eof < 0) exit
!        read (108,5100,iostat=eof) titldum
!        if (eof < 0) exit
!        read (108,5100,iostat=eof) titldum
!        if (eof < 0) exit
!        read (108,5100,iostat=eof) titldum
!        if (eof < 0) exit
!      end if
      read (108,*,iostat=eof) dep_imp(j)
      if (eof < 0) exit
      read (108,5100,iostat=eof) titldum
      if (eof < 0) exit
      read (108,5100,iostat=eof) titldum
      if (eof < 0) exit
      read (108,5100,iostat=eof) titldum
      if (eof < 0) exit
      read (108,*,iostat=eof) evpot
      if (eof < 0) exit
      read (108,*,iostat=eof) dis_stream(j)
      if (eof < 0) exit
!! armen & stefan changes for SWAT-C
      read (108,*,iostat=eof) cf(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) cfh(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) cfdec(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) sed_con(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) orgn_con(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) orgp_con(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) soln_con(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) solp_con(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_solpl
      if (eof < 0) exit
      read (108,*,iostat=eof) pot_k(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) n_reduc(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) n_lag
      if (eof < 0) exit
      read (108,*,iostat=eof) n_ln(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) n_lnco(j)
!-------------------------------------------------------Moriasi 4/8/2014
      if (eof < 0) exit
      read (108,*,iostat=eof) surlag(j)
      if (eof < 0) exit

      read (108,*,iostat=eof) r2adj
!      read (108,*,iostat=eof) r2adj(j) !Soil retention parameter D. Moriasi 4/8/2014
! not used

      if (eof < 0) exit
      read (108,*,iostat=eof) cmn(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) cdn(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) nperco(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) phoskd(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) psp(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) sdnco(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) iwetile(j)
      if (eof < 0) exit
      read (108,*,iostat=eof) iwetgw(j)
      exit
   end do

   if (iwetile(j) <= 0) iwetile(j) = 0
   if (iwetgw(j) <= 0) iwetgw(j) = 0

   if (n_reduc(j) <= 0.) n_reduc(j) = 300.
   ! if (n_lag(j) <= 0.) n_lag(j) = 0.25 ! not used
   if (n_ln(j) <= 0.) n_ln(j) = 2.0
   if (n_lnco(j) <= 0.) n_lnco(j) = 2.0

!!    compare .hru input values to .bsn input values
   if (escohru > 1.e-4) esco(j) = escohru
   if (epcohru > 1.e-4) epco(j) = epcohru

!!    set default values
   if (dep_imp(j) <= 0.) dep_imp(j) = depimp_bsn
   if (surlag(j) <= 0.) surlag(j) = surlag_bsn
   if (cdn(j) <= 0.) cdn(j) = cdn_bsn
   if (nperco(j) <= 0.) nperco(j) = nperco_bsn
   if (cmn(j) <= 0.) cmn(j) = cmn_bsn
   if (phoskd(j) <= 0.) phoskd(j) = phoskd_bsn
   if (psp(j) <= 0.) psp(j) = psp_bsn
   if (sdnco(j) <= 0.) sdnco(j) = sdnco_bsn
!New and modified parameters D. Moriasi 4/8/2014
!   if (r2adj(j) <= 0.) r2adj(j) = r2adj_bsn ! not used
!   if (r2adj(j) > 0.95) r2adj(j) = 0.95 ! not used
!! comment the following line for the hru_fraction data !!
   if (hru_fr(j) <= 0.) hru_fr(j) = .0000001
   if (slsubbsn(j) <= 0.) slsubbsn(j) = 50.0
   if (hru_slp(j) <= 0.0001) hru_slp(j) = .0001
   if (hru_slp(j) >= 1.0) hru_slp(j) = 1.0
   if (slsoil(j) <= 0.)  slsoil(j) = slsubbsn(j)
   if (esco(j) <= 0.) esco(j) = .95
!     if (dep_imp(j) <= 0.) dep_imp(j) = 6000.
!     esco(j) = 1. - esco(j)
   if (epco(j) <= 0. .or. epco(j) > 1.) epco(j) = 1.0
!   if (evpot(j) <= 0.) evpot(j) = 0.5 ! not used
   if (dis_stream(j) <= 0.) dis_stream(j) = 35.0

!! armen & stefan changes for SWAT-C
   if (cf(j) <= 0.) cf(j)= 1.0
   if (cfh(j) <= 0.) cfh(j)= 1.0
   if (cfdec(j) <= 0.) cfdec(j)= 0.055
!! armen & stefan end


!!    calculate USLE slope length factor
   xm = .6 * (1. - Exp(-35.835 * hru_slp(j)))
   sin_sl = Sin(Atan(hru_slp(j)))
   usle_ls(j) = (slsubbsn(j)/22.128)**xm * (65.41 * sin_sl *&
      &sin_sl + 4.56 * sin_sl + .065)

!!    other calculations
   hru_km(j) = sub_km(i) * hru_fr(j)
   hru_ha(j) = hru_km(j) * 100.
   lat_sed(j) = lat_sed(j) * 1.e-3     !!mg/L => g/L
   pot_vol(j) = pot_volmm
!   pot_volx(j) = pot_volxmm(j) ! not used
!   pot_tile(j) = pot_tilemm(j) ! not used

   xx = 10. * pot_volmm * hru_ha(j) / 1000000.  !! mg/L * m3 * 1000L/m3 * t/1,000,000,000   Srini pothole
   pot_sed(j) = pot_nsed * xx
   pot_san(j) = 0.
   pot_sil(j) = 0.
   pot_cla(j) = pot_nsed * xx
   pot_sag(j) = 0.
   pot_lag(j) = 0.

   close (108)
   return
5100 format (a)
end
