!> @file readlwq.f90
!> file containing the subroutine readlwq
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the lake water quality input file (.lwq).
!> This file contains data related to initial pesticide and nutrient levels
!> in the lake/reservoir and transformation processes occuring within the
!> lake/reservoir. Data in the lake water quality input file is assumed to
!> apply to all reservoirs in the watershed.
!> @param[in] ii reservoir number (none)
subroutine readlwq(ii)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii           |none          |reservoir number
!!    res_vol(:)   |m**3          |reservoir volume
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chlar(:)      |none          |chlorophyll-a production coefficient for
!!                                 |reservoir
!!    ires1(:)      |none          |beginning of mid-year nutrient settling
!!                                 |"season"
!!    ires2(:)      |none          |end of mid-year nutrient settling "season"
!!    lkpst_conc(:) |mg/m**3       |pesticide concentration in lake water
!!    lkpst_koc(:)  |m**3/g        |pesticide partition coefficient between
!!                                 |water and sediment in lake water
!!    lkpst_mix(:)  |m/day         |mixing velocity (diffusion/dispersion) in
!!                                 |lake water for pesticide
!!    lkpst_rea(:)  |1/day         |pesticide reaction coefficient in lake water
!!    lkpst_rsp(:)  |m/day         |resuspension velocity in lake water for
!!                                 |pesticide sorbed to sediment
!!    lkpst_stl(:)  |m/day         |settling velocity in lake water for
!!                                 |pesticide sorbed to sediment
!!    lkpst_vol(:)  |m/day         |pesticide volatilization coefficient in lake
!!                                 |water
!!    lkspst_act(:) |m             |depth of active sediment layer in lake for
!!                                 |for pesticide
!!    lkspst_bry(:) |m/day         |pesticide burial velocity in lake bed
!!                                 |sediment
!!    lkspst_conc(:)|mg/m**3       |pesticide concentration in lake bed sediment
!!    lkspst_rea(:) |1/day         |pesticide reaction coefficient in lake bed
!!                                 |sediment
!!    nsetlr(1,:)   |m/day         |nitrogen settling rate for mid-year
!!                                 |period (read in as m/year and converted to
!!                                 |m/day)
!!    nsetlr(2,:)   |m/day         |nitrogen settling rate for remainder of
!!                                 |year (read in as m/year and converted to
!!                                 |m/day)
!!    psetlr(1,:)   |m/day         |phosphorus settling rate for mid-year
!!                                 |period (read in as m/year and converted to
!!                                 |m/day)
!!    psetlr(2,:)   |m/day         |phosphorus settling rate for remainder of
!!                                 |year (read in as m/year and converted to
!!                                 |m/day)
!!    res_nh3(:)    |kg N          |amount of ammonia in reservoir
!!    res_no2(:)    |kg N          |amount of nitrite in reservoir
!!    res_no3(:)    |kg N          |amount of nitrate in reservoir
!!    res_orgn(:)   |kg N          |amount of organic N in reservoir
!!    res_orgp(:)   |kg P          |amount of organic P in reservoir
!!    res_solp(:)   |kg P          |amount of soluble P in reservoir
!!    seccir(:)     |none          |water clarity coefficient for reservoir
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag
!!    lkarea
!!    nh3i        |mg N/L        |initial concentration of ammonia in reservoir
!!    no2i        |mg N/L        |initial concentration of nitrite in reservoir
!!    no3i        |mg N/L        |initial concentration of nitrate in reservoir
!!    orgni       |mg N/L        |initial concentration of organic N in
!!                               |reservoir
!!    orgpi       |mg P/L        |initial concentration of organic P in
!!                               |reservoir
!!    solpi       |mg P/L        |initial concentration of soluble P in
!!                               |reservoir
!!    titldum     |NA            |title line of .lwq file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: ii
   character (len=80) :: titldum
   real*8 :: lkarea, nh3i, no2i, no3i, orgni, orgpi, solpi
   integer :: eof

   eof = 0
   titldum = ""
   orgpi = 0.
   solpi = 0.
   orgni = 0.
   no3i = 0.
   nh3i = 0.
   no2i = 0.

!!    read lake water quality data
   do
      read (106,1000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,1000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,*,iostat=eof) ires1(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) ires2(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) psetlr(1,ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) psetlr(2,ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) nsetlr(1,ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) nsetlr(2,ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) chlar(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) seccir(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) orgpi
      if (eof < 0) exit
      read (106,*,iostat=eof) solpi
      if (eof < 0) exit
      read (106,*,iostat=eof) orgni
      if (eof < 0) exit
      read (106,*,iostat=eof) no3i
      if (eof < 0) exit
      read (106,*,iostat=eof) nh3i
      if (eof < 0) exit
      read (106,*,iostat=eof) no2i
      if (eof < 0) exit
      read (106,1000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_conc(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_rea(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_vol(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_koc(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_stl(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_rsp(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkpst_mix(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkspst_conc(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkspst_rea(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkspst_bry(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) lkspst_act(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) theta_n(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) theta_p(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) con_nirr(ii)
      if (eof < 0) exit
      read (106,*,iostat=eof) con_pirr(ii)
      if (eof < 0) exit
      exit
   end do

!!    convert units
   psetlr(1,ii) = psetlr(1,ii) / 365.         !m/yr -> m/day
   psetlr(2,ii) = psetlr(2,ii) / 365.
   nsetlr(1,ii) = nsetlr(1,ii) / 365.
   nsetlr(2,ii) = nsetlr(2,ii) / 365.
!     set initial n and p concentrations --> (ppm) * (m^3) / 1000 = kg
!                                            ppm = t/m^3 * 10^6
   res_solp(ii) = solpi * res_vol(ii) / 1000.
   res_orgp(ii) = orgpi * res_vol(ii) / 1000.
   res_no3(ii) = no3i  * res_vol(ii) / 1000.
   res_no2(ii) = no2i  * res_vol(ii) / 1000.
   res_nh3(ii) = nh3i  * res_vol(ii) / 1000.
   res_orgn(ii) = orgni * res_vol(ii) / 1000.

!!    lake pesticide mass
   lkpst_mass(ii) = lkpst_conc(ii) * res_vol(ii)
   lkarea = br1(ii) * res_vol(ii) ** br2(ii)
   lkspst_mass(ii) = lkspst_conc(ii) * lkspst_act(ii) * lkarea * 10000.

   close (106)

   return
1000 format (a80)
end
