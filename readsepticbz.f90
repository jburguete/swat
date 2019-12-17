!> @file readsepticbz.f90
!> file containing the subroutine readsepticbz
!> @author
!> modified by Javier Burguete


!> this subroutine reads data from the septic input file (.sep).  This file
!> contains information related to septic tanks modeled or defined at the
!> watershed level
subroutine readsepticbz

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads data from the septic input file (.sep).  This file
!!    contains information related to septic tanks modeled or defined at the
!!    watershed level

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru             |none          |HRU number

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bz_z(:)          |mm            |Depth of biozone layer
!!    bz_thk(:)        |mm            |thickness of biozone
!!    bio_bd(:)        |kg/m^3        |density of biomass
!!    coeff_bod_conv(:)|none          |BOD to live bacteria biomass conversion factor
!!    coeff_bod_dc(:)  |m^3/day       |BOD decay rate coefficient
!!    coeff_denitr(:)  |none          |Denitrification rate coefficient
!!    coeff_fc1(:)     |none          |field capacity calibration parameter 1
!!    coeff_fc2(:)     |none          |field capacity calibration parameter 2
!!    coeff_fecal(:)   |m^3/day       |Fecal coliform bacteria decay rate coefficient
!!    coeff_mrt(:)     |none          |mortality rate coefficient
!!    coeff_nitr(:)    |none          |Nitrification rate coefficient
!!    coeff_plq(:)     |none          |Conversion factor for plaque from TDS
!!    coeff_rsp(:)     |none          |respiration rate coefficient
!!    coeff_slg1(:)    |none          |slough-off calibration parameter
!!    coeff_slg2(:)    |none          |slough-off calibration parameter
!!    sep_cap(:)       |none          |Number of permanent residents in the hourse
!!    isep_typ(:)      |none          |Septic system type
!!    isep_opt(:)      |none          |Septic system operation flag (1=active,2=failing,3=not operated)
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof              |none          |end of file flag (=-1 if eof, else =0)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   character (len=80) :: titldum
   integer :: j, eof

!!    initialize variables
   eof = 0
   j = ihru


!! read septic parameters
   do
      read (172,1000) titldum
      read (172,*,iostat=eof) isep_typ(j)
      if (eof < 0) exit
      if (isep_typ(j) <= 0) return
      read (172,*,iostat=eof) isep_iyr(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) isep_opt(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) sep_cap(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) bz_area(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) isep_tfail(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) bz_z(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) bz_thk(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) sep_strm_dist(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) sep_den(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) bio_bd(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_bod_dc(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_bod_conv(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_fc1(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_fc2(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_fecal(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_plq(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_mrt(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_rsp(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_slg1(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_slg2(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_nitr(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_denitr(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_pdistrb(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_psorpmax(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_solpslp(j)
      if (eof < 0) exit
      read (172,*,iostat=eof) coeff_solpintc(j)
      exit
   end do

   coeff_mrt(j) = 0.01 * coeff_mrt(j)
   coeff_rsp(j) = 0.01 * coeff_rsp(j)
   coeff_slg1(j) = 0.001 * coeff_slg1(j)
   coeff_nitr(j) = 0.01 * coeff_nitr(j)
   coeff_denitr(j) = 0.01 * coeff_denitr(j)

   !!Convert QSTE from volume to depth unit, mm
   qstemm(j) = sptqs(isep_typ(j)) * sep_cap(j) / bz_area(j) * 1000.

!!    set default values for undefined parameters
   if (isep_iyr(j)==0) isep_iyr(j) = iyr
   if (bz_z(j) <= 1.e-6) bz_z(j) = 500.
   if (bz_thk(j) <= 1.e-6) bz_thk(j) = 20.
   if (bio_bd(j) <= 1.e-6) bio_bd(j) = 1000.
   if (coeff_bod_dc(j) <= 1.e-6) coeff_bod_dc(j) = 9.33
   if (coeff_bod_conv(j) <= 1.e-6) coeff_bod_conv(j) = 0.42
   if (coeff_fc1(j) <= 1.e-6) coeff_fc1(j) = 30.0
   if (coeff_fc2(j) <= 1.e-6) coeff_fc2(j) = 0.7
   if (coeff_fecal(j) <= 1.e-6) coeff_fecal(j) = 0.11
   if (coeff_plq(j) <= 1.e-6) coeff_plq(j) = 0.10
   if (coeff_mrt(j) <= 1.e-6) coeff_mrt(j) = 0.025
   if (coeff_rsp(j) <= 1.e-6) coeff_rsp(j) = 0.0156
   if (coeff_slg1(j) <= 1.e-6) coeff_slg1(j) = 4.e-8
   if (coeff_slg2(j) <= 1.e-6) coeff_slg2(j) = 1.5
   if (coeff_nitr(j) <= 1.e-6) coeff_nitr(j) = 0.086
   if (coeff_denitr(j) <= 1.e-6) coeff_denitr(j) = 0.00432


   close (172)
1000 format (a)
   return
end
