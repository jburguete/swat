!> @file readrte.f90
!> file containing the subroutine readrte
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the reach (main channel) input file
!> (.rte). This file contains data related to channel attributes. Only
!> one reach file should be made for each subbasin. If multiple HRUs are
!> modeled within a subbasin, the same .rte file should be listed for all
!> HRUs in file.cio
subroutine readrte

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |reach number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    alpha_bnk(:)  |days        |alpha factor for bank storage recession curve
!!    alpha_bnke(:) |none        |Exp(-alpha_bnk(:))
!!    ch_bed_bd(:)  |(g/cc)      |bulk density of channel bed sediment (1.1-1.9)
!!    ch_bed_d50(:) |            |D50(median) particle size diameter of channel
!!                               |  bed sediment (micrometers) (0.001 - 20)
!!    ch_bed_kd(:)  |            |erodibility of channel bed sediment by jet test
!!                               | (Peter Allen needs to give more info on this)
!!    ch_bnk_bd(:)  |(g/cc)      |bulk density of channel bank sediment (1.1-1.9)
!!    ch_bnk_kd(:)  |            |erodibility of channel bank sediment by jet test
!!                               | (Peter Allen needs to give more info on this)
!!    ch_bnk_d50(:) |            |D50(median) particle size diameter of channel
!!                               |  bank sediment (0.001 - 20)
!!    ch_cov(1,:)    |none        |channel erodibility factor (0.0-1.0)
!!                               |0 non-erosive channel
!!                               |1 no resistance to erosion
!!    ch_cov(2,:)    |none        |channel cover factor (0.0-1.0)
!!                               |0 channel is completely protected from
!!                               |  erosion by cover
!!                               |1 no vegetative cover on channel
!!    ch_d(:)       |m           |average depth of main channel
!!    ch_di(:)      |m           |initial depth of main channel
!!    ch_eqn        |            |sediment routine methods (DAILY):
!!                               | 0 = original SWAT method
!!                               | 1 = Bagnold's
!!                               | 2 = Kodatie
!!                               | 3 = Molinas WU
!!                               | 4 = Yang
!!    ch_erod(:)    |none        |channel erodibility factor (0.0-1.0)
!!                               |0 non-erosive channel
!!                               |1 no resistance to erosion
!!    ch_k(2,:)      |mm/hr       |effective hydraulic conductivity of
!!                               |main channel alluvium
!!    ch_l(2,:)      |km          |length of main channel
!!    ch_li(:)      |km          |initial length of main channel
!!    ch_n(2,:)      |none        |Manning's "n" value for the main channel
!!    ch_onco(:)    |ppm         |channel organic n concentration
!!    ch_opco(:)    |ppm         |channel organic p concentration
!!    ch_s(2,:)      |m/m         |average slope of main channel
!!    ch_si(:)      |m/m         |initial slope of main channel
!!    ch_w(2,:)      |m           |average width of main channel
!!    ch_wdr(:)     |m/m         |channel width to depth ratio
!!    chside(:)     |            |change in horizontal distance per unit
!!                               |  vertical distance (0.0 - 5)
!!                               |0 = for vertical channel bank
!!                               |5 = for channel bank with gentl side slope
!!    prf(:)        |none        |Reach peak rate adjustment factor for sediment
!!                               |routing in the channel. Allows impact of
!!                               |peak flow rate on sediment routing and
!!                               |channel reshaping to be taken into account.
!!    tc_bed        |N/m2        |critical shear stress of channel bed
!!    tc_bnk        |N/m2        |critical shear stress of channel bank
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bedsize
!!    banksize
!!    eof         |none          |end of file flag
!!    j           |none          |counter
!!    mo
!!    SC
!!    sumerod
!!    titldum     |NA            |title line of .rte file (not used elsewhere)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ SUBROUTINES/FUNCTIONS ~ ~ ~ ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   character (len=80) :: titldum
   real*8 :: bedsize, bnksize, SC, sumerod
   integer :: eof, j, mo

   eof = 0
   j = irch
   do
      read (103,5000) titldum
      read (103,*) ch_w(2,j)
      read (103,*) ch_d(j)
      read (103,*) ch_s(2,j)
      read (103,*) ch_l(2,j)
      read (103,*) ch_n(2,j)
      read (103,*) ch_k(2,j)
      read (103,*) ch_cov(1,j)
      read (103,*,iostat=eof) ch_cov(2,j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_wdr(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) alpha_bnk(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) icanal(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_onco(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_opco(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) chside(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_bnk_bd(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_bed_bd(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_bnk_kd(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_bed_kd(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_bnk_d50(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_bed_d50(j)
      if (eof < 0) exit
      read (103,5000,iostat=eof) tc_bnk(j)
      if (eof < 0) exit
      read (103,5000,iostat=eof) tc_bed(j)
      if (eof < 0) exit
      read (103,5100,iostat=eof) (ch_erodmo(j,mo), mo = 1,12)
      if (eof < 0) exit
      read (103,*,iostat=eof) ch_eqn(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) prf(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) spcon(j)
      if (eof < 0) exit
      read (103,*,iostat=eof) spexp(j)
      if (eof < 0) exit
      exit
   end do

!!    set default values for parameters
   if (ch_s(2,j) <= 0.) ch_s(2,j) = .0001
   if (ch_n(2,j) <= 0.01) ch_n(2,j) = .01
   if (ch_n(2,j) >= 0.70) ch_n(2,j) = 0.70
   if (ch_l(2,j) <= 0.) ch_l(2,j) = .0010
   if (ch_wdr(j) <= 0.) ch_wdr(j) = 3.5
   if (chside(j) <= 1.e-6) chside(j) = 2.0
   if (tc_bnk(j) <= 0.) tc_bnk(j)=0. !! Critical shear stress (N.m^2)
   if (tc_bed(j) <= 0.) tc_bed(j)=0. !! Critical shear stress (N.m^2)
   if (prf(j) <= 0.) prf(j) = prf_bsn
   if (spcon(j) <= 0.) spcon(j) = spcon_bsn
   if (spexp(j) <= 0.) spexp(j) = spexp_bsn

   if (ch_eqn(j) <= 0) then
      ch_eqn(j)=0 !! SWAT Default sediment routing routine
      if (ch_cov(1,j) <= 0.0) ch_cov(1,j) = 0.0
      if (ch_cov(2,j) <= 0.0) ch_cov(2,j) = 0.0
      if (ch_cov(1,j) >= 1.0) ch_cov(1,j) = 1.0
      if (ch_cov(2,j) >= 1.0) ch_cov(2,j) = 1.0
   else
      if (ch_cov(1,j) <= 0.0) ch_cov(1,j) = 1.0
      if (ch_cov(2,j) <= 0.0) ch_cov(2,j) = 1.0
      if (ch_cov(1,j) >= 25.) ch_cov(1,j) = 25.
      if (ch_cov(2,j) >= 25.) ch_cov(2,j) = 25.
   end if


!!    Bank material is assumed to be silt type partcile if not given.
   if (ch_bnk_d50(j) <= 1.e-6) ch_bnk_d50(j) = 50. !! Units are in Micrometer
   if (ch_bnk_d50(j) > 10000) ch_bnk_d50(j) = 10000.


   bnksize = ch_bnk_d50(j)/1000.  !! Units conversion Micrometer to Millimeters
!!    Channel sediment particle size distribution
!!    Clayey bank
   if (bnksize <= 0.005) then
      ch_bnk_cla(j) = 0.65
      ch_bnk_sil(j) = 0.15
      ch_bnk_san(j) = 0.15
      ch_bnk_gra(j) = 0.05
   end if

!!    Silty bank
   if (bnksize > 0.005 .and. bnksize <= 0.05) then
      ch_bnk_sil(j) = 0.65
      ch_bnk_cla(j) = 0.15
      ch_bnk_san(j) = 0.15
      ch_bnk_gra(j) = 0.05
   end if

!!    Sandy bank
   if (bnksize > 0.05 .and. bnksize <= 2.) then
      ch_bnk_san(j) = 0.65
      ch_bnk_sil(j) = 0.15
      ch_bnk_cla(j) = 0.15
      ch_bnk_gra(j) = 0.05
   end if

!!    Gravel bank
   if (bnksize > 2.) then
      ch_bnk_gra(j) = 0.65
      ch_bnk_san(j) = 0.15
      ch_bnk_sil(j) = 0.15
      ch_bnk_cla(j) = 0.05
   end if

!!    Bed material is assumed to be sand type partcile if not given.
   if (ch_bed_d50(j) <= 1.e-6) ch_bed_d50(j) = 500 !! Units are in Micrometer
   if (ch_bed_d50(j) > 10000) ch_bed_d50(j) = 10000.

!!    Channel sediment particle size distribution
!!    Clayey bed
   bedsize = ch_bed_d50(j)/1000.  !! Units conversion Micrometer to Millimeters
   if (bedsize <= 0.005) then
      ch_bed_cla(j) = 0.65
      ch_bed_sil(j) = 0.15
      ch_bed_san(j) = 0.15
      ch_bed_gra(j) = 0.05
   end if

!!    Silty bed
   if (bedsize > 0.005 .and. bedsize <= 0.05) then
      ch_bed_sil(j) = 0.65
      ch_bed_cla(j) = 0.15
      ch_bed_san(j) = 0.15
      ch_bed_gra(j) = 0.05
   end if

!!    Sandy bed
   if (bedsize > 0.05 .and. bedsize <= 2.) then
      ch_bed_san(j) = 0.65
      ch_bed_sil(j) = 0.15
      ch_bed_cla(j) = 0.15
      ch_bed_gra(j) = 0.05
   end if

!!    Gravel bed
   if (bedsize > 2.) then
      ch_bed_gra(j) = 0.65
      ch_bed_san(j) = 0.15
      ch_bed_sil(j) = 0.15
      ch_bed_cla(j) = 0.05
   end if

!!    Bulk density of channel bank sediment
   if (ch_bnk_bd(j) <= 1.e-6) ch_bnk_bd(j) = 1.40 !! Silty loam bank

!!    Bulk density of channel bed sediment
   if (ch_bed_bd(j) <= 1.e-6) ch_bed_bd(j) = 1.50  !! Sandy loam bed


!!    An estimate of Critical shear stress if it is not given (N/m^2)
!! Critical shear stress based on silt and clay %
!! Critical Shear Stress based on Julian and Torres (2005)
!!    Units of critical shear stress (N/m^2)
   if  (tc_bnk(j) <= 1.e-6) then
      SC = (ch_bnk_sil(j) + ch_bnk_cla(j)) * 100.
      tc_bnk(j) = (0.1 + (0.1779*SC) + (0.0028*(SC)**2)&
         &- ((2.34E-05)*(SC)**3)) * ch_cov(1,j)
   end if

   if  (tc_bed(j) <= 1.e-6) then
      SC = (ch_bed_sil(j) + ch_bed_cla(j)) * 100.
      tc_bed(j) = (0.1 + (0.1779*SC) + (0.0028*(SC)**2)&
         &- ((2.34E-05)*(SC)**3)) * ch_cov(2,j)
   end if

!!  An estimate of channel bank erodibility coefficient from jet test if it is not available
!!  Units of kd is (cm^3/N/s)
!!  Base on Hanson and Simon, 2001
   if (ch_bnk_kd(j) <= 1.e-6) then
      if (tc_bnk(j) <= 1.e-6) then
         ch_bnk_kd(j) = 0.2
      else
         ch_bnk_kd(j) = 0.2 / sqrt(tc_bnk(j))
      end if
   end if

!!  An estimate of channel bed erodibility coefficient from jet test if it is not available
!!  Units of kd is (cm^3/N/s)
!!  Base on Hanson and Simon, 2001
   if (ch_bed_kd(j) <= 1.e-6) then
      if (tc_bed(j) <= 1.e-6) then
         ch_bed_kd(j) = 0.2
      else
         ch_bed_kd(j) = 0.2 / sqrt(tc_bed(j))
      end if
   end if

   sumerod = 0.
   do mo = 1, 12
      sumerod = sumerod + ch_erodmo(j,mo)
   end do

   if (sumerod < 1.e-6) then
      do mo = 1, 12
         ch_erodmo(j,mo) = ch_cov(1,j)
      end do
   end if

!!    set default values for mike van liew
   if (ch_onco(j) <= 0.) ch_onco(j) = ch_onco_bsn
   if (ch_opco(j) <= 0.) ch_opco(j) = ch_opco_bsn
!!    set default values for mike van liew


!!    initialize variables for channel degradation
   ch_di(j) = ch_d(j)
   ch_li(j) = ch_l(2,j)
   ch_si(j) = ch_s(2,j)
   ch_wi(j) = ch_w(2,j)

   close (103)
   return
5000 format (a)
5100 format (12f6.2)
end
