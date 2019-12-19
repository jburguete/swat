!> @file readgw.f90
!> file containing the suroutine readgw
!> @author
!> modified by Javier Burguete

!> this subroutine reads the parameters from the HRU/subbasin groundwater
!> input file (.gw)
subroutine readgw

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |subbasin number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    alpha_bf(:) |1/days        |alpha factor for groundwater recession curve
!!    alpha_bf_d(:) | 1/days     |alpha factor for groudwater recession curve of the deep aquifer
!!    alpha_bfe(:)|none          |Exp(-alpha_bf(:))
!!    alpha_bfe_d (:) |1/days    |Exp(-alpha_bf_d(:)) for deep aquifer
!!    ch_revap(:) |none          |revap coeff: this variable controls the amount
!!                               |of water moving from bank storage to the root
!!                               |zone as a result of soil moisture depletion
!!    deepst(i)   |mm H2O        |depth of water in deep aquifer
!!    delay(:)    |days          |groundwater delay: time required for water
!!                               |leaving the bottom of the root zone to
!!                               |reach the shallow aquifer
!!    gw_delaye(:)|none          |Exp(-1./(delay(:))
!!    gw_revap(:) |none          |revap coeff: this variable controls the amount
!!                               |of water moving from the shallow aquifer to
!!                               |the root zone as a result of soil moisture
!!                               |depletion
!!    gw_spyld(:) |m**3/m**3     |specific yield for shallow aquifer
!!    gwht(:)     |m             |groundwater height
!!    gwminp(:)   |mg P/L        |soluble P concentration in groundwater
!!                               |loading to reach
!!    gwno3(:)    |mg N/L        |nitrate-N concentration in groundwater
!!                               |loading to reach
!!    gwqmn(:)    |mm H2O        |threshold depth of water in shallow aquifer
!!                               |required before groundwater flow will occur
!!    rchrg_dp(:) |none          |recharge to deep aquifer: the fraction of
!!                               |root zone percolation that reaches the deep
!!                               |aquifer
!!    revapmn(:)  |mm H2O        |threshold depth of water in shallow aquifer
!!                               |required to allow revap to occur
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    shallst_n(:)|ppm NO3-N     |nitrate concentration in shallow aquifer
!!                               |converted to kg/ha
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag
!!    hlife_ngw
!!    titldum     |NA            |title line for .gw file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   character (len=80) :: titldum
   real*8 :: hlife_ngw
   integer :: eof, j

   eof = 0
   hlife_ngw = 0.0
   j = ihru

   do
      read (110,5000) titldum
      read (110,*) shallst(j)
      read (110,*) deepst(j)
      read (110,*) delay(j)
      read (110,*) alpha_bf(j)
      read (110,*) gwqmn(j)
      read (110,*) gw_revap(j)
      read (110,*) revapmn(j)
      read (110,*) rchrg_dp(j)
      read (110,*,iostat=eof) gwht(j)
      if (eof < 0) exit
      read (110,*,iostat=eof) gw_spyld(j)
      if (eof < 0) exit
      read (110,*,iostat=eof) shallst_n(j)
      if (eof < 0) exit
      read (110,*,iostat=eof) gwminp(j)
      if (eof < 0) exit
      read (110,*,iostat=eof) hlife_ngw
      if (eof < 0) exit
!! organic n and p in the lateral flow     - by J.Jeong BREC 2011
      read (110,*,iostat=eof) lat_orgn(j)
      if (eof < 0) exit
      read (110,*,iostat=eof) lat_orgp(j)
      if (eof < 0) exit
      read (110,*,iostat=eof) alpha_bf_d(j)
      exit
   end do

!!    set default values for mike van liew
   if (hlife_ngw <= 0.) hlife_ngw = hlife_ngw_bsn
!!    set default values for mike van liew

!!    set default values
   if (deepst(j) <= 0.) deepst(j) = 1000.
   if (hlife_ngw <= 0.) hlife_ngw = 365.
   if (lat_orgn(j) <= 1.e-6) lat_orgn(j) = 0.
   if (lat_orgp(j) <= 1.e-6) lat_orgp(j) = 0.

!!    perform additional calculations
   alpha_bfe(j) = Exp(-alpha_bf(j))
   if(delay(j) < .1) delay(j) = .1
   gw_delaye(j) = Exp(-1./(delay(j) + 1.e-6))
   shallst_n(j) = shallst_n(j) * shallst(j) / 100.
   gw_nloss(j) = Exp(-.693 / hlife_ngw)

!!    alpha baseflow factor for deep aquifer according to Yi Luo
   alpha_bfe_d(j) = Exp(-alpha_bf_d(j))


!! assign values to channels
   ch_revap(i) = gw_revap(j)

!! assign values to channels
   ch_revap(i) = gw_revap(j)

   close (110)
   return
5000 format (a)
end
