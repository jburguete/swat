!> @file pkq.f90
!> file containing the subroutine pkq
!> @author
!> modified by Javier Burguete

!> this subroutine computes the peak runoff rate for each HRU
!> and the entire subbasin using a modification of the rational formula
!> @param[in] iwave
!> flag to differentiate calculation of HRU and subbasin sediment calculation
!> (none)\n
!> iwave = 0 for HRU MUSLE(sedyld) each hru is calculated independently using
!> hru area and adjusted channel length\n
!> iwave = 1 subbasin # for subbasin MUSLE is computed for entire subbasin using
!> hru weighted KLSCP
!> @param[in] j HRU number (none)
subroutine pkq(iwave, j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    iwave        |none          |flag to differentiate calculation of HRU and
!!                                |subbasin sediment calculation
!!                                |iwave = 0 for HRU
!!                                |iwave = subbasin # for subbasin
!!    j            |none          |HRU number
!!    al5          |none          |fraction of daily rainfall that occurs
!!                                |during 0.5h highest intensity
!!    hru_km(:)    |km^2          |area of HRU in square kilometers
!!    qday         |mm H2O        |surface runoff that reaches main channel
!!                                |during day in HRU
!!    sub_km(:)    |km^2          |area of subbasin in square kilometers
!!    sub_qd(:)    |mm H2O        |surface runoff that reaches main channel
!!                                |during day in subbasin
!!    sub_tc(:)    |hr            |time of concentration for subbasin
!!    sub_tran(:)  |mm H2O        |transmission losses on day in subbasin
!!    tconc(:)     |hr            |time of concentration for HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    peakr        |m^3/s         |peak runoff rate
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    altc         |
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log
!!    SWAT: Expo

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 Expo
   integer, intent (in) :: iwave, j
   real*8 :: altc

   if (iwave > 0) then
      !! subbasin sediment calculations
      altc = 1. - Expo(2. * sub_tc(iwave) * Log(1. - al5))
      peakr = altc * (sub_qd(iwave) + sub_tran(iwave)) / sub_tc(iwave) !! mm/h
      peakr = peakr * sub_km(iwave) / 3.6                              !! m^3/s
   else
      !! HRU sediment calculations
      altc = 1. - Expo(2. * tconc(j) * Log(1. - al5))
      peakr = altc * qday / tconc(j)           !! mm/h
      peakr = peakr * hru_km(j) / 3.6          !! m^3/s
   end if

   return
end
