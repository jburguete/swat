!> @file eiusle.f90
!> file containing the subroutine eiusle
!> @author
!> modified by Javier Burguete

!> this subroutine computes the USLE erosion index (EI)
subroutine eiusle(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none        |HRU number
!!    ab          |mm H2O      |lowest value xa can have
!!    amp_r(:,:)  |none        |alpha factor for rain (mo max 0.5h rain)
!!    idg(:)      |none        |array location of random number seed
!!                             |used for a given process
!!    i_mo        |none        |month being simulated
!!    ovrlnd(:)   |mm H2O      |overland flow onto HRU from upstream
!!                             |routing unit
!!    precipday   |mm H2O      |amount of water reaching soil surface
!!    rndseed(:,:)|none        |random number generator seed
!!    snomlt      |mm H2O      |amount of snow melt in HRU on current day
!!    tconc(:)    |hr          |time of concentration
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units                  |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    usle_ei     |100(ft-tn in)/(acre-hr)|USLE rainfall erosion index
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ajp         |mm H2O      |highest value xa can have
!!    pkrf        |none        |intermediate calculation
!!    pkrf30      |mm/hr       |maximum 30-min. storm intensity
!!    preceff     |mm H2O      |amount of rainfall impacting ground surface
!!                             |on day in HRU
!!    xa          |none        |fraction of daily rainfall occuring during
!!                             |half-hour of max intensity rain
!!    xb          |none        |intermediate calculation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log, Log10
!!    SWAT: Expo, Atri

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 Expo, Atri
   integer, intent(in) :: j
   real*8 :: ajp, pkrf, pkrf30, preceff, xa, xb

   preceff = precipday   !- snomlt - ovrlnd(j)
   if (preceff > 1.e-4) then
      ajp = 1. - Expo(-125. / preceff)
      xa = Atri(ab, amp_r(i_mo,hru_sub(j)), ajp, rndseed(idg(4),j))
      xb = -2. * Log(1. - xa)
      pkrf30 = 2. * preceff * xa
      pkrf = 2. * xb * preceff
      usle_ei = preceff * (12.1 + 8.9 * (Log10(pkrf) - .4343)) *&
         &preceff / 1000.
      if (usle_ei < 1.e-4) usle_ei = 0.
      usle_eifac(j) = usle_ei
   endif

   return
end
