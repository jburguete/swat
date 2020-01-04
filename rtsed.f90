!> @file rtsed.f90
!> file containing the subroutine rtsed
!> @author
!> modified by Javier Burguete

!> this subroutine routes sediment from subbasin to basin outlets
!> deposition is based on fall velocity and degradation on stream
!> @param[in] jrch reach number
subroutine rtsed(jrch)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    jrch        |none          |reach number
!!    ch_cov(1,:)  |none          |channel erodibility factor (0.0-1.0)
!!                               |0 non-erosive channel
!!                               |1 no resistance to erosion
!!    ch_cov(2,:)  |none          |channel cover factor (0.0-1.0)
!!                               |0 channel is completely protected from
!!                               |  erosion by cover
!!                               |1 no vegetative cover on channel
!!    ch_d(:)     |m             |average depth of main channel
!!    ch_di(:)    |m             |initial depth of main channel
!!    ch_li(:)    |km            |initial length of main channel
!!    ch_n(2,:)   |none          |Manning's "n" value for the main channel
!!    ch_s(2,:)   |m/m           |average slope of main channel
!!    ch_si(:)    |m/m           |initial slope of main channel
!!    ch_w(2,:)   |m             |average width of main channel
!!    ch_wdr(:)   |m/m           |channel width to depth ratio
!!    ideg        |none          |channel degredation code
!!                               |0: do not compute channel degradation
!!                               |1: compute channel degredation (downcutting
!!                               |   and widening)
!!    inum2       |none          |inflow hydrograph storage location number
!!    phi(5,:)    |m^3/s         |flow rate when reach is at bankfull depth
!!    prf(:)      |none          |Reach peak rate adjustment factor for sediment
!!                               |routing in the channel. Allows impact of
!!                               |peak flow rate on sediment routing and
!!                               |channel reshaping to be taken into account
!!    rchdep      |m             |depth of flow on day
!!    rnum1       |none          |fraction of overland flow
!!    sdti        |m^3/s         |average flow on day in reach
!!    sedst(:)    |metric tons   |amount of sediment stored in reach
!!    spcon       |none          |linear parameter for calculating sediment
!!                               |reentrained in channel sediment routing
!!    spexp       |none          |exponent parameter for calculating sediment
!!                               |reentrained in channel sediment routing
!!    varoute(3,:)|metric tons   |sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_d(:)     |m             |average depth of main channel
!!    ch_s(2,:)   |m/m           |average slope of main channel
!!    ch_w(2,:)   |m             |average width of main channel
!!    peakr       |m^3/s         |peak runoff rate in channel
!!    sedst(:)    |metric tons   |amount of sediment stored in reach
!!    sedrch      |metric tons   |sediment transported out of channel
!!                               |during time step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cych
!!    cyin
!!    dat2        |m             |change in channel depth during time step
!!    deg         |metric tons   |sediment reentrained in water by channel
!!                               |degradation
!!    deg1
!!    deg2
!!    dep         |metric tons   |sediment deposited on river bottom
!!    depdeg      |m             |depth of degradation/deposition from original
!!    depnet      |metric tons   |
!!    outfract
!!    qdin        |m^3 H2O       |water in reach during time step
!!    sedin
!!    vc          |m/s           |flow velocity in reach
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Abs, Max
!!    SWAT: ttcoef

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: jrch
   real*8 :: cych, cyin, dat2, deg, deg1, deg2, dep, depdeg, depnet, outfract,&
      &qdin, sedin, vc

   sedin = 0.0

   if (rtwtr > 0. .and. rchdep > 0.) then

!! initialize water in reach during time step
      qdin = rtwtr + rchstor(jrch)

!! do not perform sediment routing if no water in reach
      if (qdin > 0.01) then

!! initialize sediment in reach during time step
         sedin = varoute(3,inum2) * (1. - rnum1) + sedst(jrch)
         !sedinorg = sedin ! not used
!! initialize reach peak runoff rate
         peakr = prf(jrch) * sdti

!! calculate flow velocity
         if (rchdep < .010) then
            vc = 0.01
         else
            vc = peakr / rcharea
         end if
         if (vc > 5.) vc = 5.

!! JIMMY'S NEW IMPROVED METHOD for sediment transport
         cyin = sedin / qdin
         cych = spcon(jrch) * vc ** spexp(jrch)
         depnet = qdin * (cych - cyin)
         if (Abs(depnet) < 1.e-6) depnet = 0.

!!  tbase is multiplied so that erosion is proportional to the traveltime,
!!  which is directly related to the length of the channel
!!  Otherwise for the same discharge rate and sediment deficit
!!  the model will erode more sediment per unit length of channel
!!  from a small channel than a larger channel. Modification made by Balaji Narasimhan

         if (depnet > 1.e-6) then
            deg = depnet
            !! First the deposited material will be degraded before channel bed
            if (deg >= depch(jrch)) then
               deg1 = depch(jrch)
               deg2 = (deg - deg1) * ch_erodmo(jrch,i_mo) * ch_cov(2,jrch)
            else
               deg1 = deg
               deg2 = 0.
            endif
            dep = 0.
         else
            dep = -depnet
            deg = 0.
            deg1 = 0.
            deg2 = 0.
         endif

         depch(jrch) = depch(jrch) + dep - deg1
         if (depch(jrch) < 1.e-6) depch(jrch) = 0.

         sedin = sedin + deg1 + deg2 - dep
         if (sedin < 1.e-6) sedin = 0.

         outfract = rtwtr / qdin
         if (outfract > 1.) outfract = 1.

         sedrch = sedin * outfract
         if (sedrch < 1.e-6) sedrch = 0.

         sedst(jrch) = sedin - sedrch
         if (sedst(jrch) < 1.e-6) sedst(jrch) = 0.

!!  In this default sediment routing sediment is not tracked by particle size
         rch_san = 0.
         rch_sil = sedrch  !! As particles are not tracked by size, the sediments
         rch_cla = 0.      !! in reach is assumed to be silt for mass conservation
         rch_sag = 0.
         rch_lag = 0.
         rch_gra = 0.

!!    Bank erosion
         rchdy(55,jrch) = 0.
!!    Channel Degredation
         rchdy(56,jrch) = deg2
!!    Channel Deposition
         rchdy(57,jrch) = dep
!!    Floodplain Deposition
         rchdy(58,jrch) = 0.
!!    Total suspended sediments
         rchdy(59,jrch) = sedrch / rtwtr * 1.e6

!!    Organic nitrogen and Organic Phosphorus contribution from channel erosion
         ch_orgn(jrch) = deg2 * ch_onco(jrch) / 1000.
         ch_orgp(jrch) = deg2 * ch_opco(jrch) / 1000.

!! compute changes in channel dimensions
         if (ideg == 1) then
            depdeg = ch_d(jrch) - ch_di(jrch)
            if (depdeg < ch_si(jrch) * ch_li(jrch) * 1000.) then
               if (qdin > 1400000.) then
                  dat2 =  358.6 * rchdep * ch_s(2,jrch) * ch_cov(1,jrch)
                  ch_d(jrch) = ch_d(jrch) + dat2
                  ch_w(2,jrch) = ch_wdr(jrch) * ch_d(jrch)
                  ch_s(2,jrch) = ch_s(2,jrch) - dat2 / (ch_l(2,jrch) * 1000.)
                  ch_s(2,jrch) = Max(.0001, ch_s(2,jrch))
                  call ttcoef(jrch)
               endif
            endif
         endif

      else
         sedrch = 0.
         rch_san = 0.
         rch_sil = 0.
         rch_cla = 0.
         rch_sag = 0.
         rch_lag = 0.
         rch_gra = 0.
         sedst(jrch) = sedin
!!    Bank erosion
         rchdy(55,jrch) = 0.
!!    Channel Degredation
         rchdy(56,jrch) = 0.
!!    Channel Deposition
         rchdy(57,jrch) = 0.
!!    Floodplain Deposition
         rchdy(58,jrch) = 0.
!!    Total suspended sediments
         rchdy(59,jrch) = 0.

      endif !! end of qdin > 0.01 loop

   endif  !! end of rtwtr and rchdep > 0 loop

   return
end
