!> @file route.f90
!> file containing the subroutine route
!> @author
!> modified by Javier Burguete

!> this subroutine simulates channel routing
!> @param[in] i current day in simulation--loop counter (julian date)
!> @param[in] jrch reach number (none)
!> @param[in] k inflow hydrograph storage location number (none)
subroutine route(i, jrch, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |julian date   |current day in simulation--loop counter
!!    jrch        |none          |reach number
!!    k           |none          |inflow hydrograph storage location number
!!    alpha_bnke(:)|none         |Exp(-alpha_bnk(:))
!!    bankst(:)   |m^3 H2O       |bank storage
!!    ch_eqn      |              |sediment routing methods
!!                               | 0 = original SWAT method
!!                               | 1 = Bagnold's
!!                               | 2 = Kodatie
!!                               | 3 = Molinas Wu
!!                               | 4 = Yang
!!    ch_l(2,:)   |km            |length of main channel
!!    ch_revap(:) |none          |revap coeff: this variable controls the amount
!!                               |of water moving from bank storage to the root
!!                               |zone as a result of soil moisture depletion
!!    ch_w(2,:)   |m             |average width of main channel
!!    da_ha       |ha            |area of watershed in hectares
!!    hru_sub(:)  |none          |subbasin number for HRU
!!    ievent      |none          |rainfall/runoff code
!!                               |0 daily rainfall/curve number technique
!!                               |1 sub-daily rainfall/Green&Ampt/hourly
!!                               |  routing
!!                               |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    irte        |none          |water routing method:
!!                               |0 variable storage method
!!                               |1 Muskingum method
!!    iwq         |none          |stream water quality code
!!                               |0 do not model stream water quality
!!                               |1 model stream water quality (QUAL2E)
!!    nhru        |none          |number of HRUs in watershed
!!    pet_day     |mm H2O        |potential evapotranspiration on day
!!    rchdep      |m             |depth of flow on day
!!    rnum1       |none          |fraction of overland flow
!!    rttlc       |m^3 H2O       |transmission losses from reach on day
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    sub_fr(:)   |none          |fraction of watershed area in subbasin
!!    varoute(3,:)|metric tons   |sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    revapday    |m^3 H2O       |amount of water moving from bank storage
!!                               |into the soil profile or being taken
!!                               |up by plant roots in the bank storage zone
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    sedrch      |metric tons   |sediment transported out of reach on day
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    qdbank      |m^3 H2O       |streamflow contribution from bank storage
!!    rnum1i      |none          |1 - rnum1
!!    subwtr
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min
!!    SWAT: rchinit, rtday, rtmusk, rthvsc, rthmusk, rtsed, rtsed2, rthsed,
!!          watqual2, watqual, noqual, hhwatqual, hhnoqual, rtpest, rthpest,
!!          rtbact, irr_rch, rchuse, rtout
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i, jrch, k 
   real*8 :: qdbank, rnum1i, subwtr
   integer :: ii

   !inum3 is the subbasin for stream-aquifer interaction
   !inum5 is the landscape within the subbasin
   iru = inum5

!! initialize variables for route command loop
   call rchinit(jrch, k)

   vel_chan(jrch) = 0.

!! route water through reach
   if (ievent == 0) then
      if (irte == 0) call rtday(jrch, k)
      if (irte == 1) call rtmusk(i, jrch, k)
   else
      if (irte == 0) call rthvsc(jrch, k)
      if (irte == 1) call rthmusk(i, jrch, k)
   endif

!! average daily water depth for sandi doty 09/26/07
   dep_chan(jrch) = rchdep

!! if reach is an irrigation canal, restrict outflow
   if (icanal(jrch) == 1) then
      rchstor(jrch) = rchstor(jrch) + rtwtr
      rtwtr = 0.
   end if

!! add transmission losses to bank storage/deep aquifer in subbasin

   if (rttlc > 0.) then
      bankst(jrch) = bankst(jrch) + rttlc * (1. - trnsrch)
      if (da_ha > 1.e-9) then
         subwtr = rttlc * trnsrch / (da_ha * sub_fr(jrch) * 10.)
         do ii = hru1(jrch), hru1(jrch) + hrutot(jrch) - 1
            deepst(ii) = deepst(ii) + subwtr
         end do
      end if
   end if

!! compute revap from bank storage
   revapday = ch_revap(jrch) * pet_day * ch_l(2,jrch) * ch_w(2,jrch)
   revapday = Min(revapday,bankst(jrch))
   bankst(jrch) = bankst(jrch) - revapday

!! compute contribution of water in bank storage to streamflow
   qdbank = bankst(jrch) * (1. - alpha_bnke(jrch))
   bankst(jrch) = bankst(jrch) - qdbank
   rtwtr = rtwtr + qdbank
   if (ievent > 0) then
      do ii = 1, nstep
         hrtwtr(ii) = hrtwtr(ii) + qdbank / dfloat(nstep)
      end do
   end if


!! perform in-stream sediment calculations
   sedrch = 0.
   rch_san = 0.
   rch_sil = 0.
   rch_cla = 0.
   rch_sag = 0.
   rch_lag = 0.
   rch_gra = 0.
   ch_orgn(jrch) = 0.
   ch_orgp(jrch) = 0.
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

!! do not perform sediment routing for headwater subbasins
   !! when i_subhw = 0
   if (i_subhw == 0 .and. jrch == k) then
      rnum1i = 1. - rnum1
      if (ievent == 0) then
         if (rtwtr > 0. .and. rchdep > 0.) then
            sedrch  = varoute(3,k)  * rnum1i
            rch_san = varoute(23,k) * rnum1i
            rch_sil = varoute(24,k) * rnum1i
            rch_cla = varoute(25,k) * rnum1i
            rch_sag = varoute(26,k) * rnum1i
            rch_lag = varoute(27,k) * rnum1i
            rch_gra = varoute(28,k) * rnum1i
         end if
      else
         do ii = 1, nstep
            if (hrtwtr(ii) > 0. .and. hdepth(ii) > 0.) then
               hsedyld(ii) = hhvaroute(3,k,ii) * rnum1i
               sedrch = sedrch + hsedyld(ii)
               rch_san = 0.
               rch_sil = rch_sil + hsedyld(ii)  !!All are assumed to be silt type particles
               rch_cla = 0.
               rch_sag = 0.
               rch_lag = 0.
               rch_gra = 0.
            end if
         end do
      end if
   else
      if (ievent == 0) then
         if (ch_eqn(jrch) == 0) then
            call rtsed(jrch, k)
         else
            call rtsed2(jrch, k)
         end if
      else
         call rthsed(jrch, k)
         do ii = 1, nstep
            if (hrtwtr(ii) > 0. .and. hdepth(ii) > 0.) then
               sedrch = sedrch + hsedyld(ii)
               rch_sil = rch_sil + hsedyld(ii)  !!All are assumed to be silt type particles
            end if
         end do

      end if
   end if

!! perform in-stream nutrient calculations
   if (ievent == 0) then
      select case (iwq)
       case (2)
         call watqual2(jrch, k)
       case (1)
         call watqual(i, jrch, k)
       case (0)
         call noqual(jrch, k)
      end select
   else
      select case (iwq)
       case (1)
         call hhwatqual(jrch, k)
       case (0)
         call hhnoqual(jrch, k)
      end select
   end if

!! perform in-stream pesticide calculations
   if (ievent == 0) then
      call rtpest(jrch, k)
   else
      call rthpest(jrch, k)
   end if

!! perform in-stream bacteria calculations
   call rtbact(jrch, k)

!! remove water from reach for irrigation
   call irr_rch(jrch)

!! remove water from reach for consumptive water use
   call rchuse(jrch)

!! summarize output/determine loadings to next routing unit
   call rtout(jrch, k)

   return
end
