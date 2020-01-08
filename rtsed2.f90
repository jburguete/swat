!> @file rtsed2.f90
!> file containing the subroutine rtsed2
!> @author
!> Balaji Narasimhan,\n
!> Peter Allen,\n
!> modified by Javier Burguete

!> this subroutine routes sediment from subbasin to basin outlets
!> deposition is based on fall velocity and degradation on stream.
!> Modification to the original SWAT sediment routine.
!> Bagnolds strempower, Kodatie (Modified Simons-Li associates), Molinas&Wu
!> strempower and Yang's sand-gravel equation approaches combined with
!> Einstein's deposition equation plus particle size tracking
!> @param[in] jrch reach number
!> @param[in] k inflow hydrograph storage location number (none)
subroutine rtsed2(jrch, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    jrch        |none          |reach number
!!    k           |none          |inflow hydrograph storage location number
!!    ch_cov(1,:) |none          |channel bank cover factor (0.0-1.0)
!!                               |0 channel is completely protected from
!!                               |  erosion by cover
!!                               |1 no vegetative cover on channel
!!    ch_cov(2,:) |none          |channel bed cover factor (0.0-1.0)
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
!!    adddep
!!    akod_a
!!    akod_b
!!    akod_c
!!    akod_d
!!    alog10cychppm
!!    asinea
!!    bedrt
!!    bedsize
!!    bnkcla
!!    bnkgra
!!    bnkrt
!!    bnkrte
!!    bnksan
!!    bnksil
!!    c
!!    clain
!!    cych
!!    cychppm
!!    cychv
!!    cychw
!!    cyin
!!    dat2        |m             |change in channel depth during time step
!!    deg         |metric tons   |sediment reentrained in water by channel
!!                               |degradation
!!    deg1
!!    deg1cla
!!    deg1gra
!!    deg1lag
!!    deg1sag
!!    deg1san
!!    deg1sil
!!    degcla
!!    deggra
!!    degremain
!!    degrte
!!    degsan
!!    degsil
!!    dep         |metric tons   |sediment deposited on river bottom
!!    depcla
!!    depdeg      |m             |depth of degradation/deposition from original
!!    depgra
!!    deplag
!!    depnet      |metric tons   |
!!    depsag
!!    depsan
!!    depsil
!!    effbnkbed
!!    fpratio
!!    grain
!!    lagin
!!    outfract
!!    pbank
!!    pbed
!!    pdep
!!    qcych
!!    qdin        |m^3 H2O       |water in reach during time step
!!    sagin
!!    sanin
!!    sedin
!!    SFbank
!!    silin
!!    Tbank
!!    Tbed
!!    topw
!!    Tou
!!    USpower
!!    var1
!!    var2
!!    var3
!!    var4
!!    var5
!!    var6
!!    var56
!!    vc          |m/s           |flow velocity in reach
!!    vcla
!!    vgra
!!    vlag
!!    vsag
!!    vsan
!!    vsh
!!    vsil
!!    w50
!!    watdep
!!    x
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt, Max, Log10, Min, Exp
!!    SWAT: ttcoef

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: jrch, k
!! Fall velocity Based on equation 1.36 from SWRRB manual
   real*8, parameter :: vcla = 411.0 * ((0.002)**2.) / (3600.),&
      &vgra = 411.0 * ((2.00)**2.) / (3600.),&
      &vlag = 411.0 * ((0.50)**2.) / (3600.),&
      &vsag = 411.0 * ((0.03)**2.) / (3600.),&
      &vsan = 411.0 * ((0.20)**2.) / (3600.),&
      &vsil = 411.0 * ((0.01)**2.) / (3600.)
   real*8 :: adddep, akod_a, akod_b, akod_c, akod_d, alog10cychppm, asinea,&
      &bedrt, bedsize, bnkcla, bnkgra, bnkrt, bnkrte, bnksan, bnksil, c, clain,&
      &cych, cychppm, cychv, cychw, cyin, dat2, deg, deg1, deg1cla, deg1gra,&
      &deg1lag, deg1sag, deg1san, deg1sil, degcla, deggra, degremain, degrte,&
      &degsan, degsil, dep, depcla, depdeg, depgra, deplag, depnet, depsag,&
      &depsan, depsil, effbnkbed, fpratio, grain, lagin, outfract, pbank, pbed,&
      &pdep, qcych, qdin, rh, sagin, sanin, sedin, SFbank, silin, Tbank, Tbed,&
      &topw, Tou, USpower, var1, var2, var3, var4, var5, var6, var56, vc, vsh,&
      &w50, watdep, x, xx

   if (rtwtr > 0. .and. rchdep > 0.) then

!! initialize water in reach during time step
      qdin = rtwtr + rchstor(jrch)

!! initialize sediment in reach during time step
      xx = 1. - rnum1
      sedin = varoute(3, k) * xx + sedst(jrch)
      sanin = varoute(23,k) * xx + sanst(jrch)
      silin = varoute(24,k) * xx + silst(jrch)
      clain = varoute(25,k) * xx + clast(jrch)
      sagin = varoute(26,k) * xx + sagst(jrch)
      lagin = varoute(27,k) * xx + lagst(jrch)
      grain = varoute(28,k) * xx + grast(jrch)
      !sedinorg = sedin ! not used

!! do not perform sediment routing if no water in reach
      if (qdin > 0.01) then

!! initialize reach peak runoff rate
         peakr = sdti
         if (ch_eqn(jrch) == 1) peakr = peakr * prf(jrch)

!! calculate peak flow velocity
         if (rcharea < .010) then
            vc = 0.01
         else
            vc = peakr / rcharea
         end if

         if (vc > 5.) vc = 5.

!! JIMMY'S NEW IMPROVED METHOD for sediment transport
         cyin = 0.
         cych = 0.
         depnet = 0.
         deg = 0.

         deg1 = 0.
         deg1san = 0.
         deg1sil = 0.
         deg1cla = 0.
         deg1sag = 0.
         deg1lag = 0.
         deg1gra = 0.

         degrte = 0.
         deggra = 0.
         degsan = 0.
         degsil = 0.
         degcla = 0.
         bnksan = 0.
         bnksil = 0.
         bnkcla = 0.
         bnkgra = 0.
         bnkrte = 0.
         dep = 0.
         depsan = 0.
         depsil = 0.
         depcla = 0.
         depsag = 0.
         deplag = 0.
         depgra = 0.
         bnkrt = 0.
         bedrt = 0.
         effbnkbed = 0.

         c = chside(jrch)
         pbed = phi(6,jrch)
         pbank = 2. * rchdep * Sqrt(1. + c * c)

         if (rchdep <= ch_d(jrch)) then
            topw = phi(6,jrch) + 2. * rchdep * c
            fpratio = 0.
            watdep = rchdep
         else
            topw = 5 * ch_w(2,jrch) + 2. * (rchdep - ch_d(jrch)) * 4.
            adddep = rchdep - ch_d(jrch)
            !! Area Ratio of water in flood plain to total cross sectional area
            fpratio = (rcharea - phi(1,jrch) - ch_w(2,jrch) * adddep) / rcharea
            fpratio = Max(0.,fpratio)
            watdep = ch_d(jrch)
         end if

!! Applied Bank Shear Stress
!!    Equations from Eaton and Millar (2004)
         SFbank = 10**(-1.4026 * Log10((pbed/pbank) + 1.5) + 2.247)

         Tou = 9800. * rchdep * ch_s(2,jrch)

         asinea = 1. / Sqrt((1.**2) + (c**2))

         Tbank = Tou * (SFbank/100.) * (topw + pbed) * asinea / (4. * rchdep)

         Tbed  = Tou * (1. - (SFbank/100.)) * (topw/(2.*pbed) + 0.5)

!!    Potential Bank Erosion rate in metric tons per day
!!    Assumed on an average Only one bank eroding due to meandering of channel
         bnkrte = ch_bnk_kd(jrch) * (Tbank - tc_bnk(jrch)) * 1e-06
         if (bnkrte < 0.) bnkrte = 0.
         bnkrte = bnkrte * ch_l(2,jrch) * 1000.* (watdep * Sqrt(1. + c * c))&
            &* ch_bnk_bd(jrch) * 86400.

!!    Potential Bed degradation rate in metric tons per day
         degrte = ch_bed_kd(jrch) * (Tbed - tc_bed(jrch)) * 1e-06
         if (degrte < 0.) degrte = 0.
         degrte = degrte * ch_l(2,jrch) * 1000.* phi(6,jrch)&
            &* ch_bed_bd(jrch) * 86400.

!!    Relative potential for bank/bed erosion
         if (bnkrte + degrte > 1.e-6) then
            bnkrt = bnkrte / (bnkrte + degrte)
         else
            bnkrt = 1.0
         end if
         bnkrt = Min(1.0, bnkrt)

!!    Relative potential for bed erosion
         bedrt = 1. - bnkrt

         select case (ch_eqn(jrch))
          case (1) ! Bagnold

!!         Streampower for sediment calculated based on Bagnold (1977) concept
           cych = spcon(jrch) * vc ** spexp(jrch)

          case (2) ! Kodatie

!!          Silt-bed rivers
            bedsize = ch_bed_d50(jrch)/1000.
            if (bedsize <= 0.05) then
               akod_a = 281.4
               akod_b = 2.622
               akod_c = 0.182
               akod_d = 0
            end if

!!          Very fine to fine-bed rivers
            if (bedsize > 0.05 .and. bedsize <= 0.25) then
               akod_a = 2829.6
               akod_b = 3.646
               akod_c = 0.406
               akod_d = 0.412
            end if

!!          Medium to very coarse sand-bed rivers
            if (bedsize > 0.25 .and. bedsize <= 2.) then
               akod_a = 2123.4
               akod_b = 3.300
               akod_c = 0.468
               akod_d = 0.613
            end if

!!          Gravel bed rivers
            if (bedsize > 2.) then
               akod_a = 431884.8
               akod_b = 1.000
               akod_c = 1.000
               akod_d = 2.000
            end if

!!          Maximum bedmaterial transport capacity in metric tons/m/day
            qcych =akod_a*(vc**akod_b)*(rchdep**akod_c)*(ch_s(2,jrch)**akod_d)

!!          Maximum bedmaterial transport capacity in metric tons/day
            cych = qcych/qdin * (topw + phi(6,jrch))/2.

          case (3) ! Molinas and Wu

!!          Streampower calculated based on Molinas and Wu model (2001)
!!          Settling velocity for ch_d50(jrch) particle size
            w50 = 411.0 * ((ch_bed_d50(jrch)/1000.)**2.) / (3600.)
            USpower = (vc ** 3.)/((2.65 - 1.) * 9.81 * rchdep * w50&
               &* ((Log10(rchdep/(ch_bed_d50(jrch) * 1e-006)))**2.))

!!          Concentration in ppm by weight
            cychppm = 1430. * (0.86 + Sqrt(USpower)) * ((USpower)**1.5)&
               &/ ((0.016 + USpower))

!!          Concentration by weight (g/g)
            cychw = cychppm / 1e+06

!!          Concentration by volume
            cychv = (cychw / (cychw + (1. - cychw)*2.65))

!!          concentration in mg/l or metric Tons / m^3
            cych =  cychv * 2.65

          case (4) ! Yang's sand

!!          Shear velocity
            rh = rcharea / (pbed + pbank)
            vsh = Sqrt(9.81 * rh * ch_s(2,jrch))

            var1 = vsh * ch_bed_d50(jrch) * 1e-006 / 1.16e-6

            if (var1 < 70) then
               if (var1 < 1.2) var1 = 1.2
               var2 = 2.5/(Log10(var1) - 0.06)
             else
               var2 = 2.05
            end if

!!          Settling velocity for ch_d50(jrch) particle size
            w50 = 411.0 * ((ch_bed_d50(jrch)/1000.)**2.) / (3600.)

            var3 = w50 * ch_bed_d50(jrch) * 1e-006 / 1.16e-6

            var4 = vsh / w50

            var5 = vc * ch_s(2,jrch) / w50

            var6 = var2 * ch_s(2,jrch)

            var56 = var5 - var6

!!          This needs to be checked for accuracy
            if (var56 <= 0.) var56 = 1e-06

            bedsize = ch_bed_d50(jrch)/1000.
            if (bedsize <= 2.) then
               !!  Yangs sand equation for particles less than 2 mm (2000 mircometer)
               alog10cychppm = 5.435 - 0.286*Log10(var3) - 0.457*Log10(var4)&
                  &+(1.799 - 0.409*Log10(var3) - 0.314*Log10(var4))&
                  &*Log10(var56)
            end if

            if (bedsize > 2.) then
               !!  Yangs gravel equation for particles between 2mm and 10mm
               alog10cychppm = 6.681 - 0.633*Log10(var3) - 4.816*Log10(var4)&
                  &+(2.784 - 0.305*Log10(var3) - 0.282*Log10(var4))&
                  &*Log10(var56)
            end if

!!          Concentration in ppm by weight
            cychppm = 10**alog10cychppm

!!          Concentration by weight
            cychw = cychppm / 1e+06

!!          Concentration by volume
            cychv = (cychw / (cychw + (1. - cychw)*2.65))

!!          concentration in metric Tons / m^3
            cych =  cychv * 2.65

         end select

!!    Incoming sediment concentration
         cyin = sedin / qdin

!!    Potential sediment Transport capacity
         depnet = qdin * (cych - cyin)

         if (depnet .LE. 1.e-6) then
            depnet = 0.
            bnkrte = 0.
            degrte = 0.
         else
            !! First the deposited material will be degraded before channel bed or bank erosion
            if (depnet >= depch(jrch)) then
               !! Effective erosion
               effbnkbed = depnet - depch(jrch)
               !! Effective bank erosion
               if (effbnkbed*bnkrt <= bnkrte) bnkrte = effbnkbed*bnkrt
               bnksan = bnkrte * ch_bnk_san(jrch)
               bnksil = bnkrte * ch_bnk_sil(jrch)
               bnkcla = bnkrte * ch_bnk_cla(jrch)
               bnkgra = bnkrte * ch_bnk_gra(jrch)

               !! Effective bed erosion
               if (effbnkbed*bedrt <= degrte) degrte = effbnkbed*bedrt
               degsan = degrte * ch_bed_san(jrch)
               degsil = degrte * ch_bed_sil(jrch)
               degcla = degrte * ch_bed_cla(jrch)
               deggra = degrte * ch_bed_gra(jrch)

               deg1 = depch(jrch)
               deg1san = depsanch(jrch)
               deg1sil = depsilch(jrch)
               deg1cla = depclach(jrch)
               deg1sag = depsagch(jrch)
               deg1lag = deplagch(jrch)
               deg1gra = depgrach(jrch)

               depch(jrch) = 0.
               depsanch(jrch) = 0.
               depsilch(jrch) = 0.
               depclach(jrch) = 0.
               depsagch(jrch) = 0.
               deplagch(jrch) = 0.
               depgrach(jrch) = 0.

            else

               bnkrte = 0.
               degrte = 0.
               degsan = 0.
               degsil = 0.
               degcla = 0.
               deggra = 0.
               bnksan = 0.
               bnksil = 0.
               bnkcla = 0.
               bnkgra = 0.

               depch(jrch) = depch(jrch) - depnet
               deg1 = depnet

               if (depclach(jrch) >= depnet) then
                  depclach(jrch) = depclach(jrch) - depnet
                  deg1cla = depnet
               else
                  degremain = depnet - depclach(jrch)
                  deg1cla = depclach(jrch)
                  depclach(jrch) = 0.
                  if (depsilch(jrch) >= degremain) then
                     depsilch(jrch) = depsilch(jrch) - degremain
                     deg1sil = degremain
                  else
                     degremain = degremain - depsilch(jrch)
                     deg1sil = depsilch(jrch)
                     depsilch(jrch) = 0.
                     if (depsagch(jrch) >= degremain) then
                        depsagch(jrch) = depsagch(jrch) - degremain
                        deg1sag = degremain
                     else
                        degremain = degremain - depsagch(jrch)
                        deg1sag = depsagch(jrch)
                        depsagch(jrch) = 0.
                        if (depsanch(jrch) >= degremain) then
                           depsanch(jrch) = depsanch(jrch) - degremain
                           deg1san = degremain
                        else
                           degremain = degremain - depsanch(jrch)
                           deg1san = depsanch(jrch)
                           depsanch(jrch) = 0.
                           if (deplagch(jrch) >= degremain) then
                              deplagch(jrch) = deplagch(jrch) - degremain
                              deg1lag = degremain
                           else
                              degremain = degremain - deplagch(jrch)
                              deg1lag = deplagch(jrch)
                              deplagch(jrch) = 0.
                              if (depgrach(jrch) >= degremain) then
                                 depgrach(jrch) = depgrach(jrch) - degremain
                                 deg1gra = degremain
                              else
                                 deg1gra = depgrach(jrch)
                                 depgrach(jrch) = 0.
                              endif
                           endif
                        endif
                     endif
                  endif
               endif

            endif

         end if

         if (depch(jrch) < 1.e-6) then
            depch(jrch) = 0.
            depsanch(jrch) = 0.
            depsilch(jrch) = 0.
            depclach(jrch) = 0.
            depsagch(jrch) = 0.
            deplagch(jrch) = 0.
            depgrach(jrch) = 0.
         end if

!! Deposition calculated based on Einstein Equation
         xx = 1.055 * 1000. * ch_l(2,jrch) / (vc * rchdep)

!! Gravel deposition
         x = Min(20., xx * vgra)
         pdep = Min((1. - Exp(-x)), 1.)
         depgra = grain * pdep
         dep = depgra

!! sand deposition
         x = Min(20., xx * vsan)
         pdep = Min((1. - Exp(-x)), 1.)
         depsan = sanin * pdep
         dep = dep + depsan

!! Silt deposition
         x = Min(20., xx * vsil)
         pdep = Min((1. - Exp(-x)), 1.)
         depsil = silin * pdep
         dep = dep + depsil

!! Clay deposition
         x = Min(20., xx * vcla)
         pdep = Min((1. - Exp(-x)), 1.)
         depcla = clain * pdep
         dep = dep + depcla

!! Small aggregates deposition
         x = Min(20., xx * vsag)
         pdep = Min((1. - Exp(-x)), 1.)
         depsag = sagin * pdep
         dep = dep + depsag

!! Large aggregates deposition
         x = Min(20., xx * vlag)
         pdep = Min((1. - Exp(-x)), 1.)
         deplag = lagin * pdep
         dep = dep + deplag

!!    Particles deposited on Floodplain (only silt and clay type particles)
         depfp(jrch)    = depfp(jrch) + (depsil + depcla) * fpratio
         depsilfp(jrch) = depsilfp(jrch) + depsil * fpratio
         depclafp(jrch) = depclafp(jrch) + depcla * fpratio

!!    Remaining is deposited in the channel
         depch(jrch)    = depch(jrch)    + dep - (depsil + depcla) * fpratio
         depsilch(jrch) = depsilch(jrch) + depsil * (1. - fpratio)
         depclach(jrch) = depclach(jrch) + depcla * (1. - fpratio)
         depsanch(jrch) = depsanch(jrch) + depsan
         depsagch(jrch) = depsagch(jrch) + depsag
         deplagch(jrch) = deplagch(jrch) + deplag
         depgrach(jrch) = depgrach(jrch) + depgra

         sedin  = sedin + degrte + bnkrte + deg1    - dep
         grain  = grain + deggra + bnkgra + deg1gra - depgra
         sanin  = sanin + degsan + bnksan + deg1san - depsan
         silin  = silin + degsil + bnksil + deg1sil - depsil
         clain  = clain + degcla + bnkcla + deg1cla - depcla
         sagin  = sagin + deg1sag - depsag
         lagin  = lagin + deg1lag - deplag

         if (sedin  < 1.e-6) then
            sedin = 0.
            sanin = 0.
            silin = 0.
            clain = 0.
            sagin = 0.
            lagin = 0.
            grain = 0.
         end if

         outfract = rtwtr / qdin
         if (outfract > 1.) outfract = 1.

         sedrch =  sedin * outfract
         rch_san = sanin * outfract
         rch_sil = silin * outfract
         rch_cla = clain * outfract
         rch_sag = sagin * outfract
         rch_lag = lagin * outfract
         rch_gra = grain * outfract

         if (sedrch  < 1.e-6) then
            sedrch = 0.
            rch_san = 0.
            rch_sil = 0.
            rch_cla = 0.
            rch_sag = 0.
            rch_lag = 0.
            rch_gra = 0.
         end if

         sedst(jrch) = sedin - sedrch
         sanst(jrch) = sanin - rch_san
         silst(jrch) = silin - rch_sil
         clast(jrch) = clain - rch_cla
         sagst(jrch) = sagin - rch_sag
         lagst(jrch) = lagin - rch_lag
         grast(jrch) = grain - rch_gra

         if (sedst(jrch) < 1.e-6) then
            sedst(jrch) = 0.
            sanst(jrch) = 0.
            silst(jrch) = 0.
            clast(jrch) = 0.
            sagst(jrch) = 0.
            lagst(jrch) = 0.
            grast(jrch) = 0.
         end if

!!    Bank erosion
         rchdy(55,jrch) = bnkrte
!!    Channel Degredation
         rchdy(56,jrch) = degrte
!!    Channel Deposition (Only new deposits during the current time step)
         if (depch(jrch) >= depprch(jrch)) then
            rchdy(57,jrch) = depch(jrch) - depprch(jrch)
         else
            rchdy(57,jrch) = 0.
         end if
!!    Floodplain Deposition (Only new deposits during the current time step)
         if (depfp(jrch) >= depprfp(jrch)) then
            rchdy(58,jrch) = depfp(jrch) - depprfp(jrch)
         else
            rchdy(58,jrch) = 0.
         end if
!!    Total suspended sediments (only silt and clay)
         rchdy(59,jrch) = (rch_sil + rch_cla)/rtwtr * 1.e6

!!    Deposition during the previous time step
         depprch(jrch) = depch(jrch)  !! Channel
         depprfp(jrch) = depfp(jrch)  !! Flood plain

!!    Organic nitrogen and Organic Phosphorus contribution from channel erosion
!!    Only bank erosion is assumed to contribute to channel erosion
         ch_orgn(jrch) = bnkrte * ch_onco(jrch) / 1000.
         ch_orgp(jrch) = bnkrte * ch_opco(jrch) / 1000.

!! compute changes in channel dimensions
         if (ideg == 1) then
            depdeg = ch_d(jrch) - ch_di(jrch)
            if (depdeg < ch_si(jrch) * ch_li(jrch) * 1000.) then
               if (qdin > 1400000.) then
                  dat2 = 358.6 * rchdep * ch_s(2,jrch) * ch_cov(1,jrch)
                  ch_d(jrch) = ch_d(jrch) + dat2
                  ch_w(2,jrch) = ch_wdr(jrch) * ch_d(jrch)
                  ch_s(2,jrch) = ch_s(2,jrch) - dat2 / (ch_l(2,jrch) * 1000.)
                  ch_s(2,jrch) = Max(.0001, ch_s(2,jrch))
                  call ttcoef(jrch)
               endif
            endif
         endif

      else

         sedst(jrch) = sedin
         sanst(jrch) = sanin
         silst(jrch) = silin
         clast(jrch) = clain
         sagst(jrch) = sagin
         lagst(jrch) = lagin
         grast(jrch) = grain

      end if !! end of qdin > 0.01 loop

   end if  !! end of rtwtr and rchdep > 0 loop

   return
end
