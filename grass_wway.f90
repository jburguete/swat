!> @file grass_wway.f90
!> file containing the subroutine grass_wway
!> @author
!> modified by Javier Burguete

!> this subroutine controls the grass waterways
!> @param[in] j HRU number (none)
subroutine grass_wway(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    surfq(:)        |mm H2O        |amount of water in surface runoff generated
!!    grwat_d(:)      |m             |Depth of Grassed waterway
!!    grwat_i(:)      |none          |On/off Flag for waterway simulation
!!    grwat_l(:)      |km            |Length of Grass Waterway
!!    grwat_n(:)      |none          |Mannings's n for grassed waterway
!!    grwat_s(:)      |m/m           |Slope of grass waterway
!!    grwat_spcon(:)  |none          |sediment transport coefficant defined by user
!!    grwat_w(:)      |none          |Width of grass waterway
!!    tc_gwat(:)      |none          |Time of concentration for Grassed waterway and its drainage area
!!    mhru
!!    peakr           |m^3/s         |peak runoff rate for the day
!!    rcharea         |m^2           |cross-sectional area of flow
!!    rchdep          |m             |depth of flow on day
!!    sedyld(:)       |metric tons   |daily soil loss caused by water erosion
!!    surfq(:)        |mm H2O        |surface runoff generated on day in HRU
!!    sedyld(:)       |metric tons   |daily soil loss caused by water erosion

!!    wat_phi1(:)     |m^2           |cross-sectional area of flow at bankfull
!!                                   |depth
!!    wat_phi5(:)     |m^3/s         |flow rate when reach is at bankfull depth
!!    wat_phi6(:)     |m             |bottom width of main channel
!!    wat_phi7(:)     |m             |depth of water when reach is at bankfull
!!                                   |depth
!!    wat_phi8(:)     |m/s           |average velocity when reach is at
!!                                   |bankfull depth
!!    wat_phi9(:)     |m/s           |wave celerity when reach is at
!!                                   |bankfull depth
!!    wat_phi10(:)    |hr            |storage time constant for reach at
!!                                   |bankfull depth (ratio of storage to
!!                                   |discharge)
!!    wat_phi13(:)    |m/s           |average velocity when reach is at
!!                                   |0.1 bankfull depth (low flow)
!!                                   |of storage to discharge)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sedyld(:)       |metric tons   |daily soil loss caused by water erosion
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chflow_m3   |m^3/s         |runoff in CMS
!!    chflow_day  |m^3/day       |runoff
!!    cych
!!    cyin
!!    k           |none          |total number of HRUs plus this HRU number
!!    depnet
!!    p
!!    rh
!!    sed_frac
!!    sed_remove  |%             |percent of sediment capture in VFS
!!    sedin       |mg            |sediment in waterway
!!    sedint      |mg            |sediment into waterway channel
!!    sedout      |mg            |sediment out of waterway channel
!!    sedtrap
!!    sf_area     |m^2           |area of waterway sides in sheetflow
!!    sf_depth
!!    sf_sed      |kg/m^2        |sediment loads on sides of waterway
!!    surq_frac
!!    surq_remove |%             |percent of surface runoff capture in VFS
!!    vc          |m/s           |flow velocity in reach
!!    xrem
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Sqrt, Log, Max
!!    SWAT: Qman

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 Qman
   integer, intent(in) :: j
   real*8 :: chflow_day, chflow_m3, cych, cyin, depnet, p, rh, sed_frac,&
      &sed_remove, sedin, sedint, sedout, sedtrap, sf_area, sf_depth, sf_sed,&
      &surq_frac, surq_remove, vc, xrem, xx
   integer :: k

!! do this only if there is surface runoff this day
   if (surfq(j) > 0.001) then

!!        compute channel peak rate using SCS triangular unit hydrograph
!!  Calculate average flow based on 3 hours of runoff
      chflow_day = 1000. * surfq(j) * hru_km(j)
      chflow_m3 = chflow_day/10800
      peakr = 2. * chflow_m3 / (1.5 * tc_gwat(j))

!! if peak rate is greater than bankfull discharge
      if (peakr > wat_phi5(j)) then
         rcharea = wat_phi1(j)
         rchdep = grwat_d(j)
      else
!!          find the crossectional area and depth for todays flow
!!          by iteration method at 1cm interval depth
!!          find the depth until the discharge rate is equal to volrt
         sdti = 0.
         rchdep = 0.

         Do While (sdti < peakr)
            rchdep = rchdep + 0.01
            rcharea = (wat_phi6(j) + 8 * rchdep) * rchdep
            p = wat_phi6(j) + 2. * rchdep * Sqrt(1. + 8 * 8)
            rh = rcharea / p
            sdti = Qman(rcharea, rh, grwat_n(j), grwat_s(j))
         end do
      end if

!!        Sediment yield (kg) from fraction of area drained by waterway

      sedin = sedyld(j)
!! Calculate sediment losses in sheetflow at waterway sides

!! calculate area of sheeflow in m^2 assumne *:1 side slope 8.06 = (8^2+1^2)^.5
      sf_area = (grwat_d(j) - rchdep) * 8.06 * grwat_l(j) * 1000
!! Adjust Area to account for flow nonuniformities White and Arnold 2009 found half of flow in VFS
!!handled by 10% of VFS area. Waterways likely even more concentrated Assume only 20% of sideslope acts as filters
      sf_depth = 0.
      sf_sed = 0.
      if (sf_area > 1.e-6) then
         sf_area = sf_area * 0.20
!! calculate runoff depth over sheetflow area in mm
         sf_depth=surfq(j)  * hru_km(j) * 1000000/sf_area
!! Calculate sediment load on sheetflow area kg/ha
         sf_sed = sedin * 1000 / sf_area
!! Calculate runoff and sediment losses taken from mostly from filter.f
      end if

      if (sf_area > 0.) then
!!  surq_remove = 75.8 - 10.8 * Log(sf_depth) + 25.9
!!     &    * Log(sol_k(1,j))
         !! Simpler form derived from vfsmod simulations. r2 = 0.57 Publication pending white and arnold 2008

         surq_remove = 95.6 - 10.79 * Log(sf_depth)
         if (surq_remove > 100.) surq_remove = 100.
         if (surq_remove < 0.) surq_remove = 0.

         sed_remove = 79.0 - 1.04 * sf_sed + 0.213 * surq_remove
         if (sed_remove > 100.) sed_remove = 100.
         if (sed_remove < 0.) sed_remove = 0.

      Else
         sed_remove = 0
         surq_remove = 0
      endif
      sedint = sedin * (1. - sed_remove / 100.)

!!        calculate flow velocity
      if (rcharea > 1.e-4) then
         vc = peakr / rcharea
         if (vc > wat_phi9(j)) vc = wat_phi9(j)
      else
         vc = 0.001
      end if

!!        compute deposition in the waterway
!! if there is significant flow calculate
      if (chflow_m3 > 1.e-4) then
!! Calculate sediment concentration in inflow mg/m^3
         cyin = sedint / chflow_day
!! Calculate sediment transport capacity mg/m^3
         cych = grwat_spcon(j) * vc ** 1.5
!! Calculate deposition in mg
         depnet = chflow_day * (cyin - cych)
         if (depnet < 0.) depnet = 0
         if (depnet > sedint) depnet = sedint
      else
         cyin = 0.
         cych = 0.
         depnet = 0.
      endif
!! Calculate sediment out of waterway channel
      sedout = sedint - depnet

!! Calculate total fraction of sediment and surface runoff transported
      if (sedyld(j) < .0001) sedyld(j) = .0001
      sed_frac =  sedout/sedyld(j)


      surq_frac = 1 - surq_remove/100

!! Subtract reductions from sediment, nutrients, bacteria, and pesticides NOT SURFACE RUNOFF to protect water balance
      sedtrap = sedyld(j) * (1. - sed_frac)
      sedyld(j) = sedyld(j) * sed_frac
      sedminpa(j) = sedminpa(j) * sed_frac
      sedminps(j) = sedminps(j) * sed_frac
      sedorgp(j) = sedorgp(j) * sed_frac
      surqsolp(j) = surqsolp(j) * surq_frac
      sedorgn(j) = sedorgn(j) * sed_frac
      surqno3(j) = surqno3(j) * surq_frac

      if (sedtrap <= lagyld(j)) then
         lagyld(j) = lagyld(j) - sedtrap
      else
         xrem = sedtrap - lagyld(j)
         lagyld(j) = 0.
         if (xrem <= sanyld(j)) then
            sanyld(j) = sanyld(j) - xrem
         else
            xrem = xrem - sanyld(j)
            sanyld(j) = 0.
            if (xrem <= sagyld(j)) then
               sagyld(j) = sagyld(j) - xrem
            else
               xrem = xrem - sagyld(j)
               sagyld(j) = 0.
               if (xrem <= silyld(j)) then
                  silyld(j) = silyld(j) - xrem
               else
                  xrem = xrem - silyld(j)
                  silyld(j) = 0.
                  if (xrem <= clayld(j)) then
                     clayld(j) = clayld(j) - xrem
                  else
                     clayld(j) = 0.
                  end if
               end if
            end if
         end if
      end if
      sanyld(j) = Max(0., sanyld(j))
      silyld(j) = Max(0., silyld(j))
      clayld(j) = Max(0., clayld(j))
      sagyld(j) = Max(0., sagyld(j))
      lagyld(j) = Max(0., lagyld(j))

!! Calculate pesticide removal
!! based on the sediment and runoff removal only
      if (hrupest(j) == 1) then
         xx = 1. - sed_remove / 100.
         do k = 1, npmx
            pst_surq(k,j) = pst_surq(k,j) * surq_frac
            pst_sed(k,j) = pst_sed(k,j) * xx
         end do
      end if
!! compute bacteria reductions
      bactrop = bactrop * surq_frac
      bactrolp = bactrolp * surq_frac
      bactsedp = bactsedp * sed_frac
      bactsedlp = bactsedlp * sed_frac

   Endif
   return
end
