!> @file urbanhr.f90
!> file containing the subroutine urbanhr
!> @author
!> modified by Javier Burguete

!> this subroutine computes loadings from urban areas using the
!> a build-up/wash-off algorithm at subdaily time intervals
!> @param[in] j HRU number (none)
subroutine urbanhr(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j            |none           |HRU number
!!    al5          |none           |fraction of daily rainfall that occurs
!!                                 |during 0.5h highest intensity
!!    curbden(:)   |km/ha          |curb length density in HRU
!!    dirtmx(:)    |kg/curb km     |maximum amount of solids allowed to
!!                                 |build up on impervious surfaces
!!    fimp(:)      |fraction       |fraction of HRU area that is
!!                                 |impervious (both directly and
!!                                 |indirectly connected)
!!    hru_km(:)    |km^2           |area of HRU in square kilometers
!!    iida         |julian date    |day being simulated (current julian date)
!!    isweep(:,:,:)|julian date    |date of street sweeping operation
!!    iurban(:)    |none           |urban simulation code:
!!                                 |0  no urban sections in HRU
!!                                 |1  urban sections in HRU, simulate using
!!                                 |   USGS regression equations
!!                                 |2  urban sections in HRU, simulate using
!!                                 |   build up/wash off algorithm
!!    nro(:)       |none           |sequence number of year in rotation
!!    nsweep(:)    |none           |sequence number of street sweeping operation
!!                                 |within the year
!!    peakr        |m^3/s          |peak runoff rate
!!    precipday    |mm H2O         |precipitation for the day in HRU
!!    sedorgn(:)   |kg N/ha        |amount of organic nitrogen in surface runoff
!!                                 |in HRU for the day
!!    sedorgp(:)   |kg P/ha        |amount of organic phosphorus in surface
!!                                 |runoff in HRU for the day
!!    sedyld(:)    |metric tons    |daily soil loss caused by water erosion
!!    surfq(:)     |mm H2O         |surface runoff for the day in HRU
!!    surqno3(:)   |kg N/ha        |amount of NO3-N in surface runoff in HRU for
!!                                 |the day
!!    surqsolp(:)  |kg P/ha        |amount of soluble phosphorus in surface
!!                                 |runoff in HRU for the day
!!    tconc(:)     |hr             |time of concentration
!!    thalf(:)     |days           |time for the amount of solids on
!!                                 |impervious areas to build up to 1/2
!!                                 |the maximum level
!!    tnconc(:)    |mg N/kg sed    |concentration of total nitrogen in
!!                                 |suspended solid load from impervious
!!                                 |areas
!!    tno3conc(:)  |mg NO3-N/kg sed|concentration of NO3-N in suspended
!!                                 |solid load from impervious areas
!!    tpconc(:)    |mg P/kg sed    |concentration of total phosphorus in
!!                                 |suspended solid load from impervious
!!                                 |areas
!!    twash(:)     |days           |time that solids have built-up on streets
!!    urbcoef(:)   |1/mm           |wash-off coefficient for removal of
!!                                 |constituents from an impervious surface
!!    ubntss(:)    |metric tons    |TSS loading from urban impervious cover
!!    urblu(:)     |none           |urban land type identification number from
!!                                 |urban database
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sedorgn(:)  |kg N/ha       |amount of organic nitrogen in surface runoff
!!                               |in HRU for the day
!!    sedorgp(:)  |kg P/ha       |amount of organic phosphorus in surface runoff
!!                               |in HRU for the day
!!    sedyld(:)   |metric tons   |daily soil loss caused by water erosion
!!    surqno3(:)  |kg N/ha       |amount of NO3-N in surface runoff in HRU for
!!                               |the day
!!    surqsolp(:) |kg P/ha       |amount of soluble phosphorus in surface runoff
!!                               |in HRU for the day
!!    twash(:)    |
!!    ubntss(:)   |metric tons   |TSS loading from urban impervious cover
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dirt        |kg/curb km    |amount of solids built up on impervious
!!                               |surfaces
!!    dirto       |kg/ha         |amount of solids built up on impervious
!!                               |surfaces at the beginning of time step
!!    sus_sol     |kg            |suspended solid loading in surface runoff
!!                               |from urban area
!!    qdt
!!    tn          |kg            |total nitrogen in surface runoff from
!!                               |urban area
!!    tno3
!!    tp          |kg            |total phosphorus in surface runoff from
!!                               |urban area
!!    urbk        |1/hr          |
!!    xx
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max
!!    SWAT: sweep

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: dirt, dirto, qdt, sus_sol, tn, tno3, tp, urbk, xx
   integer :: k

   do k = 1, nstep

      select case (iurban(j))

       case (1,2)       !! build-up/wash-off algorithm

         !! rainy day: no build-up, street cleaning allowed

         qdt = ubnrunoff(k) * 60./ dfloat(idt) !urban runoff in mm/hr
         if (qdt > 0.025 .and. surfq(j) > 0.1) then   ! SWMM : 0.001 in/hr (=0.0254mm/hr)

            !! calculate amount of dirt on streets prior to wash-off
            dirto = dirtmx(urblu(j)) * twash(j) / (thalf(urblu(j)) + twash(j))

            !! calculate wash-off of solids
            urbk = urbcoef(urblu(j)) * qdt

            dirt = dirto * Exp (- urbk * dfloat(idt) / 60.)
            if (dirt < 1.e-6) dirt = 0.0

            !! set time to correspond to lower amount of dirt
            twash(j) = thalf(urblu(j)) * dirt / (dirtmx(urblu(j)) - dirt)
            !! amounts are kg/ha
            sus_sol = Max(0., (dirto - dirt) * curbden(urblu(j)))
            tn = tnconc(urblu(j)) * sus_sol / 1.e6
            tp = tpconc(urblu(j)) * sus_sol / 1.e6
            tno3 = tno3conc(urblu(j)) * sus_sol / 1.e6

            ubntss(k) = (.001 * sus_sol * hru_ha(j)) * fimp(urblu(j))
            xx = 1. - fimp(urblu(j))
            surqno3(j) = tno3 * fimp(urblu(j)) + surqno3(j) * xx
            sedorgn(j) = (tn - tno3) * fimp(urblu(j)) + sedorgn(j) * xx
            sedorgp(j) = .75 * tp * fimp(urblu(j)) + sedorgp(j) * xx
            surqsolp(j) = .25 * tp * fimp(urblu(j)) + surqsolp(j) * xx
         else
            !! no surface runoff
            twash(j) = twash(j) + idt / 1440.

            !! perform street sweeping
            !! isweep is not initialized in any part, assumed to be 0
!            if (isweep(j) > 0 .and. iida >= isweep(j)) then
!               call sweep(j)
!            else if (phusw(j) > 0.0001) then
!               if (igro(j) == 0) then
!                  if (phubase(j) > phusw(j)) call sweep(j)
!               else
!                  if (phuacc(j) > phusw(j)) call sweep(j)
!               end if
!            end if
            if (phusw(j) > 0.0001) then
               if (igro(j) == 0) then
                  if (phubase(j) > phusw(j)) call sweep(j)
               else
                  if (phuacc(j) > phusw(j)) call sweep(j)
               end if
            end if

         end if
      end select
      sus_sol=0

      ! Compute evaporation of water (initial dabstraction) from impervious cover
      !init_abstrc(j) = init_abstrc(j) - etday / nstep ! not used
      !init_abstrc(j) = Max(0.,init_abstrc(j)) ! not used
   end do

   !! perform street sweeping
   if(surfq(j) < 0.1) then
      !! isweep is not initialized in any part, assumed to be 0
      !if (isweep(j) > 0 .and. iida >= isweep(j)) then
      !   call sweep(j)
      !else if (phusw(j) > 0.0001) then
      if (phusw(j) > 0.0001) then
         if (igro(j) == 0) then
            if (phubase(j) > phusw(j)) then
               call sweep(j)
            endif
         else
            if (phuacc(j) > phusw(j)) then
               call sweep(j)
            endif
         end if
      end if
   end if

   return
end
