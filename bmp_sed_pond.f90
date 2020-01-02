!> @file bmp_sed_pond.f90
!> file containing the subroutine bmp_sed_pond
!> @author
!> modified by Javier Burguete

!> this subroutine routes water and sediment through a sedimentation pond in
!> the subbasin
!> param[in] sb subbasin or reach number
!> param[in] kk filter id number in the subbasin
!> param[inout] flw stormwater runoff coming in/out of pond at a time step
!> param[inout] sed overland flow sediment coming in/out of pond at a time
!> step
subroutine bmp_sed_pond(sb, kk, flw, sed)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sb          |none          |subbasin or reach number
!!    kk          |none          |filter id number in the subbasin
!!    hru_sub(:)  |none          |subbasin in which HRU/reach is located
!!    i_mo        |none          |current month of simulation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    pnd_sed(:)  |kg/L          |ratio of sediment to water in pond
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dp
!!    hpnd
!!    hweir
!!    ii
!!    ksat
!!    ksed
!!    mxh
!!    mxvol
!!    pdia
!!    phead
!!    pipeflow
!!    qevap
!!    qin
!!    qout
!!    qpipe
!!    qpnd
!!    qpndi
!!    qtrns
!!    qweir
!!    sedconcweir
!!    sedpnde
!!    sedpndi
!!    splw
!!    spndconc
!!    sub_ha
!!    sweir
!!    td
!!    tsa
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max, Min, Exp
!!    SWAT: Pipeflow

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 Pipeflow
   integer, intent(in) :: sb, kk
   real*8, dimension(4,0:nstep), intent(inout) :: flw, sed
   real*8 :: dp, hpnd, hweir, ksat, ksed, mxh, mxvol, pdia, phead, qevap, qin,&
      &qout, qpipe, qpnd, qpndi, qtrns, qweir, sedconcweir, sedpnde, sedpndi,&
      &splw, spndconc, sub_ha, sweir, td, tsa
   integer :: ii

   sub_ha = da_ha * sub_fr(sb)
   qin = 0.; qout = 0.; sedpndi = 0.; sedpnde = 0.
   flw(2,:) = 0.; sed(2,:) = 0.

   !! Initialize parameters, coefficients, etc
   tsa = sp_sa(sb,kk)     !total surface area of pond (m^2)
   mxvol = sp_pvol(sb,kk)    !max. capacity of the basin (m^3)
   mxh = sp_pvol(sb,kk) / sp_sa(sb,kk) !max. depth of water, m
   splw = sp_bpw(sb,kk)       !spillway overflow weir width (m)
   pdia = sp_pd(sb,kk)       !outflow orifice pipe diameter (mm)
   ksat = sp_k(sb,kk)      !saturated hydraulic conductivity (mm/hr)
   dp = sp_dp(sb,kk) * 1000. !median particle size of TSS, micrometer

   !! Get initial values from previous day
   qpnd = sp_qi(sb,kk) !m^3
   spndconc = sp_sedi(sb,kk)
   qevap = 0

   do ii=1,nstep
      qweir = 0.; qtrns = 0.; qpipe = 0.
      !inflow = runoff + precipitation
      qin = flw(1,ii) * 10. * (sub_ha - tsa / 10000.) +&
         &precipdt(ii) * tsa / 1000.  !m^3
      qpndi = qpnd + qin !m^3
      hpnd = qpndi / tsa !ponding depth on the filter surface, m

      if (qin<0.001.and.qpndi<0.001)then
         !No inflow/outflow during this time step
         flw(2,ii) = 0.
      else
         !spillway overflow
         if (hpnd > mxh) then
            hweir = Max(0.,hpnd - mxh)  !water depth over weir crest, m
            !weir overflow, m^3
            qweir = 3.33 * splw * hweir ** 1.5 * idt * 60.
            hpnd = Max(0.,(qpndi - qweir) / tsa) !m
            !overflow amount is no larger than surplus water above spillway height
            if (qweir>qpndi-mxvol) then
               qweir = Max(0.,qpndi - mxvol) !weir outflow
               hpnd = mxh
            end if
            qpnd = Max(0.,qpndi - qweir)
         else
            qpnd = qpndi
         end if

         !Transmission loss through infiltration
         qtrns = ksat * tsa / 1000./ 60. * idt
         qpnd = qpnd - qtrns
         bmp_recharge(sb) = bmp_recharge(sb)&
            &+ qtrns / (sub_ha * 10000. - tsa) * 1000.

         if (qpnd<0) then
            qpnd = 0.
            qtrns = 0.
         endif

         !Evapotranspiration loss
         qevap = tsa * sub_etday(sb) / 1000. / 1440. * idt !m^3
         if(qevap<1e-6) qevap = 0.
         qpnd = qpnd - qevap
         if (qpnd<0) then
            qpnd = 0.
            qevap = 0.
         endif

         !Outflow through orifice pipe
         hpnd = qpnd / tsa  !m
         phead = Max(0.,hpnd * 1000. - pdia / 2.)  !mm
!            if (phead>pdia/2.) then
         qpipe = pipeflow(pdia,phead) * idt *60. !m^3
!            else
!               qpipe = qout *  0.9
!            endif

         !update out flow, m^3
         qout = qpipe
         qpnd = Max(0.,qpnd - qout)

         !outflow normalized to subbasin area, mm
         flw(1,ii) = qin / ((sub_ha - tsa / 10000.) *10.)
         flw(2,ii) = qout / ((sub_ha - tsa / 10000.) *10.) !mm
         flw(3,ii) = qweir / ((sub_ha - tsa / 10000.) *10.) !mm
         flw(4,ii) = qtrns / (sub_ha *10000. - tsa) * 1000.  !mm
      endif

      !---------------------------------------------------------
      !! TSS removal
      sweir = 0.

      !Sediment bypass through spillway overflow, tons
      if (qweir>0.001) then
         if(qin>0.001) then
            sedconcweir = sed(1,ii) / qin !tons/m3
         else
            sedconcweir = spndconc !tons/m3
         endif
         sweir = sedconcweir * qweir !tons
         sedpndi = sedpnde + sed(1,ii) - sweir !tons
         spndconc = sedpndi / qpndi !tons/m3
      else

         if(qpndi>0.001) then
            sedpndi = sedpnde + sed(1,ii)  !tons
            spndconc = sedpndi / qpndi ! tons/m3
         else
            sedpndi = 0
            spndconc = 0
         endif

      endif

      ! sediment conc at the beginning of the time step
      spndconc = spndconc * 1.e6 !mg/l

      !Estimate TSS removal due to sedimentation
      if (spndconc>sp_sede(sb,kk)) then
         ksed = Min(134.8,41.1 * hpnd ** -(0.999))  !decay coefficient, Huber et al. 2006
         td = qpnd / qpipe / nstep !detention time, day
         spndconc = (spndconc - sp_sede(sb,kk)) * Exp(-ksed * td) +&
            &sp_sede(sb,kk)
      endif

      !Sediment coming out of the pond
      sed(2,ii) = spndconc * qpipe * 1.e-6  !tons
      sedpnde = spndconc * qpnd * 1.e-6 !tons
      sed(3,ii) = sweir

      ! total sediment removed from the pond, tons
      sp_sed_cumul(sb,kk) = sp_sed_cumul(sb,kk) + sedpndi - sedpnde

   end do

   ! Store end-of-day values for next day
   sp_qi(sb,kk) = qpnd
   sp_sedi(sb,kk) = spndconc


   return
end subroutine bmp_sed_pond
