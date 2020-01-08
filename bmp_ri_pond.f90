!> @file bmp_ri_pond.f90
!> file containing the subroutine bmp_ri_pond
!> @author
!> modified by Javier Burguete

!> this subroutine routes water through a retention irrigation pond in the
!> subbasin
!> param[in] sb subbasin or reach number
!> param[in] kk pond id number in the subbasin
!> param[inout] riflw stormwater runoff coming in/out of pond at a time step
!> param[inout] rised overland flow sediment coming in/out of pond at a time
!> step
subroutine bmp_ri_pond(sb, kk, riflw, rised)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sb          |none          |subbasin or reach number
!!    kk          |none          |pond id number in the subbasin
!!    riflw(:,:)  |mm            |stormwater runoff coming in/out of pond at a time step
!!    rised(:,:)  |mm            |overland flow sediment coming in/out of pond at a time step
!!    ri_sed(:)   |tons          |total sediment deposited in the pond
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ri_sed(:)   |tons          |total sediment deposited in the pond
!!    riflw(:,:)  |mm            |stormwater runoff coming in/out of pond at a time step
!!    rised(:,:)  |mm            |overland flow sediment coming in/out of pond at a time step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hpnd
!!    ii          |none          |time step
!!    inflw
!!    insed
!!    ksat
!!    ksed
!!    mxh
!!    mxvol
!!    outflw
!!    outsed
!!    qet
!!    qin
!!    qout
!!    qpnd
!!    qpump
!!    qseep
!!    sedconc
!!    sedpnde
!!    sedpndi
!!    sedpump
!!    sub_ha
!!    td
!!    tsa
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max, Min, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: sb, kk
   real*8, dimension(4,0:nstep), intent(inout) :: riflw, rised
   real*8, dimension(0:nstep) :: inflw, insed, outflw, outsed
   real*8 :: hpnd, ksat, ksed, mxh, mxvol, qet, qin, qout, qpnd, qpump, qseep,&
      &sedconc, sedpnde, sedpndi, sedpump, sub_ha, td, tsa
   integer :: ii

   sub_ha = da_ha * sub_fr(sb)
   outflw = 0.; outsed = 0.
   qseep = 0.; qet = 0.;

   !! Initialize parameters, coefficients, etc
   tsa = ri_sa(sb,kk)     !total surface area of pond (m^2)
   mxvol = ri_vol(sb,kk)    !max. capacity of the basin (m^3)
   mxh = ri_dep(sb,kk)  !max. depth of water, m
   ksat = ri_k(sb,kk)      !saturated hydraulic conductivity (mm/hr)

   !! Get initial values from previous day
   qpnd = ri_qi(sb,kk) !m^3
   sedpnde =  ri_sedi(sb,kk)
   inflw(:) = riflw(1,:)
   insed(:) = rised(1,:)


   do ii=1,nstep

      qout = 0.

      !inflow = runoff + precipitation
      qin = inflw(ii) * 10. * (sub_ha - tsa / 10000.) +&
         &precipdt(ii) * tsa / 1000.  !m^3

      !update ponded water volume
      qpnd = qpnd + qin

      !bypass flow when pond is full
      if (qpnd>mxvol) then
         qout = qpnd - mxvol !m3
         outflw(ii) = qout
         outsed(ii) = insed(ii) * qout / qin !tons
         qpnd = mxvol
      end if

      !initial sediment
      sedpndi = sedpnde + insed(ii)
      if (qpnd>0) then
         sedconc = sedpndi / qpnd * 1.e6 !mg/l
      else
         sedconc = 0.
      endif

      !Transmission loss through infiltration
      qseep = ksat * tsa / 1000./ 60. * idt !m^3
      bmp_recharge(sb) = bmp_recharge(sb)&
         &+ qseep / (sub_ha * 10000. - tsa) * 1000.

      !Evapotranspiration loss
      qet = ri_evrsv(sb,kk) * tsa * pet_day / 1000. / 1440. * idt !m^3

      !water pumped for irrigation
      qpump =  Max(0.,ri_pmpvol(kk,ii))

      !mass balance
      qpnd = qpnd - qseep - qet - qpump
      if (qpnd<0) qpnd = 0
      hpnd = qpnd / (mxvol / mxh)

      !Estimate TSS removal due to sedimentation
      if (sedconc>12.) then ! assume 12mg/l as equilibrium concentration, , Huber et al. 2006
         ksed = Min(134.8,41.1 * hpnd ** (-0.999))  !decay coefficient, Huber et al. 2006
         td = 1. / nstep !detention time, day
         sedconc = (sedconc - 12.) * Exp(-ksed * td) + 12.
      endif

      !sediment pumped
      sedpump = qpump * sedconc / 1.e6 !tons

      !sediment deposition
      sedpnde = qpnd *sedconc / 1.e6 !tons

      !ri_sed(sb,kk) seems to be always 0
      !ri_sed_cumul(sb,kk) = ri_sed(sb,kk) + sedpndi - sedpnde
      !ri_sed_cumul(sb,kk) = sedpndi - sedpnde ! not used

      riflw(1,ii) = qin / (sub_ha *10000. - tsa) * 1000.  !mm
      riflw(2,ii) = outflw(ii) / (sub_ha *10000. - tsa) * 1000.
      riflw(3,ii) = qpump / (sub_ha *10000. - tsa) * 1000.
      riflw(4,ii) = qseep / (sub_ha *10000. - tsa) * 1000.
      rised(3,:) = sedpump
      rised(2,:) = outsed(:)
   End do


   ! Store end-of-day values for next day
   ri_qi(sb,kk) = qpnd
   ri_qloss(kk,1) = qet
   ri_qloss(kk,2) = qseep
   ri_sedi(sb,kk) = sedpnde

   return
end subroutine bmp_ri_pond
