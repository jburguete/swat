!> @file orgncswat2.f90
!> file containing the subroutine orgncswat2
!> @author
!> modified by Javier Burguete

!> this subroutine calculates the amount of organic nitrogen removed in
!> surface runoff - when using CSWAT==2 it
!> @param[in] iwave
!> flag to differentiate calculation of HRU and subbasin sediment calculation
!> (none)\n
!> iwave = 0 for HRU\n
!> iwave = subbasin # for subbasin
!> @param[in] j HRU number
subroutine orgncswat2(iwave, j)


!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    iwave         |none         |flag to differentiate calculation of HRU and
!!                                |subbasin sediment calculation
!!                                |iwave = 0 for HRU
!!                                |iwave = subbasin # for subbasin
!!    j             |none         |HRU number
!!    da_ha         |ha           |area of watershed in hectares
!!    enratio       |none         |enrichment ratio calculated for day in HRU
!!    erorgn(:)     |none         |organic N enrichment ratio, if left blank
!!                                |the model will calculate for every event
!!    sedyld(:)     |metric tons  |daily soil loss caused by water erosion in
!!                                |HRU
!!    sol_bd(:,:)   |Mg/m**3      |bulk density of the soil
!!    sol_z(:,:)    |mm           |depth to bottom of soil layer
!!    sub_bd(:)     |Mg/m^3       |bulk density in subbasin first soil layer
!!    sub_fr(:)     |none         |fraction of watershed area in subbasin
!!    sub_orgn(:)   |kg N/ha      |amount of nitrogen stored in all organic
!!    sedc_d(:)     |kg C/ha      |amount of C lost with sediment
!!                                |pools
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sedorgn(:)    |kg N/ha      |amount of organic nitrogen in surface runoff
!!                                |in HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    conc        |              |concentration of organic N in soil
!!    er          |none          |enrichment ratio
!!    wt1         |none          |conversion factor (mg/kg => kg/ha)
!!    xx          |kg N/ha       |amount of organic N in first soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: iwave ,j
   !YW = WIND EROSION (T/HA)
   !PMRT_44 RATIO OF SOLUBLE C CONCENTRATION IN RUNOFF TO PERCOLATE(0.1_1.)
   real*8, parameter :: PRMT_21 = 1000., PRMT_44 = .5, YW = 0.
   real*8 :: CO, conc, CS, DK, er, latc_clyr, perc_clyr, QBC, sol_latc,&
      &sol_mass, sol_percc, sol_thick, TOT, V, VBC, wt1, YBC, YEW, YOC, X1, X3,&
      &xx, xx1, Y1
   integer :: k
   perc_clyr = 0.

   !!for debug purpose by zhang
   !if (iyr == 1991 .and. i==235) then
   ! write(*,*) 'stop'
   !end if

   if (iwave <= 0) then
      !! HRU calculations
      !xx = sol_n(1,j) + sol_fon(1,j) + sol_mn(1,j)
      xx = sol_LSN(1,j)+sol_LMN(1,j)+sol_HPN(1,j)+sol_HSN(1,j) !+sol_BMN(1,j)
      !wt = sol_bd(1,j) * sol_z(1,j) * 10. (tons/ha)
      !wt1 = wt/1000
      !! conversion factor
      wt1 = sol_bd(1,j) * sol_z(1,j) / 100.

      !! enrichment ratio
      if (erorgn(j) > .001) then
         er = erorgn(j)
      else
         er = enratio
      end if

   else
      !! subbasin calculations
      xx = sub_orgn(iwave)
      wt1 = sub_bd(iwave) * sol_z(1,j) / 100.

      er = enratio
   end if

   conc = xx * er / wt1

   if (iwave <= 0) then
      !! HRU calculations
      sedorgn(j) = .001 * conc * sedyld(j) / hru_ha(j)
   else
      !! subbasin calculations
      sedorgn(j) = .001 * conc * sedyld(j) / (da_ha * sub_fr(iwave))
   end if

   !! update soil nitrogen pools only for HRU calculations
   if (iwave <= 0 .and. xx > 1.e-6) then
      xx1 = (1. - sedorgn(j) / xx)

      !!add by zhang to update soil nitrogen pools

      sol_LSN(1,j) = sol_LSN(1,j) * xx1
      sol_LMN(1,j) = sol_LMN(1,j) * xx1
      sol_HPN(1,j) = sol_HPN(1,j) * xx1
      sol_HSN(1,j) = sol_HSN(1,j) * xx1
   end if

   !Calculate runoff and leached C&N from micro-biomass
   latc_clyr = 0.
   !kg/ha
   sol_mass = (sol_z(1,j) / 1000.) * 10000. * sol_bd(1,j)* 1000. * (1- sol_rock(1,j) / 100.)


   QBC=0.    !c loss with runoff or lateral flow
   VBC=0.    !c los with vertical flow
   YBC=0.    !BMC LOSS WITH SEDIMENT
   TOT=sol_HPC(1,j)+sol_HSC(1,j)+sol_LMC(1,j)+sol_LSC(1,j) !Total organic carbon in layer 1
   ! Not sure whether should consider enrichment ratio or not!
   YEW = Min((sedyld(j)/hru_ha(j)+YW/hru_ha(j))/(sol_mass/1000.),.9) !fraction of soil erosion of total soil mass
   X1=1.-YEW
   !ER enrichment ratio
   !YSD water erosion
   !YW wind erosion
   YOC=YEW*TOT !Organic C loss with sediment
   sol_HSC(1,j)=sol_HSC(1,j)*X1
   sol_HPC(1,j)=sol_HPC(1,j)*X1
   sol_LS(1,j)=sol_LS(1,j)*X1
   sol_LM(1,j)=sol_LM(1,j)*X1
   sol_LSL(1,j)=sol_LSL(1,j)*X1
   sol_LSC(1,j)=sol_LSC(1,j)*X1
   sol_LMC(1,j)=sol_LMC(1,j)*X1
   sol_LSLC(1,j)=sol_LSLC(1,j)*X1
   sol_LSLNC(1,j)=sol_LSC(1,j)-sol_LSLC(1,j)
   if (surfq(j) > 0) then
      !write(*,*) 'stop'
   end if
   IF(sol_BMC(1,j)>.01) THEN
      !KOC FOR CARBON LOSS IN WATER AND SEDIMENT(500._1500.) KD = KOC * C
      sol_WOC(1,j) = sol_LSC(1,j)+sol_LMC(1,j)+sol_HPC(1,j)+sol_HSC(1,j)+sol_BMC(1,j)
      DK=.0001*PRMT_21*sol_WOC(1,j)
      !X1=PO(LD1)-S15(LD1)
      X1 = sol_por(1,j)*sol_z(1,j)-sol_wpmm(1,j) !mm
      IF (X1 <= 0.) THEN
         X1 = 0.01
      END IF
      XX=X1+DK
      V = surfq(j) + sol_prk(1,j) + flat(1,j)
      !QD surface runoff
      IF(V>1.E-10)THEN
         X3=sol_BMC(1,j)*(1.-Exp(-V/XX)) !loss of biomass C
         CO=X3/(sol_prk(1,j) + PRMT_44*(surfq(j)+flat(1,j))) !CS is the horizontal concentration
         CS=PRMT_44*CO                                     !CO is the vertical concentration
         VBC=CO*(sol_prk(1,j)) !!!    sol_prk(:,:) |mm H2O        |percolation from soil layer on current day
         sol_BMC(1,j)=sol_BMC(1,j)-X3
         QBC=CS*(surfq(j)+flat(1,j))
         !     COMPUTE WBMC LOSS WITH SEDIMENT
         IF(YEW>0.)THEN
            CS=DK*sol_BMC(1,j)/XX
            YBC=YEW*CS
         END IF
      END IF
   END IF

   sol_BMC(1,j)=sol_BMC(1,j)-YBC
   surfqc_d(j) = QBC*(surfq(j)/(surfq(j)+flat(1,j)+1.e-6))

   sol_latc = QBC*(flat(1,j)/(surfq(j)+flat(1,j)+1.e-6))
   !sol_percc = VBC ! not used
   sedc_d(j) = YOC + YBC

   latc_clyr = latc_clyr + sol_latc
   DO k=2,sol_nly(j)
      !if (sol_prk(k,j) > 0 .and. k == sol_nly(j)) then
      !write (*,*) 'stop'
      !end if
      sol_thick = sol_z(k,j)-sol_z(k-1,j)
      sol_WOC(k,j) = sol_LSC(k,j)+sol_LMC(k,j)+sol_HPC(k,j)+sol_HSC(k,j)
      Y1=sol_BMC(k,j)+VBC
      VBC=0.
      IF(Y1>=.01)THEN
         V=sol_prk(k,j) + flat(k,j)
         IF(V>0.)VBC=Y1*(1.-Exp(-V/(sol_por(k,j)*sol_thick-sol_wpmm(k,j)+.0001*PRMT_21*sol_WOC(k,j))))
      END IF
      sol_latc = VBC*(flat(k,j)/(sol_prk(k,j) + flat(k,j)+1.e-6))
      sol_percc = VBC-sol_latc
      sol_BMC(k,j)=Y1-VBC

      !! calculate nitrate in percolate
      perc_clyr = perc_clyr + sol_percc

      latc_clyr = latc_clyr + sol_latc
   END DO

   latc_d(j) = latc_clyr
   percc_d(j) = perc_clyr


   return
end
