!> @file carbon_zhang2.f90
!> file containing the subroutine carbon_zhang2
!> @author
!> modified by Javier Burguete

!> @param j HRU number
subroutine carbon_zhang2(j)
   use parm
   implicit none
   !!============================================
   !!Input variables
   !!    sol_bd(:,:)   |Mg/m**3       |bulk density of the soil
   !!    sol_st(:,:)   |mm H2O        |amount of water stored in the soil layer on
   !!                                 |current day
   !!    sol_fc(:,:)   |mm H2O        |amount of water available to plants in soil
   !!                                 |layer at field capacity (fc - wp),Index:(layer,HRU)
   !!    sol_wp(:,:)   |mm H20/mm soil|water content of soil at -1.5 MPa (wilting
   !!                                 |point)
   !!    sol_wpmm(:,:) |mm H20        |water content of soil at -1.5 MPa (wilting
   !!                                 |point)
   !!  sol_percc(k,j)
   !!  sol_latc(k,j)
   !!==============================================
   !!Transput variables;
   !!  sol_HSC(:,:)    : mass of C present in slow humus (kg ha-1)
   !!  sol_HSN(:,:)    : mass of N present in slow humus (kg ha-1)
   !!  sol_HPC(:,:)    : mass of C present in passive humus (kg ha-1)
   !!  sol_HPN(:,:)    : mass of N present in passive humus (kg ha-1)
   !!  sol_LM(:,:)     : mass of metabolic litter (kg ha-1)
   !!  sol_LMC(:,:     : mass of C in metabolic litter (kg ha-1)
   !!  sol_LMN(:,:)    : mass of N in metabolic litter (kg ha-1)
   !!  sol_LS(:,:)     : mass of structural litter (kg ha-1)
   !!  sol_LSC(:,:)    : mass of C in structural litter (kg ha-1)
   !!  sol_LSL(:,:)    : mass of lignin in structural litter (kg ha-1)
   !!  sol_LSN(:,:)    : mass of N in structural litter (kg ha-1)
   !!  STD(:)          : standing dead (kg ha-1)                                               (Not used)
   !!  STDL(:)         : mass of lignin in standing dead (kg ha-1)                             (Not used)
   !!  STDN(:)         : mass of N in standing dead (dead plants + sorbed from soil; kg ha-1)  (Not used)
   !!  STDNEl(:)       : standing dead N after enrichment with sorbed N in a soil layer (kg ha-1)
   !!  sol_NO3(:,:)    : weight of NO3-N in soil layer (kg ha-1)
   !!  sol_NH3(:,:)    : weight of NH3-N in soil layer (kg ha-1)

   !!==============================================
   !!read in parameters
   !HPR     : rate of transformation of passive humus under optimal conditions (subsurface
   !layers = 0.000012 day-1) (Parton et al.,1993, 1994)
   !HSR     : rate of transformation of slow humus under optimal conditions (all layers
   != 0.0005 day-1) (Parton et al., 1993, 1994; Vitousek et al., 1993)
   !KOC     : liquid–solid partition coefficient for microbial biomass (10^3 m^3 Mg-1)
   !BMR     : rate of transformation of microbial biomass and associated products under optimal
   !            conditions (surface = 0.0164 day-1; all other layers = 0.02 day-1) (Parton et al., 1993, 1994)
   !LMR     : rate of transformation of metabolic litter under optimal conditions (surface =
   !0.0405 day-1; all other layers = 0.0507 day-1) (Parton et al., 1994)
   !Sf      : fraction of mineral N sorbed to litter: 0.05 for surface litter, 0.1 for belowground litter

   !Cf      : carbon fraction of organic materials 0.42; from data of Pinck et al., 1950)
   !LSR     : rate of potential transformation of structural litter under optimal conditions
   !(surface = 0.0107 day-1; all other layers= 0.0132 day-1) (Parton et al., 1994)
   !PRMT_51 !COEF ADJUSTS MICROBIAL ACTIVITY FUNCTION IN TOP SOIL LAYER (0.1_1.),
   !PRMT_45 !COEF IN CENTURY EQ ALLOCATING SLOW TO PASSIVE HUMUS(0.001_0.05) ORIGINAL VALUE = 0.003,


   !!==============================================
   !! local variables
   !ABCO2   : allocation from biomass to CO2; 0.6 (surface litter), 0.85–0.68*(CLAF + SILF) (all other layers) (Parton et al., 1993, 1994)
   !ABL     : carbon allocation from biomass to leaching; ABL = (1-exp(-f/(0.01* SW+ 0.1*(KdBM)*DB)) (Williams, 1995)
   !ABP     : allocation from biomass to passive humus; 0 (surface litter), 0.003 + 0.032*CLAF (all other layers) (Parton et al., 1993, 1994)
   !ALMCO2  : allocation from metabolic litter to CO2; 0.6 (surface litter), 0.55 (all other layers) (Parton et al., 1993, 1994)
   !ALSLCO2 : allocation from lignin of structural litter to CO2; 0.3 (Parton et al., 1993, 1994)
   !ALSLNCO2: allocation from non-lignin of structural litter to CO2; 0.6 (surface litter), 0.55 (all other layers) (Parton et al., 1993, 1994)
   !APCO2   : allocation from passive humus to CO2; 0.55 (Parton et al., 1993, 1994)
   !ASCO2   : allocation from slow humus to CO2; 0.55 (Parton et al., 1993, 1994)
   !ASP     : allocation from slow humus to passive; 0 (surface litter), 0.003-0.009*CLAF (all other layers) (Parton et al., 1993, 1994)
   !BMC     : mass of C in soil microbial biomass and associated products (kg ha-1)
   !BMCTP   : potential transformation of C in microbial biomass (kg ha-1 day-1)
   !BMN     : mass of N in soil microbial biomass and associated products (kg ha-1)
   !BMNTP   : potential transformation of N in microbial biomass (kg ha-1 day-1)

   !CDG     : soil temperature control on biological processes
   !CNR     : C/N ratio of standing dead
   !CS      : combined factor controlling biological processes [CS = sqrt(CDG×sut)* 0.8*OX*X1), CS < 10; CS = 10, CS>=10 (Williams, 1995)]
   !DBp     : soil bulk density of plow layer (Mg m-3) (Not used)
   !HSCTP   : potential transformation of C in slow humus (kg ha-1 day-1)
   !HSNTP   : potential transformation of N in slow humus (kg ha-1 day-1)
   !HPCTP   : potential transformation of C in passive humus (kg ha-1 day-1)
   !HPNTP   : potential transformation of N in passive humus (kg ha-1 day-1)
   !LMF     : fraction of the litter that is metabolic
   !LMNF    : fraction of metabolic litter that is N (kg kg-1)
   !LMCTP   : potential transformation of C in metabolic litter (kg ha-1 day-1)
   !LMNTP   : potential transformation of N in metabolic litter (kg ha-1 day-1)
   !LSCTP   : potential transformation of C in structural litter (kg ha-1 day-1)
   !LSF     : fraction of the litter that is structural
   !LSLF    : fraction of structural litter that is lignin (kg kg-1)
   !LSNF    : fraction of structural litter that is N (kg kg-1)
   !LSLCTP  : potential transformation of C in lignin of structural litter (kg ha-1 day-1)
   !LSLNCTP : potential transformation of C in nonlignin structural litter (kg ha-1 day-1)
   !LSNTP   : potential transformation of N in structural litter (kg ha-1 day-1)
   !NCBM    : N/C ratio of biomass
   !NCHP    : N/C ratio passive humus
   !NCHS    : N/C ratio of the slow humus
   !OX      : oxygen control on biological processes with soil depth
   !sut     : soil water control on biological processes
   !X1      : tillage control on residue decomposition (Not used)
   !XBMT    : control on transformation of microbial biomass by soil texture and structure.
   !Its values: surface litter layer = 1; all other layers = 1-0.75*(SILF + CLAF) (Parton et al., 1993, 1994)
   !XLSLF   : control on potential transformation of structural litter by lignin fraction
   !of structural litter [XLSLF = exp(-3* LSLF) (Parton et al., 1993, 1994)]
   integer, intent(in) :: j
   real*8, parameter :: wdn = 0., HSR = 5.4799998E-04, HPR = 1.2000000E-05,&
      &A1CO2=.55, APCO2=.55, ASCO2=.60
   real*8 :: sol_mass, sol_min_n
   real*8 :: fc, wc, sat, void, sut, cdg, OX, CS
   real*8 :: X1,X3, XX
   real*8 :: BMR
   real*8 :: LSR, LSCTA, LSLCTA, LSLNCTA, LSNTA, LMCTA, LMNTA, BMCTA, BMNTA, HSCTA, HSNTA, HPCTA, HPNTA
   real*8 :: LSCTP, LSLCTP, LSLNCTP, LSNTP, LMR, LMCTP, LMNTP, BMCTP,HSCTP, HSNTP, HPCTP, HPNTP
   real*8 :: NCHP, NCBM, NCHS, ABCO2, ABP, ASP, A1, ASX, APX
   real*8 :: PRMT_51, PRMT_45
   real*8 :: DF1, DF2, DF3, DF4, DF5, DF6, ADD, ADF1, ADF2, ADF3, ADF4, ADF5
   real*8 :: TOT
   real*8 :: PN1, PN2, PN3, PN5, PN6, PN7, PN8, PN9
   real*8 :: SUM1, CPN1, CPN2, CPN3, CPN4, CPN5
   real*8 :: WMIN,DMDN
   real*8 :: BMNTP, decr, hmn, hmp, hmp_rate, RLR, rmn1, rmp, RTO, rwn, XBM
   integer :: k, kk


   !!    zero new carbon variables for output.hru
   cmup_kgh(j) = 0.
   cmtot_kgh(j) = 0.
   !! initilize local variables

   !!for debug purpose by zhang
   !if (iyr == 1941 .and. i==134) then
   !write(*,*) 'stop'
   !end if

   !calculate tillage factor using DSSAT
   if (tillage_switch(j) .eq. 1 .and. tillage_days(j) .le. 30) then
      tillage_factor(j) = 1.6
   else
      tillage_factor(j) = 1.0
   end if
   !calculate tillage factor using DSSAT


   !!calculate C/N dynamics for each soil layer
   !!===========================================
   do k = 1, sol_nly(j)
      if (k == 1) then
         !10 cm / 1000 = 0.01m; 1 ha = 10000 m2; ton/m3; * 1000 --> final unit is kg/ha; rock fraction is considered
         sol_mass = (10) / 1000.* 10000. * sol_bd(k,j)* 1000. *(1- sol_rock(k,j) / 100.)
      else
         sol_mass = (sol_z(k,j) - sol_z(k-1,j)) / 1000.* 10000. * sol_bd(k,j)* 1000. *(1- sol_rock(k,j) / 100.)
      end if

      !!If k = 1, then using temperature, soil moisture in layer 2 to calculate decomposition factor
      !!Not
      if (k == 1) then
         kk = 2
      else
         kk = k
      end if
      !! mineralization can occur only if temp above 0 deg
      !check sol_st soil water content in each soil ayer mm H2O
      if (sol_tmp(k,j) > 0. .AND. sol_st(k,j) > 0.) then
         !!from Armen
         !!compute soil water factor - sut
         fc = sol_fc(k,j) + sol_wpmm(k,j)        ! units mm
         wc = sol_st(k,j) + sol_wpmm(k,j)        ! units mm
         sat = sol_ul(k,j) + sol_wpmm(k,j)       ! units mm
         void = sol_por(k,j) * (1. - wc / sat)   ! fraction

         X1 = wc - sol_wpmm(k,j)
         IF(X1<0.)THEN
            sut = .1 * (sol_st(kk,j) / sol_wpmm(k,j))**2
         ELSE
            sut = .1 + .9 * Sqrt(sol_st(k,j) / sol_fc(k,j))
         END IF
         sut = Min(1., sut)
         sut = Max(.05, sut)
         !check X1, FC, S15

         !!compute tillage factor (X1)
         !use the tillfactor module from Armen
         !X1 = ftilf(tillagef(kk,j), wc, sat)
         X1 = 1.0

         !calculate tillage factor using DSSAT
         if (tillage_switch(j) .eq. 1 .and. tillage_days(j) .le. 30) then
            if (k == 1) then
               X1 = 1.6
            else
               if (sol_z(k,j) .le. tillage_depth(j)) then
                  X1 = 1.6
               elseif (sol_z(k-1,j) .lt. tillage_depth(j)) then
                  X1 = 1.0 + 0.6*(tillage_depth(j) - sol_z(k-1,j))/(sol_z(k,j) - sol_z(k-1,j))
               end if
            end if
         else
            X1 = 1.0
         end if
         !calculate tillage factor using DSSAT

         !!compute soil temperature factor
         cdg = fcgd(sol_tmp(k,j))

         !!compute oxygen (OX)
         OX=1.-0.8*((sol_z(kk,j)+sol_z(kk-1,j))/2)/&
            &(((sol_z(kk,j)+sol_z(kk-1,j))/2)&
            &+Exp(18.40961-0.023683632*((sol_z(kk,j)+sol_z(kk-1,j))/2)))
         !! compute combined factor
         CS=Min(10.,Sqrt(cdg*sut)*0.9*OX*X1)

         !! call denitrification (to use void and cdg factor)
         if (cdg > 0. .and. void <= 0.1) then
            call ndenit(k,j,cdg,wdn,void)
         end if
         wshd_dnit = wshd_dnit + wdn * hru_dafr(j)
         wdntl = wdntl + wdn



         sol_min_n = sol_no3(k,j) + sol_nh3(k,j)


         !lignin content in structural litter (fraction)
         RLR = Min(0.8,sol_LSL(k,j)/(sol_LS(k,j) + 1.E-5))


         !PRMT_51 =0.   !COEF ADJUSTS MICROBIAL ACTIVITY FUNCTION IN TOP SOIL LAYER (0.1_1.),
         PRMT_51 = 1.
         !!The following codes are clculating of the N:C ration in the newly formed SOM for each pool
         !!please note that in the surface layer, no new materials enter Passive pool, therefore, no NCHP is
         !!calculated for the first layer.
         IF(k==1)THEN
            CS=CS*PRMT_51
            ABCO2=.55
            BMR=.0164
            LMR=.0405
            LSR=.0107
            NCHP=.1
            XBM=1.
            !     COMPUTE N/C RATIOS
            X1 = 0.1*(sol_LSN(k,j)+sol_LMN(k,j))/(sol_rsd(k,j)/1000+1.E-5) !relative notrogen content in residue (%)
            IF(X1>2.)THEN
               NCBM=.1
            ELSE
               IF(X1>.01)THEN
                  NCBM=1./(20.05-5.0251*X1)
               ELSE
                  NCBM=.05
               END IF
            END IF
            NCHS=NCBM/(5.*NCBM+1.)
         ELSE
            ABCO2 = 0.17 + 0.0068 * sol_sand(k,j)
            BMR=.02
            LMR=.0507
            LSR=.0132
            XBM=.25+.0075*sol_sand(k,j)
            X1 = 1000. * sol_min_n/(sol_mass/1000)
            IF(X1>7.15)THEN
               NCBM=.33
               NCHS=.083
               NCHP=.143
            ELSE
               NCBM=1./(15.-1.678*X1)
               NCHS=1./(20.-1.119*X1)
               NCHP=1./(10.-.42*X1)
            END IF
         END IF
         ABP=.003+.00032*sol_clay(k,j)
         PRMT_45 = 0.  !COEF IN CENTURY EQ ALLOCATING SLOW TO PASSIVE HUMUS(0.001_0.05) ORIGINAL VALUE = 0.003,
         PRMT_45 = 5.0000001E-02
         ASP=Max(.001,PRMT_45-.00009*sol_clay(k,j))
         !     POTENTIAL TRANSFORMATIONS STRUCTURAL LITTER
         X1=LSR*CS*Exp(-3.*RLR)
         LSCTP=X1*sol_LSC(k,j)
         LSLCTP=LSCTP*RLR
         LSLNCTP=LSCTP*(1.-RLR)
         LSNTP=X1*sol_LSN(k,j)
         !     POTENTIAL TRANSFORMATIONS METABOLIC LITTER
         X1=LMR*CS
         LMCTP=sol_LMC(k,j)*X1
         LMNTP=sol_LMN(k,j)*X1
         !     POTENTIAL TRANSFORMATIONS MICROBIAL BIOMASS
         X1=BMR*CS*XBM
         BMCTP=sol_BMC(k,j)*X1
         BMNTP=sol_BMN(k,j)*X1
         !     POTENTIAL TRANSFORMATIONS SLOW HUMUS
         X1=HSR*CS
         HSCTP=sol_HSC(k,j)*X1
         HSNTP=sol_HSN(k,j)*X1
         !     POTENTIAL TRANSFORMATIONS PASSIVE HUMUS
         X1=CS*HPR
         HPCTP=sol_HPC(k,j)*X1
         HPNTP=sol_HPN(k,j)*X1
         !     ESTIMATE N DEMAND
         A1=1.-A1CO2
         ASX=1.-ASCO2-ASP
         APX=1.-APCO2

         PN1=LSLNCTP*A1*NCBM               !Structural Litter to Biomass
         PN2=.7*LSLCTP*NCHS               !Structural Litter to Slow
         PN3=LMCTP*A1*NCBM                 !Metabolic Litter to Biomass
         !PN4=BMCTP*ABL*NCBM                !Biomass to Leaching (calculated in NCsed_leach)
         PN5=BMCTP*ABP*NCHP                !Biomass to Passive
         PN6=BMCTP*(1.-ABP-ABCO2)*NCHS     !Biomass to Slow
         PN7=HSCTP*ASX*NCBM                !Slow to Biomass
         PN8=HSCTP*ASP*NCHP                !Slow to Passive
         PN9=HPCTP*APX*NCBM                !Passive to Biomass

         !     COMPARE SUPPLY AND DEMAND FOR N
         SUM1=0.
         CPN1=0.
         CPN2=0.
         CPN3=0.
         CPN4=0.
         CPN5=0.
         X1=PN1+PN2
         IF(LSNTP<X1)THEN
            CPN1=X1-LSNTP
         ELSE
            SUM1=SUM1+LSNTP-X1
         END IF
         IF(LMNTP<PN3)THEN
            CPN2=PN3-LMNTP
         ELSE
            SUM1=SUM1+LMNTP-PN3
         END IF
         X1=PN5+PN6
         IF(BMNTP<X1)THEN
            CPN3=X1-BMNTP
         ELSE
            SUM1=SUM1+BMNTP-X1
         END IF
         X1=PN7+PN8
         IF(HSNTP<X1)THEN
            CPN4=X1-HSNTP
         ELSE
            SUM1=SUM1+HSNTP-X1
         END IF
         IF(HPNTP<PN9)THEN
            CPN5=PN9-HPNTP
         ELSE
            SUM1=SUM1+HPNTP-PN9
         END IF
         !     WNH3(ISL)=WNH3(ISL)+SUM1
         !total available N
         WMIN=Max(1.E-5,sol_NO3(k,j) + sol_NH3(k,j)+SUM1)
         !total demand for potential tranformaiton of SOM
         DMDN=CPN1+CPN2+CPN3+CPN4+CPN5

         X3=1.
         !     REDUCE DEMAND IF SUPPLY LIMITS
         IF(WMIN<DMDN) then
            X3=WMIN/DMDN
         end if
         !     ACTUAL TRANSFORMATIONS
         IF(CPN1>0.)THEN
            LSCTA=LSCTP*X3
            LSNTA=LSNTP*X3
            LSLCTA=LSLCTP*X3
            LSLNCTA=LSLNCTP*X3
         ELSE
            LSCTA=LSCTP
            LSNTA=LSNTP
            LSLCTA=LSLCTP
            LSLNCTA=LSLNCTP
         END IF
         IF(CPN2>0.)THEN
            LMCTA=LMCTP*X3
            LMNTA=LMNTP*X3
         ELSE
            LMCTA=LMCTP
            LMNTA=LMNTP
         END IF
         IF(CPN3>0.)THEN
            BMCTA=BMCTP*X3
            BMNTA=BMNTP*X3
         ELSE
            BMCTA=BMCTP
            BMNTA=BMNTP
         END IF
         IF(CPN4>0.)THEN
            HSCTA=HSCTP*X3
            HSNTA=HSNTP*X3
         ELSE
            HSCTA=HSCTP
            HSNTA=HSNTP
         END IF
         IF(CPN5>0.)THEN
            HPCTA=HPCTP*X3
            HPNTA=HPNTP*X3
         ELSE
            HPCTA=HPCTP
            HPNTA=HPNTP
         END IF

         !Recalculate demand using actural transformations
         !revised from EPIC code by Zhang
         PN1=LSLNCTA*A1*NCBM               !Structural Litter to Biomass
         PN2=.7*LSLCTA*NCHS               !Structural Litter to Slow
         PN3=LMCTA*A1*NCBM                 !Metabolic Litter to Biomass
         !PN4=BMCTP*ABL*NCBM                !Biomass to Leaching (calculated in NCsed_leach)
         PN5=BMCTA*ABP*NCHP                !Biomass to Passive
         PN6=BMCTA*(1.-ABP-ABCO2)*NCHS     !Biomass to Slow
         PN7=HSCTA*ASX*NCBM                !Slow to Biomass
         PN8=HSCTA*ASP*NCHP                !Slow to Passive
         PN9=HPCTA*APX*NCBM                !Passive to Biomass
         !     COMPARE SUPPLY AND DEMAND FOR N
         SUM1=0.
         CPN1=0.
         CPN2=0.
         CPN3=0.
         CPN4=0.
         CPN5=0.
         X1=PN1+PN2
         IF(LSNTA<X1)THEN
            CPN1=X1-LSNTA
         ELSE
            SUM1=SUM1+LSNTA-X1
         END IF
         IF(LMNTA<PN3)THEN
            CPN2=PN3-LMNTA
         ELSE
            SUM1=SUM1+LMNTA-PN3
         END IF
         X1=PN5+PN6
         IF(BMNTA<X1)THEN
            CPN3=X1-BMNTA
         ELSE
            SUM1=SUM1+BMNTA-X1
         END IF
         X1=PN7+PN8
         IF(HSNTA<X1)THEN
            CPN4=X1-HSNTA
         ELSE
            SUM1=SUM1+HSNTA-X1
         END IF
         IF(HPNTA<PN9)THEN
            CPN5=PN9-HPNTA
         ELSE
            SUM1=SUM1+HPNTA-PN9
         END IF
         !total available N
         WMIN=Max(1.E-5,sol_NO3(k,j) + sol_NH3(k,j)+SUM1)
         !total demand for potential tranformaiton of SOM
         DMDN=CPN1+CPN2+CPN3+CPN4+CPN5


         !supply - demand
         sol_RNMN(k,j)=SUM1-DMDN
         !     UPDATE
         IF(sol_RNMN(k,j)>0.)THEN
            sol_NH3(k,j)=sol_NH3(k,j)+sol_RNMN(k,j)
         ELSE
            X1=sol_NO3(k,j)+sol_RNMN(k,j)
            IF(X1<0.)THEN
               sol_RNMN(k,j)=-sol_NO3(k,j)
               sol_NO3(k,j)=1.E-10
            ELSE
               sol_NO3(k,j)=X1
            END IF
         END IF
         DF1=LSNTA

         DF2=LMNTA
         !!DF represents Demand from

         !calculate P flows
         !! compute humus mineralization on active organic p
         hmp_rate = 1.4* (HSNTA + HPNTA)/(sol_HSN(k,j) + sol_HPN(k,j) + 1.e-6)
         hmp = hmp_rate*sol_orgp(k,j)
         hmp = Min(hmp, sol_orgp(k,j))
         sol_orgp(k,j) = sol_orgp(k,j) - hmp
         sol_solp(k,j) = sol_solp(k,j) + hmp

         !! compute residue decomp and mineralization of
         !! fresh organic n and p (upper two layers only)
         decr = (LSCTA + LMCTA)/(sol_LSC(k,j) + sol_LMC(k,j) + 1.e-6)
         decr = Min(1., decr)
         rmp = decr * sol_fop(k,j)

         sol_fop(k,j) = sol_fop(k,j) - rmp
         sol_solp(k,j) = sol_solp(k,j) + .8 * rmp
         sol_orgp(k,j) = sol_orgp(k,j) + .2 * rmp
         !calculate P flows


         LSCTA = Min(sol_LSC(k,j),LSCTA)
         sol_LSC(k,j)=Max(1.E-10,sol_LSC(k,j)-LSCTA)
         LSLCTA = Min(sol_LSLC(k,j),LSLCTA)
         sol_LSLC(k,j)=Max(1.E-10,sol_LSLC(k,j)-LSLCTA)
         sol_LSLNC(k,j)=Max(1.E-10,sol_LSLNC(k,j)-LSLNCTA)
         LMCTA=Min(sol_LMC(k,j),LMCTA)
         IF (sol_LM(k,j) > 0.) THEN
            RTO = Max(0.42,sol_LMC(k,j)/sol_LM(k,j))
            sol_LM(k,j) = sol_LM(k,j) - LMCTA/RTO
            sol_LMC(k,j) = sol_LMC(k,j) - LMCTA
         END IF
         sol_LSL(k,j)=Max(1.E-10,sol_LSL(k,j)-LSLCTA/.42)
         sol_LS(k,j)=Max(1.E-10,sol_LS(k,j)-LSCTA/.42)

         X3=APX*HPCTA+ASX*HSCTA+A1*(LMCTA+LSLNCTA)
         sol_BMC(k,j)=sol_BMC(k,j)-BMCTA+X3
         DF3=BMNTA-NCBM*X3
         !!DF3 is the supply of BMNTA - demand of N to meet the Passive, Slow, Metabolic, and Non-lignin Structural
         !! C pools transformaitons into microbiomass pool
         X1=.7*LSLCTA+BMCTA*(1.-ABP-ABCO2)
         sol_HSC(k,j)=sol_HSC(k,j)-HSCTA+X1
         DF4=HSNTA-NCHS*X1
         !!DF4 Slow pool supply of N - N demand for microbiomass C transformed into slow pool
         X1=HSCTA*ASP+BMCTA*ABP
         sol_HPC(k,j)=sol_HPC(k,j)-HPCTA+X1
         DF5=HPNTA-NCHP*X1
         !!DF5 Passive pool demand of N - N demand for microbiomass C transformed into passive pool
         DF6=sol_min_n-sol_NO3(k,j)-sol_NH3(k,j)
         !!DF6 Supply of mineral N - available mineral N = N demanded from mineral pool
         !SMS(10,ISL)=SMS(10,ISL)-DF6
         ADD=DF1+DF2+DF3+DF4+DF5+DF6
         ADF1=Abs(DF1)
         ADF2=Abs(DF2)
         ADF3=Abs(DF3)
         ADF4=Abs(DF4)
         ADF5=Abs(DF5)
         TOT=ADF1+ADF2+ADF3+ADF4+ADF5
         XX=ADD/(TOT+1.E-10)
         sol_LSN(k,j)=Max(.001,sol_LSN(k,j)-DF1+XX*ADF1)
         sol_LMN(k,j)=Max(.001,sol_LMN(k,j)-DF2+XX*ADF2)
         sol_BMN(k,j)=sol_BMN(k,j)-DF3+XX*ADF3
         sol_HSN(k,j)=sol_HSN(k,j)-DF4+XX*ADF4
         sol_HPN(k,j)=sol_HPN(k,j)-DF5+XX*ADF5
         sol_RSPC(k,j)=.3*LSLCTA+A1CO2*(LSLNCTA+LMCTA)+ABCO2*BMCTA+ASCO2*HSCTA+APCO2*HPCTA
         rspc_d(j) = rspc_d(j) +  sol_RSPC(k,j)
         sol_rsd(k,j)= sol_LS(k,j)+sol_LM(k,j)
         sol_orgn(k,j) = sol_HPN(k,j)
         sol_aorgn(k,j) = sol_HSN(k,j)
         sol_fon(k,j) = sol_LMN(k,j) + sol_LSN(k,j)
         sol_cbn(k,j) = 100*(sol_LSC(k,j)+sol_LMC(k,j) +sol_HSC(k,j) + sol_HPC(k,j) + sol_BMC(k,j))/sol_mass

         !!    carbon outputs for .hru file
         if (k == 1) cmup_kgh(j) = sol_cbn(k,j) * sol_mass / 100.
         cmtot_kgh(j) = cmtot_kgh(j) + sol_cbn(k,j) * sol_mass / 100.


         !! summary calculations
         !! calculations are based on century model, and not alighned with SWAT old algorithm yet.
         if (curyr > nyskip) then
            hmn = sol_RNMN(k,j)
            wshd_hmn = wshd_hmn + hmn * hru_dafr(j)
            rwn = HSNTA
            wshd_rwn = wshd_rwn + rwn * hru_dafr(j)

            wshd_hmp = wshd_hmp + hmp * hru_dafr(j)
            rmn1 = (LSNTA+LMNTA)
            wshd_rmn = wshd_rmn + rmn1 * hru_dafr(j)
            wshd_rmp = wshd_rmp + rmp * hru_dafr(j)
            wshd_dnit = wshd_dnit + wdn * hru_dafr(j)
            hmntl = hmntl + hmn
            rwntl = rwntl + rwn
            hmptl = hmptl + hmp
            rmn2tl = rmn2tl + rmn1
            rmptl = rmptl + rmp
            wdntl = wdntl + wdn
         end if




      end if

   end do

   !write (*,*) iyr,i, DeltaBMC
   !write (*,*) iyr,i, wdntl
   return
end
