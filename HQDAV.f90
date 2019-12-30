!> @file HQDAV.f90
!> file containing the subroutine HQDAV
!> @author
!> Jaehak Jeong,\n
!> modified by Javier Burguete

!> this subprogram computes flow area and depth given rate in a reach.
!> Adopted from APEX1501 by Jaehak Jeong 2017
SUBROUTINE HQDAV(A,CBW,QQ,SSS,ZCH,ZX,CHW,FPW,jrch)

   use parm
   implicit none

   real*8, intent (inout) :: A, ZX, CHW, FPW
   real*8, intent (in) :: CBW, QQ, SSS, ZCH
   integer, intent (in) :: jrch

   real*8 :: RFPW, RFPX, ACH, AFP, DFDZ, DFQ, FU, FU1, P, Q, QCH, QFP, R, X1, X6, ZX1, ZZ
   integer :: IT, NBCF, NBCX
   ZX=.5*ZCH
   RFPW = ch_w(2,jrch) * 4. !width of floodplain
   RFPX = SQRT(ch_s(2,jrch)) * RFPW / ch_n(1,jrch)

   DO IT=1,10
      IF(QQ>QCAP(jrch))THEN
         ZX=MAX(ZX,ZCH)
         ZZ=ZX-ZCH
         !COMPUTE CH FLOW ABOVE QCAP
         ACH=CHXA(jrch)+ZZ*ch_w(2,jrch)
         R=ACH/CHXP(jrch)
         QCH=ACH*R**.66667*RCHX(jrch)
         CHW=ch_w(2,jrch)
         !COMPUTE FP FLOW
         AFP=ZZ*(RFPW-ch_w(2,jrch))
         QFP=AFP*ZZ**.66667*RFPX/RFPW
         Q=QCH+QFP
         A=ACH+AFP
         FPW=RFPW
         NBCF=1
      ELSE
         X1=ZX*RCSS(jrch)
         A=ZX*(CBW+X1)
         P=CBW+2.*SSS*ZX
         Q=A**1.66667*RCHX(jrch)/P**.66667
         CHW=CBW+2.*X1
         FPW=0.
         NBCX=1
      END IF
      FU=Q-QQ
      X6=MAX(1.,QQ)
      IF(ABS(FU/X6)<.001)EXIT
      IF(IT==1)THEN
         DFQ=-.1*ZX
      ELSE
         DFDZ=(FU-FU1)/(ZX-ZX1)
         DFQ=FU/DFDZ
      END IF
      FU1=FU
      ZX1=ZX
      ZX=ZX-DFQ
   END DO
   RETURN
END
