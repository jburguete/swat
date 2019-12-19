!> @file estimate_ksat.f90
!> file containing the subroutine estimate_ksat
!> @author
!> modified by Javier Burguete

!> This subroutine calculates ksat value for a soil layer
!> given the % of clay in the soil layer

!> Background: published work of Walter Rawls. Calculated
!> ksat values based on soil texture (sand, silt and clay).
!> Idea: there exists a relationship between % clay and Ksat.
!> Equations used in this subroutine are based on the above
!> idea (Jimmy Willimas)
!> @param[in] perc_clay clay percentage (%)
!> @param[out] esti_ksat estimated ksat
subroutine estimate_ksat(perc_clay,esti_ksat)

   implicit none

!  real*8, parameter :: exksat = 5.0
   real*8, intent(in) :: perc_clay
   real*8, intent(out) :: esti_ksat
   real*8 :: xc

! print *,"Enter the % clay in the soil layer"
! read *,perc_clay

   xc = 100.0 - perc_clay
   esti_ksat = 12.7 * xc / (xc + exp(11.45 - 0.097 * xc)) + 1.0

! print *,"The estimated ksat value is ",min(esti_ksat,exksat)

   return

end subroutine estimate_ksat
