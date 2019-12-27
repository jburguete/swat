!> @file ndenit.f90
!> file containing the subroutine ndenit
!> @author
!> modified by Javier Burguete

!> this subroutine computes denitrification
!> @param[in] k
!> @param[in] j
!> @param[in] cdg
!> @param[in] wdn
!> @param[out] void
subroutine ndenit(k, j, cdg, wdn, void)

   use parm
   implicit none
   integer, intent(in) :: k, j
   real*8, intent(in) :: cdg, void
   real*8, intent(out) :: wdn
   real*8 vof

   vof = 1. / (1. + (void / 0.04)**5)
   wdn = sol_no3(k,j) * (1. - Exp(-cdn(j) * cdg * vof * sol_cbn(k,j)))
   sol_no3(k,j) = sol_no3(k,j) - wdn

   return
end
