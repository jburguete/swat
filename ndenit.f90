subroutine ndenit(k,j,cdg,wdn,void)
!!    this subroutine computes denitrification

   use parm, only: cdn, sol_cbn, sol_no3
   implicit none
   integer, intent(in) :: k,j
   real*8, intent(in) :: cdg, void
   real*8, intent(out) :: wdn
   real*8 vof

   wdn = 0.
   vof = 1. / (1. + (void/0.04)**5)
   wdn = sol_no3(k,j) * (1. - Exp(-cdn(j) * cdg * vof * sol_cbn(k,j)))
   sol_no3(k,j) = sol_no3(k,j) - wdn

   return
end
