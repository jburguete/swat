!> @file rootfr.f90
!> file containing the subroutine rootfr
!> @author
!> Armen R. Kemanian,\n
!> modified by Javier Burguete

!> this subroutine distributes dead root mass through the soil profile
!> @param[in] j HRU number
subroutine rootfr(j)
   !! code developed by Armen R. Kemanian in 2008
   !! March, 2009 further adjustments expected

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    a
!!    b
!!    c
!!    d
!!    cum_d
!!    cum_rd
!!    cum_rf
!!    k
!!    l
!!    sol_thick
!!    x1
!!    x2
!!    xx
!!    xx1
!!    xx2
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   ! Normalized Root Density = 1.15*exp[-11.7*NRD] + 0.022, where NRD = normalized rooting depth
   ! Parameters of Normalized Root Density Function from Dwyer et al 19xx
   real*8, parameter :: a = 1.15, b = 11.7, c = 0.022, d = 0.12029 ! Integral of Normalized Root Distribution Function
   ! from 0 to 1 (normalized depth) = 0.12029
   real*8 :: sol_thick(sol_nly(ihru))
   real*8 :: cum_d, cum_rd, cum_rf, x1, x2, xx, xx1, xx2
   integer :: k, l

   if (stsol_rd(j) < 1.e-6) then
      rtfr(1) = 1
      return
   endif


   cum_d = 0.
   cum_rf = 0.
   sol_thick(:) = 0.
   rtfr = 0.

   k = 0
   do l=1, sol_nly(j)
      if (l == 1) then
         sol_thick(l) = sol_z(l,j)
      else
         sol_thick(l) = sol_z(l,j) - sol_z(l-1,j)
      end if

      cum_d = cum_d + sol_thick(l)
      if (cum_d >= stsol_rd(j)) cum_rd = stsol_rd(j)
      if (cum_d < stsol_rd(j)) cum_rd = cum_d
      x1 = (cum_rd - sol_thick(l)) / stsol_rd(j)
      x2 = cum_rd / stsol_rd(j)
      xx1 = -b * x1
      if (xx1 > 20.) xx1 = 20.
      xx2 = -b * x2
      if (xx2 > 20.) xx2 = 20.
      rtfr(l) = (a / b * (Exp(xx1) - Exp(xx2)) + c *(x2 - x1)) / d
      xx = cum_rf
      cum_rf = cum_rf + rtfr(l)
      if (cum_rf > 1.) then
         rtfr(l) = 1. - xx
         cum_rf = 1.0
      end if
      k = l
      if (cum_rd >= stsol_rd(j)) Exit

   end do

   !!  ensures that cumulative fractional root distribution = 1
   if (k > 0) then
      do l=1,k ! exits loop on the same layer as the previous loop
         rtfr(l) = rtfr(l) / cum_rf
      end do
   end if

end subroutine
