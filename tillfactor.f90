!> @file tillfactor.f90
!> file containing the subroutine tillfactor
!> @author
!> modified by Javier Burguete

!> this procedure increases tillage factor (tillagef(l,j) per layer for each
!> operation. The tillage factor settling will depend of soil moisture
!> (tentatively) and must be called every day. For simplicity the settling is
!> calculated now at the soil carbon sub because soil water content is
!> available.
!> @param[in] j HRU number (none)
!> @param[in] bmix biological mixing efficiency: this number is zero for tillage
!> operations (none)
!> @param[inout] emix mixing efficiency (none)
!> @param[in] dtil depth of mixing (mm)
!> @param[in] sol_thick

!> The tillage factor depends on the cumulative soil disturbance rating = csdr
!> For simplicity, csdr is a function of emix.
!> First step is to calculate "current" csdr by inverting tillage factor
!> function.
!> The effect of texture on tillage factor (ZZ) is removed first (and recovered
!> at the end of the procedure).
!> \f\[YY = tillagef(l,j) / ZZ\f\]
!> Since the tillage factor function is non linear, iterations are needed.
!> \f$XX=0.5\f$ is the initial value that works OK for the range of values
!> observed.
!> If a layer is only partially tilled then emix is corrected accordingly
subroutine tillfactor(j,bmix,emix,dtil,sol_thick)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j             |none          |HRU number
!!    bmix          |none          |biological mixing efficiency: this
!!                                 |number is zero for tillage operations
!!    emix          |none          |mixing efficiency
!!    dtil          |mm            |depth of mixing
!!    sol_thick
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    emix          |none          |mixing efficiency
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    csdr
!!    l           |none          |counter
!!    m1
!!    m2
!!    xx
!!    xx1
!!    xx2
!!    yy
!!    zz
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: j
   real*8, intent (in) :: bmix
   real*8, intent (inout) :: emix
   real*8, intent (in) :: dtil
   real*8, intent (in) :: sol_thick(sol_nly(j))
   real*8, parameter :: m1 = 1, m2 = 2
   real*8 :: csdr, xx, xx1, xx2, yy, zz
   integer :: l

   emix = emix - bmix ! this is to avoid affecting tillage factor with biological mixing

   if (emix > 0.) then

      do l=1, sol_nly(j)

         if (sol_z(l,j) <= dtil) then
            !emix = emix !redundant
         else if (sol_z(l,j) > dtil .AND. sol_z(l-1,j) < dtil) then
            emix = emix * (dtil - sol_z(l-1,j)) / sol_thick(l)
         else
            emix = 0.
         end if

         ! to save computation time if emix = 0 here then the other layers can be avoided
         ! tillage always proceeds from top to bottom
         if (emix <= 0.) exit

         xx = 0.
         zz = 3. + (8. - 3.) * Exp(-5.5 * sol_clay(l,j)/100.)
         yy = tillagef(l,j) / zz

         ! empirical solution for x when y is known and y=x/(x+Exp(m1-m2*x))
         if (yy > 0.01) then
            xx1 = yy ** Exp(-0.13 + 1.06 * yy)
            xx2 = Exp(0.64 + 0.64 * yy ** 100.)
            xx = xx1 * xx2
         end if

         csdr = xx + emix
         tillagef(l,j) = zz * (csdr / (csdr + Exp(m1 - m2 * csdr)))

      end do

   end if

   return
end subroutine
