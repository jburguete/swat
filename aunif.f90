!> @file aunif.f90
!> file containing the function aunif
!> @author
!> modified by Javier Burguete

!> This function generates random numbers ranging from 0.0 to 1.0.
!> In the process of calculating the random number, the seed (x1) is
!> set to a new value.
!> This function implements the prime-modulus generator
!> \f\[xi=16807\,xi\,\textrm{mod}\left(2^{31}-1\right)\f\]
!> using code which ensures that no intermediate result uses more than
!> 31 bits.
!> The theory behind the code is summarized in \cite Bratley83
!> @param[inout] x1
!> random number generator seed (integer) where \f$0 < x1 < 2147483647\f$
!> @return random number ranging from 0.0 to 1.0
real*8 function aunif (x1) result (unif)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    x1         |none           |random number generator seed (integer)
!!                               |where  0 < x1 < 2147483647
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name       |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    unif       |none           |random number between 0.0 and 1.0
!!    x1         |none           |random number generator seed (integer)
!!                               |set to new value where 0 < x1 < 2147483647
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name       |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    x2         |none           |variable to hold calculation results
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!     ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
   implicit none

   integer, intent (inout) :: x1
   integer :: x2

   x2 = x1 / 127773
   x1 = 16807 * (x1 - x2 * 127773) - x2 * 2836
   if (x1 < 0) x1 = x1 + 2147483647
   unif = x1 * 4.656612875d-10

   return
end
