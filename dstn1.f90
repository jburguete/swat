!> @file dstn1.f90
!> file containing the function dstn1
!> @author
!> modified by Javier Burguete

!> this function computes the distance from the mean of a normal
!> distribution with mean = 0 and standard deviation = 1, given two
!> random numbers
!> @param[in] rn1 first random number
!> @param[in] rn2 second random number
!> @return distance from the mean
real*8 function dstn1(rn1,rn2) result (r_dstn1)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    rn1        |none          |first random number
!!    rn2        |none          |second random number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name       |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt, Log, Cos

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   implicit none
   real*8, intent (in) :: rn1, rn2

   r_dstn1 = Sqrt(-2. * Log(rn1)) * Cos(6.283185 * rn2)

   return
end
