!> @file log_normal.f90
!> file containing the function log_normal
!> @author
!> modified by Javier Burguete

!> this function generates a random number from a lognormal distribution curve
!> for estimating constituent concentration in the effluent of urban bmps
!> given mean and standard deviation values.
!> Jaehak Jeong, 2017
!> @param[in] mu mean value
!> @param[in] standard deviation
!> @return value generated for distribution
real*8 function log_normal(mu, sig)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mu          |mg/l          |mean value
!!    sig         |mg/l          |standard deviation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    c
!!    PI          |none          |pi number
!!    r
!!    theta
!!    tmp
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: random_number, Sqrt, Log, Sin, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   implicit none

   real*8, intent(in) :: mu, sig
   real*8, parameter :: PI=3.141592653589793238462
   real*8 :: temp(2), c, r, theta
   call random_number(temp)
   r = Sqrt(-2.0d0 * Log(temp(1)))
   theta = 2.0d0 * PI * temp(2)
   c = mu + sig * r * Sin(theta)
   log_normal = Exp(c)
end
