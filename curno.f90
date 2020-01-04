!> @file curno.f90
!> file containing the subroutine curno
!> @author
!> modified by Javier Burguete

!> this subroutine determines the curve numbers for moisture conditions
!> I and III and calculates coefficents and shape parameters for the
!> water retention curve.
!> The coefficents and shape parameters are calculated by one of two methods:\n
!> the default method is to make them a function of soil water,\n
!> the alternative method (labeled new) is to make them a function of
!> accumulated PET, precipitation and surface runoff
!> @param[in] cnn SCS runoff curve number for moisture condition II
!> @param[in] h HRU number
subroutine curno(cnn,h)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cnn         |none          |SCS runoff curve number for moisture
!!                               |condition II
!!    h           |none          |HRU number
!!    sol_sumfc(:)|mm H2O        |amount of water held in soil profile at
!!                               |field capacity
!!    sol_sumul(:)|mm H2O        |amount of water held in soil profile at
!!                               |saturation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cn1(:)      |none          |SCS runoff curve number for moisture
!!                               |condition I
!!    cn3(:)      |none          |SCS runoff curve number for moisture
!!                               |condition III
!!    sci(:)      |none          |retention coefficient for cn method based on
!!                               |plant ET
!!    smx(:)      |none          |retention coefficient for cn method based on
!!                               |soil moisture
!!    wrt1(:)    |none          |1st shape parameter for calculation of
!!                               |water retention
!!    wrt2(:)    |none          |2nd shape parameter for calculation of
!!                               |water retention
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    c2          |none          |variable used to hold calculated value
!!    rto3        |none          |fraction difference between CN3 and CN1
!!                               |retention parameters
!!    rtos        |none          |fraction difference between CN=99 and CN1
!!                               |retention parameters
!!    s3          |none          |retention parameter for CN3
!!    sumfc_ul
!!    smxold
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max
!!    SWAT: ascrv

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8, intent (in) :: cnn
   integer, intent (in) :: h
   real*8 :: c2, rto3, rtos, s3, sumfc_ul, smxold

   cn2(h) = cnn
   smxold = 0.
   if (cn1(h) > 1.e-6) smxold = 254.* (100. / cn1(h) - 1.)

!! calculate moisture condition I and III curve numbers
   c2 = 100. - cnn
   cn1(h) = cnn - 20. * c2 / (c2 + Exp(2.533 - 0.0636 * c2))
   cn1(h) = Max(cn1(h), .4 * cnn)
   cn3(h) = cnn * Exp(.006729 * c2)

!! calculate maximum retention parameter value
   smx(h) = 254. * (100. / cn1(h) - 1.)

!! calculate retention parameter value for CN3
   s3 = 254. * (100. / cn3(h) - 1.)

!! calculate fraction difference in retention parameters
   rto3 = 1. - s3 / smx(h)
   rtos = 1. - 2.54 / smx(h)

   sumfc_ul = sol_sumfc(h)
!! calculate shape parameters
   call ascrv(rto3,rtos,sumfc_ul,sol_sumul(h),wrt1(h),wrt2(h))

   if (curyr == 0) then
      sci(h) = 0.9 * smx(h)
   else
      sci(h) = (1. - ((smxold - sci(h)) / smxold)) * smx(h)      !! plant ET
   end if

   return
end
