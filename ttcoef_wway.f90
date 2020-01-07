!> @file ttcoef_wway.f90
!> file containing the subroutine ttcoef_wway
!> @author
!> modified by Javier Burguete

!> this subroutine computes travel time coefficients for routing
!> along the main channel - grassed waterways
subroutine ttcoef_wway(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    grwat_d(:)  |m             |average depth of main channel
!!    grwat_l(:)  |km            |length of main channel
!!    grwat_n(:)  |none          |Manning's "n" value for the main channel
!!    grwat_s(:)  |m/m           |average slope of main channel
!!    grwat_w(:)  |m             |average width of main channel
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name           |units     |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    wat_phi1(:)    |m^2       |cross-sectional area of flow at bankfull
!!                              |depth
!!    wat_phi2(:)    |none      |
!!    wat_phi3(:)    |none      |
!!    wat_phi4(:)    |none      |
!!    wat_phi5(:)    |m^3/s     |flow rate when reach is at bankfull depth
!!    wat_phi6(:)    |m         |bottom width of main channel
!!    wat_phi7(:)    |m         |depth of water when reach is at bankfull
!!    wat_phi8(:)    |m/s       |average velocity when reach is at
!!                              |bankfull depth
!!    wat_phi9(:)    |m/s       |wave celerity when reach is at
!!                              |bankfull depth
!!    wat_phi10(:)   |hr        |storage time constant for reach at
!!                              |bankfull depth (ratio of storage to
!!                              |discharge)
!!    wat_phi11(:)   |m/s       |average velocity when reach is at
!!                              |0.1 bankfull depth (low flow)
!!    wat_phi12(:)   |m/s       |wave celerity when reach is at
!!                              |0.1 bankfull depth (low flow)
!!    wat_phi13(:)   |hr        |storage time constant for reach at
!!                              |0.1 bankfull depth (low flow) (ratio
!!                              |of storage to discharge)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    a           |m^2           |cross-sectional area of channel
!!    aa          |none          |area/area=1 (used to calculate velocity with
!!                               |Manning's equation)
!!    b           |m             |bottom width of channel
!!    chsslope(:) |none          |change in horizontal distance per unit
!!                               |change in vertical distance on channel side
!!                               |slopes; always set to 2 (slope=1/2)
!!    d           |m             |depth of flow
!!    fps         |none          |change in horizontal distance per unit
!!                               |change in vertical distance on floodplain side
!!                               |slopes; always set to 4 (slope=1/4)
!!    p           |m             |wetting perimeter
!!    qq1         |m^3/s         |flow rate for a specified depth
!!    rh          |m             |hydraulic radius of channel
!!    tt1         |km s/m        |time coefficient for specified depth
!!    tt2         |km s/m        |time coefficient for bankfull depth
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt
!!    SWAT: Qman

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
   use parm
   implicit none

   real*8 Qman
   integer, intent(in) :: j
   real*8, parameter :: aa = 1.! , fps = 4. ! not used
   real*8 :: a, b, chsslope, d, p, rh

!!    If side slope is not set in .rte file then assume this default
!!    If it is main reach default side slope to 2:1 if it is a waterway default to 8:1
   chsslope = 8.

   d = grwat_d(j)
   b = grwat_w(j) - 2. * d * chsslope


!!    check if bottom width (b) is < 0
   if (b <= 0.) then
      b = .5 * grwat_w(j)
      chsslope = (grwat_w(j) - b) / (2. * d)
   end if
   wat_phi6(j) = b
   ! wat_phi7(j) = d not used

!!    compute flow and travel time at bankfull depth
   p = b + 2. * d * Sqrt(chsslope * chsslope + 1.)
   a = b * d + chsslope * d * d
   rh = a / p
   wat_phi1(j) = a
   wat_phi5(j) = Qman(a, rh, grwat_n(j), grwat_s(j))
   wat_phi9(j) = Qman(aa, rh, grwat_n(j), grwat_s(j)) * 5. / 3.
   ! wat_phi8(j) = Qman(aa, rh, grwat_n(j), grwat_s(j)) ! not used
   ! wat_phi9(j) = wat_phi8(j) * 5. / 3. ! not used
   ! wat_phi10(j) = grwat_l(j) / wat_phi9(j) / 3.6 ! not used
   ! tt2 = grwat_l(j) * a / wat_phi5(j) ! not used

!!    compute flow and travel time at 1.2 bankfull depth
   ! d = 1.2 * grwat_d(j) ! not used
   ! a = a + (grwat_w(j) * grwat_d(j) + fps * (d - grwat_d(j)) ** 2) ! not used
   ! p = p + 4.*grwat_w(j) + (0.4 * grwat_d(j) * Sqrt(fps * fps + 1.)) ! not used
   ! rh = a / p ! not used
   ! qq1 = Qman(a, rh, grwat_n(j), grwat_s(j)) ! not used
   ! tt1 = grwat_l(j) * a / qq1 ! not used

!!    compute flow and travel time at 0.1 bankfull depth
   ! d = 0.1 * grwat_d(j) ! not used
   ! p = b + 2. * d * Sqrt(chsslope * chsslope + 1.) ! not used
   ! a = b * d + chsslope * d * d ! not used
   ! rh = a / p ! not used
   ! qq1 = Qman(a, rh, grwat_n(j), grwat_s(j)) ! not used
   ! tt1 = grwat_l(j) * a / qq1 ! not used
   ! wat_phi11(j) = Qman(aa, rh, grwat_n(j), grwat_s(j)) ! not used
   ! wat_phi12(j) = wat_phi11(j) * 5. / 3. not used
   ! wat_phi13(j) = grwat_l(j) / wat_phi12(j) / 3.6 ! not used

   return
end
