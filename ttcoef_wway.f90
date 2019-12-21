!> @file ttcoef_wway.f90
!> file containing the subroutine ttcoef_wway
!> @author
!> modified by Javier Burguete

!> this subroutine computes travel time coefficients for routing
!> along the main channel - grassed waterways
subroutine ttcoef_wway

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
!!    k           |none          |dummy argument (HRU number)
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


   real*8, parameter :: aa = 1., fps = 4.
   real*8 :: a, b, chsslope, d, p, qq1, rh, tt1, tt2
   integer :: k

   k = ihru

!!    If side slope is not set in .rte file then assume this default
!!    If it is main reach default side slope to 2:1 if it is a waterway default to 8:1
   !if (chside(k) <= 1.e-6) then
   chsslope = 8.
   !else
   !   chsslope = chside(k)
   !end if

   d = grwat_d(k)
   b = grwat_w(k) - 2. * d * chsslope


!!    check if bottom width (b) is < 0
   if (b <= 0.) then
      b = .5 * grwat_w(k)
      chsslope = (grwat_w(k) - b) / (2. * d)
   end if
   wat_phi6(k) = b
   wat_phi7(k) = d

!!    compute flow and travel time at bankfull depth
   p = b + 2. * d * Sqrt(chsslope * chsslope + 1.)
   a = b * d + chsslope * d * d
   rh = a / p
   wat_phi1(k) = a
   wat_phi5(k) = Qman(a, rh, grwat_n(k), grwat_s(k))
   wat_phi8(k) = Qman(aa, rh, grwat_n(k), grwat_s(k))
   wat_phi9(k) = wat_phi8(k) * 5. / 3.
   wat_phi10(k) = grwat_l(k) / wat_phi9(k) / 3.6
   tt2 = grwat_l(k) * a / wat_phi5(k)

!!    compute flow and travel time at 1.2 bankfull depth
   d = 1.2 * grwat_d(k)
   a = a + (grwat_w(k) * grwat_d(k) + fps * (d - grwat_d(k)) ** 2)
   p=p + 4.*grwat_w(k) + (0.4 * grwat_d(k) * Sqrt(fps * fps + 1.))
   rh = a / p
   qq1 = Qman(a, rh, grwat_n(k), grwat_s(k))
   tt1 = grwat_l(k) * a / qq1

!!    compute flow and travel time at 0.1 bankfull depth
   d = 0.1 * grwat_d(k)
   p = b + 2. * d * Sqrt(chsslope * chsslope + 1.)
   a = b * d + chsslope * d * d
   rh = a / p
   qq1 = Qman(a, rh, grwat_n(k), grwat_s(k))
   tt1 = grwat_l(k) * a / qq1
   wat_phi11(k) = Qman(aa, rh, grwat_n(k), grwat_s(k))
   wat_phi12(k) = wat_phi11(k) * 5. / 3.
   wat_phi13(k) = grwat_l(k) / wat_phi12(k) / 3.6

   return
end
