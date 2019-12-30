!> @file subwq.f90
!> file containing the subroutine subwq
!> @author
!> modified by Javier Burguete

!> this subroutine computes HRU loadings of chlorophyll-a, CBOD,
!> and dissolved oxygen to the main channel
!> @param[in] j HRU number (none)
subroutine subwq(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    enratio     |none          |enrichment ratio calculated for day in HRU
!!    hru_km(:)   |km^2          |area of HRU in square kilometers
!!    qtile       |mm H2O        |drainage tile flow in soil layer for the
!!                               |day
!!    sedorgn(:)  |kg N/ha       |amount of organic nitrogen in surface runoff
!!                               |in HRU for the day
!!    sedorgp(:)  |kg P/ha       |amount of organic phosphorus in surface runoff
!!                               |in HRU for the day
!!    sedyld(:)   |metric tons   |daily soil loss caused by water erosion in
!!                               |HRU
!!    sol_cbn(:,:)|%             |percent organic carbon in soil layer
!!    surfq(:)    |mm H2O        |surface runoff generated on day in HRU
!!    surqno3(:)  |kg N/ha       |amount of NO3-N in surface runoff in HRU for
!!                               |the day
!!    surqsolp(:) |kg P/ha       |amount of soluble phosphorus in surface runoff
!!                               |in HRU for the day
!!    t_ov(:)     |hr            |time for flow from farthest point in subbasin
!!                               |to enter a channel
!!    tmpav(:)    |deg C         |average air temperature on current day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cbodu(:)    |mg/L          |carbonaceous biological oxygen demand of
!!                               |surface runoff on current day in HRU
!!    chl_a(:)    |microgram/L   |chlorophyll-a concentration in water yield
!!                               |on current day in HRU
!!    doxq(:)     |mg/L          |dissolved oxygen concentration in the surface
!!                               |runoff on current day in HRU
!!    soxy        |mg/L          |dissolved oxygen saturation concentration
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    wtmp        |deg K         |temperature of surface runoff
!!    v1          |none          |variable to hold intermediate calculation
!!                               |result
!!    v2          |none          |variable to hold intermediate calculation
!!                               |result
!!    xx          |none          |variable to hold intermediate calculation
!!                               |result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max, Exp, Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: v1, v2, wtmp, xx

   if (qdr(j) > 1.e-4) then

      !! calculcate water temperature
      !! Stefan and Preudhomme. 1993.  Stream temperature estimation
      !!from air temperature.  Water Res. Bull. p. 27-45
      !! SWAT manual 2.3.13
      wtmp = Max(0.1, 5.0 + 0.75 * tmpav(j))
      wtmp = wtmp + 273.15    !! deg C to deg K

      xx = 100. * (sedorgn(j) + surqno3(j)) / qdr(j)   !100*kg/ha/mm = ppm
      chl_a(j) = chla_subco * xx

      !! calculate organic carbon loading to main channel
      if (cswat == 2) then
         xx = sedc_d(j) * hru_ha(j)
      else
         xx = sol_cbn(1,j) * enratio * sedyld(j) * 10.
      end if


      !! calculate carbonaceous biological oxygen demand (CBOD)
      cbodu(j) = cbodu(j) + 2.7 * xx / (qdr(j) * hru_km(j)) !jaehak 2016

      !! calculate dissolved oxygen saturation concentration
      !! QUAL2E equation III-29
      v1 = 1. / wtmp
      v2 = v1
      xx = -139.34410 + 1.575701E05 * v2
      v2 = v2 * v1
      xx = xx - 6.642308E07 * v2
      v2 = v2 * v1
      xx = xx + 1.243800E10 * v2
      v2 = v2 * v1
      xx = xx - 8.621949E11 * v2
      soxy = Exp(xx)

      !! calculate actual dissolved oxygen concentration
      doxq(j) = Min(soxy, soxy * Exp(-0.1 * cbodu(j)))
   else
      chl_a(j) = 0.
      cbodu(j) = 0.
      doxq(j) = 0.
   end if

   return
end
