!> @file recmon.f90
!> file containing the subroutine recmon
!> @author
!> modified by Javier Burguete

!> this subroutine inputs measured loadings to the stream network
!> for routing through the watershed where the records are summarized
!> on a monthly basis
!> @param[in] k file number (none)
subroutine recmon(k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k             |none         |file number
!!    bactlpmon(:,:,:)|# cfu/100ml|average amount of less persistent bacteria
!!                                |loaded to stream on a given day in the
!!                                |month
!!    bactpmon(:,:,:)|# cfu/100ml |average amount of persistent bacteria
!!                                |loaded to stream on a given day in the
!!                                |month
!!    cbodmon(:,:,:) |kg/day      |average amount of CBOD loaded to stream
!!                                |on a given day in the month
!!    chlamon(:,:,:) |kg/day      |average amount of chlorophyll a loaded
!!                                |to stream on a given day in the month
!!    cmtlmon(1,:,:,:)|kg/day      |average amount of conservative metal #1
!!                                |loaded to stream on a given day in the
!!                                |month
!!    cmtlmon(2,:,:,:)|kg/day      |average amount of conservative metal #2
!!                                |loaded to stream on a given day in the
!!                                |month
!!    cmtlmon(3,:,:,:)|kg/day      |average amount of conservative metal #3
!!                                |loaded to stream on a given day in the
!!                                |month
!!    curyr         |none         |year of simulation
!!    disoxmon(:,:,:)|kg/day      |average amount of dissolved oxygen loaded to
!!                                |stream on a given day in the month
!!    ievent        |none         |rainfall/runoff code
!!                                |0 daily rainfall/curve number technique
!!                                |1 sub-daily rainfall/Green&Ampt/hourly
!!                                |  routing
!!                                |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ihout         |none         |hydrograph storage location number
!!    flomon(:,:,:) |m^3/day      |average amount of water loaded to stream
!!                                |on a given day in the month
!!    minpmon(:,:,:)|kg P/day     |average amount of soluble P loaded to
!!                                |stream on a given day in the month
!!    i_mo          |none         |month of simulation
!!    mvaro         |none         |max number of variables routed through the
!!                                |reach
!!    nh3mon(:,:,:) |kg N/day     |average amount of NH3-N loaded to
!!                                |stream on a given day in the month
!!    no2mon(:,:,:) |kg N/day     |average amount of NO2-N loaded to
!!                                |stream on a given day in the month
!!    no3mon(:,:,:) |kg N/day     |average amount of NO3-N loaded to
!!                                |stream on a given day in the month
!!    orgnmon(:,:,:)|kg N/day     |average amount of organic N loaded to
!!                                |stream on a given day in the month
!!    orgpmon(:,:,:)|kg P/day     |average amount of organic P loaded to
!!                                |stream on a given day in the month
!!    sedmon(:,:,:) |metric tons/d|average amount of sediment loaded to
!!                                |stream on a given day in the month
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name             |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhvaroute(2,:,:) |m^3          |volume of water
!!    hhvaroute(3,:,:) |metric tons  |sediment
!!    hhvaroute(4,:,:) |kg N         |organic N
!!    hhvaroute(5,:,:) |kg P         |organic P
!!    hhvaroute(6,:,:) |kg N         |NO3-N
!!    hhvaroute(7,:,:) |kg P         |mineral (soluble) P
!!    hhvaroute(13,:,:)|kg           |chlorophyll-a
!!    hhvaroute(14,:,:)|kg N         |NH3
!!    hhvaroute(15,:,:)|kg N         |NO2
!!    hhvaroute(16,:,:)|kg           |carbonaceous biological oxygen demand
!!    hhvaroute(17,:,:)|kg           |dissolved oxygen
!!    hhvaroute(18,:,:)|# cfu/100ml  |persistent bacteria
!!    hhvaroute(19,:,:)|# cfu/100ml  |less persistent bacteria
!!    hhvaroute(20,:,:)|kg           |conservative metal #1
!!    hhvaroute(21,:,:)|kg           |conservative metal #2
!!    hhvaroute(22,:,:)|kg           |conservative metal #3
!!    varoute(2,:)     |m^3          |volume of water
!!    varoute(3,:)     |metric tons  |sediment
!!    varoute(4,:)     |kg N         |organic N
!!    varoute(5,:)     |kg P         |organic P
!!    varoute(6,:)     |kg N         |NO3-N
!!    varoute(7,:)     |kg P         |mineral (soluble) P
!!    varoute(13,:)    |kg           |chlorophyll-a
!!    varoute(14,:)    |kg N         |NH3
!!    varoute(15,:)    |kg N         |NO2
!!    varoute(16,:)    |kg           |CBOD
!!    varoute(17,:)    |kg           |dissolved oxygen
!!    varoute(18,:)    |# cfu/100ml  |persistent bacteria
!!    varoute(19,:)    |# cfu/100ml  |less persistent bacteria
!!    varoute(20,:)    |kg           |conservative metal #1
!!    varoute(21,:)    |kg           |conservative metal #2
!!    varoute(22,:)    |kg           |conservative metal #3
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name       |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    fn         |none          |number of time steps (float)
!!    ii         |none          |counter
!!    j          |none          |counter
!!    xx         |none          |auxiliar variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Dfloat

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: k
   real*8 :: fn, xx
   integer :: ii, j

!! zero flow out variables
   do j = 1, mvaro
      varoute(j,ihout) = 0.
      do ii = 1, nstep
         hhvaroute(j,ihout,ii) = 0.
      end do
   end do

   varoute(2,ihout) = flomon(k,curyr,i_mo)
   varoute(3,ihout) = sedmon(k,curyr,i_mo)
   varoute(4,ihout) = orgnmon(k,curyr,i_mo)
   varoute(5,ihout) = orgpmon(k,curyr,i_mo)
   varoute(6,ihout) = no3mon(k,curyr,i_mo)
   varoute(7,ihout) = minpmon(k,curyr,i_mo)
   varoute(11,ihout) = solpstmon(k,curyr,i_mo)
   varoute(12,ihout) = srbpstmon(k,curyr,i_mo)
   varoute(13,ihout) = chlamon(k,curyr,i_mo)
   varoute(14,ihout) = nh3mon(k,curyr,i_mo)
   varoute(15,ihout) = no2mon(k,curyr,i_mo)
   varoute(16,ihout) = cbodmon(k,curyr,i_mo)
   varoute(17,ihout) = disoxmon(k,curyr,i_mo)
   varoute(18,ihout) = bactpmon(k,curyr,i_mo)
   varoute(19,ihout) = bactlpmon(k,curyr,i_mo)
   varoute(20,ihout) = cmtlmon(1,k,curyr,i_mo)
   varoute(21,ihout) = cmtlmon(2,k,curyr,i_mo)
   varoute(22,ihout) = cmtlmon(3,k,curyr,i_mo)

   !! Assumed equal distribution of sediment
   xx = sedmon(k,curyr,i_mo)
   varoute(23,ihout) = xx * 0.   ! sand
   varoute(24,ihout) = xx * 1.   ! silt
   varoute(25,ihout) = xx * 0.   ! cla
   varoute(26,ihout) = xx * 0.   ! sag
   varoute(27,ihout) = xx * 0.   ! lag
   varoute(28,ihout) = 0.                    ! gravel

   if (ievent > 0) then
      fn = Dfloat(nstep)
      do ii = 1, nstep
         hhvaroute(2,ihout,ii) = flomon(k,curyr,i_mo) / fn
         hhvaroute(3,ihout,ii) = xx / fn
         hhvaroute(4,ihout,ii) = orgnmon(k,curyr,i_mo) / fn
         hhvaroute(5,ihout,ii) = orgpmon(k,curyr,i_mo) / fn
         hhvaroute(6,ihout,ii) = no3mon(k,curyr,i_mo) / fn
         hhvaroute(7,ihout,ii) = minpmon(k,curyr,i_mo) / fn
         hhvaroute(11,ihout,ii) = solpstmon(k,curyr,i_mo) / fn
         hhvaroute(12,ihout,ii) = srbpstmon(k,curyr,i_mo) / fn
         hhvaroute(13,ihout,ii) = chlamon(k,curyr,i_mo) / fn
         hhvaroute(14,ihout,ii) = nh3mon(k,curyr,i_mo) / fn
         hhvaroute(15,ihout,ii) = no2mon(k,curyr,i_mo) / fn
         hhvaroute(16,ihout,ii) = cbodmon(k,curyr,i_mo) / fn
         hhvaroute(17,ihout,ii) = disoxmon(k,curyr,i_mo) / fn
         hhvaroute(18,ihout,ii) = bactpmon(k,curyr,i_mo) / fn
         hhvaroute(19,ihout,ii) = bactlpmon(k,curyr,i_mo) / fn
         hhvaroute(20,ihout,ii) = cmtlmon(1,k,curyr,i_mo) / fn
         hhvaroute(21,ihout,ii) = cmtlmon(2,k,curyr,i_mo) / fn
         hhvaroute(22,ihout,ii) = cmtlmon(3,k,curyr,i_mo) / fn
      end do
   end if

   return
end
