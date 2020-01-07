!> @file reccnst.f90
!> file containing the subroutine reccnst
!> @author
!> modified by Javier Burguete

!> this subroutine inputs measured loadings to the stream network
!> for routing through the watershed where the records are averaged
!> over the entire period of record
!> @param[in] k file number (none)
subroutine reccnst(k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k           |none          |file number
!!    bactlpcnst(:)|# cfu/100ml  |average daily less persistent bacteria
!!                               |loading to reach
!!    bactpcnst(:)|# cfu/100ml   |average daily persistent bacteria loading
!!                               |to reach
!!    cbodcnst(:) |kg/day        |average daily CBOD loading to reach
!!    chlacnst(:) |kg/day        |average daily chlorophyll-a loading to reach
!!    cmtlcnst(1,:)|kg/day        |average daily conservative metal #1 loading
!!    cmtlcnst(2,:)|kg/day        |average daily conservative metal #2 loading
!!    cmtlcnst(3,:)|kg/day        |average daily conservative metal #3 loading
!!    disoxcnst(:)|kg/day        |average daily dissolved oxygen loading to reach
!!    flocnst(:)  |m^3 H2O/day   |average daily water loading to reach
!!    ievent      |none          |rainfall/runoff code
!!                               |0 daily rainfall/curve number technique
!!                               |1 sub-daily rainfall/Green&Ampt/hourly
!!                               |  routing
!!                               |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ihout       |none          |hydrograph storage location number
!!    mvaro       |none          |max number of variables routed through the
!!                               |reach
!!    minpcnst(:) |kg P/day      |average daily soluble P loading to reach
!!    nh3cnst(:)  |kg N/day      |average daily ammonia loading to reach
!!    no2cnst(:)  |kg N/day      |average daily nitrite loading to reach
!!    no3cnst(:)  |kg N/day      |average daily nitrate loading to reach
!!    orgncnst(:) |kg N/day      |average daily organic N loading to reach
!!    orgpcnst(:) |kg P/day      |average daily organic P loading to reach
!!    sedcnst(:)  |metric tons/d |average daily sediment loading to reach
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
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    fn          |none          |number of time steps (float)
!!    ii          |none          |counter
!!    j           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Dfloat

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: k
   real*8 :: fn
   integer :: ii, j

!! zero flow out variables
   do j = 1, mvaro
      varoute(j,ihout) = 0.
      if (ievent > 0) then
         do ii = 1, nstep
            hhvaroute(j,ihout,ii) = 0.
         end do
      endif
   end do

   varoute(2,ihout) = flocnst(k)
   varoute(3,ihout) = sedcnst(k)
   varoute(4,ihout) = orgncnst(k)
   varoute(5,ihout) = orgpcnst(k)
   varoute(6,ihout) = no3cnst(k)
   varoute(7,ihout) = minpcnst(k)
   varoute(11,ihout) = solpstcnst(k)
   varoute(12,ihout) = srbpstcnst(k)
   varoute(13,ihout) = chlacnst(k)
   varoute(14,ihout) = nh3cnst(k)
   varoute(15,ihout) = no2cnst(k)
   varoute(16,ihout) = cbodcnst(k)
   varoute(17,ihout) = disoxcnst(k)
   varoute(18,ihout) = bactpcnst(k)
   varoute(19,ihout) = bactlpcnst(k)
   varoute(20,ihout) = cmtlcnst(1,k)
   varoute(21,ihout) = cmtlcnst(2,k)
   varoute(22,ihout) = cmtlcnst(3,k)

   !! Assumed equal distribution of sediment
   varoute(23,ihout) = sedcnst(k) * 0.   ! sand
   varoute(24,ihout) = sedcnst(k) * 1.   ! silt
   varoute(25,ihout) = sedcnst(k) * 0.   ! cla
   varoute(26,ihout) = sedcnst(k) * 0.   ! sag
   varoute(27,ihout) = sedcnst(k) * 0.   ! lag
   varoute(28,ihout) = 0.                    ! gravel

   if (ievent > 0) then
      fn = Dfloat(nstep)
      do ii = 1,nstep
         hhvaroute(2,ihout,ii) = flocnst(k) / fn
         hhvaroute(3,ihout,ii) = sedcnst(k) / fn
         hhvaroute(4,ihout,ii) = orgncnst(k) / fn
         hhvaroute(5,ihout,ii) = orgpcnst(k) / fn
         hhvaroute(6,ihout,ii) = no3cnst(k) / fn
         hhvaroute(7,ihout,ii) = minpcnst(k) / fn
         hhvaroute(11,ihout,ii) = solpstcnst(k) / fn
         hhvaroute(12,ihout,ii) = srbpstcnst(k) / fn
         hhvaroute(13,ihout,ii) = chlacnst(k) / fn
         hhvaroute(14,ihout,ii) = nh3cnst(k) / fn
         hhvaroute(15,ihout,ii) = no2cnst(k) / fn
         hhvaroute(16,ihout,ii) = cbodcnst(k) / fn
         hhvaroute(17,ihout,ii) = disoxcnst(k) / fn
         hhvaroute(18,ihout,ii) = bactpcnst(k) / fn
         hhvaroute(19,ihout,ii) = bactlpcnst(k) / fn
         hhvaroute(20,ihout,ii) = cmtlcnst(1,k) / fn
         hhvaroute(21,ihout,ii) = cmtlcnst(2,k) / fn
         hhvaroute(22,ihout,ii) = cmtlcnst(3,k) / fn
      end do
   end if


   return
end
