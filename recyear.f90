!> @file recyear.f90
!> file containing the subroutine recyear
!> @author
!> modified by Javier Burguete

!> this subroutine inputs measured loadings to the stream network
!> for routing through the watershed where the records are summarized
!> on an annual basis
!> @param[in] k file number (none)
subroutine recyear(k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k           |none          |file number
!!    bactlpyr(:,:)|# cfu/100ml  |average daily loading of less persistent
!!                               |bacteria for year
!!    bactpyr(:,:)|# cfu/100ml   |average daily loading of persistent bacteria
!!                               |for year
!!    cbodyr(:,:) |kg/day        |average daily loading of CBOD for year
!!    chlayr(:,:) |kg/day        |average daily loading of chlorophyll-a for year
!!    cmtlyr(1,:,:)|kg/day        |average daily loading of conservative metal #1
!!                               |for year
!!    cmtlyr(2,:,:)|kg/day        |average daily loading of conservative metal #2
!!                               |for year
!!    cmtlyr(3,:,:)|kg/day        |average daily loading of conservative metal #3
!!                               |for year
!!    curyr       |none          |year of simulation
!!    disoxyr(:,:)|kg/day        |average daily loading of dissolved oxygen for
!!                               |year
!!    floyr(:,:)  |m**3/d        |average daily water loading for year
!!    ihout       |none          |hydrograph storage location number
!!    minpyr(:,:) |kg P/day      |average daily mineral P loading for year
!!    mvaro       |none          |max number of variables routed through the
!!                               |reach
!!    nh3yr(:,:)  |kg N/day      |average daily NH3-N loading for year
!!    no2yr(:,:)  |kg N/day      |average daily NO2-N loading for year
!!    no3yr(:,:)  |kg N/day      |average daily NO3-N loading for year
!!    orgnyr(:,:) |kg N/day      |average daily organic N loading for year
!!    orgpyr(:,:) |kg P/day      |average daily organic P loading for year
!!    sedyr(:,:)  |metric tons/d |average daily sediment loading for year
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name             |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    fn          |none          |number of time steps (float)
!!    ii          |none          |counter
!!    j           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

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
      do ii = 1, nstep
         hhvaroute(j,ihout,ii) = 0.
      end do
   end do


   varoute(2,ihout) = floyr(k,curyr)
   varoute(3,ihout) = sedyr(k,curyr)
   varoute(4,ihout) = orgnyr(k,curyr)
   varoute(5,ihout) = orgpyr(k,curyr)
   varoute(6,ihout) = no3yr(k,curyr)
   varoute(7,ihout) = minpyr(k,curyr)
   varoute(11,ihout) = solpstyr(k,curyr)
   varoute(12,ihout) = srbpstyr(k,curyr)
   varoute(13,ihout) = chlayr(k,curyr)
   varoute(14,ihout) = nh3yr(k,curyr)
   varoute(15,ihout) = no2yr(k,curyr)
   varoute(16,ihout) = cbodyr(k,curyr)
   varoute(17,ihout) = disoxyr(k,curyr)
   varoute(18,ihout) = bactpyr(k,curyr)
   varoute(19,ihout) = bactlpyr(k,curyr)
   varoute(20,ihout) = cmtlyr(1,k,curyr)
   varoute(21,ihout) = cmtlyr(2,k,curyr)
   varoute(22,ihout) = cmtlyr(3,k,curyr)

   !! Assumed equal distribution of sediment
   varoute(23,ihout) = sedyr(k,curyr) * 0.   ! sand
   varoute(24,ihout) = sedyr(k,curyr) * 1.   ! silt
   varoute(25,ihout) = sedyr(k,curyr) * 0.   ! cla
   varoute(26,ihout) = sedyr(k,curyr) * 0.   ! sag
   varoute(27,ihout) = sedyr(k,curyr) * 0.   ! lag
   varoute(28,ihout) = 0.                        ! gravel

   if (ievent > 0) then
      fn = Dfloat(nstep)
      do ii = 1, nstep
         hhvaroute(2,ihout,ii) = floyr(k,curyr) / fn
         hhvaroute(3,ihout,ii) = sedyr(k,curyr) / fn
         hhvaroute(4,ihout,ii) = orgnyr(k,curyr) / fn
         hhvaroute(5,ihout,ii) = orgpyr(k,curyr) / fn
         hhvaroute(6,ihout,ii) = no3yr(k,curyr) / fn
         hhvaroute(7,ihout,ii) = minpyr(k,curyr) / fn
         hhvaroute(11,ihout,ii) = solpstyr(k,curyr) / fn
         hhvaroute(12,ihout,ii) = srbpstyr(k,curyr) / fn
         hhvaroute(13,ihout,ii) = chlayr(k,curyr) / fn
         hhvaroute(14,ihout,ii) = nh3yr(k,curyr) / fn
         hhvaroute(15,ihout,ii) = no2yr(k,curyr) / fn
         hhvaroute(16,ihout,ii) = cbodyr(k,curyr) / fn
         hhvaroute(17,ihout,ii) = disoxyr(k,curyr) / fn
         hhvaroute(18,ihout,ii) = bactpyr(k,curyr) / fn
         hhvaroute(19,ihout,ii) = bactlpyr(k,curyr) / fn
         hhvaroute(20,ihout,ii) = cmtlyr(1,k,curyr) / fn
         hhvaroute(21,ihout,ii) = cmtlyr(2,k,curyr) / fn
         hhvaroute(22,ihout,ii) = cmtlyr(3,k,curyr) / fn
      end do
   end if

   return
end
