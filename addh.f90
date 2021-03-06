!> @file addh.f90
!> file containing the subroutine addh
!> @author
!> modified by Javier Burguete

!> this subroutine adds loadings from two sources for routing
!> @param[in] j
!> hydrograph storage location number of first dataset to be added (none)
!> @param[in] k
!> inflow hydrograph storage location number of second dataset to be added
!> (none)
subroutine addh(j, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |hydrograph storage location number of
!!                               |first dataset to be added
!!    k           |none          |hydrograph storage location number of
!!                               |second dataset to be added
!!    hhvaroute(:,:,:)|varies    |routing storage array for hourly time step
!!    ievent      |none          |rainfall/runoff code
!!                               |0 daily rainfall/curve number technique
!!                               |1 sub-daily rainfall/Green&Ampt/hourly
!!                               |  routing
!!                               |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ihout       |none          |outflow hydrograph storage location number
!!    mvaro       |none          |max number of variables routed through the
!!                               |reach
!!    varoute(:,:)|varies        |routing storage array
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhvaroute(:,:,:)|varies        |routing storage array for hourly time step
!!    varoute(:,:)    |varies        |routing storage array for daily time step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name       |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii         |none          |counter
!!    kk         |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j, k
   integer :: ii, kk

!! zero flow out variables
   do ii = 1, mvaro
      varoute(ii,ihout) = 0.
   end do

!! add loadings and store in new hydrograph location
   if (varoute(2,j) + varoute(2,k) > 0.1) then
      varoute(1,ihout) = (varoute(1,j) * varoute(2,j) +&
         &varoute(1,k) * varoute(2,k)) / (varoute(2,j) + varoute(2,k))
      varoute(18,ihout) = (varoute(18,j) * varoute(2,j) +&
         &varoute(18,k) * varoute(2,k)) / (varoute(2,j) + varoute(2,k))
      varoute(19,ihout) = (varoute(19,j) * varoute(2,j) +&
         &varoute(19,k) * varoute(2,k)) / (varoute(2,j) + varoute(2,k))
   end if
   do ii = 2, 17
      varoute(ii,ihout) = varoute(ii,j) + varoute(ii,k)
   end do
   do ii = 20, mvaro
      varoute(ii,ihout) = varoute(ii,j) + varoute(ii,k)
   end do
   if (varoute(2,ihout) < 0.) varoute(2,ihout) = 0.


!! add hydrograph points (hourly time step)
   if (ievent > 0) then
      do kk = 1, nstep  ! modified for urban modeling by J.Jeong 4/15/2008
         if (hhvaroute(2,j,kk) + hhvaroute(2,k,kk) > 0.1) then
            hhvaroute(1,ihout,kk) = (hhvaroute(1,j,kk) *&
               &hhvaroute(2,j,kk) + hhvaroute(1,k,kk) *&
               &hhvaroute(2,k,kk)) / (hhvaroute(2,j,kk) + hhvaroute(2,k,kk))
         end if
      end do
      do ii = 2, mvaro
         do kk = 1, nstep  ! modified for urban modeling by J.Jeong 4/15/2008
            hhvaroute(ii,ihout,kk) = hhvaroute(ii,j,kk) + hhvaroute(ii,k,kk)
         end do
      end do

      DO kk = 1, nstep
         QHY(kk,ihout,IHX(1))=QHY(kk,j,IHX(1))+QHY(kk,k,IHX(1)) !flood routing jaehak 2017
      END DO

   endif

   do ii = 29, mvaro
      varoute(ii,j) = 0.
      varoute(ii,k) = 0.
   end do

   return
end
