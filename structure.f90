!> @file structure.f90
!> file containing the subroutine structure
!> @author
!> A. Van Griensven, Hydrology-Vrije Universiteit Brussel, Belgium.\n
!> Modified by Javier Burguete

!> this subroutine adjusts dissolved oxygen content for aeration at
!> structures.
!> @param[in] k reach number (none)
subroutine structure(k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k            |none          |reach number
!!    hhvaroute(:,:)|varies       |hourly routing storage array
!!    ievent       |none          |rainfall/runoff code
!!                                |0 daily rainfall/curve number technique
!!                                |1 sub-daily rainfall/Green&Ampt/hourly
!!                                |  routing
!!                                |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ihout        |none          |hydrograph storage location number for
!!                                |output
!!    mvaro        |none          |max number of variables routed through the
!!                                |reach
!!    rnum1        |none          |aeration coefficient
!!    varoute(:,:) |varies        |daily routing storage array
!!    varoute(1,:) |deg C         |water temperature
!!    varoute(2,:) |m^3 H2O       |water
!!    varoute(17,:)|kg O2         |dissolved oxygen
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhvaroute(:,:)|varies      |hourly routing storage array
!!    soxy        |mg O2/L       |saturation concentration of dissolved oxygen
!!    varoute(:,:)|varies        |daily routing storage array
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    disoxin     |mg O2/L       |dissolved oxygen concentration
!!    ii          |none          |counter
!!    jj          |none          |counter
!!    reak        |none          |aeration coefficient
!!    wtmp        |deg C         |water temperature
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: Oxygen_saturation

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 Oxygen_saturation
   integer, intent(in) :: k
   real*8 :: disoxin, reak, wtmp
   integer :: ii, jj

!! initialize variables
   reak = rnum1
   if (reak <= 1.) reak = 1.
   soxy = 0.

!! daily array
   do ii = 1, mvaro
      varoute(ii,ihout) = varoute(ii,k)
   end do
   if (varoute(2,k) > 0.001) then
      wtmp = varoute(1,k)
      !! calculate saturation concentration for dissolved oxygen
      !! QUAL2E section 3.6.1 equation III-29
      soxy = Oxygen_saturation(wtmp)
      disoxin = varoute(17,k) * 1000. / varoute(2,k)
      disoxin = soxy - ((soxy - disoxin) / reak)
      if (disoxin < 0.) disoxin = 0.
      varoute(17,ihout) = disoxin * varoute(2,k) / 1000.
   else
      varoute(17,ihout)=0.
   end if

!! subdaily array
   if (ievent > 0) then
      do ii = 1, nstep
         do jj = 1, mvaro
            hhvaroute(jj,ihout,ii) = hhvaroute(jj,k,ii)
         end do
         soxy = 0.
         if (hhvaroute(2,k,ii) > 0.0001) then
            wtmp = hhvaroute(1,k,ii)
            !! calculate saturation concentration for dissolved oxygen
            !! QUAL2E section 3.6.1 equation III-29
            soxy = Oxygen_saturation(wtmp)
            disoxin = hhvaroute(17,k,ii) * 1000. / hhvaroute(2,k,ii)
            disoxin = soxy - ((soxy - disoxin) / reak)
            if (disoxin < 0.) disoxin = 0.
            hhvaroute(17,ihout,ii) = disoxin * hhvaroute(2,k,ii) / 1000.
         end if
      end do
   end if

   return
end
