!> @file apex_day.f90
!> file containing the subroutine apex_day
!> @author
!> modified by Javier Burguete

!> this subroutine inputs measured loadings to the stream network for
!> routing through the watershed where the records are summarized on a
!> daily basis
!> @param[in] i current day in simulation--loop counter (julian date)
!> @param[in] k reach number or file number (none)
subroutine apex_day(i, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |julian date   |current day in simulation--loop counter
!!    k           |none          |reach number or file number
!!    id1         |julian date   |first day of simulation in year
!!    ievent      |none          |rainfall/runoff code
!!                               |0 daily rainfall/curve number technique
!!                               |1 sub-daily rainfall/Green&Ampt/hourly
!!                               |  routing
!!                               |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    ifirstr(:)  |none          |measured data search code
!!                               |0 first day of measured data located in file
!!                               |1 first day of measured data not located in
!!                               |file
!!    ihout       |none          |hydrograph storage location number
!!    iyr         |year          |current year of simulation (actual year)
!!    mvaro       |none          |max number of variables routed through the
!!                               |reach
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
!!    ifirstr(:)       |none         |measured data search code
!!                                   |0 first day of measured data located in
!!                                   |file
!!                                   |1 first day of measured data not located in
!!                                   |file
!!    varoute(2,:)     |m^3          |volume of water
!!    varoute(3,:)     |metric tons  |sediment
!!    varoute(4,:)     |kg N         |organic N
!!    varoute(5,:)     |kg P         |organic P
!!    varoute(6,:)     |kg N         |NO3-N
!!    varoute(7,:)     |kg P         |mineral (soluble) P
!!    varoute(13,:)    |kg           |chlorophyll-a
!!    varoute(14,:)    |kg N         |NH3
!!    varoute(15,:)    |kg N         |NO2
!!    varoute(16,:)    |kg           |carbonaceous biological oxygen demand
!!    varoute(17,:)    |kg           |dissolved oxygen
!!    varoute(18,:)    |# bact       |persistent bacteria
!!    varoute(19,:)    |# bact       |less persistent bacteria
!!    varoute(20,:)    |kg           |conservative metal #1
!!    varoute(21,:)    |kg           |conservative metal #2
!!    varoute(22,:)    |kg           |conservative metal #3
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    flodaya     |m^3 H2O       |water loading to reach on day
!!    fn          |none          |number of time steps (float)
!!    idapa       |julian date   |julian date of record
!!    ii          |none          |counter
!!    iypa        |year          |year of record
!!    j           |none          |counter
!!    minpdaya    |kg N          |nitrite loading to reach on day
!!    no3daya     |kg N          |nitrate loading to reach on day
!!    orgndaya    |kg N          |organic N loading to reach on day
!!    orgpdaya    |kg P          |organic P loading to reach on day
!!    seddaya     |metric tons   |sediment loading to reach on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i, k
   real*8, dimension (mapex) :: flodaya, minpdaya, no3daya, orgndaya, orgpdaya,&
      &seddaya
   integer, dimension (mapex) :: idapa, iypa 
   real*8 :: fn
   integer :: ii, j

   flodaya = 0.
   minpdaya = 0.
   no3daya = 0.
   orgndaya = 0.
   orgpdaya = 0.
   seddaya = 0.
   idapa = 0
   iypa = 0

   do j = 1, mvaro
      varoute(j,ihout) = 0.
      if (ievent > 0) then
         do ii = 1, nstep
            hhvaroute(j,ihout,ii) = 0.
         end do
      endif
   end do

!!   read from apex measured file
!!   check if idap/iyp are same

!!    changes for Mike Winchell if apex.swt files has a different start date than the SWAT run
   if (ifirsta(k) == 1) then
      do
         read (112+k,*) idapa(k), iypa(k), flodaya(k),&
            &seddaya(k), orgndaya(k), orgpdaya(k), no3daya(k),&
            &minpdaya(k)
         if(idapa(k) == id1 .and. iypa(k) == iyr) exit
      end do
      ifirsta(k) = 0
   endif

   if (iypa(k) == iyr .and. idapa(k) == i) then
      varoute(2,ihout) = flodaya(k)
      varoute(3,ihout) = seddaya(k)
      varoute(4,ihout) = orgndaya(k)
      varoute(5,ihout) = orgpdaya(k)
      varoute(6,ihout) = no3daya(k)
      varoute(7,ihout) = minpdaya(k)
      varoute(14,ihout) = 0.0
      varoute(15,ihout) = 0.0
      varoute(11,ihout) = 0.0
      varoute(12,ihout) = 0.0
      varoute(13,ihout) = 0.0
      varoute(14,ihout) = 0.0
      varoute(15,ihout) = 0.0
      varoute(16,ihout) = 0.0
      varoute(17,ihout) = 0.0
      varoute(18,ihout) = 0.0
      varoute(19,ihout) = 0.0
      varoute(20,ihout) = 0.0
      varoute(21,ihout) = 0.0
      varoute(22,ihout) = 0.0
      if (curyr /= nbyr .and. iida /= idal) then
         read (112+k,*) idapa(k), iypa(k), flodaya(k),&
            &seddaya(k), orgndaya(k), orgpdaya(k), no3daya(k),&
            &minpdaya(k)
      endif
   else
      varoute(2,ihout) = 0.0
      varoute(3,ihout) = 0.0
      varoute(4,ihout) = 0.0
      varoute(5,ihout) = 0.0
      varoute(6,ihout) = 0.0
      varoute(7,ihout) = 0.0
      varoute(14,ihout) = 0.0
      varoute(15,ihout) = 0.0
      varoute(11,ihout) = 0.0
      varoute(12,ihout) = 0.0
      varoute(13,ihout) = 0.0
      varoute(14,ihout) = 0.0
      varoute(15,ihout) = 0.0
      varoute(16,ihout) = 0.0
      varoute(17,ihout) = 0.0
      varoute(18,ihout) = 0.0
      varoute(19,ihout) = 0.0
      varoute(20,ihout) = 0.0
      varoute(21,ihout) = 0.0
      varoute(22,ihout) = 0.0
   endif

   if (ievent > 0) then
      fn = Dfloat(nstep)
      do ii = 1, nstep
         hhvaroute(2,ihout,ii) = flodaya(k) / fn
         hhvaroute(3,ihout,ii) = seddaya(k) / fn
         hhvaroute(4,ihout,ii) = orgndaya(k) / fn
         hhvaroute(5,ihout,ii) = orgpdaya(k) / fn
         hhvaroute(6,ihout,ii) = no3daya(k) / fn
         hhvaroute(7,ihout,ii) = minpdaya(k) / fn

      end do
   end if

   return
end
