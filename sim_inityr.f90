!> @file sim_inityr.f90
!> file containing the subroutine sim_inityr
!> @author
!> modified by Javier Burguete

!> this subroutine initializes variables at the beginning of the year
subroutine sim_inityr

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    curyr       |none          |current year in simulation (sequence)
!!    iskip       |none          |flag for calculations performed only for
!!                               |the first year of simulation
!!    nyskip      |none          |number of years to not print output
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    anano3(:)   |kg N/ha       |total amount of NO3 applied during the year
!!                               |in auto-fertilization
!!    bio_aams(:) |metric tons   |average annual biomass in the HRU
!!    bio_yrms(:) |metric tons/ha|annual biomass (dry weight) in the HRU
!!    hrumono(:,:)|varies        |HRU monthly output data array
!!    lai_aamx(:) |none          |maximum leaf area index for the entire period
!!                               |of simulation in the HRU
!!    lai_yrmx(:) |none          |maximum leaf area index for the year in the
!!                               |HRU
!!    nafert(:)   |none          |sequence number of auto-fert application
!!                               |within the year
!!    nair(:)     |none          |sequence number of auto-irrigation application
!!                               |within the year
!!    ncf(:)      |none          |sequence number of continuous fertilization
!!                               |operation within the year
!!    nfert(:)    |none          |sequence number of fertilizer application
!!                               |within the year
!!    ngr(:)      |none          |sequence number of grazing operation
!!                               |within the year
!!    nirr(:)     |none          |sequence number of irrigation application
!!                               |within the year
!!    nrelease(:) |none          |sequence number of impound/release operation
!!                               |within the year
!!    nsweep(:)   |none          |sequence number of street sweeping operation
!!                               |within the year
!!    phubase(:)  |heat units    |base zero total heat units (used when no land
!!                               |cover is growing)
!!    rchmono(:,:)|varies        |reach monthly output array
!!    resoutm(:,:)|varies        |reservoir monthly output array
!!    submono(:,:)|varies        |subbasin monthly output array
!!    wtrmon(:,:) |varies        |HRU monthly output data array for impoundments
!!    yldaa(:)    |metric tons   |average annual yield in the HRU
!!    yldanu(:)   |metric tons/ha|annual yield (dry weight) in the HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   !! initialize variables/arrays at beginning of every year
   !ncpest = 1 ! not used
   !nfert = 1 ! not used
   !nsweep= 1 ! not used
   !nirr = 1 ! not used
   !ngr = 1 ! not used
   !ncf = 1 ! not used
   yldanu = 0.
   bio_yrms = 0.
   lai_yrmx = 0.
   anano3 = 0.
   !tauton = 0. ! not used
   !tautop = 0. ! not used
   !tcfrtn = 0. ! not used
   !tcfrtp = 0. ! not used
   !tfertn = 0. ! not used
   !tfertp = 0. ! not used
   !tgrazn = 0. ! not used
   !tgrazp = 0. ! not used


   if (iskip == 0) then
      !! initialize arrays which store printout data
      if (curyr > nyskip) then
         hrumono = 0.
         wtrmon = 0.
         submono = 0.
         rchmono = 0.
         resoutm = 0.
         yldaa = 0.
         bio_aams = 0.
         lai_aamx = 0.
         iskip = 1
      endif
   endif

   return
end
