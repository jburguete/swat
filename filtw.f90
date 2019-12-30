!> @file filtw.f90
!> file containing the subroutine filtw
!> @author
!> modified by Javier Burguete

!> this subroutine calculates the reduction of pollutants in surface runoff
!> due to an edge of field filter or buffer strip
!> @param[in] j HRU number (none)
subroutine filtw(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    bactrolp    |# colonies/ha |less persistent bacteria transported to main
!!                               |channel with surface runoff
!!    bactrop     |# colonies/ha |persistent bacteria transported to main
!!                               |channel with surface runoff
!!    bactsedlp   |# colonies/ha |less persistent bacteria transported with
!!                               |sediment in surface runoff
!!    bactsedp    |# colonies/ha |persistent bacteria transported with
!!                               |sediment in surface runoff
!!    curyr       |none          |current year of simulation
!!    fsred(:)    |none          |reduction in bacteria loading from filter
!!                               |strip
!!    hru_dafr(:) |none          |fraction of watershed area in HRU
!!    hrupest(:)  |none          |pesticide use flag:
!!                               | 0: no pesticides used in HRU
!!                               | 1: pesticides used in HRU
!!    npmx        |none          |number of different pesticides used in
!!                               |the simulation
!!    nyskip      |none          |number of years to skip output summarization
!!                               |and printing
!!    pst_sed(:,:)|kg/ha         |pesticide loading from HRU sorbed onto
!!                               |sediment
!!    pst_surq(:,:)|kg/ha        |amount of pesticide type lost in surface
!!                               |runoff on current day in HRU
!!    sbactrolp   |# colonies/ha |average annual number of less persistent
!!                               |bacteria transported to main channel
!!                               |with surface runoff in solution
!!    sbactrop    |# colonies/ha |average annual number of persistent bacteria
!!                               |transported to main channel with surface
!!                               |runoff in solution
!!    sbactsedlp  |# colonies/ha |average annual number of less persistent
!!                               |bacteria transported with sediment in
!!                               |surface runoff
!!    sbactsedp   |# colonies/ha |average annual number of persistent bacteria
!!                               |transported with sediment in surface runoff
!!    sedminpa(:) |kg P/ha       |amount of active mineral phosphorus sorbed to
!!                               |sediment in surface runoff in HRU for day
!!    sedminps(:) |kg P/ha       |amount of stable mineral phosphorus sorbed to
!!                               |sediment in surface runoff in HRU for day
!!    sedorgn(:)  |kg N/ha       |amount of organic nitrogen in surface runoff
!!                               |in HRU for the day
!!    sedorgp(:)  |kg P/ha       |amount of organic phosphorus in surface runoff
!!                               |in HRU for the day
!!    sedyld(:)   |metric tons   |daily soil loss caused by water erosion
!!    surqno3(:)  |kg N/ha       |amount of NO3-N in surface runoff in HRU for
!!                               |the day
!!    surqsolp(:) |kg P/ha       |amount of soluble phosphorus in surface runoff
!!                               |in HRU for the day
!!    trapeff(:)  |none          |filter strip trapping efficiency (used for
!!                               |everything but bacteria)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bactrolp    |# colonies/ha |less persistent bacteria transported to main
!!                               |channel with surface runoff
!!    bactrop     |# colonies/ha |persistent bacteria transported to main
!!                               |channel with surface runoff
!!    bactsedlp   |# colonies/ha |less persistent bacteria transported with
!!                               |sediment in surface runoff
!!    bactsedp    |# colonies/ha |persistent bacteria transported with
!!                               |sediment in surface runoff
!!    pst_sed(:,:)|kg/ha         |pesticide loading from HRU sorbed onto
!!                               |sediment
!!    pst_surq(:,:)|kg/ha        |amount of pesticide type lost in surface
!!                               |runoff on current day in HRU
!!    sbactrolp   |# colonies/ha |average annual number of less persistent
!!                               |bacteria transported to main channel
!!                               |with surface runoff in solution
!!    sbactrop    |# colonies/ha |average annual number of persistent bacteria
!!                               |transported to main channel with surface
!!                               |runoff in solution
!!    sbactsedlp  |# colonies/ha |average annual number of less persistent
!!                               |bacteria transported with sediment in
!!                               |surface runoff
!!    sbactsedp   |# colonies/ha |average annual number of persistent bacteria
!!                               |transported with sediment in surface runoff
!!    sedminpa(:) |kg P/ha       |amount of active mineral phosphorus sorbed to
!!                               |sediment in surface runoff in HRU for day
!!    sedminps(:) |kg P/ha       |amount of stable mineral phosphorus sorbed to
!!                               |sediment in surface runoff in HRU for day
!!    sedorgn(:)  |kg N/ha       |amount of organic nitrogen in surface runoff
!!                               |in HRU for the day
!!    sedorgp(:)  |kg P/ha       |amount of organic phosphorus in surface runoff
!!                               |in HRU for the day
!!    sedyld(:)   |metric tons   |daily soil loss caused by water erosion
!!    surqno3(:)  |kg N/ha       |amount of NO3-N in surface runoff in HRU for
!!                               |the day
!!    surqsolp(:) |kg P/ha       |amount of soluble phosphorus in surface runoff
!!                               |in HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k           |none          |counter
!!    xx          |none          |reduction factor
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: xx
   integer :: k

!! compute filter strip reduction
   xx = fsred(j)
   bactrop = bactrop * xx
   bactrolp = bactrolp * xx
   bactsedp = bactsedp * xx
   bactsedlp = bactsedlp * xx
   xx = 1. - trapeff(j)
   sedorgn(j) = sedorgn(j) * xx
   surqno3(j) = surqno3(j) * xx
   sedorgp(j) = sedorgp(j) * xx
   sedminpa(j) = sedminpa(j) * xx
   sedminps(j) = sedminps(j) * xx
   surqsolp(j) = surqsolp(j) * xx
   sedyld(j) = sedyld(j) * xx
   sanyld(j) = sanyld(j) * xx
   silyld(j) = silyld(j) * xx
   clayld(j) = clayld(j) * xx
   sagyld(j) = sagyld(j) * xx
   lagyld(j) = lagyld(j) * xx
   if (hrupest(j) == 1) then
      do k = 1, npmx
         pst_surq(k,j) = pst_surq(k,j) * xx
         pst_sed(k,j) = pst_sed(k,j) * xx
      end do
   end if

!! summary calculations
   if (curyr > nyskip) then
      xx = hru_dafr(j)
      sbactrop = sbactrop + bactrop * xx
      sbactrolp = sbactrolp + bactrolp * xx
      sbactsedp = sbactsedp + bactsedp * xx
      sbactsedlp = sbactsedlp + bactsedlp * xx
   end if

   return
end
