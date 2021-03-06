!> @file apply.f90
!> file containing the subroutine apply
!> @author
!> modified by Javier Burguete

!> this subroutine applies pesticide
!> @param[in] j HRU number
subroutine apply(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j            |none             |HRU number
!!    ap_ef(:)     |none             |application efficiency (0-1)
!!    curyr        |none             |current year of simulation
!!    drift(:)     |kg               |amount of pesticide drifting onto main
!!                                   |channel in subbasin
!!    driftco(:)   |none             |coefficient for pesticide drift directly
!!                                   |onto stream
!!    hru_dafr(:)  |km**2/km**2      |fraction of watershed area in HRU
!!    hru_km(:)    |km**2            |area of HRU in square kilometers
!!    ipest        |none             |pesticide identification number from
!!                                   |pest.dat
!!    irtpest      |none             |the sequence number of the pesticide type
!!                                   |in NPNO(:) which is to be routed through
!!                                   |the watershed
!!    laiday(:)    |none             |leaf area index
!!    nope(:)      |none             |sequence number of pesticide in NPNO(:)
!!    nro(:)       |none             |sequence number of year in rotation
!!    nyskip       |none             |number of years to skip output
!!                                   |summarization/printing
!!    plt_pst(:,:) |kg/ha            |pesticide on plant foliage
!!    pst_dep      |kg/ha            |depth of pesticide in soil
!!    pst_kg       |kg/ha            |amount of pesticide applied to HRU
!!    sol_pst(:,:,1)|kg/ha           |pesticide in first layer of soil
!!    wshd_pstap(:)|kg/ha            |total amount of pesticide type applied in
!!                                   |watershed during simulation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    drift(:)    |kg            |amount of pesticide drifting onto main
!!                               |channel in subbasin
!!    plt_pst(:,:)|kg/ha         |pesticide on plant foliage
!!    sol_pst(:,:,1)|kg/ha       |pesticide in first layer of soil
!!    wshd_pstap(:)|kg/ha        |total amount of pesticide type applied in
!!                               |watershed during simulation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    gc          |none          |fraction of ground covered by plant foliage
!!    k           |none          |sequence number of pesticide in NPNO(:)
!!    kk          |none          |pesticide identification number from
!!                               |pest.dat
!!    xx          |kg/ha         |amount of pesticide applied to HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: Erfc

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: gc, xx
   integer :: kk, k, nly

!! initialize local variables

   kk = ipest

   k = nope(kk)

   xx = pst_kg

!! calculate amount of pesticide drifting onto main channel in subbasin
   xx = xx * ap_ef(kk)

!   added for pesticide incorporation 3/31/08 gsm
   if (pst_dep > 1.e-6) then
      do nly = 1,sol_nly(j)
         if (nly == 1) then
            if (pst_dep < sol_z(nly,j)) then
               sol_pst(k,j,1) = sol_pst(k,j,1) + xx
               exit
            endif
         else
            if (pst_dep > sol_z((nly-1),j) .and.&
               &pst_dep < sol_z(nly,j)) then
               sol_pst(k,j,nly) = sol_pst(k,j,nly) + xx
               exit
            endif
         endif
      enddo
   else
!   added above for pesticide incorporation 3/31/08 gsm

!! calculate ground cover
      gc = (1.99532 - Erfc(1.333 * laiday(j) - 2.)) / 2.1
      if (gc < 0.) gc = 0.

!! update pesticide levels on ground and foliage
      plt_pst(k,j) = plt_pst(k,j) + gc * xx
      sol_pst(k,j,1) = sol_pst(k,j,1) + (1. - gc) * xx
!! added endif for pesticide incorporation 3/31/08 gsm
   endif

!! summary calculations
   if (curyr > nyskip) then
      wshd_pstap(k) = wshd_pstap(k) + pst_kg *&
         &ap_ef(kk) * hru_dafr(j)
   end if

   return
end
