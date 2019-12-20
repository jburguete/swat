!> @file rteinit.f90
!> file containing the subroutine rteinit
!> @author
!> modified by Javier Burguete

!> This subroutine reads in the areas associated with files processed with
!> the recday, recepic, recmon and recyear commands, calculates subbasin
!> areas, calculates reach and hydrograph node drainage areas.
subroutine rteinit

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    da_km       |km2           |area of the watershed in square kilometers
!!    icodes(:)   |none          |routing command code:
!!                               |0 = finish       9 = save
!!                               |1 = subbasin    10 = recday
!!                               |2 = route       11 = reccnst
!!                               |3 = routres     12 = structure
!!                               |4 = transfer    13 =
!!                               |5 = add         14 = saveconc
!!                               |6 = rechour     15 =
!!                               |7 = recmon
!!                               |8 = recyear
!!    ihouts(:)   |none          |For ICODES equal to
!!                               |0: not used
!!                               |1,2,3,5,7,8,10,11: hydrograph storage
!!                               |                     location number
!!                               |4: departure type
!!                               |   (1=reach)
!!                               |   (2=reservoir)
!!                               |9: hydrograph storage location of data to
!!                               |   be printed to event file
!!                               |14:hydrograph storage location of data to
!!                               |   be printed to saveconc file
!!    inum1s(:)   |none          |For ICODES equal to
!!                               |0: not used
!!                               |1: subbasin number
!!                               |2: reach number
!!                               |3: reservoir number
!!                               |4: reach or res # flow is diverted from
!!                               |5: hydrograph storage location of 1st
!!                               |   dataset to be added
!!                               |7,8,9,10,11,14: file number
!!    inum2s(:)   |none          |For ICODES equal to
!!                               |0,1,7,8,10,11: not used
!!                               |2,3: inflow hydrograph storage location
!!                               |4: destination type
!!                               |   (1=reach)
!!                               |   (2=reservoir)
!!                               |5: hydrograph storage location of 2nd
!!                               |   dataset to be added
!!                               |9,14:print frequency
!!                               |   (0=daily)
!!                               |   (1=hourly)
!!    inum3s(:)   |none          |For ICODES equal to
!!                               |0,1,2,3,5,7,8,10,11: not used
!!                               |4: destination number. Subbasin, reach or
!!                               |   reservoir receiving water
!!                               |9: print format
!!                               |   (0=normal, fixed format)
!!                               |   (1=txt format for AV interface,recday)
!!    inum4s(:)   |none          |For ICODES equal to
!!                               |0,2,3,5,7,8,9,10,11: not used
!!                               |1: GIS code printed to output file
!!                               |   (optional)
!!                               |4: rule code governing transfer of water
!!                               |   (1=fraction transferred out)
!!                               |   (2=min flow left)
!!                               |   (3=exact amount transferred)
!!    mhyd        |none          |maximum number of hydrograph nodes
!!    rnum1s(:)   |none          |For ICODES equal to
!!                               |0,1,3,5,9: not used
!!                               |2: Fraction of flow in channel
!!                               |4: amount of water transferred (as
!!                               |   defined by INUM4S)
!!                               |7,8,10,11: drainage area in square kilometers
!!                               |   associated with the record file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hyd_dakm(:) |km**2         |total drainage area of hydrograph in square
!!                               |kilometers
!!    rch_dakm(:) |km**2         |total drainage area contributing to flow at
!!                               |the outlet (pour point) of the reach in
!!                               |square kilometers
!!    sub_fr(:)   |km2/km2       |fraction of total watershed area contained
!!                               |in subbasin
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    icd         |none          |variable to hold value for icodes(:)
!!    iht         |none          |variable to hold value for ihouts(:)
!!    inm1        |none          |variable to hold value for inum1s(:)
!!    inm2        |none          |variable to hold value for inum2s(:)
!!    inm3        |none          |variable to hold value for inum3s(:)
!!    inm4        |none          |variable to hold value for inum4s(:)
!!    isb
!!    j
!!    j1
!!    jseq
!!    k
!!    rnm1        |none          |variable to hold value for rnum1s(:)
!!    xx          |none          |auxiliar variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   real*8 :: rnm1, xx
   integer :: icd, iht, inm1, inm2, inm3, inm4, isb, j, j1, jseq, k

!!    calculate watershed area in hectares
   da_ha = da_km * 100.
!!    calculate subbasin fraction of watershed area
   do k = 1, subtot
      sub_fr(k) = sub_km(k) / da_km
   end do
!!    calculate HRU fraction of watershed area
   do k = 1, nhru
      hru_dafr(k) = hru_km(k) / da_km
      isb = hru_sub(k)
      if (idplt(k) > 0) then
         if (bio_e(idplt(k)) > 1.e-6) then
            subfr_nowtr(isb) = subfr_nowtr(isb) + hru_fr(k)
         end if
      end if
   end do

!!    read in areas associated with .fig record files
   do k = 1, mhyd

      icd = icodes(k)
      rnm1 = rnum1s(k)
      iht = ihouts(k)

      select case (icd)
       case (0)   !! FINISH command
         exit
       !! RECHOUR, RECMON, RECYEAR, RECDAY or RECCNST command
       case (6, 7, 8, 10, 11)
         hyd_dakm(iht) = rnm1
      end select
   end do

!!    calculate subbasin areas
!!    calculate reach and hydrograph drainage areas
   do k = 1, mhyd

      icd = icodes(k)
      iht = ihouts(k)
      inm1 = inum1s(k)
      inm2 = inum2s(k)
      inm3 = inum3s(k)
      rnm1 = rnum1s(k)
      inm4 = inum4s(k)

      select case (icd)
       case (0)   !! FINISH command
         exit

       case (1)   !! SUBBASIN command
         hyd_dakm(iht) = da_km * sub_fr(inm1)
         rch_dakm(inm1) = hyd_dakm(iht)
         if (iroutunit == 1) then
            j1 = hru1(inm1) - 1
            do j = 1, hrutot(inm1)
               jseq = j1 + j
               hru_fr(jseq) = 0.
               hru_dafr(jseq) = 0.
            end do
         end if

       case (2)   !! ROUTE command
         hyd_dakm(iht) = hyd_dakm(inm2)
         rch_dakm(inm1) = hyd_dakm(inm2)

       case (3, 18)   !! ROUTRES or ROUTEOVER command
         hyd_dakm(iht) = hyd_dakm(inm2)

       case (5)   !! ADD command
         hyd_dakm(iht) = hyd_dakm(inm1) + hyd_dakm(inm2)

       case (17)   !! ROUTEUNIT command
         hyd_dakm(iht) = rnm1
         j1 = hru1(inm2) - 1
         do j = 1, hrutot(inm2)
            jseq = j1 + j
            xx = hru_rufr(inm1,j) * daru_km(inm2,inm1)
            hru_fr(jseq) = hru_fr(jseq) + xx / sub_km(inm2)
            hru_dafr(jseq) = hru_dafr(jseq) + xx / da_km
         end do

      end select
   end do

   return
end
