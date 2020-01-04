!> @file command.f90
!> file containing the subroutine command
!> @author
!> modified by Javier Burguete

!> for every day of simulation, this subroutine steps through the command
!> lines in the watershed configuration (.fig) file. Depending on the
!> command code on the .fig file line, a command loop is accessed
!> @param[in] i current day in simulation--loop counter (julian date)
subroutine command(i)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |julian date   |current day in simulation--loop counter
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
!!                               |4: water source type
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
!!                               |0,1,5,7,8,10,11: not used
!!                               |2,3: subbasin number
!!                               |4: destination number. Reach or
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
!!                               |   (2=min volume or flow left)
!!                               |   (3=exact amount transferred)
!!    mhyd        |none          |maximum number of hydrographs
!!    rnum1s(:)   |none          |For ICODES equal to
!!                               |0,1,3,5,9: not used
!!                               |2: fraction of overland flow
!!                               |4: amount of water transferred (as
!!                               |   defined by INUM4S)
!!                               |7,8,10,11: drainage area in square kilometers
!!                               |   associated with the record file
!!                               |12: rearation coefficient
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihout       |none          |variable to hold value for ihouts(:)
!!    inum1       |none          |variable to hold value for inum1s(:)
!!    inum2       |none          |variable to hold value for inum2s(:)
!!    inum3       |none          |variable to hold value for inum3s(:)
!!    inum4       |none          |variable to hold value for inum4s(:)
!!    rnum1       |none          |variable to hold value for rnum1s(:)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    icode    |none          |variable to hold value for icodes(:)
!!    ii       |none          |counter
!!    j        |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: subbasin, print_hyd, route, bmp_det_pond, bmp_wet_pond, sumhyd,
!!          routres, transfer, addh, rechour, recmon, recyear, save, recday,
!!          reccnst, structure, apex_day, saveconc, routeunit, routels

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i
   integer :: icode, ii, inum1, j

   j = 0

   do ii = 1, mhyd_bsn
      j = j + 1
      icode = icodes(j)
      ihout = ihouts(j)
      inum1 = inum1s(j)
      inum2 = inum2s(j)
      inum3 = inum3s(j)
      rnum1 = rnum1s(j)
      inum4 = inum4s(j)
      inum5 = inum5s(j)
      inum6 = inum6s(j)
      inum7 = inum7s(j)
      inum8 = inum8s(j)

      select case (icode)
       case (0)
         return
       case (1)
         call subbasin(i, inum1)
         call print_hyd(i)
       case (2)
         call route(i, inum1)
         if (dtp_onoff(inum1)==1) call bmp_det_pond(inum1) !route detention pond J.Jeong feb 2010
         if (wtp_onoff(inum1)==1) call bmp_wet_pond(inum1) !route wetention pond J.Jeong june 2010
         call sumhyd
         call print_hyd(i)
       case (3)
         call routres(inum1)
         call sumhyd
       case (4)
         call transfer(inum1)
       case (5)
         call addh(inum1)
         call sumhyd
         call print_hyd(i)
       case (6)
         call rechour(inum1)
         call sumhyd
         call print_hyd(i)
       case (7)
         call recmon(inum1)
         call sumhyd
       case (8)
         call recyear(inum1)
         call sumhyd
       case (9)
         call save(inum1)
       case (10)
         call recday(inum1)
         call sumhyd
       case (11)
         call reccnst(inum1)
         call sumhyd
       case (12)
         call structure(inum1)
       case (13)
         call apex_day(i, inum1)
       case (14)
         call saveconc(inum1)
       case (17)
         call routeunit(inum1)
         call sumhyd
       case (18)
         inum8 = 1
         call routels(inum1, inum1)
         call sumhyd
      end select

   end do

   return
end
