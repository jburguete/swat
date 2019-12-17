!> @file readchm.f90
!> file containing the subroutine readchm
!> @author
!> modified by Javier Burguete

!> This subroutine reads data from the HRU/subbasin soil chemical input
!> file (.chm). This file contains initial amounts of pesticides/nutrients
!> in the first soil layer. (Specifics about the first soil layer are given
!> in the .sol file.) All data in the .chm file is optional input.
subroutine readchm

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru        |none          |HRU number
!!    mlyr        |none          |maximum number of soil layers
!!    mpst        |none          |maximum number of pesticides used in
!!                               |watershed
!!    nope(:)     |none          |sequence number of pesticide in NPNO(:)
!!    npmx        |none          |number of different pesticides used in
!!                               |the simulation
!!    npno(:)     |none          |array of unique pesticides used in
!!                               |watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hrupest(:)    |none          |pesticide use flag:
!!                                 | 0: no pesticides used in HRU
!!                                 | 1: pesticides used in HRU
!!    nope(:)       |none          |sequence number of pesticide in NPNO(:)
!!    npmx          |none          |number of different pesticides used in
!!                                 |the simulation
!!    npno(:)       |none          |array of unique pesticides used in
!!                                 |watershed
!!    plt_pst(:,:)  |kg/ha         |pesticide on plant foliage
!!    sol_pst(:,:,1)|mg/kg         |pesticide concentration in soil
!!    pst_enr(:,:)  |none          |pesticide enrichment ratio
!!    sol_no3(:,:)  |mg N/kg       |concentration of nitrate in soil layer
!!    sol_orgn(1,:) |mg N/kg soil  |organic N concentration in top soil layer
!!    sol_orgp(1,:) |mg P/kg soil  |organic P concentration in top soil layer
!!    sol_solp(1,:) |mg P/kg soil  |soluble P concentration in top soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag
!!    j           |none          |counter
!!    k           |none          |counter
!!    l           |none          |HRU number
!!    newpest     |none          |new pesticide flag
!!    pltpst      |kg/ha         |pesticide on plant foliage
!!    pstenr      |none          |pesticide enrichment ratio
!!    pstnum      |none          |pesticide number
!!    solpst      |mg/kg         |pesticide concentration in soil
!!    titldum     |NA            |title line for .chm file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   character (len=80) :: titldum
   real*8 :: pltpst, solpst, pstenr
   integer :: j, eof, k, l, newpest, pstnum

   eof = 0
   l = ihru


   do
      read (106,5000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,5000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,5000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,5100,iostat=eof) (sol_no3(j,l), j = 1, mlyr)
      if (eof < 0) exit
      read (106,5100,iostat=eof) (sol_orgn(j,l), j = 1, mlyr)
      if (eof < 0) exit
      read (106,5100,iostat=eof) (sol_solp(j,l), j = 1, mlyr)
      if (eof < 0) exit
      read (106,5100,iostat=eof) (sol_orgp(j,l), j = 1, mlyr)
      if (eof < 0) exit
      read (106,5100,iostat=eof) (pperco_sub(j,l), j = 1, mlyr)
      if (eof < 0) exit
      read (106,5000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,5000,iostat=eof) titldum
      if (eof < 0) exit
      read (106,5000,iostat=eof) titldum
      if (eof < 0) exit
!! end do

      do j = 1, mpst
         pstnum = 0
         pltpst = 0.
         solpst = 0.
         pstenr = 0.
         read (106,*,iostat=eof) pstnum, pltpst, solpst, pstenr
         if (pstnum > 0) then
            hrupest(l) = 1
            newpest = 0
            do k = 1, npmx
               if (pstnum == npno(k)) then
                  newpest = 1
                  exit
               endif
            end do
            if (newpest == 0) then
               npno(npmx) = pstnum
               nope(pstnum) = npmx
               npmx = npmx + 1
            end if

            k = nope(pstnum)
            plt_pst(k,l) = pltpst
            sol_pst(k,l,1) = solpst
            pst_enr(k,l) = pstenr
         end if

         if (eof < 0) exit
      end do
      exit


   end do

   close (106)

   do j = 1, mlyr
      if (pperco_sub(j,l) <= 1.e-6) pperco_sub(j,l) = pperco_bsn
   end do

   return
5000 format (a)
5100 format (27x,10f12.2)
end
