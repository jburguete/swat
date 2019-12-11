subroutine hruallo

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!   This subroutine calculates the number of management operation types, etc.
!!   used in the simulation. These values are used to allocate array sizes for
!!   processes occurring in the HRU.

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mapp        |none        |max number of applications
!!    mcr         |none        |max number of crops grown per year
!!    mcut        |none        |max number of cuttings per year
!!    mgr         |none        |max number of grazings per year
!!    mlyr        |none        |max number of soil layers
!!    mnr         |none        |max number of years of rotation
!!    pstflg(:)   |none        |flag for types of pesticide used in watershed
!!                             |array location is pesticide ID number
!!                             |0: pesticide not used
!!                             |1: pesticide used
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chmfile     |NA          |HRU soil chemical data file name (.chm)
!!    depth(:)    |mm          |depth to bottom of soil layer
!!    eof         |none        |end of file flag (=-1 if eof, else =0)
!!    hru         |none        |number of HRUs in subbasin
!!    hrufile     |NA          |name of HRU general data file name (.hru)
!!    j           |none        |counter
!!    k           |none        |counter
!!    lyrtot      |none        |total number of layers in profile
!!    mgt_op      |none        |manangement operation code
!!    mgt1i       |none        |sixth parameter in mgt file operation line
!!    mgtfile     |NA          |HRU management data file name (.mgt)
!!    pstnum      |none        |pesticide ID number from database file
!!    solfile     |NA          |HRU soil data file name (.sol)
!!    titldum     |NA          |input lines in .sub that are not processed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max
!!    SWAT: caps

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   character (len=13) :: hrufile, mgtfile, solfile, chmfile
   character (len=80) ::  titldum
   integer :: eof, j, k, lyrtot
   integer :: mgt1i, pstnum
   integer :: mgt_op, mcri
   integer :: iopera_sub, ilnds
   real*8 :: depth(25)

   do j= mhru1, mhru
      mgtfile = ""
      solfile = ""
      chmfile = ""
      read (25,5300,iostat=eof)hrufile, mgtfile, solfile, chmfile, ilnds
      if (eof < 0) return
      if (ilnds > 0) then
         ils_nofig = 1
      end if
      call caps(mgtfile)
      call caps(solfile)
      call caps(chmfile)
      open (9,file=solfile,recl=350)
      !! calculate # of soil layers in profile
      depth = 0.
      lyrtot = 0
      read (9,6000) titldum
      read (9,6000) titldum
      read (9,6000) titldum
      read (9,6000) titldum
      read (9,6000) titldum
      read (9,6000) titldum
      read (9,6000) titldum
      read (9,6100) (depth(k), k = 1, 25)
      do k = 1, 25
         if (depth(k) <= 0.001) lyrtot = k - 1
         if (depth(k) <= 0.001) exit
      end do
      mlyr = Max(mlyr,lyrtot)
      close (9)
      open (10,file=mgtfile)

!!  calculate max number of operations per hru
      iopera_sub = 1
      mcri = 0
      do k = 1, 30
         read (10,6000) titldum
      end do

      do k = 1, 1000
         read (10,6300,iostat=eof) mgt_op, mgt1i
         if (eof < 0) exit
         if (mgt_op == 1) then
            mcri = mcri + 1
         end if
         if (mgt_op == 4 .and. mgt1i > 0) pstflg(mgt1i) = 1
         iopera_sub = iopera_sub + 1
      end do
      iopera = Max(iopera,iopera_sub)
      mcr = Max(mcr,mcri)

      close (10)            !!   nubz test

      open (11,file=chmfile)
      eof = 0
      do
         do k = 1, 11
            read (11,6000,iostat=eof) titldum
            if (eof < 0) exit
         end do
         if (eof < 0) exit
         do
            pstnum = 0
            read (11,*,iostat=eof) pstnum
            if (eof < 0) exit
            if (pstnum > 0) pstflg(pstnum) = 1
         end do
         if (eof < 0) exit
      end do
      close (11)
   end do    ! hru loop

   return
5300 format (4a13,52x,i6)
6000 format (a80)
6100 format (27x,25f12.2)
6300 format (16x,i2,1x,i4)
end
