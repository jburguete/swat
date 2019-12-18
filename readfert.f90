!> @file readfert.f90
!> file containing the subroutine readfert
!> @author
!> modified by Javier Burguete

!> this subroutine reads input parameters from the fertilizer/manure
!> (i.e. nutrient) database (fert.dat)
subroutine readfert

!!     ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!     name        |units           |definition
!!     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!     mfdb        |none            |maximum number of fertilizers in
!!                                  |database
!!     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!     ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!     name        |units           |definition
!!     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!     bactkddb(:) |none            |bacteria partition coefficient:
!!                                  |1: all bacteria in solution
!!                                  |0: all bacteria sorbed to soil particles
!!     bactlpdb(:) |# cfu/g manure  |concentration of less persistent
!!                                  |bacteria in manure(fertilizer)
!!     bactpdb(:)  |# cfu/g manure  |concentration of persistent bacteria
!!                                  |in manure(fertilizer)
!!     fertnm(:)   |NA              |name of fertilizer
!!     fminn(:)    |kg minN/kg fert |fraction of mineral N (NO3 + NH3)
!!     fminp(:)    |kg minP/kg fert |fraction of mineral P
!!     fnh3n(:)    |kg NH3-N/kg minN|fraction of NH3-N in mineral N
!!     forgn(:)    |kg orgN/kg fert |fraction of organic N
!!     forgp(:)    |kg orgP/kg fert |fraction of organic P
!!     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!     ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!     name        |units           |definition
!!     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!     bctpdb
!!     bctlpdb
!!     bctkddb
!!     eof         |none            |end of file flag
!!     ffminn
!!     ffminp
!!     ffnh3n
!!     fnm
!!     fforgn
!!     fforgp
!!     it          |none            |counter which represents the array
!!                                  |storage number of the pesticide data
!!                                  |the array storage number is used by the
!!                                  |model to access data for a specific
!!                                  |fertilizer
!!     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!     ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8 :: bctpdb, bctlpdb, bctkddb, ffminn, ffminp, ffnh3n, fforgn, fforgp
   character (len=8) :: fnm
   integer :: eof, it

   eof = 0

   do
      !! initialize local variables
      bctkddb = 0.
      bctlpdb = 0.
      bctpdb = 0.
      ffminn = 0.
      ffminp = 0.
      ffnh3n = 0.
      fforgn = 0.
      fforgp = 0.
      fnm = ""
      it = 0

      read (107,5000,iostat=eof) it, fnm, ffminn, ffminp, fforgn,&
         &fforgp, ffnh3n, bctpdb, bctlpdb, bctkddb

      if (eof < 0 .or. it == 0) exit

      fertnm(it) = fnm
      fminn(it) = ffminn
      fminp(it) = ffminp
      forgn(it) = fforgn
      forgp(it) = fforgp
      fnh3n(it) = ffnh3n
      bactpdb(it) = bctpdb
      bactlpdb(it) = bctlpdb
      bactkddb(it) = bctkddb
   end do

   close (107)
   return
!! format changed for prem - last two bacteria inputs ok for all
5000 format (i4,1x,a8,6f8.3,2f11.0)
end
