!> @file readtill.f90
!> file containing the subroutine readtill
!> @author
!> modified by Javier Burguete

!> this subroutine reads input data from tillage database (till.dat)
subroutine readtill

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mtil        |none          |maximum number of tillage operations in
!!                               |till.dat
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    deptil(:)   |mm            |depth of mixing caused by operation
!!    effmix(:)   |none          |mixing efficiency of operation
!! drainmod tile equations   06/2006
!!    ranrns      |mm            |random roughness of a given tillage operation
!! drainmod tile equations   06/2006
!!    tillnm(:)   |NA            |8-character name for the tillage
!!                               |operation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dtil
!!    emix
!!    eof         |none          |end of file flag
!!    it          |none          |counter which represents the array
!!                               |storage number of the tillage data
!!                               |the array storage number is used by the
!!                               |model to access data for a specific
!!                               |tillage operation
!!    rrns
!!    tlnm
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

!! drainmod tile equations  - addition random roughness 06/2006
   real*8 :: dtil, emix, rrns
!! drainmod tile equations  - addition random roughness 06/2006
   character (len=8) :: tlnm
   integer :: eof, it

   eof = 0

   do
      dtil = 0.
      emix = 0.
      it = 0
!! drainmod tile equations   06/2006
      rrns = 0.
!! drainmod tile equations   06/2006
      tlnm = ""

!! drainmod tile equations   06/2006
      read (105,5000,iostat=eof) it, tlnm, emix, dtil, rrns
!! drainmod tile equations   06/2006

      if (eof < 0) exit

      if (it == 0) exit

      tillnm(it) = tlnm
      effmix(it) = emix
      deptil(it) = dtil
!! drainmod tile equations   06/2006
      ranrns(it) = rrns
!! drainmod tile equations   06/2006

   end do

   close (105)
   return
!! drainmod tile equations  - addition random roughness 06/2006
5000 format (i4,3x,a8,8x,f8.3,8x,f8.3,8x,f8.3)
!! drainmod tile equations  - addition random roughness 06/2006
end
