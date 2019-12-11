subroutine readsdr

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads data from the HRU/subbasin management input file
!!    (.mgt). This file contains data related to management practices used in
!!    the HRU/subbasin.

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name      |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!                                  |urban.dat
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name      |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!              |daily
!!    drain_co(:)|mm/day        |drainage coefficient
!!    latksatf(:)|none          |multiplication factor to determine conk(j1,j) from sol_k(j1,j) for HRU
!!    pc(:)     |mm/hr          |pump capacity (default pump capacity = 1.042mm/hr or 25mm/day)
!!    re(:)     |mm             |effective radius of drains
!!    sdrain(:) |mm             |distance between two drain tubes or tiles
!!    sstmaxd(:)|mm             |static maximum depressional storage; read from .sdr
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    eof         |none          |end of file flag (=-1 if eof, else = 0)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mgt3i       |none          |third management parameter out of .mgt
!!                               |file (definition changes depending on
!!                               |mgt_op)
!!    titldum     |NA            |title line from input dataset
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: Jdt

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   character (len=80) :: titldum
   integer :: eof

   do
      read (112,5000,iostat=eof) titldum
      if (eof < 0) exit
!!      read scheduled operations
      read (112,*,iostat=eof) re(ihru)
      if (eof < 0) exit
      read (112,*,iostat=eof) sdrain(ihru)
      if (eof < 0) exit
      read (112,*,iostat=eof) drain_co(ihru)
      if (eof < 0) exit
      read (112,*,iostat=eof) pc(ihru)
      if (eof < 0) exit
      read (112,*,iostat=eof) latksatf(ihru)
      if (eof < 0) exit
      read (112,*,iostat=eof) sstmaxd(ihru)
      if (eof < 0) exit
   end do

   close (112)

   return
5000 format (a)

end
