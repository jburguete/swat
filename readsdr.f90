!> @file readsdr.f90
!> file containing the subroutine readsdr
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU/subbasin management input file
!> (.mgt). This file contains data related to management practices used in
!> the HRU/subbasin.
!> @param[in] j HRU number (none)
subroutine readsdr(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name      |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j         |none             |HRU number
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
!!    eof         |none          |end of file flag (=-1 if eof, else = 0)
!!    titldum     |NA            |title line from input dataset
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   character (len=80) :: titldum
   integer :: j, eof

   do
      read (112,5000,iostat=eof) titldum
      if (eof < 0) exit
!!      read scheduled operations
      read (112,*,iostat=eof) re(j)
      if (eof < 0) exit
      read (112,*,iostat=eof) sdrain(j)
      if (eof < 0) exit
      read (112,*,iostat=eof) drain_co(j)
      if (eof < 0) exit
      read (112,*,iostat=eof) pc(j)
      if (eof < 0) exit
      read (112,*,iostat=eof) latksatf(j)
      if (eof < 0) exit
      read (112,*,iostat=eof) sstmaxd(j)
      if (eof < 0) exit
   end do

   close (112)

   return
5000 format (a)

end
