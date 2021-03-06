!> @file readwus.f90
!> file containing the subroutine readwus
!> @author
!> modified by Javier Burguete

!> This subroutine reads data from the HRU/subbasin water use input file
!> (.wus). The water use file extracts water from the subbasin and it is
!> considered to be lost from the watershed. These variables should be used
!> to remove water transported outside the watershed.
!> @param[in] i subbasin number
subroutine readwus(i)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |subbasin number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    wudeep(:,:) |10^4 m^3/day  |average daily water removal from the deep
!!                               |aquifer for the month
!!    wupnd(:,:)  |10^4 m^3/day  |average daily water removal from the pond
!!                               |for the month
!!    wurch(:,:)  |10^4 m^3/day  |average daily water removal from the reach
!!                               |for the month
!!    wushal(:,:) |10^4 m^3/day  |average daily water removal from the shallow
!!                               |aquifer for the month
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag
!!    j           |none          |counter
!!    l           |none          |HRU number
!!    mon         |none          |counter
!!    swudp
!!    swupnd
!!    swush
!!    titldum     |NA            |title line of .wus file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i
   character (len=80) :: titldum
   integer :: eof, j, l, mon
   real*8 :: swudp(12), swupnd(12), swush(12)

   eof = 0
   swupnd = 0.
   swush = 0.
   swudp = 0.

   do
      read (105,5300,iostat=eof) titldum
      if (eof < 0) exit
      read (105,5300,iostat=eof) titldum
      if (eof < 0) exit
      read (105,5300,iostat=eof) titldum
      if (eof < 0) exit
      read (105,5100,iostat=eof) (swupnd(mon),mon = 1,6)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (swupnd(mon),mon = 7,12)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (wurch(mon,i),mon = 1,6)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (wurch(mon,i),mon = 7,12)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (swush(mon),mon = 1,6)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (swush(mon),mon = 7,12)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (swudp(mon),mon = 1,6)
      if (eof < 0) exit
      read (105,5100,iostat=eof) (swudp(mon),mon = 7,12)
      exit
   end do

   do j = 1, hrutot(i)
      l = nhru + j
      do mon = 1, 12
         wupnd(mon,l) = swupnd(mon)
         wushal(mon,l) = swush(mon)
         wudeep(mon,l) = swudp(mon)
      end do
   end do

   close (105)

   return
5100 format (6f10.1)
5300 format (a80)
end
