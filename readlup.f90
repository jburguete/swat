!> @file readlup.f90
!> file containing the subroutine readlup
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU/subbasin management input file
!> (.mgt). This file contains data related to management practices used in
!> the HRU/subbasin.
subroutine readlup

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!                                 |urban.dat
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    fname
!!    ida_lup
!!    iyr_lup
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag (=-1 if eof, else = 0)
!!    fnam
!!    iday
!!    iyear
!!    mo
!!    no_up
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: Jdt

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer :: Jdt
   character*13 fnam
   integer :: eof, iday, iyear, mo, no_up

!!    read lup file
   do
      read (122,100,iostat=eof) no_up, mo, iday, iyear, fnam
      if (eof < 0 .or. no_up == 0) exit

!!     compute julian days
      ida_lup(no_up) = Jdt(ndays,iday,mo)
      iyr_lup(no_up) = iyear

      fname(no_up) = fnam
   end do

100 format (4i5,1x,a)

   return
end
