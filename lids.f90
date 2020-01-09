!> @file lids.f90
!> file containing the subroutine lids
!> @author
!> modified by Javier Burguete

!> call subroutines to simulate green roof, rain garden, cistern and porous
!> pavement processes
!> @param[in] sb subbasin number (none)
!> @param[in] j HRU number (none)
!> @param[in] k subdaily time index (none)
!> @param[in] lid_prec
!> precipitation depth a LID receives in a simulation time interval (mm)
subroutine lids(sb,j,k,lid_prec)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~ ~ ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sb               |none          |Subbasin number
!!    j                |none          |HRU number
!!    k                |none          |Subdaily time index
!!    lid_prec         |mm            |Precipitation depth a LID receives in a simulation time interval
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    fr
!!    jj               |none          |counter

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: lid_greenroof, lid_raingarden, lid_cistern, lid_porpavement

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: sb, j, k
   real*8, intent(in) :: lid_prec
   real*8 :: fr
   integer :: jj

   jj = urblu(j)

   lid_farea = 0.
   lid_sw_add = 0.
   fr = 0.

   if (gr_onoff(sb,jj)==1) then
      select case (jj)
       case (14, 18)
         fr = 0.80
       case (19, 20, 21, 24)
         fr = 0.40
      end select
      lid_farea(1,j) = gr_farea(sb,jj) * fr ! assuming fractions of roof area to impervious area is 0.80 and 0.40 for regidential and commercial/industrial, respectively
      call lid_greenroof(sb,j,k,lid_prec)
   end if

   if (rg_onoff(sb,jj)==1) then
      lid_farea(2,j) = rg_farea(sb,jj)
      call lid_raingarden(sb,j,k,lid_prec)
   end if

   if (cs_onoff(sb,jj)==1) then
      select case (jj)
       case (14, 18)
         fr = 0.80
       case (19, 20, 21, 24)
         fr = 0.40
      end select
      if (cs_grcon(sb,jj)==0) then
         lid_farea(3,j) = cs_farea(sb,jj)
         call lid_cistern(sb,j,k,lid_prec)
      else
         lid_farea(3,j) = gr_farea(sb,jj) * fr ! assuming fractions of roof area to impervious area is 0.80 and 0.40 for regidential and commercial/industrial, respectively
         call lid_cistern(sb,j,k,lid_qsurf(j,1))
         lid_qsurf(j,1) = 0.
      end if
   end if

   if (pv_onoff(sb,jj)==1) then
      select case (jj)
       case (14, 18)
         fr = 0.05
       case (19, 20, 21, 24, 26, 30, 31)
         fr = 0.40
      end select
      lid_farea(4,j) = pv_farea(sb,jj) * fr
      call lid_porpavement(sb,j,k,lid_prec)
   end if

   return
end subroutine
