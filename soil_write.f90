!> @file soil_write.f90
!> file containing the subroutine soil_write
!> @author
!> modified by Javier Burguete

!> this subroutine writes output to the output.sol file
!> @param[in] i current day in simulation - loop counter (julian date)
subroutine soil_write(i)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |julian date   |current day in simulation--loop counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |counter
!!    l           |none          |counter
!!    solno3_t
!!    solorgn_t
!!    solorgp_t
!!    solp_t
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: i
   real*8 :: solno3_t, solorgn_t, solorgp_t, solp_t
   integer :: j, l

   do j = 1,nhru
      solp_t = 0.
      solno3_t = 0.
      solorgn_t = 0.
      solorgp_t = 0.
      do l = 1,sol_nly(j)
         solp_t = solp_t + sol_solp(l,j)
         solno3_t = solno3_t + sol_no3(l,j)

         select case (cswat)
          case (0)
            solorgn_t = solorgn_t + sol_orgn(l,j)
          case (1)
            solorgn_t = solorgn_t + sol_n(l,j)
          case (2)
            solorgn_t = solorgn_t + sol_HSN(l,j) + sol_HPN(l,j)
         end select

         solorgp_t = solorgp_t + sol_orgp(l,j)
      end do
      write (121,1000) i, subnum(j), hruno(j), sol_rsd(1,j), solp_t,&
         &solno3_t, solorgn_t, solorgp_t, cnday(j)
   end do

   return
1000 format ('SNU   ',i4,1x,a5,a4,1x,6f10.2)
end
