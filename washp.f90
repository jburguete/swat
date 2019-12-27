!> @file modparm.f90
!> file containing the module parm
!> @author
!> modified by Javier Burguete Tolosa

!> this subroutine calculates the amount of pesticide washed off the plant
!> foliage and onto the soil
!> @param[in] j HRU number
subroutine washp(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j             |none          |HRU number
!!    hrupest(:)    |none          |pesticide use flag:
!!                                 | 0: no pesticides used in HRU
!!                                 | 1: pesticides used in HRU
!!    npmx          |none          |number of different pesticides used in
!!                                 |the simulation
!!    npno(:)       |none          |array of unique pesticides used in watershed
!!    plt_pst(:,:)  |kg/ha         |pesticide on plant foliage
!!    pst_wof(:)    |none          |fraction of pesticide on foliage which
!!                                 |is washed-off by a rainfall event
!!    sol_pst(:,:,1)|kg/ha         |pesticide in first layer of soil
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    plt_pst(:,:)  |kg/ha         |pesticide on plant foliage
!!    sol_pst(:,:,1)|kg/ha         |pesticide in first layer of soil
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k           |none          |counter
!!    kk          |none          |pesticide number from pest.dat
!!    xx          |kg/ha         |amount of pesticide washed off foliage
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: xx
   integer :: k, kk

   do k = 1, npmx
      kk = npno(k)
      if (kk > 0 .and. plt_pst(k,j) >= 0.0001) then
         xx = plt_pst(k,j)
         if (pst_wof(kk) < 1.) xx = xx * pst_wof(kk)
         sol_pst(k,j,1) = sol_pst(k,j,1) + xx
         plt_pst(k,j) = plt_pst(k,j) - xx
      end if
   end do

   return
end
