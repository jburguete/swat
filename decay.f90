!> @file decay.f90
!> file containing the subroutine decay
!> @author
!> modified by Javier Burguete

!> this subroutine calculates degradation of pesticide in the soil and on
!> the plants
!> @param[in] j HRU number
subroutine decay(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j             |none          |HRU number
!!    decay_f(:)    |none          |exponential of the rate constant for
!!                                 |degradation of the pesticide on foliage
!!    decay_s(:)    |none          |exponential of the rate constant for
!!                                 |degradation of the pesticide in soil
!!    hru_dafr(:)   |none          |fraction of watershed area in HRU
!!    hrupest(:)    |none          |pesticide use flag:
!!                                 | 0: no pesticides used in HRU
!!                                 | 1: pesticides used in HRU
!!    npmx          |none          |number of different pesticides used in
!!                                 |the simulation
!!    npno(:)       |none          |array of unique pesticides used in watershed
!!    plt_pst(:,:)  |kg/ha         |pesticide on plant foliage
!!    sol_nly(:)    |none          |number of layers in soil profile
!!    sol_pst(:,:,:)|kg/ha         |pesticide in soil layer
!!    wshd_pstdg(:) |kg pst/ha     |amount of pesticide lost through degradation
!!                                 |in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    plt_pst(:,:)  |kg/ha         |pesticide on plant foliage
!!    sol_pst(:,:,:)|kg/ha         |pesticide in soil layer
!!    wshd_pstdg(:) |kg pst/ha     |amount of pesticide lost through degradation
!!                                 |in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k           |none          |counter
!!    kk          |none          |pesticide number from pest.dat
!!    l           |none          |counter (soil layers)
!!    x1          |kg/ha         |amount of pesticide present at beginning of
!!                               |day
!!    xx          |kg/ha         |amount of pesticide present at end of day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   integer :: k, kk, l
   real*8 :: x1, xx

   do k = 1, npmx
      kk = npno(k)
      if (kk > 0) then

         !! calculate degradation in soil
         do l = 1, sol_nly(j)
            x1 = sol_pst(k,j,l)
            if (x1 >= 0.0001) then
               xx = x1 * decay_s(kk)
               wshd_pstdg(k) = wshd_pstdg(k) + (x1 - xx) * hru_dafr(j)
               sol_pst(k,j,l) = xx
            end if
         end do

         !! calculate degradation off plant foliage
         x1 = plt_pst(k,j)
         if (x1 >= 0.0001) then
            xx = x1 * decay_f(kk)
            wshd_pstdg(k) = wshd_pstdg(k) + (x1 - xx) * hru_dafr(j)
            plt_pst(k,j) = xx
         end if
      end if
   end do

   return
end
