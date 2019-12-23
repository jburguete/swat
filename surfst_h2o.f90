!> @file surfst_h2o.f90
!> file containing the subroutine surfst_h2o
!> @author
!> modified by Javier Burguete

!> this subroutine determines the net surface runoff reaching the
!> main channel on a given day. The net amount of water reaching
!> the main channel can include water in surface runoff from the
!> previous day and will exclude surface runoff generated on the
!> current day which takes longer than one day to reach the main
!> channel
subroutine surfst_h2o

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    brt(:)      |none          |fraction of surface runoff that takes
!!                               |one day or less to reach the subbasin
!!                               |outlet
!!    hhqday(:)   |mm H2O        |surface runoff generated in HRU on the
!!                               |current hour at current day
!!    ihru        |none          |HRU number
!!    surf_bs1(:)|mm H2O        |amount of surface runoff lagged over one
!!                               |day
!!    surfq(:)    |mm H2O        |surface runoff generated in HRU on the
!!                               |current day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bsprev      |mm H2O        |surface runoff lagged from prior day
!!    hhqday(:)   |mm H2O        |surface runoff generated in HRU on the
!!                               |current hour at current day
!!    qday        |mm H2O        |surface runoff loading to main channel
!!                               |from HRU on current day
!!    surf_bs1(:)|mm H2O        |amount of surface runoff lagged over one
!!                               |day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    k           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max, Sum

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer :: j, k

   j = ihru

   if (ievent == 0) then

      bsprev = surf_bs1(j)
      surf_bs1(j) = Max(1.e-6, surf_bs1(j) + surfq(j))
      qday = surf_bs1(j) * brt(j)
      surf_bs1(j) = surf_bs1(j) - qday

   else
      !subdaily runoff lag (applies only to runoff from pervious surface)
      bsprev = hhsurf_bs1(j,nstep)  ! lag from previous day J.Jeong 4/06/2009

      do k=1,nstep

         !! Left-over (previous timestep) + inflow (current  timestep)
         hhsurf_bs1(j,k) = Max(0., bsprev + hhqday(k))

         !! new estimation of runoff and sediment reaching the main channel
         hhqday(k) = hhsurf_bs1(j,k) * brt(j)
         hhsurf_bs1(j,k) = hhsurf_bs1(j,k) - hhqday(k)

         !! lagged at the end of time step
         bsprev = hhsurf_bs1(j,k)
      end do

      !! daily total yield from the HRU
      qday = Sum(hhqday) + Sum(ubnrunoff)
   end if

   return
end
