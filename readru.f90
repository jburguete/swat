!> @file readru.f90
!> file containing the subroutine readru
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the sub input file (.sub).
!> This file contains data related to routing
!> @param[in] i subbasin number (none)
!> @param[in] k subbasin number (none)
subroutine readru(i, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definitionov
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |subbasin number
!!    k           |none          |subbasin number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chk_ru      |(mm/hr)       |eff hydr cond
!!    chl_ru      |(km)          |channel length
!!    chn_ru      |              |Manning's N tributary channels
!!    chs_ru      |(m/m)         |ave slope
!!    chw_ru      |(mm/km)       |ave width
!!    da_ru       |ha            |area of routing unit
!!    eof
!!    ix
!!    j
!!    ovn_ru      |              |Manning's N value overland flow
!!    ovs         |(m)           |average slope steepness
!!    ovsl        |(m)           |average slope length
!!    sumk
!!    tck
!!    titldum

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: i, k
   character (len=80) :: titldum
   real*8 :: chk_ru, chl_ru, chn_ru, chs_ru, chw_ru, da_ru, ix,&
      &ovn_ru, ovs, ovsl, sumk, tck
   integer :: eof, j

   eof = 0
   tck = 0.
   do
      read (113,5000,iostat=eof) titldum
      if (eof < 0) exit
      !       read (113,*,iostat=eof) tck
      !       if (eof < 0) exit
      read (113,*,iostat=eof) da_ru
      if (eof < 0) exit
      read (113,*,iostat=eof) ovsl
      if (eof < 0) exit
      read (113,*,iostat=eof) ovs
      if (eof < 0) exit
      read (113,*,iostat=eof) ovn_ru
      if (eof < 0) exit
      read (113,*,iostat=eof) chl_ru
      if (eof < 0) exit
      read (113,*,iostat=eof) chs_ru
      if (eof < 0) exit
      read (113,*,iostat=eof) chw_ru
      if (eof < 0) exit
      read (113,*,iostat=eof) chk_ru
      if (eof < 0) exit
      read (113,*,iostat=eof) chn_ru
      if (eof < 0) exit
      exit
   end do

   if (ovsl < 1.e-6) ovsl = 50.

   do j = 1, hrutot(i)
      read (113,*) ix, hru_rufr(iru,j)
   end do

   !! compute weighted K factor for sediment transport capacity
   sumk = 0.
   do j = 1, hrutot(i)
      sumk = sumk + usle_k(j) * hru_rufr(iru,j)
   end do
   ru_k(k,iru) = sumk
   ru_ovsl(k,iru) = ovsl
   ru_ovs(k,iru) = ovs
   ru_ktc(k,iru) = tck ! tck is not readed
   !daru_km(k,iru) = da_ru

5000 format (a)
   return
end
