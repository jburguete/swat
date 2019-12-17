subroutine substor

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine stores and lags lateral soil flow and nitrate

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bss1(:)      |mm H2O       |amount of lateral flow lagged
!!    bss2(:)      |kg N/ha      |amount of nitrate in lateral flow lagged

!!    bss3(:)      |mm           |amount of tile flow lagged
!!    bss4(:)      |kg N/ha      |amount of nitrate in tile flow lagged
!!    hrupest(:)    |none         |pesticide use flag:
!!                                | 0: no pesticides used in HRU
!!                                | 1: pesticides used in HRU
!!    ihru          |none         |HRU number
!!    lat_pst(:)    |kg pst/ha    |amount of pesticide in lateral flow in HRU
!!                                |for the day
!!    lat_ttime(:)  |none         |Exponential of the lateral flow travel time
!!    latno3(:)     |kg N/ha      |amount of NO3-N in lateral flow in HRU for
!!                                |the day
!!    latq(:)       |mm H2O       |amount of water in lateral flow in HRU for
!!                                |the day
!!    qtile(:)      |mm H2O       |amount of water in tile flow in HRU for the day
!!    tile_ttime(:) |none         |Exponential of tile flow travel time
!!    pst_lag(:,3,:)|kg pst/ha    |amount of pesticide lagged
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bss1(:)      |mm H2O        |amount of lateral flow lagged
!!    bss2(:)      |kg N/ha       |amount of nitrate in lateral flow lagged
!!    bss3(:)      |mm            |amount of tile flow lagged
!!    bss4(:)      |kg N/ha       |amount of nitrate in tile flow lagged
!!    bssprev       |mm H2O        |lateral flow lagged from prior day of
!!                                 |simulation
!!    lat_pst(:)    |kg pst/ha     |amount of pesticide in lateral flow in HRU
!!                                 |for the day
!!    latno3(:)     |kg N/ha       |amount of NO3-N in lateral flow in HRU for
!!                                 |the day
!!    latq(:)       |mm H2O        |amount of water in lateral flow in HRU for
!!                                 |the day
!!    pst_lag(:,3,:)|kg pst/ha     |amount of pesticide lagged
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    k           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer :: j, k

   j = ihru

   bssprev = bss1(j)
   bss1(j) = bss1(j) + latq(j)
   bss2(j) = bss2(j) + latno3(j)
   bss3(j) = bss3(j) + qtile
   bss4(j) = bss4(j) + tileno3(j)
   if (bss1(j) < 1.e-6) bss1(j) = 0.0
   if (bss2(j) < 1.e-6) bss2(j) = 0.0
   if (bss3(j) < 1.e-6) bss3(j) = 0.0
   if (bss4(j) < 1.e-6) bss4(j) = 0.0

   if (hrupest(j) == 1) then
      do k = 1, npmx
         if (pst_lag(k,3,j) < 1.e-6) pst_lag(k,3,j) = 0.0
         !MFW, 3/3/12: Modified lagged pesticide to include decay in lag
         pst_lag(k,3,j) = (pst_lag(k,3,j) * decay_s(npno(k))) + lat_pst(k)
         ! pst_lag(k,3,j) = pst_lag(k,3,j) + lat_pst(k)
      end do
   end if

   latq(j) = bss1(j) * lat_ttime(j)
   latno3(j) = bss2(j) * lat_ttime(j)
   qtile = bss3(j) * tile_ttime(j)
   tileno3(j) = bss4(j) * tile_ttime(j)
   if (latq(j) < 1.e-6) latq(j) = 0.
   if (latno3(j) < 1.e-6) latno3(j) = 0.
   if (qtile < 1.e-6) qtile = 0.
   if (tileno3(j) < 1.e-6) tileno3(j) = 0.
   if (hrupest(j) == 1) then
      do k = 1, npmx
         lat_pst(k) = pst_lag(k,3,j) * lat_ttime(j)
      end do
   end if

   bss1(j) = bss1(j) - latq(j)
   bss2(j) = bss2(j) - latno3(j)
   bss3(j) = bss3(j) - qtile
   bss4(j) = bss4(j) - tileno3(j)
   if (hrupest(j) == 1) then
      do k = 1, npmx
         pst_lag(k,3,j) = pst_lag(k,3,j) - lat_pst(k)
      end do
   end if

   return
end
