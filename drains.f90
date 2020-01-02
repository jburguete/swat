!> @file drains.f90
!> file containing the subroutine drains
!> @author
!> modified by Javier Burguete

!> this subroutine finds the effective lateral hydraulic conductivity
!> and computes drainage or subirrigation flux
!> @param[in] j HRU number
subroutine drains(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    conk(:,:)   |mm/hr         |lateral saturated hydraulic conductivity for each profile
!!                               |layer in a give HRU. For example (conk(2,1) is conductivity
!!                               |of layer from sol_z(1,1) to sol_z(2,1) in HRU1
!!    curyr       |none          |current year in simulation (sequence)
!!    drain_co(:) |mm/day        |drainage coefficient
!!    ddrain(:)   |mm            |depth of drain tube from the soil surface
!!    hru_slp(:)  |m/m           |average slope steepness in HRU
!!    id1         |julian date   |first day of simulation in current year
!!    latksatf(:) |none          |multiplication factor to determine conk(j1,j) from sol_k(j1,j) for HRU
!!    pc(:)       |mm/hr         |pump capacity (default pump capacity = 1.042mm/hr or 25mm/day)
!!    sdrain(:)   |mm            |distance between two drain tubes or tiles
!!    sstmaxd(:)  |mm            |static maximum depressional storage; read from .sdr
!!    sol_k(:,:)  |mm/hr         |saturated hydraulic conductivity of soil
!!                               |layer
!!    sol_nly(:)  |none          |number of layers in soil profile
!!    sol_z(:,:)  |mm            |depth to bottom of each profile layer in a given HRU
!!    stmaxd(:)   |mm            |maximum surface depressional storage for the day in a given HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    qtile       |mm H2O        |drainage tile flow in soil profile for the day

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    above       |mm            |depth of top layer considered
!!    ad
!!    ap
!!    depth       |mm            |actual depth from surface to impermeable layer
!!    cone        |mm/hr         |effective saturated lateral conductivity - based
!!                               |on water table depth and conk/sol_k of layers
!!    ddranp      |mm            |a variable used to indicate distance slightly less
!!                               |than ddrain. Used to prevent calculating subirrigation
!!                               |when water table is below drain bottom or when it is empty
!!    deep        |mm            |total thickness of saturated zone
!!    depth       |mm            |effective depth to impermeable layer from soil surface
!!                               |effective depth may be smaller than actual depth to account
!!                               |for convergence near drain tubes
!!    dflux       |mm/hr         |drainage flux
!!    dg          |mm            |depth of soil layer
!!    dot         |mm            |actual depth from impermeable layer to water level
!!                               |above drain during subsurface irrigation
!!    em          |mm            |distance from water level in the drains to water table
!!                               |at midpoint: em is negative during subirrigation
!!    gee         |none          |factor -g- in Kirkham equation
!!    gee1
!!    gee2
!!    gee3
!!    hdrain      |mm            |equivalent depth from water surface in drain tube to
!!                               |impermeable layer
!!    j1          |none          |counter
!!    k2
!!    k3
!!    k4
!!    k5
!!    k6
!!    m
!!    nlayer      |none          |number of layers to be used to determine cone
!!    pi
!!    stor        |mm            |surface storage for the day in a given HRU
!!    storro	  |mm            |surface storage that must b
!!                               |can move to the tile drain tube
!!    sum1
!!    x
!!    xx
!!    y1          |mm            |dummy variable for dtwt
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log, Tan, Cosh, Cos, Abs
!!    SWAT: depstor

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8, parameter :: pi = 22./7.
   real*8 :: above, ad, ap, cone, ddranp, deep, depth, dflux, dg, dot, em, gee,&
      &gee1, gee2, gee3, hdmin, hdrain, k2, k3, k4, k5, k6, stor, storro, sum1,&
      &x, xx, y1
   integer :: j1, m, nlayer

   !! initialize variables

   wnan = 0
   y1 = dep_imp(j) - wt_shall
   if (y1 > dep_imp(j)) y1 = dep_imp(j)
   above = 0.
   gee1 = 0.

!! find number of soil layers
   nlayer = 1
   do j1 = 1, mlyr
      if(sol_z(j1,j) > 0.) nlayer = j1
   end do

!! find effective lateral hydraulic conductivity for the profile in hru j
   do j1 = 1, nlayer
      if(y1 > sol_z(j1,j)) then
         wnan(j1) = 0.
      else
         wnan(j1) = sol_z(j1,j) - y1
         x = sol_z(j1,j) -  above
         if(wnan(j1) > x) wnan(j1) = x
      end if
      above = sol_z(j1,j)
   end do
   sum1 = 0.
   deep = 0.
   do j1=1,nlayer
      conk(j1,j) = sol_k(j1,j) * latksatf(j) !Daniel 2/26/09
      sum1 = sum1 + wnan(j1) * conk(j1,j)
      deep = deep + wnan(j1)
   end do
   if ((deep <= 0.001).or.(sum1 <= 0.001)) then
      sum1 = 0.
      deep = 0.001
      do j1=1,nlayer
         !! Compute layer depth ! Daniel 10/05/07
         dg = 0.
         if(j1 == 1) then
            dg = sol_z(j1,j)
         else
            dg = sol_z(j1,j) - sol_z(j1-1,j)
         end if
         !! Compute layer depth ! Daniel 10/05/07
         sum1 = sum1 + conk(j1,j) * dg !Daniel 10/09/07
         deep = deep + dg   !Daniel 10/09/07
      end do
   end if
   cone = sum1 / deep

   !! calculate parameters hdrain and gee1
   ad = dep_imp(j) - ddrain(j)
   ap = 3.55 - ((1.6 * ad) / sdrain(j)) + 2 * ((2 / sdrain(j))**2)
   if (ad/sdrain(j) < 0.3) then
      hdrain= ad / (1 + ((ad / sdrain(j)) * (((8 / pi) *&
         &Log(ad / re(j)) - ap))))
   else
      hdrain = ad
   end if
   !! calculate Kirkham G-Factor, gee
   k2 = Tan((pi * ((2. * ad) - re(j))) / (4. * dep_imp(j)))
   k3 = Tan((pi * re(j)) / (4. * dep_imp(j)))
   do m=1,2
      k4 = (pi * m * sdrain(j)) / (2. * dep_imp(j))
      k5 = (pi * re(j)) / (2. * dep_imp(j))
      k6 = (pi * (2. * ad - re(j))) / (2. * dep_imp(j))
      gee2 = (Cosh(k4) + Cos(k5)) / (Cosh(k4) - Cos(k5))
      gee3 = (Cosh(k4) - Cos(k6)) / (Cosh(k4) + Cos(k6))
      gee1 = gee1 + Log(gee2 * gee3)
   end do
   xx = k2 / k3
   if (xx < 1.) then
      gee = 1.
   else
      gee = 2 * Log(k2 / k3) + 2 * gee1
   end if
   if (gee < 1.) gee = 1.
   if (gee > 12.) gee = 12.

   !! calculate drainage and subirrigation flux section
   ! drainage flux for ponded surface
   depth = ddrain(j) + hdrain
   hdmin = depth - ddrain(j)
   if (ismax == 1) then
      call depstor(j) ! dynamic stmaxd(j): compute current HRU stmaxd based
      ! on cumulative rainfall and cum. intensity
   else
      stmaxd(j) = sstmaxd(j)
   end if
   storro = 0.2 * stmaxd(j) !surface storage that must be filled before surface
   !water can move to the tile drain tube
   !! Determine surface storage for the day in a given HRU (stor)
   !initialize stor on the beginning day of simulation, Daniel 9/20/2007
   if (curyr == 1 .and. iida == id1) then
      stor= 0.
   end if
   if (potsa(j) <= 0.) then ! determine stor
      stor = precipday - inflpcp - etday !Daniel 10/05/07
      if(surfq(j) > 0.0) stor=stmaxd(j)
   else
      stor = pot_vol(j)/(potsa(j)*1000)
   endif
   if (hdrain < hdmin) hdrain = hdmin
   if ((stor > storro) .and. (y1 < 5.0)) then
      dflux= (12.56637*24.0*cone*(depth-hdrain+stor))/&
         &(gee*sdrain(j)) !eq.10
      if(dflux > drain_co(j)) dflux = drain_co(j) !eq.11
   else
! subirrigation flux
      em = depth -y1 - hdrain
      if (em < -1.0) then
         ddranp = ddrain(j) - 1.0
         dot = hdrain + dep_imp(j) - depth
         dflux=4.0*24.0*cone*em*hdrain*(2.0+em/dot)/sdrain(j)**2
         if ((depth-hdrain) >= ddranp) dflux=0.
         if (Abs(dflux) > pc(j)) then
            dflux = -pc(j)*24.0
         end if
! drainage flux - for WT below the surface and for ponded depths < storro (S1)
      else
         dflux=4.0*24.0*cone*em*(2.0*hdrain+em)/sdrain(j)**2 !eq.5
         if (dflux > drain_co(j)) dflux=drain_co(j) !eq.11
         if (dflux < 0.) dflux=0.
         if (em < 0.) dflux=0.
      end if
   end if
   qtile=dflux

!     write(222,222) curyr, iida, hdrain, gee1, gee  !Daniel 3/1/09
!222   format(1x,4x,i4,4x,i3,4x,3f12.3)
   return
end
