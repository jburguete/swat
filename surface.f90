!> @file surface.f90
!> file containing the subroutine surface
!> @author
!> modified by Javier Burguete

!> this subroutine models surface hydrology at any desired time step
!> @param[in] i current day in simulation--loop counter (julian date)
!> @param[in] j HRU number (none)
!> @param[in] sb subbasin number (none)
subroutine surface(i, j, sb)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |julian date   |current day in simulation--loop counter
!!    j           |none          |HRU number
!!    ovrlnd(:)   |mm H2O        |overland flow onto HRU from upstream
!!                               |routing unit
!!    peakr       |mm/hr         |peak runoff rate
!!    precipday   |mm H2O        |effective precipitation for the day in HRU
!!    qday        |mm H2O        |surface runoff loading to main channel
!!                               |for day
!!    surfq(:)    |mm H2O        |surface runoff generated in HRU during
!!                               |the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    precipday   |mm H2O        |effective precipitation for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hruvirr
!!    ii          |none          |counter
!!    irfr
!!    irmmdt
!!    kk          |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: canopyint, snom, crackvol, dailycn, volq, crackflow, surfst_h2o,
!!    SWAT: alph, pkq, tran, eiusle, ovr_sed, cfactor, ysed

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: i, j, sb
   real*8, dimension(nstep) :: irmmdt
   real*8 :: hruvirr, irfr
   integer :: ii, kk

   irmmdt = 0.

   !! compute canopy interception
   if (idplt(j) > 0) then
      call canopyint(j)
   end if

   !! compute snow melt
   call snom(j)

   !! output by elevation band to output.snw
   if (isnow == 1) then
      write(115,1010) i, iyr, subnum(j), hruno(j), (snoeb(ii,j), ii = 1,10)
   end if

   if (isnow ==1) then
      write (116,1010) i, iyr, subnum(j), hruno(j), (tavband(ii,j), ii = 1, 10)
   end if

   !! compute crack volume
   if (icrk == 1) call crackvol(j)

   !! add overland flow from upstream routing unit
   precipday = precipday + ovrlnd(j)
   ! ovrlnd_dt seems to be always 0
   !if (nstep > 0) then
   !   do ii = 1, nstep
   !      precipdt(ii+1) = precipdt(ii+1) + ovrlnd_dt(j,ii)
   !   end do
   !end if

   !! add irrigation from retention-irrigation ponds to soil water
   if (ri_luflg(j)==1) then
      kk = hru_sub(j)
      irfr = hru_km(j) * (1. - fimp(urblu(j))) / ri_subkm(kk)
      do ii=1,nstep
         !amount irrigated in hru
         hruvirr = ri_totpvol(ii) * irfr !m3
         irmmdt(ii) = hruvirr / (hru_km(j)&
            &* (1.- fimp(urblu(j))) * 1000.) !mm/dt

         !add irrigated water to soil water content
         do kk=1,sol_nly(j)
            if(irmmdt(ii)<sol_ul(kk,j)-sol_st(kk,j)) then
               sol_st(kk,j) = sol_st(kk,j) + irmmdt(ii)
               exit
            else
               sol_st(kk,j) = sol_ul(kk,j)
               irmmdt(ii) = irmmdt(ii) - (sol_ul(kk,j)-sol_st(kk,j))
            end if
         end do

      end do
   end if

   !!calculate subdaily curve number value
   call dailycn(j)

   !! compute runoff - surfq in mm H2O
   if (precipday > 0.1) then
      call volq(j)
      !! bmp adjustment
      surfq(j) = surfq(j) * bmp_flo(j)
      !! adjust runoff for loss into crack volume
      if (surfq(j) > 0. .and. icrk == 1) call crackflow(j)
   end if

   surfq(j) = surfq(j) + qird(j)
   qird(j) = 0.

   !! calculate amount of surface runoff reaching main channel during day
   !! (qday) and store the remainder
   call surfst_h2o(j)

   !! calculate half-hour rainfall
   if (precipday > 0.01) call alph(0, j)

   if (qday > 0.0001) then
      !! compute peak rate - peakr in m3/s
      call pkq(0, j)
   end if

   if (qday > 0.0001 .and. peakr > 0.) then
      !! compute transmission losses for non-HUMUS datasets
      call tran(j)
      call eiusle(j)

      !! calculate sediment erosion by rainfall and overland flow
      call ovr_sed(j, sb)
   end if

   call cfactor(j)
   if (surfq(j) > 1.e-6 .and. peakr > 1.e-6) call ysed(0, j)

   if (qday < 0.) qday = 0.

1010 format (2(i4,1x),a5,a4,1x,10f8.1)
   return
end
