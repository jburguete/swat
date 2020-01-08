!> @file snom.f90
!> file containing the subroutine snom
!> @author
!> modified by Javier Burguete

!> this subroutine predicts daily snom melt when the average air
!> temperature exceeds 0 degrees Celsius
!> @param[in] j HRU number
subroutine snom(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j            |none          |HRU number
!!    elevb(:,:)   |m             |elevation at center of band
!!    elevb_fr(:,:)|none          |fraction of subbasin area within elevation
!!                                |band
!!    iida         |julian date   |day being simulated (current julian date)
!!    pcpband(:,:) |mm H2O        |precipitation for the day in band in HRU
!!    precipday    |mm H2O        |precipitation on the current day in the HRU
!!    sub_sftmp    |deg C         |Snowfall temperature
!!                                |Mean air temperature at which precipitation
!!                                |is equally likely to be rain as snow/freezing
!!                                |rain.
!!    sub_smfmn    |mm/deg C/day  |Minimum melt rate for snow during year (Dec.
!!                                |21) where deg C refers to the air temperature
!!    sub_smfmx    |mm/deg C/day  |Maximum melt rate for snow during year (June
!!                                |21) where deg C refers to the air temperature
!!                                |SMFMX and SMFMN allow the rate of snow melt
!!                                |to vary through the year. These parameters
!!                                |are accounting for the impact of soil
!!                                |temperature on snow melt.
!!    sub_smtmp    |deg C         |Snow melt base temperature
!!                                |Mean air temperature at which snow melt will
!!                                |occur.
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day
!!    snocov1      |none          |1st shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 50% snow cover
!!    snocov2      |none          |2nd shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 95% snow cover
!!    snocovmx     |mm H2O        |Minimum snow water content that corresponds
!!                                |to 100% snow cover. If the snow water content
!!                                |is less than SNOCOVMX, then a certain
!!                                |percentage of the ground will be bare.
!!    snoeb(:,:)   |mm H2O        |snow water content in elevation band on
!!                                |current day
!!    snotmp(:)    |deg C         |temperature of snow pack in HRU
!!    snotmpeb(:,:)|deg C         |temperature of snow pack in elevation band
!!    sub_timp     |none          |Snow pack temperature lag factor (0-1)
!!                                |1 = no lag (snow pack temp=current day air
!!                                |temp) as the lag factor goes to zero, the
!!                                |snow pack's temperature will be less
!!                                |influenced by the current day's air
!!                                |temperature
!!    tavband(:,:) |deg C         |average temperature for the day in band in HRU
!!    tmpav(:)     |deg C         |average air temperature on current day for
!!                                |HRU
!!    tmx(:)       |deg C         |maximum air temperature on current day for
!!                                |HRU
!!    tmxband(:,:) |deg C         |maximum temperature for the day in band in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    precipday    |mm H2O        |amount of water in effective precipitation
!!                                |in HRU
!!    precipdt(:)  |mm H2O        |precipitation for the time step during day
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day
!!    snoeb(:,:)   |mm H2O        |snow water content in elevation band on
!!                                |current day
!!    snofall      |mm H2O        |amount of precipitation falling as freezing
!!                                |rain/snow on day in HRU
!!    snomlt       |mm H2O        |amount of water in snow melt for the day in
!!                                |HRU
!!    snotmp(:)    |deg C         |temperature of snow pack in HRU
!!    snotmpeb(:,:)|deg C         |temperature of snow pack in elevation band
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    ratio
!!    smfac       |
!!    smleb       |mm H2O        |amount of snow melt in elevation band on
!!                               |current day
!!    smp         |mm H2O        |precipitation on current day for HRU
!!    snocov      |none          |fraction of HRU area covered with snow
!!    sum1        |mm H2O        |snow water content in HRU on current day
!!    xx          |none          |ratio of amount of current day's snow water
!!                               |content to the minimum amount needed to
!!                               |cover ground completely
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sin, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: ratio, smfac, smleb, smp, snocov, sum1, xx
   integer :: ii, k

   sum1 =0.
   smp =0.
   k = hru_sub(j)

   if (elevb(1,k) > 0. .and. elevb_fr(1,k) > 0.) then
!! elevation bands
      !! compute snow fall and melt for each elevation band
      do ii = 1, 10
         if (subp(j) > 0.) then
            ratio = precipday / subp(j)
            pcpband(ii,j) = ratio * pcpband(ii,j)
         else
            pcpband(ii,j) = 0.
         end if
         if (elevb_fr(ii,k) < 0.) exit
         snotmpeb(ii,j) = snotmpeb(ii,j) * (1.-sub_timp(ii,k)) +&
            &tavband(ii,j) * sub_timp(ii,k)

         if (tavband(ii,j) < sub_sftmp(ii,k)) then

            !! compute snow fall if temperature is below sftmp
            snoeb(ii,j) = snoeb(ii,j) + pcpband(ii,j)
            snofall = snofall + pcpband(ii,j) * elevb_fr(ii,k)

         else

            !! compute snow melt if temperature is above smtmp
            if (tmxband(ii,j) > sub_smtmp(ii,k)) then
               smfac = (sub_smfmx(ii,k) + sub_smfmn(ii,k)) / 2. +&
                  &Sin((iida - 81) / 58.09) *&
                  &(sub_smfmx(ii,k) - sub_smfmn(ii,k)) / 2.    !! 365/2pi = 58.09
               smleb = smfac * (((snotmpeb(ii,j) + tmxband(ii,j)) / 2.)&
                  &- sub_smtmp(ii,k))

               !! adjust for areal extent of snow cover
               if (snoeb(ii,j) < snocovmx) then
                  xx = snoeb(ii,j) / snocovmx
                  snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
               else
                  snocov = 1.
               endif
               smleb = smleb * snocov
               if (smleb < 0.) smleb = 0.
               if (smleb > snoeb(ii,j)) smleb = snoeb(ii,j)
               snoeb(ii,j) = snoeb(ii,j) - smleb
               snomlt = snomlt + smleb * elevb_fr(ii,k)
            endif
         endif
         sum1 = sum1 + snoeb(ii,j) * elevb_fr(ii,k)
         smp = smp + pcpband(ii,j) * elevb_fr(ii,k)

      end do



      !! add/sub aggregate snow fall and melt from effective precip
      !! and snow cover
      precipday = smp + snomlt - snofall
      if (precipday < 0.) precipday = 0.
      if (nstep > 0) then
         do ii = 2, nstep+1
            precipdt(ii) = precipdt(ii) + (snomlt - snofall) / nstep
            if (precipdt(ii) < 0.) precipdt(ii) = 0.
         end do
      end if
      sno_hru(j) = sum1

   else
!! no elevation bands

      ii = 1

      !! estimate snow pack temperature
      snotmp(j)=snotmp(j) * (1. - sub_timp(ii,k)) + tmpav(j) *&
         &sub_timp(ii,k)

      if (tmpav(j) <= sub_sftmp(ii,k)) then
         !! calculate snow fall
         sno_hru(j) = sno_hru(j) + precipday
         snofall = precipday
         precipday = 0.
         precipdt = 0.
      endif

      if (tmx(j) > sub_smtmp(ii,k) .and. sno_hru(j) > 0.) then
         !! adjust melt factor for time of year
         smfac = (sub_smfmx(ii,k) + sub_smfmn(ii,k)) / 2. +&
            &Sin((iida - 81) / 58.09) *&
            &(sub_smfmx(ii,k) - sub_smfmn(ii,k)) / 2.    !! 365/2pi = 58.09
         snomlt = smfac * (((snotmp(j) + tmx(j)) / 2.) - sub_smtmp(ii,k))

         !! adjust for areal extent of snow cover
         if (sno_hru(j) < snocovmx) then
            xx = sno_hru(j) / snocovmx
            snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
         else
            snocov = 1.
         endif
         snomlt = snomlt * snocov
         if (snomlt < 0.) snomlt = 0.
         if (snomlt > sno_hru(j)) snomlt = sno_hru(j)
         sno_hru(j) = sno_hru(j) - snomlt
         precipday = precipday + snomlt
         if (nstep > 0) then
            do ii = 2, nstep+1
               precipdt(ii) = precipdt(ii) + snomlt / nstep
            end do
         end if
         if (precipday < 0.) precipday = 0.
      else
         snomlt = 0.
      end if
   end if
   return
end
