!> @file hmeas.f90
!> file containing the subroutine hmeas
!> @author
!> modified by Javier Burguete

!> this subroutine reads in relative humidity data from file and
!> assigns the data to the HRUs
subroutine hmeas

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru_sub(:)  |none          |subbasin in which HRU is located
!!    id1         |julian date   |first day of simulation in current year
!!    ifirsth     |none          |relative humidity data search code
!!                               |0 first day of relative humidity data located
!!                               |  in file
!!                               |1 first day of relative humidity data not
!!                               |  located in file
!!    ihgage(:)   |none          |HRU relative humidity data code (gage # for
!!                               |relative humidity data used in HRU)
!!    iyr         |none          |beginning year of simulation
!!    mrg         |none          |maximum number of rainfall/temp gages
!!    nhru        |none          |number of HRUs in watershed
!!    nhtot       |none          |total number of relative humidity records
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ifirsth     |none          |relative humidity data search code
!!                               |0 first day of relative humidity data located
!!                               |  in file
!!                               |1 first day of relative humidity data not
!!                               |  located in file
!!    rhd(:)      |none          |relative humidity for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    idap        |julian date   |julian date of measured weather data
!!    ih          |              |counter
!!    inum3sprev  |none          |subbasin number of previous HRU
!!    iyp         |none          |last 2 digits of year measured weather data
!!    k           |none          |counter
!!    l           |none          |counter
!!    rhdbsb      |none          |generated relative humidity for subbasin
!!    rhmeas(:)   |none          |relative humidity read in from file (fraction)
!!    tmpmean     |deg C         |average temperature for the month in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: rhgen

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8, dimension (mrg) :: rhmeas
   real*8 :: rhdbsb, tmpmean
   integer :: idap, ih, inum3sprev, iyp, k, l

   !! initialize variables for the day
   rhmeas = 0.

   !! read relative humidity data from file
   if (ifirsth == 0) then
      read (138,5200) (rhmeas(l), l = 1, nhtot)
   else
      ifirsth = 0
      do
         iyp = 0
         idap = 0
         read (138,5300) iyp, idap, (rhmeas(l), l = 1, nhtot)
         if (iyp + idap <= 0) exit
         if (iyp == iyr .and. idap == id1) exit
      end do
   end if


   !! assign relative humidity data to HRUs
   inum3sprev = 0
   rhdbsb = 0.
   do k = 1, nhru
      l = hru_sub(k)
      ih = ihgage(l)
      !! generate values to replace missing data
      if (rhmeas(ih) <  -97.) then
         !! use same generated data for all HRUs in a subbasin
         if (l == inum3sprev .and. l /= 0) then
            rhd(k) = rhdbsb
         else
            call rhgen(k)
            !! set subbasin generated values
            inum3sprev = l
            rhdbsb = rhd(k)
         end if
      else
         if (rhmeas(ih) < 1. .and. rhmeas(ih) > 0.) then
            rhd(k) = rhmeas(ih)
         else
            tmpmean=(tmpmx(i_mo,l)+tmpmn(i_mo,l))/2.
            rhd(k) = Ee(rhmeas(ih)) / Ee(tmpmean)
         endif
      end if
   end do

   return
5200 format (7x,1800f8.3)
5300 format (i4,i3,1800f8.3)
end
