!> @file pesty.f90
!> file containing the subroutine pesty
!> @author
!> modified by Javier Burguete

!> @param[in] iwave
!> flag to differentiate calculation of HRU and subbasin sediment calculation
!> (none)\n
!> iwave = 0 for HRU\n
!> iwave = subbasin # for subbasin
!> @param[in] j HRU number
subroutine pesty(iwave, j)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates pesticide transported with suspended sediment

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    iwave         |none         |flag to differentiate calculation of HRU and
!!                                |subbasin sediment calculation
!!                                |iwave = 0 for HRU
!!                                |iwave = subbasin # for subbasin
!!    j             |none         |HRU number
!!    enratio       |none         |enrichment ratio calculated for day in HRU
!!    hrupest(:)    |none         |pesticide use flag:
!!                                | 0: no pesticides used in HRU
!!                                | 1: pesticides used in HRU
!!    npmx          |none         |number of different pesticides used in
!!                                |the simulation
!!    npno(:)       |none         |array of unique pesticides used in watershed
!!    pst_enr(:,:)  |none         |pesticide enrichment ratio
!!    sedyld(:)     |metric tons  |daily soil loss caused by water erosion in
!!                                |HRU
!!    sol_kp(:,:,:)|(mg/kg)/(mg/L)|pesticide sorption coefficient, Kp; the
!!                                |ratio of the concentration in the solid
!!                                |phase to the concentration in solution
!!    sol_pst(:,:,:)|kg/ha        |amount of pesticide in layer in HRU
!!    sub_pst(:,:)  |kg/ha        |amount of pesticide in layer in subbasin
!!    zdb(:,:)      |mm           |division term from net pesticide equation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    pst_sed(:,:)  |kg/ha        |pesticide loading from HRU sorbed onto
!!                                |sediment
!!    sol_pst(:,:,:)|kg/ha        |amount of pesticide in layer in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    conc        |              |concentration of pesticide in soil
!!    er          |none          |enrichment ratio for pesticides
!!    k           |none          |counter
!!    kk          |none          |pesticide number from database
!!    xx          |kg/ha         |amount of pesticide in soil
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent (in) :: iwave, j
   real*8 :: conc, er, xx
   integer :: k, kk

   do k = 1, npmx
      kk = npno(k)
      if (kk > 0) then
         if (iwave <= 0) then
            xx = sol_pst(k,j,1)
         else
            ! xx = sub_pst(kk,iwave) ! sub_pst is 0
            xx = 0.
         end if

         if (xx >= .0001) then
            conc = 100. * sol_kp(k,j,1) * xx / (zdb(k,j)+1.e-10)
            if (pst_enr(k,j) > 0.) then
               er = pst_enr(k,j)
            else
               er = enratio
            end if

            pst_sed(k,j) = .001 * sedyld(j) * conc * er / hru_ha(j)
            if (pst_sed(k,j) < 0.) pst_sed(k,j) = 0.
            if (pst_sed(k,j) > xx) pst_sed(k,j) = xx
            sol_pst(k,j,1) = xx - pst_sed(k,j)
         end if
      end if
   end do

   return
end
