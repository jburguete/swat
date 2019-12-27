!> @file nminrl.f90
!> file containing the subroutine nminrl
!> @author
!> modified by Javier Burguete

!> this subroutine estimates daily nitrogen and phosphorus
!> mineralization and immobilization considering fresh organic
!> material (plant residue) and active and stable humus material
!> @param[in] j HRU number
subroutine nminrl(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j             |none          |HRU number
!!    cmn           |none          |rate factor for humus mineralization on
!!                                 |active organic N
!!    curyr         |none          |current year of simulation
!!    hru_dafr(:)   |km**2/km**2   |fraction of watershed area in HRU
!!    icr(:)        |none          |sequence number of crop grown within the
!!                                 |current year
!!    idplt(:)      |none          |land cover code from crop.dat
!!    nactfr        |none          |nitrogen active pool fraction. The fraction
!!                                 |of organic nitrogen in the active pool.
!!    nro(:)        |none          |sequence number of year in rotation
!!    nyskip        |none          |number of years to skip output
!!                                 |summarization and printing
!!    rsdco_pl(:)   |none          |plant residue decomposition coefficient. The
!!                                 |fraction of residue which will decompose in
!!                                 |a day assuming optimal moisture,
!!                                 |temperature, C:N ratio, and C:P ratio
!!    sol_aorgn(:,:)|kg N/ha       |amount of nitrogen stored in the active
!!                                 |organic (humic) nitrogen pool
!!    sol_cbn(:,:)  |%             |percent organic carbon in soil layer
!!    sol_fon(:,:)  |kg N/ha       |amount of nitrogen stored in the fresh
!!                                 |organic (residue) pool
!!    sol_fop(:,:)  |kg P/ha       |amount of phosphorus stored in the fresh
!!                                 |organic (residue) pool
!!    sol_nly(:)    |none          |number of layers in soil profile
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool in soil layer
!!    sol_orgn(:,:) |kg N/ha       |amount of nitrogen stored in the stable
!!                                 |organic N pool
!!    sol_orgp(:,:) |kg P/ha       |amount of phosphorus stored in the organic
!!                                 |P pool in soil layer
!!    sol_rsd(:,:)  |kg/ha         |amount of organic matter in the soil
!!                                 |classified as residue
!!    sol_solp(:,:) |kg P/ha       |amount of phosohorus stored in solution
!!    sol_st(:,:)   |mm H2O        |amount of water stored in the soil layer on
!!                                 |current day
!!    sol_tmp(:,:)  |deg C         |daily average temperature of soil layer
!!    sol_ul(:,:)   |mm H2O        |amount of water held in the soil layer at
!!                                 |saturation
!!    wshd_dnit     |kg N/ha       |average annual amount of nitrogen lost from
!!                                 |nitrate pool due to denitrification in
!!                                 |watershed
!!    wshd_hmn      |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from active organic to nitrate pool in
!!                                 |watershed
!!    wshd_hmp      |kg P/ha       |average annual amount of phosphorus moving
!!                                 |from organic to labile pool in watershed
!!    wshd_rmn      |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from fresh organic (residue) to nitrate
!!                                 |and active organic pools in watershed
!!    wshd_rmp      |kg P/ha       |average annual amount of phosphorus moving
!!                                 |from fresh organic (residue) to labile
!!                                 |and organic pools in watershed
!!    wshd_rwn      |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from active organic to stable organic pool
!!                                 |in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hmntl         |kg N/ha       |amount of nitrogen moving from active
!!                                 |organic to nitrate pool in soil profile
!!                                 |on current day in HRU
!!    hmptl         |kg P/ha       |amount of phosphorus moving from the
!!                                 |organic to labile pool in soil profile
!!                                 |on current day in HRU
!!    rmn2tl        |kg N/ha       |amount of nitrogen moving from the fresh
!!                                 |organic (residue) to the nitrate(80%) and
!!                                 |active organic(20%) pools in soil profile
!!                                 |on current day in HRU
!!    rmptl         |kg P/ha       |amount of phosphorus moving from the
!!                                 |fresh organic (residue) to the labile(80%)
!!                                 |and organic(20%) pools in soil profile
!!                                 |on current day in HRU
!!    rwntl         |kg N/ha       |amount of nitrogen moving from active
!!                                 |organic to stable organic pool in soil
!!                                 |profile on current day in HRU
!!    sol_aorgn(:,:)|kg N/ha       |amount of nitrogen stored in the active
!!                                 |organic (humic) nitrogen pool
!!    sol_fon(:,:)  |kg N/ha       |amount of nitrogen stored in the fresh
!!                                 |organic (residue) pool
!!    sol_fop(:,:)  |kg P/ha       |amount of phosphorus stored in the fresh
!!                                 |organic (residue) pool
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool in soil layer
!!    sol_orgn(:,:) |kg N/ha       |amount of nitrogen stored in the stable
!!                                 |organic N pool
!!    sol_orgp(:,:) |kg P/ha       |amount of phosphorus stored in the organic
!!                                 |P pool in soil layer
!!    sol_rsd(:,:)  |kg/ha         |amount of organic matter in the soil
!!                                 |classified as residue
!!    sol_solp(:,:) |kg P/ha       |amount of phosohorus stored in solution
!!    wdntl         |kg N/ha       |amount of nitrogen lost from nitrate pool
!!                                 |by denitrification in soil profile on
!!                                 |current day in HRU
!!    wshd_dnit     |kg N/ha       |average annual amount of nitrogen lost from
!!                                 |nitrate pool due to denitrification in
!!                                 |watershed
!!    wshd_hmn      |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from active organic to nitrate pool in
!!                                 |watershed
!!    wshd_hmp      |kg P/ha       |average annual amount of phosphorus moving
!!                                 |from organic to labile pool in watershed
!!    wshd_rmn      |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from fresh organic (residue) to nitrate
!!                                 |and active organic pools in watershed
!!    wshd_rmp      |kg P/ha       |average annual amount of phosphorus moving
!!                                 |from fresh organic (residue) to labile
!!                                 |and organic pools in watershed
!!    wshd_rwn      |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from active organic to stable organic pool
!!                                 |in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ca          |
!!    cdg         |none          |soil temperature factor
!!    cdn         |
!!    cnr         |
!!    cnrf        |
!!    cpr         |
!!    cprf        |
!!    csf         |none          |combined temperature/soil water factor
!!    decr        |
!!    hmn         |kg N/ha       |amount of nitrogen moving from active organic
!!                               |nitrogen pool to nitrate pool in layer
!!    hmp         |kg P/ha       |amount of phosphorus moving from the organic
!!                               |pool to the labile pool in layer
!!    k           |none          |counter (soil layer)
!!    kk          |none          |soil layer used to compute soil water and
!!                               |soil temperature factors
!!    r4          |
!!    rdc         |
!!    rmn1        |kg N/ha       |amount of nitrogen moving from fresh organic
!!                               |to nitrate(80%) and active organic(20%)
!!                               |pools in layer
!!    rmp         |kg P/ha       |amount of phosphorus moving from fresh organic
!!                               |to labile(80%) and organic(20%) pools in layer
!!    rwn         |kg N/ha       |amount of nitrogen moving from active organic
!!                               |to stable organic pool in layer
!!    sdnco       |none          |denitrification threshold: fraction of field
!!                               | capacity
!!    sut         |none          |soil water factor
!!    wdn         |kg N/ha       |amount of nitrogen lost from nitrate pool in
!!                               |layer due to denitrification
!!    xx          |varies        |variable to hold intermediate calculation
!!                               |result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt, Max, Exp, Min, Abs, Dmax1

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: ca, cdg, cnr, cnrf, cpr, cprf, csf, decr, hmn, hmp, r4, rdc, rmn1,&
      &rmp, rwn, sut, wdn, xx
   integer :: k, kk

   do k = 1, sol_nly(j)

      if (k == 1) then
         kk = 2
      else
         kk = k
      end if

      !! mineralization can occur only if temp above 0 deg
      if (sol_tmp(kk,j) > 0.) then
         !! compute soil water factor
         if (sol_st(kk,j) < 0.) sol_st(kk,j) = .0000001
         sut = .1 + .9 * Sqrt(sol_st(kk,j) / sol_fc(kk,j))
         sut = Max(.05, sut)

         !!compute soil temperature factor
         xx = sol_tmp(kk,j)
         cdg = .9 * xx / (xx + Exp(9.93 - .312 * xx)) + .1
         cdg = Max(.1, cdg)

         !! compute combined factor
         xx = cdg * sut
         if (xx < 0.) xx = 0.
         if (xx > 1.e6) xx = 1.e6
         csf = Sqrt(xx)

         !! compute flow from active to stable pools
         rwn = .1e-4 * (sol_aorgn(k,j) * (1. / nactfr - 1.) - sol_orgn(k,j))
         if (rwn > 0.) then
            rwn = Min(rwn, sol_aorgn(k,j))
         else
            rwn = -(Min(Abs(rwn), sol_orgn(k,j)))
         endif
         sol_orgn(k,j) = Max(1.e-6, sol_orgn(k,j) + rwn)
         sol_aorgn(k,j) = Max(1.e-6, sol_aorgn(k,j) - rwn)

         !! compute humus mineralization on active organic n
         hmn = cmn(j) * csf * sol_aorgn(k,j)
         hmn = Min(hmn, sol_aorgn(k,j))
         !! compute humus mineralization on active organic p
         xx = sol_orgn(k,j) + sol_aorgn(k,j)
         if (xx > 1.e-6) then
            hmp = 1.4 * hmn * sol_orgp(k,j) / xx
         else
            hmp = 0.
         end if
         hmp = Min(hmp, sol_orgp(k,j))
         !! move mineralized nutrients between pools
         sol_aorgn(k,j) = Max(1.e-6, sol_aorgn(k,j) - hmn)
         sol_no3(k,j) = sol_no3(k,j) + hmn
         sol_orgp(k,j) = sol_orgp(k,j) - hmp
         sol_solp(k,j) = sol_solp(k,j) + hmp

         !! compute residue decomp and mineralization of
         !! fresh organic n and p (upper two layers only)
         rmn1 = 0.
         rmp = 0.
         if (k <= 2) then
            r4 = .58 * sol_rsd(k,j)

            if (sol_fon(k,j) + sol_no3(k,j) > 1.e-4) then
               cnr = r4 / (sol_fon(k,j) + sol_no3(k,j))
               if (cnr > 500.) cnr = 500.
               cnrf = Exp(-.693 * (cnr - 25.) / 25.)
            else
               cnrf = 1.
            end if

            if (sol_fop(k,j) + sol_solp(k,j) > 1.e-4) then
               cpr = r4 / (sol_fop(k,j) + sol_solp(k,j))
               if (cpr > 5000.) cpr = 5000.
               cprf = Exp(-.693 * (cpr - 200.) / 200.)
            else
               cprf = 1.
            end if

            ca = Min(cnrf, cprf, 1.)
            if (idplt(j) > 0) then
               decr = rsdco_pl(idplt(j)) * ca * csf
            else
               decr = 0.05
            end if
            decr = Max(decr_min, decr)
            decr = Min(decr, 1.)
            sol_rsd(k,j) = Dmax1(1.e-6,sol_rsd(k,j))
            rdc = decr * sol_rsd(k,j)
            sol_rsd(k,j) = sol_rsd(k,j) - rdc
            if (sol_rsd(k,j) < 0.) sol_rsd(k,j) = 0.
            rmn1 = decr * sol_fon(k,j)
            sol_fop(k,j) = Dmax1(1.e-6,sol_fop(k,j))
            rmp = decr * sol_fop(k,j)

            sol_fop(k,j) = sol_fop(k,j) - rmp
            sol_fon(k,j) = Dmax1(1.e-6,sol_fon(k,j))
            sol_fon(k,j) = sol_fon(k,j) - rmn1
            sol_no3(k,j) = sol_no3(k,j) + .8 * rmn1
            sol_aorgn(k,j) = sol_aorgn(k,j) + .2 * rmn1
            sol_solp(k,j) = sol_solp(k,j) + .8 * rmp
            sol_orgp(k,j) = sol_orgp(k,j) + .2 * rmp
         end if
!!  compute denitrification
         wdn = 0.
         if (i_sep(j) /= k .or. isep_opt(j) /= 1) then
            if (sut >= sdnco(j)) then
               wdn = sol_no3(k,j) * (1. - Exp(-cdn(j) * cdg * sol_cbn(k,j)))
            else
               wdn = 0.
            endif
            sol_no3(k,j) = sol_no3(k,j) - wdn
         end if

         !! summary calculations
         if (curyr > nyskip) then
            wshd_hmn = wshd_hmn + hmn * hru_dafr(j)
            wshd_rwn = wshd_rwn + rwn * hru_dafr(j)
            wshd_hmp = wshd_hmp + hmp * hru_dafr(j)
            wshd_rmn = wshd_rmn + rmn1 * hru_dafr(j)
            wshd_rmp = wshd_rmp + rmp * hru_dafr(j)
            wshd_dnit = wshd_dnit + wdn * hru_dafr(j)
            hmntl = hmntl + hmn
            rwntl = rwntl + rwn
            hmptl = hmptl + hmp
            rmn2tl = rmn2tl + rmn1
            rmptl = rmptl + rmp
            wdntl = wdntl + wdn
         end if
      end if
   end do


   return
end
