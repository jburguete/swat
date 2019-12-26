!> @file soil_phys.f90
!> file containing the subroutine soil_phys
!> @author
!> modified by Javier Burguete

!> this subroutine initializes soil physical properties
!> @param[in] ii HRU number
subroutine soil_phys(ii)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cn2(:)        |none          |SCS runoff curve number for moisture
!!                                 |condition II
!!    ddrain(:)     |mm            |depth to the sub-surface drain
!!    ffc(:)        |none          |initial HRU soil water content
!!                                 |expressed as fraction of field capacity
!!    hru_dafr(:)   |km2/km2       |fraction of total watershed area contained
!!                                 |in HRU
!!    ievent        |none          |rainfall/runoff code
!!                                 |0 daily rainfall/curve number technique
!!                                 |1 sub-daily rainfall/Green&Ampt/hourly
!!                                 |  routing
!!                                 |3 sub-daily rainfall/Green&Ampt/hourly routing
!!    sol_awc(:,:)  |mm H20/mm soil|available water capacity of soil layer
!!    sol_bd(:,:)   |Mg/m**3       |bulk density of the soil
!!    sol_clay(:,:) |%             |percent clay content in soil material
!!    sol_crk(:)    |none          |crack volume potential of soil
!!    sno_hru(:)    |mm H2O        |amount of water stored as snow
!!    sol_k(:,:)    |mm/hr         |saturated hydraulic conductivity of soil
!!                                 |layer
!!    sol_nly(:)    |none          |number of soil layers
!!    sol_rock(:)   |%             |percent of rock fragments in soil layer
!!    sol_silt(:,:) |%             |percent silt content in soil material
!!    sol_z(:,:)    |mm            |depth to bottom of soil layer
!!    usle_k(:)     |none          |USLE equation soil erodibility (K) factor
!!    usle_ls(:)    |none          |USLE equation length slope (LS) factor
!!    usle_p(:)     |none          |USLE equation support practice (P) factor
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    crdep(:,:)    |mm            |maximum or potential crack volume
!!    ldrain(:)     |none          |soil layer where drainage tile is located
!!    sol_avbd(:)   |Mg/m^3        |average bulk density for soil profile
!!    sol_avpor(:)  |none          |average porosity for entire soil profile
!!    sol_fc(:,:)   |mm H2O        |amount of water available to plants in soil
!!                                 |layer at field capacity (fc - wp)
!!    sol_hk(:,:)   |none          |beta coefficent to calculate hydraulic
!!                                 |conductivity
!!    sol_por(:,:)  |none          |total porosity of soil layer expressed as a
!!                                 |fraction of the total volume
!!    sol_rock(:)   |none          |exponential value that is a function of
!!                                 |percent rock
!!    sol_sand(:,:) |%             |percent sand content in soil material
!!    sol_st(:,:)   |mm H2O        |amount of water stored in the soil layer
!!                                 |on any given day (less wp water)
!!    sol_sumfc(:)  |mm H2O        |amount of water held in soil profile at
!!                                 |field capacity
!!    sol_sumul(:)  |mm H2O        |amount of water held in soil profile at
!!                                 |saturation
!!    sol_sw(:)     |mm H2O        |amount of water stored in soil profile on
!!                                 |any given day
!!    sol_ul(:,:)   |mm H2O        |amount of water held in the soil layer at
!!                                 |saturation (sat - wp water)
!!    sol_up(:,:)   |mm H2O/mm soil|water content of soil at -0.033 MPa (field
!!                                 |capacity)
!!    sol_wp(:,:)   |mm H20/mm soil|water content of soil at -1.5 MPa (wilting
!!                                 |point)
!!    sol_wpmm(:,:) |mm H20        |water content of soil at -1.5 MPa (wilting
!!                                 |point)
!!    usle_mult(:)  |none          |product of USLE K,P,LS,exp(rock)
!!    volcr(:,:)    |mm            |crack volume for soil layer
!!    wfsh(:)       |mm            |wetting front matric potential
!!    wshd_snob     |mm H20        |average amount of water stored in snow
!!                                 |at the beginning of the simulation for the
!!                                 |entire watershed
!!    wshd_sw       |mm H2O        |average amount of water stored in soil
!!                                 |for the entire watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cl
!!    dg          |mm            |depth of layer
!!    drpor
!!    j           |none          |counter
!!    nly         |none          |number of soil layers
!!    pormm       |mm            |porosity in mm depth
!!    sa
!!    si
!!    sumpor      |mm            |porosity of profile
!!    xx          |none          |variable to hold value
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Sqrt
!!    SWAT: curno

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: ii
   real*8 :: cl, dg, drpor, pormm, sa, si, sumpor, xx
   integer :: j, nly

   nly = sol_nly(ii)

!!    calculate composite usle value
   sol_rock(1,ii) = Exp(-.053 * sol_rock(1,ii))
   usle_mult(ii) = sol_rock(1,ii) * usle_k(ii) * usle_p(ii) * usle_ls(ii) * 11.8


!!    calculate water content of soil at -1.5 MPa and -0.033 MPa
   do j = 1, nly
      sol_wp(j,ii) = 0.4 * sol_clay(j,ii) * sol_bd(j,ii) / 100.
      if (sol_wp(j,ii) <= 0.) sol_wp(j,ii) = .005
      sol_up(j,ii) = sol_wp(j,ii) + sol_awc(j,ii)
      sol_por(j,ii) = 1. - sol_bd(j,ii) / 2.65
      if (sol_up(j,ii) >= sol_por(j,ii)) then
         sol_up(j,ii) = sol_por(j,ii) - .05
         sol_wp(j,ii) = sol_up(j,ii) - sol_awc(j,ii)
         if (sol_wp(j,ii) <= 0.) then
            sol_up(j,ii) = sol_por(j,ii) * .75
            sol_wp(j,ii) = sol_por(j,ii) * .25
         end if
      end if
      !! compute drainable porosity and variable water table factor - Daniel
      drpor = sol_por(j,ii) - sol_up(j,ii)
      vwt(j,ii)= ((437.13 * drpor**2)-(95.08 * drpor)+8.257)
   end do

   sa = sol_sand(1,ii) / 100.
   cl = sol_clay(1,ii) / 100.
   si = sol_silt(1,ii) / 100.
!!    determine detached sediment size distribution
!!    typical for mid-western soils in USA (Foster et al., 1980)
!!    Based on SWRRB
   det_san(ii) = 2.49 * sa * (1. - cl)   !! Sand fraction
   det_sil(ii) = 0.13 * si               !! Silt fraction
   det_cla(ii) = 0.20 * cl               !! Clay fraction
   if (cl < .25) then
      det_sag(ii) = 2.0 * cl              !! Small aggregate fraction
   else if (cl > .5) then
      det_sag(ii) = .57
   else
      det_sag(ii) = .28 * (cl - .25) + .5
   end if

   det_lag(ii) = 1. - det_san(ii) - det_sil(ii) - det_cla(ii) - det_sag(ii) !! Large Aggregate fraction

!! Error check. May happen for soils with more sand
!!    Soil not typical of mid-western USA
!!    The fraction wont add upto 1.0
   if (det_lag(ii) < 0.) then
      xx = 1. - det_lag(ii)
      det_san(ii) = det_san(ii) * xx
      det_sil(ii) = det_sil(ii) * xx
      det_cla(ii) = det_cla(ii) * xx
      det_sag(ii) = det_sag(ii) * xx
      det_lag(ii) = 0.
   end if


!!    initialize water/drainage coefs for each soil layer
   xx = 0.
   sumpor = 0.
   do j = 1, nly
      dg = sol_z(j,ii) - xx
      pormm = sol_por(j,ii) * dg
      sumpor = sumpor + pormm
      sol_ul(j,ii) = (sol_por(j,ii) - sol_wp(j,ii)) * dg
      sol_sumul(ii) = sol_sumul(ii) + sol_ul(j,ii)
      sol_fc(j,ii) = dg * (sol_up(j,ii) - sol_wp(j,ii))
      sol_sumfc(ii) = sol_sumfc(ii) + sol_fc(j,ii)
      sol_st(j,ii) = sol_fc(j,ii) * ffc(ii)
      sol_hk(j,ii) = (sol_ul(j,ii) - sol_fc(j,ii)) / sol_k(j,ii)
      if (sol_hk(j,ii) < 1.) sol_hk(j,ii) = 1.
      sol_sw(ii) = sol_sw(ii) + sol_st(j,ii)
      sol_wpmm(j,ii) = sol_wp(j,ii) * dg
      sol_sumwp(ii) = sol_sumwp(ii) + sol_wpmm(j,ii)
      crdep(j,ii) = sol_crk(ii) * 0.916 * Exp(-.0012 * sol_z(j,ii)) * dg
      volcr(j,ii) = crdep(j,ii) * (sol_fc(j,ii) - sol_st(j,ii)) / (sol_fc(j,ii))
      xx = sol_z(j,ii)
   end do
   !! initialize water table depth and soil water for Daniel
!      sol_swpwt(ii) = sol_sw(ii)
!      if (ffc(ii) > 1.) then
!        wat_tbl(ii) = (sol_sumul(ii) - ffc(ii) * sol_sumfc(ii)) /
!     &                                                      sol_z(nly,ii)
!      else
!        wat_tbl(ii) = 0.
!      end if
   !!Initializing water table depth and soil water revised by D. Moriasi 4/8/2014
   do j = 1, nly
      sol_stpwt(j,ii) = sol_st(j,ii)
   end do
   sol_swpwt(ii) = sol_sw(ii)
   wat_tbl(ii) = dep_imp(ii)- (shallst(ii)/sol_por(nly,ii))

   !!Initializing water table depth and soil water revised by D. Moriasi 4/8/2014
   !! initialize water table depth and soil water for Daniel
   sol_avpor(ii) = sumpor / sol_z(nly,ii)
   sol_avbd(ii) = 2.65 * (1. - sol_avpor(ii))


!!    define soil layer that the drainage tile is in
   if (ddrain(ii) > 0.) then
      do j = 1, nly
         if (ddrain(ii) < sol_z(j,ii)) ldrain(ii) = j
         if (ddrain(ii) < sol_z(j,ii)) exit
      end do
   else
      ldrain(ii) = 0
   endif

!!    calculate infiltration parameters for subdaily time step
   if (ievent > 0) then
      sol_sand = 0.
      sol_sand(1,ii) = 100. - sol_clay(1,ii) - sol_silt(1,ii)
      wfsh(ii) = 10. * Exp(6.5309 - 7.32561 * sol_por(1,ii) +&
         &3.809479 * sol_por(1,ii) ** 2 + 0.001583 * sol_clay(1,ii) ** 2 +&
         &0.000344 * sol_sand(1,ii) * sol_clay(1,ii) - 0.049837 *&
         &sol_por(1,ii) * sol_sand(1,ii)&
         &+ 0.001608 * sol_por(1,ii) ** 2 * sol_sand(1,ii) ** 2 +&
         &0.001602 * sol_por(1,ii) ** 2 * sol_clay(1,ii) ** 2 -&
         &0.0000136 * sol_sand(1,ii) ** 2 * sol_clay(1,ii) -&
         &0.003479 * sol_clay(1,ii) ** 2 * sol_por(1,ii) -&
         &0.000799 * sol_sand(1,ii) ** 2 * sol_por(1,ii))
   end if


!!    initialize watershed water parameters
   wshd_sw = wshd_sw + sol_sw(ii) * hru_dafr(ii)
   wshd_snob = wshd_snob + sno_hru(ii) * hru_dafr(ii)


   call curno(cn2(ii),ii) !! J.Jeong 4/18/2008
!      call curno_subd(cn2(ii),ii)  !! changed for URBAN

   return
end
