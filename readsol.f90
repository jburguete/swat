!> @file readsol.f90
!> file containing the subroutine readsol
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU/subbasin soil properties file
!> (.sol). This file contains data related to soil physical properties and
!> general chemical properties
!> @param[in] k HRU number
subroutine readsol(k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    k             |none          |HRU number
!!    mlyr          |none          |maximum number of soil layers
!!    idplt(:)      |none          |land cover/crop identification code for
!!                                 |first crop grown in HRU (the only crop if
!!                                 |there is no rotation)
!!    rdmx(:)       |m             |maximum root depth of plant
!!    rsdin(:)      |kg/ha         |initial residue cover
!!    sol_no3(:,:)  |mg N/kg       |concentration of nitrate in soil layer
!!    sol_orgn(1,:) |mg N/kg soil  |organic N concentration in top soil layer
!!    sol_orgp(1,:) |mg P/kg soil  |organic P concentration in top soil layer
!!    sol_solp(1,:) |mg P/kg soil  |soluble P concentration in top soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    anion_excl(:) |none          |fraction of porosity from which anions
!!                                 |are excluded
!!    sol_clay(:,:) |%             |percent clay content in soil material
!!    sol_rock(:,:) |%             |percent of rock fragments in soil layer
!!    sol_silt(:,:) |%             |percent silt content in soil material
!!    snam(:)       |NA            |soil series name
!!    sol_alb(:)    |none          |albedo when soil is moist
!!    sol_awc(:,:)  |mm H20/mm soil|available water capacity of soil layer
!!    sol_bd(:,:)   |Mg/m**3       |bulk density of the soil
!!    sol_cbn(:,:)  |%             |percent organic carbon in soil layer
!!    sol_crk(:)    |none          |crack volume potential of soil
!!    sol_ec(:)     |dS/m          |electrical conductivity of soil layer
!!    sol_k(:,:)    |mm/hr         |saturated hydraulic conductivity of soil
!!                                 |layer
!!    sol_nly(:)    |none          |number of soil layers
!!    sol_no3(:,:)  |mg N/kg       |concentration of nitrate in soil layer
!!    sol_orgn(1,:) |mg N/kg soil  |organic N concentration in top soil layer
!!    sol_orgp(1,:) |mg P/kg soil  |organic P concentration in top soil layer
!!    sol_rsd(:,:)  |kg/ha         |amount of organic matter in the soil layer
!!                                 |classified as residue
!!    sol_sand(:,:) |%             |percent sand content of soil material
!!    sol_solp(1,:) |mg P/kg soil  |soluble P concentration in top soil layer
!!    sol_stap(:,:) |kg P/ha       |amount of phosphorus in the soil layer
!!                                 |stored in the stable mineral phosphorus
!!                                 |pool
!!    sol_z(:,:)    |mm            |depth to bottom of soil layer
!!    sol_zmx(:)    |mm            |maximum rooting depth
!!    usle_k(:)     |none          |USLE equation soil erodibility (K) factor
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    a
!!    b
!!    c
!!    d
!!    dep_new
!!    eof
!!    j           |none          |counter
!!    nly         |none          |number of soil layers
!!    nota
!!    plt_zmx     |mm            |rooting depth of plant
!!    titldum     |NA            |title line/skipped line in .sol file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Exp, Abs
!!    SWAT: estimate_ksat, layersplit

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: k
   real*8, parameter :: a = 50.0, b = 20.0, c = 5.0, d = 2.0
   integer, parameter :: nota = 10
   character (len=80) :: titldum
   real*8 :: dep_new, plt_zmx
   integer :: eof, j, nly

!!    initialize local variables
   nly = 0
   plt_zmx = 0.

   read (107,5500) titldum
   read (107,5100) snam(k)
   read (107,5200) hydgrp(k)
   read (107,5300) sol_zmx(k)
   read (107,5400) anion_excl(k)
   read (107,5600) sol_crk(k)
   read (107,5500) titldum
   read (107,5000) (sol_z(j,k), j = 1, mlyr)


   !! calculate number of soil layers in HRU soil series
   do j = 1, mlyr
!!    khan soils
!      sol_z(j,k) = sol_z(j,k) / 5.0
      if (sol_z(j,k) <= 0.001) sol_nly(k) = j - 1
      if (sol_z(j,k) <= 0.001) exit
   enddo
   if (sol_nly(k) == 0) sol_nly(k) = 10
   nly = sol_nly(k)

   eof = 0
   do
      read (107,5000) (sol_bd(j,k), j = 1, nly)
      read (107,5000) (sol_awc(j,k), j = 1, nly)
      read (107,5000) (sol_k(j,k), j = 1, nly)
      read (107,5000) (sol_cbn(j,k), j = 1, nly)
      read (107,5000) (sol_clay(j,k), j = 1, nly)
      read (107,5000) (sol_silt(j,k), j = 1, nly)
      read (107,5000) (sol_sand(j,k), j = 1, nly)

      read (107,5000) (sol_rock(j,k), j = 1, nly)
      read (107,5000) sol_alb(k)
      read (107,5000) usle_k(k)
!    change double subscripted sol_ec statement 1/27/09 when making septic changes
      read (107,5000,iostat=eof) (sol_ec(j,k), j = 1, nly)
      if (eof < 0) exit
!    change double subscripted sol_ec statement 1/27/09 when making septic changes

      !! MJW added rev 490
      !! PH-H20
      read (107,5000,iostat=eof) (sol_ph(j,k), j = 1, nly)
      if (eof < 0) exit

      !!CaCo3 content (%)
      read (107,5000,iostat=eof) (sol_cal(j,k), j = 1, nly)
      if (eof < 0) exit

      exit
   end do

   !!Armen January 2009
   do j=1, nly
      if (sol_rock(j,k) > 98.0) sol_rock(j,k) = 98.0
      if (sol_awc(j,k) <= .01) sol_awc(j,k) = .01
      if (sol_awc(j,k) >= .80) sol_awc(j,k) = .80
      if (sol_cbn(j,k) < 1.e-6) sol_cbn(j,k) = .10
      sol_n(j,k) = sol_cbn(j,k) / 11.0
   end do
   !!Armen January 2009 end

!!    add 10mm layer at surface of soil
   if (sol_z(1,k) > 10.1) then
      sol_nly(k) = sol_nly(k) + 1
      nly = nly + 1
      do j = nly, 2, -1
         sol_z(j,k) = sol_z(j-1,k)
         sol_bd(j,k) = sol_bd(j-1,k)
         sol_awc(j,k) = sol_awc(j-1,k)
         sol_k(j,k) = sol_k(j-1,k)
         sol_cbn(j,k) = sol_cbn(j-1,k)
         !!Armen January 2009
         sol_n(j,k) = sol_n(j-1,k)
!                 sol_mc(j,k) = sol_mc(j-1,k)
!                 sol_mn(j,k) = sol_mn(j-1,k)
!                 sol_mp(j,k) = sol_mp(j-1,k)
         sol_rock(j,k) = sol_rock(j-1,k) !!! Armen 13 Jan 2008
         sol_clay(j,k) = sol_clay(j-1,k)
         sol_sand(j,k) = sol_sand(j-1,k) !!! Claire 2 Dec 2009
         sol_silt(j,k) = sol_silt(j-1,k) !!! Claire 2 Dec 2009
         sol_ph(j,k) = sol_ph(j-1,k) !! mjw rev 490
         sol_cal(j,k) = sol_cal(j-1,k) !! mjw rev 490
         !!Armen January 2009 end
!    change below double subscripted sol_ec statement 1/27/09 when making septic changes
         sol_ec(j,k) = sol_ec(j-1,k)
!    change below double subscripted sol_ec statement 1/27/09 when making septic changes
         sol_no3(j,k) = sol_no3(j-1,k)
         sol_orgn(j,k) = sol_orgn(j-1,k)
         sol_orgp(j,k) = sol_orgp(j-1,k)
         sol_solp(j,k) = sol_solp(j-1,k)
      end do
      sol_z(1,k) = 10.
   endif

   if (isproj == 2) then
      call estimate_ksat(sol_clay(j,k),sol_k(j,k))  !!  NK June 28, 2006
   endif


!!    compare maximum rooting depth in soil to maximum rooting depth of
!!    plant
   if (sol_zmx(k) <= 0.001) sol_zmx(k) = sol_z(nly,k)
   plt_zmx = 0.
   if (idplt(k) > 0) then
      if (idc(idplt(k)) > 0) then
         plt_zmx = 1000. * rdmx(idplt(k))
      end if
   end if
   if (sol_zmx(k) > 1. .and. plt_zmx > 1.) then
      sol_zmx(k) = Min(sol_zmx(k),plt_zmx)
   else
      !! if one value is missing it will set to the one available
      sol_zmx(k) = Max(sol_zmx(k),plt_zmx)
   end if

!! create a bizone layer in septic HRUs
   if (isep_opt(k) /= 0) then
      if (bz_z(k)+bz_thk(k) > sol_z(nly,k)) then
         if (sol_z(nly,k)>bz_thk(k)+10.) then !min. soil thickness for biozone layer (10mm top+biozone layer thickness)
            bz_z(k) = sol_z(nly,k) - bz_thk(k)
         else
            bz_z(k) = sol_z(nly,k)
            sol_z(nly,k) = sol_z(nly,k) + bz_thk(k)
         endif
      endif
      if (bz_z(k) > 0.) then
         call layersplit (bz_z(k), k)
         dep_new = bz_z(k) + bz_thk(k)
         call layersplit (dep_new, k)
         i_sep(k) = iseptic
      endif
   endif

   nly = sol_nly(k)

!!    set default values/initialize variables
   if (sol_alb(k) < 0.1) sol_alb(k) = 0.1
   if (anion_excl(k) <= 1.e-6) anion_excl(k) = anion_excl_bsn
   if (anion_excl(k) >= 1.) anion_excl(k) = 0.99
   if (rsdin(k) > 0.) sol_rsd(1,k) = rsdin(k)
   do j = 1, nly
      if (sol_k(j,k) <= 0.0) then
         if (hydgrp(k) == "A") then
            sol_k(j,k) = a
         else
            if (hydgrp(k) == "B") then
               sol_k(j,k) = b
            else
               if (hydgrp(k) == "C") then
                  sol_k(j,k) = c
               else
                  if (hydgrp(k) == "D") then
!            sol_k(j,k) = c
                     sol_k(j,k) = d          !Claire 12/2/09
                  else
                     sol_k(j,k) = nota
                  endif
               endif
            endif
         endif
      endif
      if (sol_bd(j,k) <= 1.e-6) sol_bd(j,k) = 1.3
      if (sol_bd(j,k) > 2.) sol_bd(j,k) = 2.0
      if (sol_awc(j,k) <= 0.) sol_awc(j,k) = .005
      !! Defaults for ph and calcium mjw average of 20,000 SSURGO soils mjw rev 490
      if (sol_cal(j,k)<= 1.e-6) sol_cal(j,k) = 2.8
      if (sol_ph(j,k)<= 1.e-6) sol_ph(j,k) = 6.5
   end do


   close (107)
   return
5000 format (27x,15f12.2)
5100 format (12x,a16)
5200 format (24x,a1)
5300 format (28x,f12.2)
5400 format (51x,f5.3)
5500 format (a80)
5600 format (33x,f5.3)
end
