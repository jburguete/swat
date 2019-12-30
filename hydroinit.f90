!> @file hydroinit.f90
!> file containing the subroutine hydroinit
!> @author
!> modified by Javier Burguete

!> This subroutine computes variables related to the watershed hydrology:
!> the time of concentration for the subbasins, lagged surface runoff,
!> the coefficient for the peak runoff rate equation, and lateral flow travel
!> time.
subroutine hydroinit

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~1
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_l(1,:)    |km            |longest tributary channel length in subbasin
!!    ch_l(2,:)    |km            |main channel length in subbasin
!!    ch_n(1,:)    |none          |Manning's "n" value for the tributary channels
!!    ch_s(1,:)    |m/m           |average slope of tributary channels
!!    da_km       |km2           |area of the watershed in square kilometers
!!    gdrain(:)   |hrs           |drain tile lag time: the amount of time
!!                               |between the transfer of water from the soil
!!                               |to the drain tile and the release of the
!!                               |water from the drain tile to the reach.
!!    hru_dafr(:) |km2/km2       |fraction of total watershed area contained
!!                               |in HRU
!!    hru_km(:)   |km2           |area of HRU in square kilometers
!!    hru_slp(:)  |m/m           |average slope steepness
!!    hru_sub(:)  |none          |subbasin in which HRU is located
!!    lat_ttime(:)|days          |lateral flow travel time
!!   tile_ttime(:)|none          |Exponential of the tile flow travel time
!!    ldrain(:)   |none          |soil layer where drainage tile is located
!!    nhru        |none          |number of HRUs in watershed
!!    ov_n(:)     |none          |Manning's "n" value for overland flow
!!    slsoil(:)   |m             |slope length for lateral subsurface flow
!!    slsubbsn(:) |m             |average slope length for subbasin
!!    sol_k(:,:)  |mm/hr         |saturated hydraulic conductivity of soil
!!                               |layer
!!    sol_nly(:)  |none          |number of layers in soil profile
!!    sub_fr(:)   |none          |fraction of total watershed area contained in
!!                               |subbasin
!!    surlag      |days          |Surface runoff lag time.
!!                               |This parameter is needed in subbasins where
!!                               |the time of concentration is greater than 1
!!                               |day. SURLAG is used to create a "storage" for
!!                               |surface runoff to allow the runoff to take
!!                               |longer than 1 day to reach the subbasin outlet
!!    tconc(:)    |hr            |time of concentration
!!    uhalpha     |              |alpha coefficient for estimating unit hydrograph
!!                               |using a gamma function (*.bsn)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    brt(:)      |none          |fraction of surface runoff within the subbasin
!!                               |which takes 1 day or less to reach the
!!                               |subbasin outlet
!!    lat_ttime(:)|none          |Exponential of the lateral flow travel time
!!   tile_ttime(:)|none          |Exponential of the tile flow travel time
!!    sub_tc(:)   |hr            |time of concentration for subbasin
!!    t_ov(:)     |hr            |time for flow from farthest point in subbasin
!!                               |to enter a channel
!!    tconc(:)    |hr            |time of concentration for HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    a
!!    b
!!    c
!!    isb         |none          |counter
!!    j           |none          |counter
!!    l           |none          |counter
!!    q
!!    ql
!!    scmx        |mm/hr         |maximum soil hydraulic conductivity
!!    sumq
!!    tb
!!    tp
!!    xx          |none          |variable to hold calculation result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Dmin1, Sqrt, Exp, Ceiling, Dfloat, Int, Float, Max
!!    SWAT: ttcoef

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   real*8 :: a, b, c, q, ql, scmx, sumq, tb, tp, xx
   integer :: isb, j, l

   do j = 1, nhru

!! subbasin !!
!!    compute time of concentration (sum of overland and channel times)
      t_ov(j) = .0556 * (slsubbsn(j)*ov_n(j)) ** .6 / hru_slp(j) ** .3
      xx = .62 * ch_l(1,j) * ch_n(1,hru_sub(j)) ** .75 /&
         &((da_km * sub_fr(hru_sub(j)))**.125 *&
         &ch_s(1,hru_sub(j))**.375)
      sub_tc(hru_sub(j)) = t_ov(j) + xx
!! end subbasin !!


!! HRU !!
!!    compute time of concentration (sum of overland and channel times)
      ch_l(1,j) = ch_l(1,j) * hru_dafr(j) / sub_fr(hru_sub(j))
      t_ov(j) = .0556 * (slsubbsn(j)*ov_n(j)) ** .6 / hru_slp(j) ** .3
      xx = .62 * ch_l(1,j) * ch_n(1,hru_sub(j)) ** .75 /&
         &((da_km*hru_dafr(j))**.125*ch_s(1,hru_sub(j))**.375)
      tconc(j) = t_ov(j) + xx

!!    compute delivery ratio
      xx = tconc(j) / sub_tc(hru_sub(j))
      dr_sub(j) = Dmin1(.95, Sqrt(xx))


!!    compute fraction of surface runoff that is reaching the main channel
      if (ievent > 0) then
         brt(j) = 1. - Exp(-surlag(j) / (tconc(j) / (idt / 60.))) !! urban modeling by J.Jeong
      else
         brt(j) = 1. - Exp(-surlag(j) / tconc(j))
      endif
      if (isproj == 2) brt(j) = 1.


!!    compute lateral flow travel time
      if (lat_ttime(j) <= 0.) then
         scmx = 0.
         do l = 1, sol_nly(j)
            if (sol_k(l,j) > scmx) then
               scmx = sol_k(l,j)
            endif
         end do
         !! unit conversion:
         !! xx = m/(mm/h) * 1000.(mm/m)/24.(h/d) / 4.
         xx = 10.4 * slsoil(j) / scmx
         lat_ttime(j) = 1. - Exp(-1./xx)
      else
         lat_ttime(j) = 1. - Exp(-1./lat_ttime(j))
      end if

      if (ldrain(j) > 0 .and. gdrain(j) > 0.01) then
         tile_ttime(j) = 1. - Exp(-24. / gdrain(j))
      else
         tile_ttime(j) = 0.
      end if

!!    compute routing coefficients for main channel
      if (ch_l(2,hru_sub(j)) > 0.) call ttcoef(hru_sub(j))
      if (j == hru1(hru_sub(j))) then
         if (alpha_bnk(hru_sub(j)) <= 0.) then
            alpha_bnk(hru_sub(j)) = alpha_bf(j)
         end if
         alpha_bnke(hru_sub(j)) = Exp(-alpha_bnk(hru_sub(j)))
      end if

   end do

   if (ievent > 0) then
!!    compute unit hydrograph for computing subbasin hydrograph from direct runoff

      do isb = 1, msub
         ql = 0.
         sumq = 0.

         tb = .5 + .6 * sub_tc(isb) + tb_adj  !baseflow time, hr

         if (tb > 96.) tb = 96.      !maximum 48hrs
         tp = .375 * tb      ! time to peak flow

         !! convert to time step (from hr), J.Jeong March 2009
         tb = Ceiling(tb * 60./ Dfloat(idt))
         tp = Int(tp * 60./ Dfloat(idt))

         if(tp==0) tp = 1
         if(tb==tp) tb = tb + 1
         itb(isb) = Int(tb)

         ! Triangular Unit Hydrograph
         if (iuh==1) then
            do l = 1, itb(isb)
               xx = Float(l)
               if (xx < tp) then           !! rising limb of hydrograph
                  q = xx / tp
               else                        !! falling limb of hydrograph
                  q = (tb - xx) / (tb - tp)
               end if
               q = Max(0.,q)
               uh(isb,l) = (q + ql) / 2.
               ql = q
               sumq = sumq + uh(isb,l)
            end do

            do l = 1, itb(isb)
               uh(isb,l) = uh(isb,l) / sumq
            end do

            ! Gamma Function Unit Hydrograph
         elseif (iuh==2) then
            l = 1; q=1.
            do while (q>0.0001)
               xx = Float(l) / tp
               q = xx ** uhalpha * Exp((1. - xx) * uhalpha)
               q = Max(0.,q)
               uh(isb,l) = (q + ql) / 2.
               ql = q
               sumq = sumq + uh(isb,l)
               l = l + 1
               if (l>3*nstep) exit
            end do
            itb(isb) = l - 1
            do l = 1, itb(isb)
               uh(isb,l) = uh(isb,l) / sumq
            end do
         endif

         !compute the number of hydrograph points for flood routing Jaehak 2017
         b = itb(isb)
         c = phi(13,isb) / dthy
         do j = 1, nhru
            a = sub_tc(hru_sub(j)) / dthy !# of timesteps
            NHY(isb) = Max(4 * nstep, Ceiling(a), Ceiling(b), Ceiling(c),&
               &NHY(isb))
         end do
         RCSS(isb) = .5 * (ch_w(2,isb) - phi(6,isb)) / ch_d(isb)
         RCHX(isb) = Sqrt(ch_s(2,isb)) / ch_n(2,isb)
         CHXA(isb) = phi(7,isb) * (phi(6,isb) + phi(7,isb) * RCSS(isb))
         CHXP(isb) = phi(6,isb) + 2. * phi(7,isb) *&
            &Sqrt(RCSS(isb) * RCSS(isb) + 1.)
         QCAP(isb) = CHXA(isb) * RCHX(isb) * (CHXA(isb) / CHXP(isb))**.66667 !bankfull flow, m3/s
      end do
   end if

   return
end
