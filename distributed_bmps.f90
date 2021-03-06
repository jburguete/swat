!> @file distributed_bmps.f90
!> file containing the subroutine distributed_bmps
!> @author
!> modified by Javier Burguete

!> this subroutine calls routines for urban BMPs in the subbasin
!> param[in] sb subbasin or reach number
subroutine distributed_bmps(sb)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff generated each hour of day
!!                               |in HRU
!!    hru_ha(:)   |ha            |area of HRU in hectares
!!    ubntss(:)   |metric tons   |TSS loading from urban impervious cover
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ubntss(:)   |metric tons   |TSS loading from urban impervious cover
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bmpfr
!!    ftqm3
!!    ftsed
!!    ii
!!    kk
!!    ri_totalflw
!!    ri_totaltss
!!    riflw
!!    riqm3
!!    rised
!!    sb
!!    sf_totalflw
!!    sf_totaltss
!!    sfflw
!!    sfsed
!!    spqm3
!!    spsed
!!    sub_ha
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Sum
!!    SWAT: bmp_sand_filter, bmp_sed_pond, bmp_ri_pond

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   !dimensions: 1=inflow/outflow, 2=pond id, 3=time step
   integer, intent(in) :: sb
   real*8, dimension(4,0:nstep) :: ftqm3, ftsed, ri_totalflw, ri_totaltss,&
      &riflw, riqm3, rised, sf_totalflw, sf_totaltss, sfflw, sfsed, spqm3, spsed
   real*8 :: bmpfr, sub_ha
   integer :: ii, kk

   sub_ha = da_ha * sub_fr(sb)
   sf_totalflw = 0.; sf_totaltss = 0.
   ri_totalflw = 0.; ri_totaltss = 0.
   sfflw = 0.; sfsed = 0.; riflw = 0.; rised = 0.;bmpfr=0.
   spqm3 = 0.; spsed=0.; ftqm3=0.; ftsed=0.; riqm3=0.

   !initialize daily recharge from distributed BMPs
   bmp_recharge(sb) = 0.

   !---------------------------------
   ! sedimentation-filtration basin
   if(num_sf(sb)>=1.and.hrnopcp(sb,nstep)<96) then

      do kk=1,num_sf(sb)
         !fraction urban runoff to the sed-fil
         sfflw(1,1:nstep) = sub_ubnrunoff(sb,1:nstep) * sf_fr(kk,sb) !mm
         sfsed(1,1:nstep) = sub_ubntss(sb,1:nstep) * sf_fr(kk,sb) !tons
         !total inflow to sedfils
         sf_totalflw(1,:) = sf_totalflw(1,:) + sfflw(1,:) !mm
         sf_totaltss(1,:) = sf_totaltss(1,:) + sfsed(1,:) !tons

         if (iyr>sf_iy(kk,sb) .or.&
            &(iyr==sf_iy(kk,sb) .and. i_mo>=sf_im(kk,sb))) then
            if (sf_typ(kk,sb)==2) then !partial scale
               call bmp_sand_filter(sb, kk, sfflw, sfsed)
               spqm3(:,:) = 0.
               spsed(:,:) = 0.
               ftqm3(:,:) = sfflw(:,:) * ((sub_ha - ft_sa(kk,sb)&  !m3
                  &/ 10000.) *10.)
               ftsed(:,:) = sfsed(:,:) !tons

               !total (aggregated) outflow from sedfils
               sf_totalflw(2,:) = sf_totalflw(2,:) + sfflw(2,:)&
                  &+ sfflw(3,:) !mm
               sf_totaltss(2,:) = sf_totaltss(2,:) + ftsed(2,:)&
                  &+ ftsed(3,:) !tons
            else if (sf_typ(kk,sb)==1) then !full scale
               !first route through sedimentation pond
               call bmp_sed_pond(sb, kk, sfflw, sfsed)

               spqm3(:,:) = sfflw(:,:) * ((sub_ha - sp_sa(kk,sb)&
                  &/ 10000.) *10.)
               spsed(:,:) = sfsed(:,:)

               ! add bypass flow to total (aggregated) outflow
               sf_totalflw(2,:) = sf_totalflw(2,:) + sfflw(3,:)
               sf_totaltss(2,:) = sf_totaltss(2,:) + sfsed(3,:)

               ! outflow from the sedimentation basin goest to sand filter
               sfflw(1,:) = sfflw(2,:)
               sfsed(1,:) = sfsed(2,:)

               ! then the outflow from sed pond goes to sand filter
               call bmp_sand_filter(sb, kk, sfflw, sfsed)

               ftqm3(:,:) = sfflw(:,:) *  ((sub_ha - ft_sa(kk,sb)&
                  &/ 10000.) *10.) !m3
               ftsed(:,:) = sfsed(:,:)  !tons

               !total (aggregated) outflow from sedfils
               sf_totalflw(2,:) = sf_totalflw(2,:) + sfflw(3,:)&
                  &+ sfflw(2,:) !mm
               sf_totaltss(2,:) = sf_totaltss(2,:) + sfsed(3,:)&
                  &+ sfsed(2,:) !tons

            else !sedimentation pond only
               call bmp_sed_pond(sb, kk, sfflw, sfsed)

               ftqm3(:,:) = 0.
               ftsed(:,:)=0.
               spqm3(:,:) = sfflw(:,:) * ((sub_ha - sp_sa(kk,sb)&
                  &/ 10000.) *10.)
               spsed(:,:) = sfsed(:,:)
               !total (aggregated) outflow from sedfils
               sf_totalflw(2,:) = sf_totalflw(2,:) + sfflw(3,:)&
                  &+ sfflw(2,:)
               sf_totaltss(2,:) = sf_totaltss(2,:) + sfsed(3,:)&
                  &+ sfsed(2,:)
            endif


         else
            ! skip bmp simulation before it's constructed
            !total (aggregated) outflow from sedfils
            sf_totalflw(2,:) = sf_totalflw(2,:) + sfflw(1,:)
            sf_totaltss(2,:) = sf_totaltss(2,:) + sfsed(1,:)
         endif

         !print out bmp result in bmp-sedfil.out
         write(77778,'(5i6,30f12.3)') iyr,iida,ii,sb,kk,Sum(spqm3(1,:)),&
            &Sum(spqm3(2,:)),Sum(spqm3(3,:)),Sum(spqm3(4,:)),&
            &Sum(spsed(1,:))*1000.,Sum(spsed(2,:))*1000.,Sum(spsed(3,:))*1000.,&
            &Sum(ftqm3(1,:)),Sum(ftqm3(2,:)),Sum(ftqm3(3,:)),Sum(ftqm3(4,:)),&
            &Sum(ftsed(1,:))*1000.,Sum(ftsed(2,:))*1000.,Sum(ftsed(3,:))*1000.

      end do
   endif

   !---------------------------------
   ! retention-irrigation (RI)
   if (num_ri(sb)>=1 .and. hrnopcp(sb,nstep)<96) then !72 hours draw-down plus another day (24hrs)

      do kk=1,num_ri(sb)

         ! skip the pond that has zero inflow
         if (ri_fr(sb,kk)==0) cycle

         !fraction urban runoff to the RI
         riflw(1,1:nstep) = sub_ubnrunoff(sb,1:nstep) *&
            &ri_fr(sb,kk)
         rised(1,1:nstep) = sub_ubntss(sb,1:nstep) * ri_fr(sb,kk)

         ! total inflow to RIs
         ri_totalflw(1,1:nstep) = ri_totalflw(1,1:nstep)&
            &+ riflw(1,1:nstep)
         ri_totaltss(1,1:nstep) = ri_totaltss(1,1:nstep)&
            &+ rised(1,1:nstep)

         ! skip bmp simulation before it's constructed
         if (iyr>ri_iy(sb,kk) .or.&
            &(iyr==ri_iy(sb,kk) .and. i_mo>=ri_im(sb,kk))) then

            call bmp_ri_pond(sb, kk, riflw, rised)

            riqm3(:,:) = riflw(:,:)* ((sub_ha - sp_sa(kk,sb)&
               &/ 10000.) *10.)
         else
            riflw(2,:) = riflw(1,:)
            rised(2,:) = rised(1,:)
         end if

         do ii=1,nstep
            !aggregate individual pond effluents in subbasin
            ri_totalflw(2,ii) = ri_totalflw(2,ii) + riflw(2,ii)
            ri_totaltss(2,ii) = ri_totaltss(2,ii) + rised(2,ii)
         end do
         write(77779,'(5i6,30f12.3)') iyr,iida,ii,sb,kk,Sum(riqm3(1,:)),&
            &Sum(riqm3(2,:)),Sum(riqm3(3,:)),Sum(riqm3(4,:)),&
            &Sum(rised(1,:))*1000.,Sum(rised(2,:))*1000.,Sum(rised(3,:))*1000.

      end do

   endif

   ! allocate bmp inflow/outflow to subbasin surface runoff volume
   sub_ubnrunoff(sb,1:nstep) = sub_ubnrunoff(sb,1:nstep) -&
      &sf_totalflw(1,1:nstep) - ri_totalflw(1,1:nstep) +&
      &sf_totalflw(2,1:nstep) + ri_totalflw(2,1:nstep)
   sub_ubntss(sb,1:nstep) = sub_ubntss(sb,1:nstep) -&
      &sf_totaltss(1,1:nstep) - ri_totaltss(1,1:nstep) +&
      &sf_totaltss(2,1:nstep) + ri_totaltss(2,1:nstep)

   return
end
