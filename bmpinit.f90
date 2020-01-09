!> @file bmpinit.f90
!> file containing the subroutine bmpinit
!> @author
!> modified by Javier Burguete

!> this subroutine sets default values for urban bmp parameters
!> @param[in] ii subbasin number
subroutine bmpinit(ii)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~ ~ ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii             |none          |subbasin number
!!    dtp_addon(:,:)  |m            |The distance between spillway levels
!!    dtp_cdis(:,:)  |none          |Discharge coeffieicne for weir/orifice flow
!!    dtp_coef(1,:)   |none          |Coefficient of 3rd degree in the polynomial equation
!!    dtp_coef(2,:)   |none          |Coefficient of 2nd degree in the polynomial equation
!!    dtp_coef(3,:)   |none          |Coefficient of 1st degree in the polynomial equation
!!    dtp_depweir(:,:)|m            |Depth of rectangular weir at different stages
!!    dtp_diaweir(:,:)|m            |Diameter of orifice hole at different stages
!!    dtp_evrsv      |none          |detention pond evaporation coefficient
!!    dtp_expont(:)  |none          |Exponent used in the exponential equation
!!    dtp_flowrate(:,:)|m3/sec      |Maximum discharge from each stage of the weir/hole
!!    dtp_imo(:)     |none          |month the reservoir becomes operational
!!    dtp_intcept(:) |none          |Intercept used in regression equations
!!    dtp_iyr(:)     |none          |year of the simulation that the reservoir
!!                                  |becomes operational
!!    dtp_lwratio(:) |none          |Ratio of length to width of water back up
!!    dtp_onoff(:)   |none          |sub-basin detention pond is associated with
!!    dtp_numweir(:) |none          |Total number of weirs in the BMP
!!    dtp_numstage(:)|none          |Total number of stages in the weir
!!    dtp_pcpret(:,:)|mm            |precipitation for different return periods (not used)
!!    dtp_reltype(:) |none          |Equations for Stage-Discharge relationship,1=exponential
!!                                  |function, 2=linear, 3=logarithmic, 4=cubic, 5=power
!!    dtp_retperd(:,:)|years        |Return period at different stages
!!    dtp_stagdis(:) |none          |0=use weir/orifice discharge equation to calculate
!!                                  |outflow, 1=use stage-dicharge relationship
!!    dtp_totwrwid(:)|m             |Total constructed width of the detention wall across
!!                                  |the creek
!!    dtp_wdratio(:,:)|none         |Width depth ratio of rectangular weirs
!!    dtp_weirdim(:,:)|none         |Weir dimensions, 1=read user input, 0=use model calculation
!!    dtp_weirtype(:,:)|none        |Type of weir: 1=rectangular and 2=circular
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    res_out(:,:,:) |m**3/day      |measured average daily outflow from the
!!                                  |reservoir for the month (needed if IRESCO=1)
!!                                  |(read in as m**3/s and converted to m**3/day)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bmpfr_ri
!!    bmpfr_sf
!!    hstg
!!    hwq
!!    j           |none          |counter
!!    k           |none          |counter
!!    qstg
!!    ri_ndt
!!    sub_ha
!!    wqv
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Sum, Max, Sqrt

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none
   integer, intent(in) :: ii
   integer :: j, k
   real*8 :: bmpfr_ri, bmpfr_sf, hstg, hwq, qstg, ri_ndt, sub_ha, wqv

   bmpfr_sf=0.; bmpfr_ri=0.; hstg = 0.; qstg=0.
   sub_ha = sub_km(ii) * 100.

   !! Detention pond
   !!----------------
   if (dtp_onoff(ii)==1) then
      if (dtp_imo(ii)<=0)      dtp_imo(ii) = 1
      if (dtp_iyr(ii)<=1000)   dtp_iyr(ii) = iyr
      if (dtp_evrsv(ii)<=0)    dtp_evrsv(ii) = 0.1
      if (dtp_lwratio(ii)>1)   dtp_lwratio(ii) = 1
      if (dtp_numweir(ii)<=0)  dtp_numweir(ii) = 1
      if (dtp_numstage(ii)<=0) dtp_numstage(ii) = 1
      if (dtp_numstage(ii)>1) then
         do k=2,dtp_numstage(ii)
            if (dtp_weirtype(k,ii)==1) dtp_addon(k,ii) = 0.
         end do
      endif

      !! Estimating design flow rate if not entered by user
      do k=1,dtp_numstage(ii)
         if (dtp_flowrate(k,ii)<=0.0) then
            dtp_flowrate(k,ii) = 0.5 * 1000.0 * dtp_pcpret(k,ii)&
               &* subdr_km(ii) / (idt*60.)
         end if

         if (dtp_flowrate(k,ii)<0) then
            print *,"Error.. Could not estimate emergency spillway volume"
            print *,"Please enter necessary data in *.pnd input file"
            print *,"for subbasin : ",ii
            !!            stop
         end if
      end do

      !! Separate cumulative flow information to individual weir
      do k=2,dtp_numstage(ii)
         dtp_flowrate(k,ii) = dtp_flowrate(k,ii) / dtp_numweir(ii)
      end do

      !!Estimate weir dimensions based on existing data
      do k=1,dtp_numstage(ii)
         if (dtp_weirdim(k,ii)==0) then  !! Estimating weir dimensions
            if (dtp_weirtype(k,ii)==2) then !! choosing weir type
               dtp_diaweir(k,ii)=(0.479081 * dtp_flowrate(k,ii)&
                  &/ dtp_cdis(k,ii))**0.4
            else !rectangular weir
               if (k==1) then
                  hstg = dtp_depweir(k,ii)
                  dtp_wrwid(k,ii) = dtp_flowrate(k,ii)&
                     &/ (1.84*dtp_cdis(k,ii)*hstg**1.5)
               else
                  qstg = 0.
                  do j=1,k-1
                     hstg = Sum(dtp_depweir(j:k,ii))
                     qstg = dtp_cdis(k,ii) * 1.84&
                        &* dtp_wrwid(ii,j) * hstg ** 1.5 !m3/s
                     dtp_flowrate(k,ii) = Max(0.,dtp_flowrate(k,ii)-qstg)

                  end do
                  dtp_wrwid(k,ii) = dtp_flowrate(k,ii)&
                     &/ (1.84*dtp_cdis(k,ii)*dtp_depweir(k,ii)**1.5)
               endif
            end if
         else  !! read user-entered data
            if (dtp_weirtype(k,ii)==1) then
               dtp_wrwid(k,ii) = dtp_wdratio(k,ii) * dtp_depweir(k,ii)
            end if
         end if
      end do
      !! divide rectangular weirs into multiple single stage weirs
      do k = 2, dtp_numstage(ii)
         dtp_addon(k,ii) = dtp_addon(k-1,ii) + dtp_depweir(k-1,ii)
      end do

      ! weir depth from the top to the bottom of each stage
      do k = dtp_numstage(ii), 2, -1
         dtp_depweir(k-1,ii) = dtp_depweir(k,ii) + dtp_depweir(k-1,ii)
      end do
   end if

   !! Wet pond
   !!----------
   if (wtp_onoff(ii)==1) then
      if (wtp_imo(ii)<=0)      wtp_imo(ii) = 1
      if (wtp_iyr(ii)<=1000)   wtp_iyr(ii) = iyr
      if (wtp_k(ii) <=0)       wtp_k(ii) = 0.01
      if (wtp_evrsv(ii)<=0)    wtp_evrsv(ii) = 0.01
      if (wtp_evrsv(ii)>1)     wtp_evrsv(ii) = 0.99
      if (wtp_sdslope(ii)<3)   wtp_sdslope(ii) = 3. !COA Manual 1.6.6.C
      if (wtp_lenwdth(ii)<2)   wtp_lenwdth(ii) = 2. !COA Manual 1.6.6.C
      if (wtp_pdia(ii)<=0)     wtp_pdia(ii) = 0.1524 !meter(=6inches), COA Manual 1.6.6
      if (wtp_plen(ii)<=2)     wtp_plen(ii) = 2 !meter
      if (wtp_pmann(ii)<=0)    wtp_pmann(ii) = 0.012 ! concrete surface
      if (wtp_ploss(ii)<=0)    wtp_ploss(ii) = 0
   end if

   !fraction runoff that directly enters the channel
   if (num_sf(ii)>=1) then
      do k=1,num_sf(ii)
         bmpfr_sf = bmpfr_sf + sf_fr(k,ii)
      end do
   endif
   if (num_ri(ii)>=1) then
      do k=1,num_ri(ii)
         bmpfr_ri = bmpfr_ri + ri_fr(ii,k)
      end do
   endif
   if (bmpfr_sf>1.or.bmpfr_ri>1) then
      write (*,*) " "
      write (*,*) "Urban BMP Warning!!"
      write (*,*) "In subbasin ", ii
      write (*,*) "The fraction runoff draining to urban BMPs"
      write (*,*) " are larger than one, so the fraction values"
      write (*,*) " were automatically reassigned"
      write (*,*) " "
      do k=1,num_sf(ii)
         sf_fr(k,ii) = sf_fr(k,ii) / bmpfr_sf
      end do
      do k=1,num_sf(ii)
         ri_fr(ii,k) = ri_fr(ii,k) / bmpfr_ri
      end do
      bmpfr_sf = 1.; bmpfr_ri=1.
   endif

   !!Retention-Irrigation
   !!---------------------
   do k=1,num_ri(ii)
      ! skip the pond that has zero inflow
      if (ri_fr(ii,k)==0) cycle

      ! determine water quality volume for defult pond sizes
      !City of Austin Design Guideline 1.6.9.2 Table 1-12
      !Retention-Irrigation for Barton Springs Zone

      hwq = (1.8 * sub_ha_imp(ii) / sub_ha_urb(ii) + 0.6) !inches
      wqv = hwq / 12. * sub_ha_urb(ii) * ri_fr(ii,k) * 107639.104167 !ft3

      if (ri_dim(ii,k)==0) then
         !Determine pond size automatically based on City of Austin's Design Guideline 1.6
         ri_vol(ii,k) = wqv * 0.028317 !m3
         ri_dep(ii,k)=1.5 !assume 1.5m as default retention pond depth
         ri_sa(ii,k) = ri_vol(ii,k) / ri_dep(ii,k)
         ri_dd (ii,k)=60.0 !drawdown time, hr
         ri_k(ii,k)=2.5
         ri_evrsv(ii,k)=0.6
      else
         !Evaluate unit variables provided by user
         if (ri_sa(ii,k)<1.or.ri_vol(ii,k)<1) then
            !City of Austin Design Guideline 1.6.5
            ri_vol(ii,k) = wqv * 0.028317 !m3
            ri_sa(ii,k) = ri_vol(ii,k) / ri_dep(ii,k)
         end if

         ri_dep(ii,k) = ri_vol(ii,k) / ri_sa(ii,k)
         if (ri_dd(ii,k)<=0) ri_dd (ii,k)=72.0
         if (ri_k(ii,k) <=0 ) ri_k(ii,k)=2.5
         if (ri_evrsv(ii,k) <=0 ) ri_evrsv(ii,k)=0.1
         if (ri_dep(ii,k)<=1) ri_dep(ii,k)=1.
      end if

      ! draw down time [number of time step]
      ri_ndt = (ri_dd(ii,k) - 12.) * 60 / idt  !minus the first dry 12 hours no-pumping

      ! pumping rate that empties the basin in 72 hours with initial 12 hour of no operatrion
      ri_pumpv(ii,k) = ri_vol(ii,k) / ri_ndt !m3/dt

      if (ri_im(ii,k)<0.or.ri_im(ii,k)>12) ri_im(ii,k) = 0
      if (ri_iy(ii,k)<1000.and.ri_iy(ii,k)>0) ri_iy(ii,k) = 0
      if (ri_iy(ii,k)==0) ri_iy(ii,k) = iyr
      if (ri_im(ii,k)==0) ri_im(ii,k) = 1

      write(77779,'(a11,i5)') 'Subbasin #:', ii  ! bmp_sedfil.out
      write(77779,'(a46)') ''
      write(77779,'(a10,i5)') 'RI #:', K   ! bmp_sedfil.out
      write(77779,'(a17,f10.1,a4)') 'Total volume =', ri_vol(ii,k),'m^3'
      write(77779,'(a17,f10.1,a4)') 'Surface area =', ri_sa(ii,k),'m^2'
      write(77779,'(a17,f10.1,a3)') 'Drawdown time =', ri_dd (ii,k),'hr'
      write(77779,'(a17)') ''

   end do

   !!Sedimentation-Filtration
   !!---------------------

   do k=1,num_sf(ii)
      write(77778,'(a11,i5)') 'Subbasin #:', ii   ! bmp_sedfil.out
      write(77778,'(a46)') ''
      write(77778,'(a10,i5)') 'SED-FIL #:', k   ! bmp_sedfil.out
      !determine water quality volume for defult pond sizes
      !City of Austin Design Guideline 1.6.2
      hwq = (0.5 + sub_ha_imp(ii) / sub_ha_urb(ii) - 0.2) !inches
      wqv = hwq / 12. * sub_ha_urb(ii) * sf_fr(k,ii) * 107639.104167 !ft3

      if (sf_dim(k,ii)==0) then
         write(77778,'(a46)') 'This SED-FIL size is automatically'&
            &// ' estimated based on WQV.'
         !Determine pond size automatically based on City of Austin's Design Guideline 1.6
         if (sf_typ(k,ii)==1.or.sf_typ(k,ii)==3) then
            ! full scale or sedimentation basin only
            sp_pvol(k,ii) = wqv * 0.028317 !m3
            sp_sa(k,ii) = sp_pvol(k,ii) / 1.5 !assume 1.5m depth for sed pond
            sp_pd(k,ii) = Sqrt(4. * 2. * sp_sa(k,ii) * 1.5**0.5&
               &/ (0.6 * 172800. * 19.6**0.5) / 3.14159) * 1000. !mm
            ft_sa(k,ii) = wqv/(7.+2.33*4.) * 0.093 !m2
            ft_fsa(k,ii) = 1.

         else
            ! partial scale
            ft_sa(k,ii) = wqv * 0.028317 !m3
!               wqv/(4.+1.33*4.) * 0.093

         end if
         sp_bpw(k,ii) = 10. !m overflow weir width
         ft_pd(k,ii) = 1524. !mm
         ft_dep(k,ii) = 420. !mm
         ft_h(k,ii) = 1200. !mm

      else
         !Evaluate unit variables given by users
         if (sf_typ(k,ii)>3) sf_typ(k,ii) = 1
         if (sp_sa(k,ii)<1) then
            !City of Austin Design Guideline 1.6.5
            sp_pvol(k,ii) = wqv * 0.028317 !m3
            sp_sa(k,ii) = sp_pvol(k,ii) / 1.5 !assume 1.5m depth for sed pond
         end if
         if (sp_pd(k,ii)<0.1) sp_pd(k,ii) = 152.4 !mm diameter, equivalent to 6in
         if (ft_h(k,ii)<5) ft_h(k,ii) = 1200. !mm
         if (ft_sa(k,ii)<1) then
            if (sf_typ(k,ii)==1) then
               ft_sa(k,ii) = wqv/(7.+2.33*ft_h(k,ii)/304.8) * 0.093 !m2
               ft_fsa(k,ii) = 1.
            else
               ft_sa(k,ii) = wqv * 0.028317 !m3
               ft_fsa(k,ii) = 1. /(4.+1.33*ft_h(k,ii)/304.8)
            end if
         end if
      end if

      !Outflow control
      if (sp_qfg(k,ii)==0) then
         sp_pd(k,ii) = Sqrt(4. * 2. * sp_sa(k,ii) * 1.5**0.5&
            &/ (0.6 * 172800. * 19.6**0.5) / 3.14159) * 1000. !mm
      end if
      if (ft_qfg(k,ii)==0) then
         ft_pd(k,ii) = Sqrt(4. * 2. * ft_sa(k,ii) * 1.5**0.5&
            &/ (0.6 * 172800. * 19.6**0.5) / 3.14159) * 1000. !mm
      end if

      !Orifice pipe for sand filter should be equal or larger than
      !sedimentation pond outlet pipe for full-type SedFils
      if (ft_pd(k,ii)<sp_pd(k,ii)) ft_pd(k,ii) = sp_pd(k,ii)

      if (ft_dep(k,ii)<100) ft_dep(k,ii) = 100.
      if (sf_ptp(k,ii)>1) sf_ptp(k,ii) = 1
      !      if (sf_typ(k,ii)==1) sf_ptp(k,ii) = 0 removed by Jaehak 2014
      !      if (sp_pd(k,ii)>254) sp_pd(k,ii) = 254. ! max 10inches dia
      if (sp_pd(k,ii)<10) sp_pd(k,ii) = 10. ! min 10mm dia
      !     if (ft_pd(k,ii)>254) ft_pd(k,ii) = 254. ! max 10inches
      if (ft_k(k,ii)<1) ft_k(k,ii) = 1.
      if (ft_dp(k,ii)<0.0001) ft_dp(k,ii) = 0.02
      if (ft_dc(k,ii)<0.01) ft_dc(k,ii) = 0.762
      if (ft_por(k,ii)<0.1) ft_por(k,ii) = 0.45
      if (tss_den(k,ii)<0.5) tss_den(k,ii) = 0.5
      if (ft_alp(k,ii)<0.1) ft_alp(k,ii) = 0.1
      if (sf_im(k,ii)<0.or.sf_im(k,ii)>12) sf_im(k,ii) = 0
      if (sf_iy(k,ii)<1000.and.sf_iy(k,ii)>0) sf_iy(k,ii) = 0
      if (sf_iy(k,ii)==0) sf_iy(k,ii) = iyr
      if (sf_im(k,ii)==0) sf_im(k,ii) = 1
      if (ft_fsa(k,ii)==0) ft_fsa(k,ii) = 0.85
      if (sf_typ(k,ii)==1) ft_fsa(k,ii) = 1


      if (sf_typ(k,ii)==1) then
         write(77778,'(a37)') 'Full type sed-fil selected'
      elseif (sf_typ(k,ii)==2) then
         write(77778,'(a40)') 'Partial type sed-fil selected'
      else
         write(77778,'(a43)') 'Sedimentation pond only selected'
      endif
      if (sf_typ(k,ii)==1.or.sf_typ(k,ii)==3) then
         write(77778,'(a18)') 'Sedimentation pond'
         write(77778,'(a17,f10.1,a4)') 'Total volume =', sp_pvol(k,ii),'m^3'
         write(77778,'(a17,f10.1,a4)') 'Surface area =', sp_sa(k,ii),'m^2'
         write(77778,'(a17,f10.1,a3)') 'Drain Pipe Dia =', sp_pd(k,ii),'mm'
         write(77778,'(a17)') ''
      endif
      if (sf_typ(k,ii)==1.or.sf_typ(k,ii)==2) then
         write(77778,'(a11)') 'Sand Filter'
         write(77778,'(a17,f10.1,a4)') 'Surface area =', ft_sa(k,ii),'m^2'
         write(77778,'(a17,f10.1,a3)') 'Max ponding =', ft_h(k,ii),'mm'
         write(77778,'(a17,f10.1,a3)') 'Filter depth =', ft_dep(k,ii),'mm'
         write(77778,'(a17,f10.1,a3)') 'Drain Pipe Dia =', ft_pd(k,ii),'mm'
         write(77778,'(a17)') ''
      end if
   end do


   return
end
