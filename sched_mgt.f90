!> @file sched_mgt.f90
!> file containing the subroutine sched_mgt
!> @author
!> modified by Javier Burguete

!> this subroutine performs all management operations
!> @param[in] j HRU number
subroutine sched_mgt(j)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j          |none             |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    biomass
!!    husc
!!    ifrt
!!    irrsalt
!!    n
!!    ncrp
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Int
!!    SWAT: plantop, irrsub, fert, apply, harvkillop, newtillmix, harvestop,
!!          killop, graze, burnop

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: biomass, husc, irrsalt
   integer :: ifrt, n, ncrp

   n = nop(j)

   select case (mgtop(n,j))

    case (1)  !! plant operation
      igro(j) = 1
      lai_init = mgt5op(n,j)
      bio_init = mgt6op(n,j)
      hi_targ(j) = mgt7op(n,j)
      bio_targ(j) = mgt8op(n,j) * 1000.
      cnop = mgt9op(n,j)
      curyr_mat(j) = mgt3iop(n,j)
      if (curyr_mat(j) == 0) igrotree(j) = 1

      idplt(j) = mgt1iop(n,j)

      if (mgt4op(n,j) < 700.) mgt4op(n,j) = 1700.
!            if (mgt4op(n,j) > 5000.) mgt4op(n,j) = 5000.
      phu_plt(j) = mgt4op(n,j)

      call plantop(j)

      if (imgt == 1) then
         write (143, 1000) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j),cpnm(idplt(j))," PLANT", phubase(j), phuacc(j),&
            &sol_sw(j),bio_ms(j), sol_rsd(1,j),sol_sumno3(j),&
            &sol_sumsolp(j)
      end if


    case (2)  !! irrigation operation
      irr_sc(j) = mgt2iop(n,j)     !!NUBZ
      irr_no(j) = mgt10iop(n,j)
      irramt(j) = mgt4op(n,j)
      irrsalt = mgt5op(n,j) !! not used
      irrefm(j) = mgt6op(n,j)
      irrsq(j) = mgt7op(n,j)
      irr_flag(j) = 1

      if (irrefm(j) < 1.e-6) irrefm(j)=1.0
      if (irr_sc(j) <= 0) irr_sc(j) = irrsc(j)
      if (irr_no(j) <= 0) irr_no(j) = irrno(j)
      if (irr_no(j) <= 0) irr_no(j) = hru_sub(j)
      if (irr_sc(j) > 2) then    !! reach and res flag ??
         call irrsub(j)
      endif

      if (imgt ==1) then
         write (143, 1002) subnum(j), hruno(j), iyr, i_mo,&
            &iida, hru_km(j), "        ",&
            &"IRRIGATE", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j), sol_sumno3(j),sol_sumsolp(j),irramt(j),&
            &irr_sc(j), irr_no(j)
1002     format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,7f10.2,10x,f10.2,70x,2i7)

      end if


    case (3)   !! fertilizer operation
      ifrt = mgt1iop(n,j)
      frt_kg = mgt4op(n,j)
      frt_surface = mgt5op(n,j)
      if (frt_surface <= 1.e-6) frt_surface = 0.2

      call fert(j, ifrt)

      if (imgt ==1) then
         write (143, 1004) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), fertnm(ifrt),&
            &"   FERT", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j), sol_sumno3(j),sol_sumsolp(j),frt_kg,&
            &fertno3, fertnh3, fertorgn, fertsolp, fertorgp
1004     format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,7f10.2,20x,f10.2,10x,5f10.2)
      endif


    case (4)   !! pesticide operation
      hrupest(j) = 1
      ipest = mgt1iop(n,j)
      pst_kg = mgt4op(n,j)
      pst_dep = mgt5op(n,j)

      call apply(j)

      if (imgt ==1) then
         write (143, 1004) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), pname(ipest),&
            &"   PEST", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j), sol_sumno3(j), sol_sumsolp(j), pst_kg
      endif

    case (5)   !! harvest and kill operation
      cnop = mgt4op(n,j)
      hi_ovr = mgt5op(n,j)
      frac_harvk = mgt6op(n,j)
      biomass = bio_ms(j)

      call harvkillop(j)

      if (imgt ==1) then
         write (143, 1001) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), cpnm(idplt(j)),&
            &"HARV/KILL", phubase(j), phuacc(j), sol_sw(j),biomass,&
            &sol_rsd(1,j), sol_sumno3(j),sol_sumsolp(j),yield,&
            &strsn_sum(j), strsp_sum(j), strstmp_sum(j), strsw_sum(j),&
            &strsa_sum(j)
!!1001  format (a5,1x,a4,3i6,2a15,8f10.2,30x,11f10.2)
1001     format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,8f10.2,30x,5f10.2, 14x,6f10.2)
      end if

      phubase(j) = 0.
      phuacc(j) = 0.

    case (6)   !! tillage operation
      idtill = mgt1iop(n,j)
      cnop = mgt4op(n,j)

      call newtillmix(j,0.0D+00)

      if (imgt ==1) then
         write (143, 1003) subnum(j), hruno(j),iyr, i_mo, iida,&
            &hru_km(j), tillnm(idtill),&
            &"TILLAGE", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j),sol_sumno3(j),sol_sumsolp(j), effmix(idtill)
1003     format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,7f10.2,30x,f10.2)
      end if

    case (7)  !! harvest only operation
      hi_bms = mgt5op(n,j)
      hi_rsd = mgt6op(n,j)
      harveff = mgt4op(n,j)
      if (harveff <= 0.) harveff = 1.0
      call harvestop(j)

      if (imgt == 1) then
         write (143, 1001) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), cpnm(idplt(j)),&
            &"HARVEST ONLY", phubase(j), phuacc(j),sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j), sol_sumno3(j), sol_sumsolp(j), yield,&
            &strsn_sum(j), strsp_sum(j), strstmp_sum(j), strsw_sum(j),&
            &strsa_sum(j), yieldgrn, yieldbms, yieldtbr, yieldrsd,&
            &yieldn, yieldp
      end if

    case (8)   !! kill operation
      call killop(j)

      if (imgt == 1) then
         write (143, 1000) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), "         ",&
            &"    KILL", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j),sol_sumno3(j),sol_sumsolp(j)
      end if

      phubase(j) = 0.
      phuacc(j) = 0.

    case (9)    !! grazing operation
      manure_id(j) = mgt2iop(n,j)
      grz_days(j) = mgt1iop(n,j)
      bio_eat(j) = mgt4op(n,j)
      bio_trmp(j) = mgt5op(n,j)
      manure_kg(j) = mgt6op(n,j)
      ndeat(j) = 0
      igrz(j) = 1

      if (manure_kg(j) <= 0.) then
         manure_kg(j) = 0.95 * mgt4op(n,j)
      end if
      call graze(j)

      if (imgt == 1) then
         write (143, 1005) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), "         ",&
            &"   GRAZE", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j),sol_sumno3(j),sol_sumsolp(j),manure_kg(j)
1005     format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,7f10.2,20x,f10.2)
      end if

    case (10)   !! auto irrigation operation
      wstrs_id(j) = mgt1iop(n,j)
      auto_wstr(j) = mgt4op(n,j)
      irr_eff(j) = mgt5op(n,j)
      irr_mx(j) = mgt6op(n,j)
      irr_asq(j) = mgt7op(n,j)
      irr_sca(j) = mgt2iop(n,j)
      irr_noa(j) = mgt10iop(n,j)
      if (irr_noa(j) <= 0) irr_noa(j) = irrno(j)
      if (irr_noa(j) <= 0) irr_noa(j) = hru_sub(j)
      if (wstrs_id(j) <= 0) wstrs_id(j) = 1
      if (irr_eff(j) > 1.) irr_eff(j) = 0.
      if (irr_eff(j) == 0.) irr_eff(j) = 1.
      if (irr_mx(j) < 1.e-6) irr_mx(j) = 25.4
      if (irr_sca(j) <= 0) irr_sca(j) = irrsc(j)
      if (imgt ==1) then
         write (143, 1010) subnum(j), hruno(j), iyr, i_mo,&
            &iida, hru_km(j), "        ",&
            &"SCHED AUTORR", phubase(j), phuacc(j), sol_sw(j), bio_ms(j),&
            &sol_rsd(1,j), sol_sumno3(j),sol_sumsolp(j)
1010     format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,7f10.2)
      end if


    case (11)   !! auto fertilizer operation
      iafrttyp(j) = mgt1iop(n,j)
      nstress(j) = mgt2iop(n,j)
      auto_nstrs(j) = mgt4op(n,j)
      auto_napp(j) = mgt5op(n,j)
      if (auto_napp(j) < 1.e-6) auto_napp(j) = 250.
      auto_nyr(j) = mgt6op(n,j)
      if (auto_nyr(j) < 1.e-6) auto_nyr(j) = 350.
      auto_eff(j) = mgt7op(n,j)
      if (auto_eff(j) <= 0.) auto_eff(j) = 1.3
      afrt_surface(j) = mgt8op(n,j)
      if (afrt_surface(j) <= 1.e-6) afrt_surface(j) = .8
      !! calculate tnylda for autofertilization
      ncrp = idplt(j)
      if (tnylda(j) < 1.e-6)tnylda(j)=150.*cnyld(ncrp)*bio_e(ncrp)
      !      if (tnylda(j) < 1.e-6)tnylda(j)=350.*cnyld(ncrp)*bio_e(ncrp)
      !         tnylda(j) = 350. * cnyld(ncrp) * bio_e(ncrp)
      !        tnylda(j) = 350. * cnyld(ncrp) * bio_e(ncrp)
      !       else
      !         tnylda(j) = 1000. * cnyld(ncrp) * bio_e(ncrp)
      !    endif

    case (12)   !! street sweeping (only if iurban=2)

      ! husc was not defined
      husc = phu_op(n,j) !?
      if (husc > 0.) then
         ! igrow is not defined
         !if (igrow == 1) then
         if (igro(j) == 1) then
            phusw(j) = husc
         endif
      endif
      sweepeff = mgt4op(n,j)
      fr_curb = mgt5op(n,j)

      if (imgt == 1) then
         write (143, 1000) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), "         ",&
            &"STREET SWEEP",phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j),sol_sumno3(j),sol_sumsolp(j)
      end if

    case (13)    !! release/impound water in rice fields
      imp_trig(j) = mgt1iop(n,j)

      if (imgt == 1) then
         write (143, 1000) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), "         ","RELEASE/IMPOUND", phubase(j),&
            &phuacc(j),sol_sw(j),bio_ms(j),sol_rsd(1,j),sol_sumno3(j),&
            &sol_sumsolp(j)
      end if

    case (14)    !! continuous fertilization operation
      fert_days(j) = mgt1iop(n,j)
      cfrt_id(j) = mgt2iop(n,j)
      ifrt_freq(j) = mgt3iop(n,j)
      cfrt_kg(j) = mgt4op(n,j)
      icfrt(j) = 1
      ndcfrt(j) = 1
      iday_fert(j) = ifrt_freq(j)

    case (15)    !! continuous pesticide operation
      cpst_id(j) = mgt1iop(n,j)
      pest_days(j) = mgt2iop(n,j)
      ipst_freq(j) = mgt3iop(n,j)
      cpst_kg(j) = mgt4op(n,j)
      icpst(j) = 1
      ndcpst(j) = 0
      iday_pest(j) = ipst_freq(j)

    case (16)   !! burning
      burn_frlb = mgt4op(n,j)
      call burnop(j)
      if (imgt == 1) then
         write (143, 1000) subnum(j), hruno(j), iyr, i_mo, iida,&
            &hru_km(j), "         ",&
            &"      BURN", phubase(j), phuacc(j), sol_sw(j),bio_ms(j),&
            &sol_rsd(1,j),sol_sumno3(j),sol_sumsolp(j)
      end if

    case (17)    !! skip a year
      yr_skip(j) = 1

   end select

   if (mgtop(n,j) /= 17) then
      nop(j) = nop(j) + 1
   end if

   if (nop(j) > nopmx(j)) then
      nop(j) = 1
   end if

1000 format (a5,1x,a4,3i6,1x,e10.5,1x,2a15,19f10.2)
   return

end
