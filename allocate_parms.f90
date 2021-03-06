!> @file allocate_parms.f90
!> file containing the subroutine allocate_parms
!> @author 
!> modified by Javier Burguete

!> this subroutine allocates array sizes
subroutine allocate_parms

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mch         |none          |max number of channels
!!    mcr         |none          |max number of crops grown per year
!!    mcrdb       |none          |max nunber of crops in crop.dat
!!    mfdb        |none          |max number of fertilizers in fert.dat
!!    mgr         |none          |max number of grazings per year
!!    mhru        |none          |max number of HRUs
!!    mhyd        |none          |max number of hydrographs
!!    mlyr        |none          |max number of soil layers
!!    mnr         |none          |max number of years of rotation
!!    mpst        |none          |max number of pesticides used in wshed
!!    mpdb        |none          |max number of pesticides in pest.dat
!!    mrecc       |none          |max number of reccnst files
!!    mrecd       |none          |max number of recday files
!!    mrech       |none          |max number of rechour files
!!    mrecm       |none          |max number of recmon files
!!    mrecy       |none          |max number of recyear files
!!    mres        |none          |max number of reservoirs
!!    mrg         |none          |max number of rainfall/temp gages
!!    nstep       |none          |max number of time steps per day
!!    msub        |none          |max number of subbasins
!!    mtil        |none          |max number of tillage types in till.dat
!!    mudb        |none          |max number of urban land types in urban.dat
!!    myr         |none          |max number of years of simulation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    imho
!!    itempa
!!    mxsubch
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: zero0, zero1, zero2, zeroini, zero_urbn

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   integer :: imho, itempa, mxsubch

   allocate (surlag(mhru))
   allocate (cdn(mhru))
   allocate (cmn(mhru))
   allocate (nperco(mhru))
   allocate (phoskd(mhru))
   allocate (psp(mhru))
   allocate (sdnco(mhru))

   mxsubch = Max(msub+1,mch+1)
   itempa = Max(mhru,mch)

!!    new arrays for routing units
   allocate (hru_rufr(mru,mhru))
   allocate (daru_km(msub,mru))
   allocate (ru_k(msub,mru))
   allocate (ru_c(msub,mru))
   allocate (ru_eiq(msub,mru))
   allocate (ru_ovs(msub,mru))
   allocate (ru_ovsl(msub,mru))
   allocate (ru_a(msub,mru))
   allocate (ru_ktc(msub,mru))
   allocate (gwq_ru(mhru))
   allocate (ils2(mhru))
   allocate (ils2flag(msub))
   allocate (ifirsthr(mrech))

!!    arrays which contain data related to the number of recday commands
   allocate (ifirstr(mrecd))

!!    arrays which contain data related to the date
!   allocate (values(8))

!!    arrays which contain data related to rainfall/temperature gages
!!     test for JRW
   allocate (elevp(mrg))
   allocate (elevt(mrg))
   allocate (ifirstpcp(mrg))
   allocate (ifirstt(mrg))

!!    apex/command variables
   allocate (ifirsta(mapex))

!! septic inputs
   allocate (isep_hru(mhru))
   allocate (sptqs(msdb+1))
   allocate (sptbodconcs(msdb+1))
   allocate (spttssconcs(msdb+1))
!   allocate (spttnconcs(msdb+1)) ! not used
   allocate (sptnh4concs(msdb+1))
   allocate (sptno3concs(msdb+1))
   allocate (sptno2concs(msdb+1))
   allocate (sptorgnconcs(msdb+1))
!   allocate (spttpconcs(msdb+1)) ! not used
   allocate (sptminps(msdb+1))
   allocate (sptorgps(msdb+1))
   allocate (sptfcolis(msdb+1))
!! pothole changes for srini
   allocate (pot_seep(mhru))
   allocate (pot_solp(mhru))
   allocate (pot_orgp(mhru))
   allocate (pot_orgn(mhru))
   allocate (pot_mps(mhru))
   allocate (pot_mpa(mhru))
   allocate (tile_solpo(mhru))
!! septic changes added 1/28/09 gsm
   allocate (i_sep(mhru))
   allocate (sep_tsincefail(mhru))
   allocate (isep_tfail(mhru))
   allocate (failyr(mhru))
   allocate (qstemm(mhru))
   allocate (bz_area(mhru))
   allocate (bio_bod(mhru))
   allocate (biom(mhru))
   allocate (rbiom(mhru))
   allocate (fcoli(mhru))
   allocate (bz_perc(mhru))
   allocate (bz_z(mhru))
   allocate (bz_thk(mhru))
   allocate (bio_bd(mhru))
!!    carbon outputs for .hru file
   allocate (cmup_kgh(mhru))
   allocate (cmtot_kgh(mhru))
!!    carbon outputs for .hru file
   allocate (coeff_bod_dc(mhru))
   allocate (coeff_bod_conv(mhru))
   allocate (coeff_fc1(mhru))
   allocate (coeff_fc2(mhru))
   allocate (coeff_fecal(mhru))
   allocate (coeff_plq(mhru))
   allocate (coeff_mrt(mhru))
   allocate (coeff_rsp(mhru))
   allocate (coeff_slg1(mhru))
   allocate (coeff_slg2(mhru))
   allocate (coeff_nitr(mhru))
   allocate (coeff_denitr(mhru))
   allocate (isep_typ(mhru))
   allocate (isep_opt(mhru))
   allocate (plqm(mhru))
   allocate (coeff_pdistrb(mhru))
   allocate (coeff_psorpmax(mhru))
   allocate (coeff_solpslp(mhru))
   allocate (coeff_solpintc(mhru))
   allocate (isep_iyr(mhru))

!! septic changes added 1/28/09 gsm
   allocate (qird(mhru))

!!    arrays which contain data related to channels
   allocate (algae(mch))
   allocate (alpha_bnk(mxsubch))
   allocate (alpha_bnke(mxsubch))
   allocate (ammonian(mch))
   allocate (bankst(mch))
   allocate (bc(4,mch))
   allocate (ch_bnk_bd(mch))
   allocate (ch_cov(2,mch))
   allocate (ch_bnk_d50(mch))
   allocate (ch_eqn(mch))
   allocate (ch_li(mch))
   allocate (ch_onco(mch))
   allocate (ch_opco(mch))
   allocate (ch_orgn(mch))
   allocate (ch_orgp(mch))
   allocate (ch_si(mch))
   allocate (ch_wdr(mch))
   allocate (ch_wi(mch))
   allocate (ch_erodmo(12,mch))
   allocate (chlora(mch))
   allocate (chpst_conc(mch))
   allocate (chpst_koc(mch))
   allocate (chpst_mix(mch))
   allocate (chpst_rea(mch))
   allocate (chpst_rsp(mch))
   allocate (chpst_stl(mch))
   allocate (chpst_vol(mch))
   allocate (dep_chan(mch))
   allocate (disolvp(mch))
   allocate (flwin(mch))
   allocate (flwout(mch))
   allocate (icanal(mch))
   allocate (nitraten(mch))
   allocate (nitriten(mch))
   allocate (organicn(mch))
   allocate (organicp(mch))
   allocate (orig_sedpstconc(mch))
   allocate (rch_bactlp(mch))
   allocate (rch_bactp(mch))
   allocate (rch_cbod(mch))
   allocate (rch_dox(mch))
   allocate (rchstor(mch))
   allocate (rk(6,mch))
   allocate (rs(7,mch))
   allocate (sedpst_act(mch))
   allocate (sedpst_bry(mch))
   allocate (sedpst_conc(mch))
   allocate (sedpst_rea(mch))
   allocate (sedst(mch))
   allocate (vel_chan(mch))

   allocate (wurch(12,mxsubch))

!!    arrays for channel added by Balaji for the new routines
   allocate (ch_bnk_san(mch))
   allocate (ch_bnk_sil(mch))
   allocate (ch_bnk_cla(mch))
   allocate (ch_bnk_gra(mch))
   allocate (ch_bed_san(mch))
   allocate (ch_bed_sil(mch))
   allocate (ch_bed_cla(mch))
   allocate (ch_bed_gra(mch))
   allocate (ch_bed_bd(mch))
   allocate (ch_bed_d50(mch))
   allocate (ch_bnk_kd(mch))
   allocate (ch_bed_kd(mch))
   allocate (tc_bnk(mch))
   allocate (tc_bed(mch))
   allocate (depfp(mch))
   allocate (depprfp(mch))
   allocate (depsilfp(mch))
   allocate (depclafp(mch))
   allocate (depch(mch))
   allocate (depprch(mch))
   allocate (depsanch(mch))
   allocate (depsilch(mch))
   allocate (depclach(mch))
   allocate (depsagch(mch))
   allocate (deplagch(mch))
   allocate (depgrach(mch))
   allocate (sanst(mch))
   allocate (silst(mch))
   allocate (clast(mch))
   allocate (sagst(mch))
   allocate (lagst(mch))
   allocate (grast(mch))

!!    arrays which contain data related to reach output
   allocate (rchaao(mrcho,mxsubch))
   allocate (rchdy(mrcho,mxsubch))
   allocate (rchmono(mrcho,mxsubch))
   allocate (rchyro(mrcho,mxsubch))

!!    arrays which contain data related to subbasins
   allocate (ch_revap(mxsubch))
   allocate (cncoef_sub(msub))
   allocate (co2(msub))
   allocate (daylmn(msub))
   allocate (drydep_no3(msub))
   allocate (drydep_nh4(msub))

!!!   atmospheric deposition by month
   allocate (rcn_mo(motot,msub))
   allocate (rammo_mo(motot,msub))
   allocate (drydep_no3_mo(motot,msub))
   allocate (drydep_nh4_mo(motot,msub))
!!!   atmospheric deposition by day
   allocate (rcn_d(msub))
   allocate (rammo_d(msub))
   allocate (drydep_no3_d(msub))
   allocate (drydep_nh4_d(msub))

   allocate (fcst_reg(msub))
   allocate (harg_petco(msub))
   allocate (hqdsave(msub,nstep*4))  !! was 49, changed for urban -> changed to 2d array J.Jeong 4/17/2009
   allocate (hsdsave(msub,nstep*4))  !! J.Jeong 4/22/2009
   allocate (hru1(msub))
   allocate (hrutot(msub))
   allocate (ihgage(msub))
   allocate (ireg(msub))
   allocate (irelh(msub))
   allocate (irgage(msub))
   allocate (isgage(msub))
   allocate (itb(msub))
   allocate (itgage(msub))
   allocate (iwgage(msub))
   allocate (latcos(msub))
   allocate (latsin(msub))
   allocate (pcpdays(msub))
   allocate (phutot(msub))
   allocate (plaps(msub))
   allocate (rammo_sub(msub))
   allocate (rcn_sub(msub))
   allocate (sub_bactlp(msub))
   allocate (sub_bactp(msub))
   allocate (sub_bd(msub))
   allocate (sub_cbod(msub))
   allocate (sub_chl(msub))
   allocate (sub_dsan(msub))
   allocate (sub_dsil(msub))
   allocate (sub_dcla(msub))
   allocate (sub_dsag(msub))
   allocate (sub_dlag(msub))
   allocate (sub_dox(msub))
   allocate (sub_elev(msub))
   allocate (sub_etday(msub))
   allocate (sub_fr(msub))
   allocate (sub_gwno3(msub))
   allocate (sub_gwsolp(msub))
   allocate (sub_gwq(msub))
   allocate (sub_gwq_d(msub))
   allocate (sub_km(msub))
   allocate (sub_lat(msub))
   allocate (sub_latq(msub))
   allocate (sub_tileq(msub))
   allocate (sub_vaptile(msub))
   allocate (sub_latno3(msub))
   allocate (sub_minpa(msub))
   allocate (sub_minps(msub))
   allocate (sub_no3(msub))
   allocate (sub_orgn(msub))
   allocate (sub_orgp(msub))

   allocate (sub_precip(msub))
   allocate (sub_qd(msub))
   allocate (sub_sedpa(msub))
   allocate (sub_sedps(msub))
   allocate (sub_sedy(msub))
   allocate (sub_sep(msub))
   allocate (sub_tileno3(msub))
   allocate (sub_snom(msub))
   allocate (sub_solp(msub))
   allocate (sub_solpst(msub))
   allocate (sub_sorpst(msub))
   allocate (sub_subp(msub))
   allocate (sub_sumfc(msub))
   allocate (sub_surfq(msub))
   allocate (sub_sw(msub))
   allocate (sub_tc(msub))
   allocate (sub_tran(msub))
   allocate (sub_wtmp(msub))
   allocate (sub_wyld(msub))
   allocate (sub_yorgn(msub))
   allocate (sub_yorgp(msub))
   allocate (subfr_nowtr(msub))
   allocate (subgis(msub))
   allocate (tlaps(msub))
   allocate (tmp_an(msub))
   allocate (wcklsp(msub))
   allocate (welev(msub))

   allocate (huminc(12,msub))
   allocate (radinc(12,msub))
   allocate (rfinc(12,msub))
   allocate (tmpinc(12,msub))

   allocate (sub_sftmp(10,msub))
   allocate (sub_smtmp(10,msub))
   allocate (sub_smfmx(10,msub))
   allocate (sub_smfmn(10,msub))
   allocate (sub_timp(10,msub))
   allocate (sub_hhqd(msub,nstep))  ! 24 changed to nstep 4 urban modeling  Oct. 19,2007
   allocate (sub_hhwtmp(msub,nstep))   ! 24 changed to nstep 4 urban modeling  Oct. 19,2007
   allocate (uh(msub,nstep*3+1))      !! was 49 changed to nstep  OCt22, 2007

   allocate (ch_k(2,mxsubch))
   allocate (ch_n(2,mxsubch))
   allocate (ch_s(2,mxsubch))
   allocate (ch_w(2,mxsubch))
   allocate (ch_l(2,mxsubch))
   allocate (ch_d(mxsubch))
   allocate (chside(mxsubch))
   allocate (ch_di(mxsubch))
   allocate (sub_pet(mxsubch))
   allocate (elevb(10,msub))
   allocate (elevb_fr(10,msub))
   allocate (amp_r(12,msub))
   allocate (dewpt(12,msub))
   allocate (pcf(12,msub))
   allocate (solarav(12,msub))
   allocate (tmpmn(12,msub))
   allocate (tmpmx(12,msub))
   allocate (tmpstdmn(12,msub))
   allocate (tmpstdmx(12,msub))
   allocate (wndav(12,msub))

   allocate (pcp_stat(12,3,msub))
   allocate (pr_w(3,12,msub))

!!    arrays which contain data related to forecast parameters
   allocate (ftmpmn(12,msub))
   allocate (ftmpmx(12,msub))
   allocate (ftmpstdmn(12,msub))
   allocate (ftmpstdmx(12,msub))
   allocate (fpcp_stat(12,3,msub))
   allocate (fpr_w(2,12,msub))
   allocate (otmpmn(12,msub))
   allocate (otmpmx(12,msub))
   allocate (otmpstdmn(12,msub))
   allocate (otmpstdmx(12,msub))
   allocate (opcp_stat(12,3,msub))
   allocate (opr_w(3,12,msub))

!!    arrays which contain data related to subbasin output
   allocate (submono(msubo,msub))
   allocate (subaao(msubo,msub))
   allocate (subyro(msubo,msub))

!!    arrays which contain data related to soil layers, HRUs
!    Drainmod tile equations  01/2006
   allocate (vwt(mlyr,mhru))
   allocate (wat_tbl(mhru))
   allocate (sol_stpwt(mlyr,mhru)) !Moriasi 4/8/2014
!    Drainmod tile equations  01/2006
   allocate (conv_wt(mlyr,mhru))
   allocate (crdep(mlyr,mhru))
   allocate (flat(mlyr,mhru))
   allocate (orig_solactp(mlyr,mhru))
   allocate (orig_solaorgn(mlyr,mhru))
   allocate (orig_solfon(mlyr,mhru))
   allocate (orig_solfop(mlyr,mhru))
   allocate (orig_solno3(mlyr,mhru))
   allocate (orig_solorgn(mlyr,mhru))
   allocate (orig_solorgp(mlyr,mhru))
   allocate (orig_solrsd(mlyr,mhru))
   allocate (orig_solsolp(mlyr,mhru))
   allocate (orig_solst(mlyr,mhru))
   allocate (orig_solstap(mlyr,mhru))
   allocate (orig_soltmp(mlyr,mhru))
   allocate (orig_volcr(mlyr,mhru))
   allocate (pperco_sub(mlyr,mhru))
   allocate (sol_actp(mlyr,mhru))
   allocate (sol_aorgn(mlyr,mhru))
   allocate (sol_awc(mlyr,mhru))
   allocate (sol_bd(mlyr,mhru))
   allocate (sol_cbn(mlyr,mhru))
   allocate (sol_clay(mlyr,mhru))
!  added 1/27/09 when making septic changes
   allocate (sol_ec(mlyr,mhru))
!  added 1/27/09 when making septic changes
   allocate (sol_fc(mlyr,mhru))
   allocate (sol_fon(mlyr,mhru))
   allocate (sol_fop(mlyr,mhru))
   allocate (sol_hk(mlyr,mhru))
   allocate (sol_k(mlyr,mhru))
   allocate (sol_nh3(mlyr,mhru))
   allocate (sol_no3(mlyr,mhru))
   allocate (sol_orgn(mlyr,mhru))
   allocate (sol_orgp(mlyr,mhru))
   allocate (sol_por(mlyr,mhru))
   allocate (sol_prk(mlyr,mhru))
   allocate (sol_rock(mlyr,mhru))
   allocate (sol_rsd(mlyr,mhru))
   allocate (sol_sand(mlyr,mhru))
   allocate (sol_silt(mlyr,mhru))
   allocate (sol_solp(mlyr,mhru))
   allocate (sol_st(mlyr,mhru))
   allocate (sol_stap(mlyr,mhru))
   allocate (sol_tmp(mlyr,mhru))
   allocate (sol_ul(mlyr,mhru))
   allocate (sol_up(mlyr,mhru))
   allocate (sol_wp(mlyr,mhru))
   allocate (sol_wpmm(mlyr,mhru))
   allocate (sol_z(mlyr,mhru))
   allocate (volcr(mlyr,mhru))

!!    arrays which contain data related to soil layers, HRUs, pesticides
   allocate (orig_solpst(mpst,mhru,mlyr))
   allocate (sol_kp(mpst,mhru,mlyr))
   allocate (sol_pst(mpst,mhru,mlyr))

!!    arrays which contain data related to reservoirs
   allocate (br(2,mres))
   allocate (chlar(mres))
   allocate (evrsv(mres))
   allocate (iflodr(2,mres))
   allocate (ires(2,mres))
   allocate (iresco(mres))
   allocate (iyres(mres))
   allocate (lkpst_conc(mres))
   allocate (lkpst_koc(mres))
   allocate (lkpst_mix(mres))
   allocate (lkpst_rea(mres))
   allocate (lkpst_rsp(mres))
   allocate (lkpst_stl(mres))
   allocate (lkpst_vol(mres))
   allocate (lkspst_act(mres))
   allocate (lkspst_bry(mres))
   allocate (lkspst_conc(mres))
   allocate (lkspst_rea(mres))
   allocate (theta_n(mres))
   allocate (theta_p(mres))
   allocate (con_nirr(mres))
   allocate (con_pirr(mres))
   allocate (mores(mres))
   allocate (ndtargr(mres))
   allocate (oflowmn_fps(mres))
   allocate (orig_lkpstconc(mres))
   allocate (orig_lkspstconc(mres))
   allocate (orig_resnh3(mres))
   allocate (orig_resno2(mres))
   allocate (orig_resno3(mres))
   allocate (orig_resorgn(mres))
   allocate (orig_resorgp(mres))
   allocate (orig_ressed(mres))
   allocate (orig_ressolp(mres))
   allocate (orig_resvol(mres))
   allocate (res_evol(mres))
   allocate (res_k(mres))
   allocate (res_nh3(mres))
   allocate (res_no2(mres))
   allocate (res_no3(mres))
   allocate (res_nsed(mres))
   allocate (res_orgn(mres))
   allocate (res_orgp(mres))
   allocate (res_psa(mres))
   allocate (res_pvol(mres))
   allocate (res_rr(mres))
   allocate (res_seci(mres))
   allocate (res_sed(mres))

   allocate (res_san(mres))
   allocate (res_sil(mres))
   allocate (res_cla(mres))
   allocate (res_sag(mres))
   allocate (res_lag(mres))
   allocate (res_gra(mres))

   allocate (res_solp(mres))
   allocate (res_sub(mres))
   allocate (res_vol(mres))
   allocate (seccir(mres))
   allocate (sed_stlr(mres))
   allocate (starg_fps(mres))
   allocate (weirc(mres))
   allocate (weirk(mres))
   allocate (weirw(mres))
   allocate (acoef(mres))
   allocate (bcoef(mres))
   allocate (ccoef(mres))
   allocate (wurtnf(mres))
   allocate (lkpst_mass(mres))
   allocate (lkspst_mass(mres))

   allocate (nsetlr(2,mres))
   allocate (psetlr(2,mres))
   allocate (oflowmx(12,mres))
   allocate (oflowmn(12,mres))
   allocate (starg(12,mres))
   allocate (wuresn(12,mres))

   !!  added per JGA for Srini by gsm 9/8/2011
   !! arrays for mangement output (output.mgt)
   allocate (sol_sumno3(mhru))
   allocate (sol_sumsolp(mhru))
   allocate (strsw_sum(mhru))
   allocate (strstmp_sum(mhru))
   allocate (strsn_sum(mhru))
   allocate (strsp_sum(mhru))
   allocate (strsa_sum(mhru))
   allocate (velsetlr(mres))

!! arrays for reservoir output
   allocate (resoutm(41,mres))
   allocate (resouta(41,mres))
   allocate (resouty(41,mres))

!!    arrays which contain data related to reservoirs, year
   allocate (res_out(mres,12,myr))

!!    arrays which contain data related to pesticides in database
   allocate (ap_ef(mpdb))
   allocate (decay_f(mpdb))
   allocate (decay_s(mpdb))
   allocate (nope(mpdb))
   allocate (pst_wof(mpdb))
   allocate (pst_wsol(mpdb))
   allocate (skoc(mpdb))

!!    arrays which contain data related to landcover/landuse in database
   allocate (alai_min(mcrdb))
   allocate (bio_e(mcrdb))
   allocate (bio_leaf(mcrdb))
   allocate (bio_n(2,mcrdb))
   allocate (bio_p(2,mcrdb))
   allocate (blai(mcrdb))
   allocate (bm_dieoff(mcrdb))
   allocate (bmx_trees(mcrdb))
   allocate (chtmx(mcrdb))
   allocate (cnyld(mcrdb))
   allocate (cpyld(mcrdb))
   allocate (cvm(mcrdb))
   allocate (dlai(mcrdb))
   allocate (ext_coef(mcrdb))
   allocate (gsi(mcrdb))
   allocate (hvsti(mcrdb))
   allocate (idc(mcrdb))
   allocate (leaf(2,mcrdb))
   allocate (mat_yrs(mcrdb))
   allocate (rdmx(mcrdb))
   allocate (rsdco_pl(mcrdb))
   allocate (rsr(2,mcrdb))
   allocate (t_base(mcrdb))
   allocate (t_opt(mcrdb))
   allocate (vpd2(mcrdb))
   allocate (wac2(2,mcrdb))
   allocate (wavp(mcrdb))
   allocate (wsyf(mcrdb))

   allocate (pltnfr1(mcrdb))
   allocate (pltnfr3(mcrdb))
   allocate (pltpfr1(mcrdb))
   allocate (pltpfr3(mcrdb))

!!    arrays which contain data related to fertilizers in database
   allocate (bactkddb(mfdb))
   allocate (bactlpdb(mfdb))
   allocate (bactpdb(mfdb))
   allocate (fminn(mfdb))
   allocate (fminp(mfdb))
   allocate (fnh3n(mfdb))
   allocate (forgn(mfdb))
   allocate (forgp(mfdb))

!!    arrays which contain data related to urban land types in database
   allocate (curbden(mudb))
   allocate (dirtmx(mudb))
   allocate (fcimp(mudb))
   allocate (fimp(mudb))
   allocate (thalf(mudb))
   allocate (tnconc(mudb))
   allocate (tno3conc(mudb))
   allocate (tpconc(mudb))
   allocate (urbcoef(mudb))
   allocate (urbcn2(mudb))

!!    arrays which contain data related to years of rotation,
!!    applications, and HRUs
   allocate (auto_wstr(mhru))
!! burn 3/5/09

   allocate (cfrt_id(mhru))
   allocate (cfrt_kg(mhru))
   allocate (cpst_id(mhru))
   allocate (cpst_kg(mhru))

   allocate (wstrs_id(mhru))
   allocate (ifrt_freq(mhru))
   allocate (ipst_freq(mhru))

   allocate (imp_trig(mhru))
   allocate (irr_asq(mhru))
   allocate (irr_mx(mhru))
   allocate (irrsq(mhru))
   allocate (irr_eff(mhru))
   allocate (irrefm(mhru))

   allocate (irr_sc(mhru))
   allocate (irr_no(mhru))
   allocate (irr_sca(mhru))
   allocate (irr_noa(mhru))
   allocate (fert_days(mhru))
   allocate (pest_days(mhru))

   !!   burn 3/5/09









!!    changes pesticide incorporation in soil 3/31/08 gsm



!!    arrays which contain data related to years of rotation,
!!    crops grown per year, and HRUs

   allocate (bio_aahv(mcr,mhru))
   allocate (bio_hv(mcr,mhru))


   allocate (hi_targ(mhru))
   allocate (idplt(mhru))
   allocate (mcrhru(mhru))
   allocate (idplrot(mcr,mhru))
   allocate (ncrops(mcr,mhru))
   allocate (orig_phu(mhru))

   allocate (phu_plt(mhru))
   allocate (nstress(mhru))
   allocate (igrotree(mhru))
   allocate (tnyld(mhru))
   allocate (tnylda(mhru))
   allocate (yldkg(mcr,mhru))
   allocate (yldn(mcr,mhru))

!!    arrays which contain data related to years of rotation,
!!    grazings per year, and HRUs
   allocate (bio_eat(mhru))
   allocate (bio_trmp(mhru))
   allocate (grz_days(mhru))
   allocate (manure_id(mhru))
   allocate (manure_kg(mhru))

!!    arrays which contain data related to years of rotation,
!!    cuttings per year, and HRUs


!!    arrays which contain data related to tillages in the database
   allocate (deptil(mtil))
   allocate (effmix(mtil))
!! drainmod tile equations   06/2006
   allocate (ranrns(mtil))
   allocate (ranrns_hru(mhru))
!! drainmod tile equations   06/2006

!!    arrays which contain data related to hydrograph nodes
   allocate (hyd_dakm(mhyd))
   allocate (icodes(mhyd))
   allocate (ihouts(mhyd))
   allocate (inum1s(mhyd))
   allocate (inum2s(mhyd))
   allocate (inum3s(mhyd))
   allocate (inum4s(mhyd))
   allocate (inum5s(mhyd))
   allocate (inum6s(mhyd))
   allocate (inum7s(mhyd))
   allocate (inum8s(mhyd))
   allocate (reccnstps(mhyd))
   allocate (recmonps(mhyd))
   allocate (rnum1s(mhyd))
   allocate (subed(mhyd))
   allocate (subnum(mhru))
   allocate (hruno(mhru))

   allocate (shyd(8,mhyd))
   allocate (varoute(mvaro,mhyd))
   allocate (vartran(mvaro,mhyd))
   allocate (hhvaroute(mvaro,mhyd,nstep))  !! from 24 to nstep for urban

!!    arrays which contain data related to HRUs
   allocate (aairr(mhru))
   allocate (afrt_surface(mhru))
   allocate (aird(mhru))
   allocate (alpha_bf(mhru))
   allocate (alpha_bfe(mhru))
   allocate (alpha_bfe_d(mhru))
   allocate (anano3(mhru))
   allocate (anion_excl(mhru))
   allocate (auto_eff(mhru))
   allocate (auto_nyr(mhru))
   allocate (auto_napp(mhru))
   allocate (auto_nstrs(mhru))
   allocate (bactlp_plt(mhru))
   allocate (bactlpq(mhru))
   allocate (bactlps(mhru))
   allocate (bactp_plt(mhru))
   allocate (bactpq(mhru))
   allocate (bactps(mhru))
   allocate (bio_aams(mhru))
   allocate (bio_min(mhru))
   allocate (bio_ms(mhru))
   allocate (bio_yrms(mhru))
   allocate (biomix(mhru))
   allocate (bp(2,mhru))
   allocate (brt(mhru))
   allocate (bw(2,mhru))
   allocate (canmx(mhru))
   allocate (canstor(mhru))
   allocate (cbodu(mhru))
   allocate (chl_a(mhru))
   allocate (chlap(mhru))
   allocate (chlaw(mhru))
   allocate (cht(mhru))
   allocate (cklsp(mhru))
   allocate (cn(3,mhru))
   allocate (cnday(mhru))
   allocate (cont_cn(20,mhru))
   allocate (cont_p(20,mhru))
   allocate (cropno_upd(20,mhru))
   allocate (curyr_mat(mhru))
   allocate (dayl(mhru))
   allocate (det_san(mhru))
   allocate (det_sil(mhru))
   allocate (det_cla(mhru))
   allocate (det_sag(mhru))
   allocate (det_lag(mhru))
   allocate (drain_co(mhru))
   allocate (ddrain(mhru))
   allocate (deepirr(mhru))
   allocate (deepst(mhru))
   allocate (delay(mhru))
   allocate (dep_imp(mhru))
   allocate (dis_stream(mhru))
   allocate (divmax(mhru))
   allocate (dormhr(mhru))
   allocate (doxq(mhru))
   allocate (drain_d(20,mhru))
   allocate (drain_idep(20,mhru))
   allocate (drain_t(20,mhru))
   allocate (drain_g(20,mhru))
   allocate (dr_sub(mhru))
   allocate (epco(mhru))
   allocate (esco(mhru))
   allocate (erorgn(mhru))
   allocate (erorgp(mhru))
   allocate (evpnd(mhru))
   allocate (evwet(mhru))
   allocate (ffc(mhru))
   allocate (filterw(mhru))
   allocate (fire_cn(20,mhru))
   allocate (flowfr(mhru))
   allocate (flowmin(mhru))
   allocate (fsred(mhru))
   allocate (gdrain(mhru))
   allocate (grwat_n(mhru))
   allocate (grwat_i(mhru))
   allocate (grwat_l(mhru))
   allocate (grwat_w(mhru))
   allocate (grwat_d(mhru))
   allocate (grwat_s(mhru))
   allocate (grwat_spcon(mhru))
   allocate (gwati(20,mhru))
   allocate (gwatn(20,mhru))
   allocate (gwatl(20,mhru))
   allocate (gwatw(20,mhru))
   allocate (gwatd(20,mhru))
   allocate (gwats(20,mhru))
   allocate (gwatspcon(20,mhru))
   allocate (gw_delaye(mhru))
   allocate (gw_nloss(mhru))
   allocate (gw_q(mhru))
   allocate (gw_qdeep(mhru))
   allocate (gw_revap(mhru))
   allocate (gw_spyld(mhru))
   allocate (gwht(mhru))
   allocate (gwminp(mhru))
   allocate (gwqmn(mhru))
   allocate (hi_upd(20,mhru))
   allocate (hru_dafr(mhru))
   allocate (hru_fr(mhru))
   allocate (hru_ha(mhru))
   allocate (hru_km(mhru))
   allocate (hru_ra(mhru))
   allocate (hru_rmx(mhru))
   allocate (hru_slp(mhru))
   allocate (hru_sub(mhru))
   allocate (hru_seq(mhru))
   allocate (hrugis(mhru))
   allocate (hrupest(mhru))
   allocate (hvstiadj(mhru))
   allocate (iafrttyp(mhru))
   allocate (icfrt(mhru))
   allocate (icpst(mhru))
   allocate (icr(mhru))
   allocate (icrmx(mhru))
   allocate (iday_fert(mhru))
   allocate (iday_pest(mhru))
   allocate (idorm(mhru))
   allocate (iflod(2,mhru))
   allocate (igro(mhru))
   allocate (igrz(mhru))
   allocate (iopday(iopera,mhru))
   allocate (iopyr(iopera,mhru))
   allocate (mgt_ops(iopera,mhru))
   allocate (ioper(mhru))
!      allocate (mcri(mhru))
   allocate (mgtop(iopera,mhru))
   allocate (idop(iopera,mhru))
   allocate (phu_op(iopera,mhru))
   allocate (mgt1iop(iopera,mhru))
   allocate (mgt2iop(iopera,mhru))
   allocate (mgt3iop(iopera,mhru))
   allocate (mgt4op(iopera,mhru))
   allocate (mgt5op(iopera,mhru))
   allocate (mgt6op(iopera,mhru))
   allocate (mgt7op(iopera,mhru))
   allocate (mgt8op(iopera,mhru))
   allocate (mgt9op(iopera,mhru))
   allocate (mgt10iop(iopera,mhru))
   allocate (nopmx(mhru))
   allocate (irramt(mhru))
   allocate (yr_skip(mhru))
   allocate (phusw(mhru))
   allocate (bio_targ(mhru))
   allocate (irr_flag(mhru))
   imho = max(mhru,20)
   allocate (ipdhru(imho))
   allocate (ipnd(2,mhru))
   allocate (irn(mhru))
   allocate (irrno(mhru))
   allocate (irrsc(mhru))
   allocate (iurban(mhru))
   allocate (lai_aamx(mhru))
   allocate (lai_yrmx(mhru))
   allocate (laiday(mhru))
   allocate (laimxfr(mhru))
   allocate (latksatf(mhru))
   allocate (laimx_upd(20,mhru))
   allocate (lat_sed(mhru))
   allocate (lat_ttime(mhru))
   allocate (latno3(mhru))
   allocate (latq(mhru))
   allocate (ldrain(mhru))
   allocate (minpgw(mhru))
   allocate (ndeat(mhru))
   allocate (ndcfrt(mhru))
   allocate (ndcpst(mhru))
   allocate (ndtarg(mhru))
   allocate (newrti(mhru))
   allocate (nmgt(mhru))
   allocate (nop(mhru))
   allocate (no3gw(mhru))
   allocate (npcp(mhru))
   allocate (nplnt(mhru))
   allocate (olai(mhru))
   allocate (orgn_con(mhru))
   allocate (orgp_con(mhru))
   allocate (orig_alai(mhru))
   allocate (orig_bioms(mhru))
   allocate (orig_deepst(mhru))
   allocate (orig_igro(mhru))
   allocate (orig_phuacc(mhru))
   allocate (orig_pndno3(mhru))
   allocate (orig_pndorgn(mhru))
   allocate (orig_pndorgp(mhru))
   allocate (orig_pndsed(mhru))
   allocate (orig_pndsolp(mhru))
   allocate (orig_pndvol(mhru))
   allocate (orig_potno3(mhru))
   allocate (orig_potsed(mhru))
   allocate (orig_potvol(mhru))
   allocate (orig_shallst(mhru))
   allocate (orig_snohru(mhru))
   allocate (orig_solcov(mhru))
   allocate (orig_solsw(mhru))
   allocate (orig_sumix(mhru))
   allocate (orig_wetno3(mhru))
   allocate (orig_wetorgn(mhru))
   allocate (orig_wetorgp(mhru))
   allocate (orig_wetsed(mhru))
   allocate (orig_wetsolp(mhru))
   allocate (orig_wetvol(mhru))
   allocate (ov_n(mhru))
   allocate (ovrlnd(mhru))
!    Drainmod tile equations  01/2006
   allocate (pc(mhru))
!    Drainmod tile equations  01/2006
   allocate (percn(mhru))
   allocate (phuacc(mhru))
   allocate (phubase(mhru))
   allocate (plantn(mhru))
   allocate (plantp(mhru))
   allocate (plt_et(mhru))
   allocate (plt_pet(mhru))
   allocate (pltfr_n(mhru))
   allocate (pltfr_p(mhru))
   allocate (pnd_chla(mhru))
   allocate (pnd_esa(mhru))
   allocate (pnd_evol(mhru))
   allocate (pnd_fr(mhru))
   allocate (pnd_k(mhru))
   allocate (pnd_no3(mhru))
   allocate (pnd_no3g(mhru))
   allocate (pnd_no3s(mhru))
   allocate (pnd_nsed(mhru))
   allocate (pnd_orgn(mhru))
   allocate (pnd_orgp(mhru))
   allocate (pnd_psa(mhru))
   allocate (pnd_psed(mhru))
   allocate (pnd_pvol(mhru))
   allocate (pnd_seci(mhru))
   allocate (pnd_sed(mhru))

   allocate (pnd_san(mhru))
   allocate (pnd_sil(mhru))
   allocate (pnd_cla(mhru))
   allocate (pnd_sag(mhru))
   allocate (pnd_lag(mhru))

   allocate (twlpnd(mhru))     !!srini pond/wet infiltration to shallow gw storage
   allocate (twlwet(mhru))     !!srini pond/wet infiltration to shallow gw storage

   allocate (pnd_solp(mhru))
   allocate (pnd_solpg(mhru))
   allocate (pnd_vol(mhru))
   allocate (pot_fr(mhru))
   allocate (pot_no3(mhru))
   allocate (pot_k(mhru))
   allocate (pot_sed(mhru))
   allocate (pot_san(mhru))
   allocate (pot_sil(mhru))
   allocate (pot_cla(mhru))
   allocate (pot_sag(mhru))
   allocate (pot_lag(mhru))
   allocate (n_reduc(mhru))
   allocate (n_ln(mhru))
   allocate (n_lnco(mhru))

   allocate (pot_vol(mhru))
   allocate (pot_tilemm(mhru))     !!NUBZ
   allocate (pot_volxmm(mhru))
   allocate (potflwi(mhru))
   allocate (potsa(mhru))
   allocate (potsedi(mhru))

   allocate (pplnt(mhru))
   allocate (prf(mch))  !Moriasi 4/8/14
   allocate (spcon(mch))
   allocate (spexp(mch))
   allocate (qdr(mhru))
   allocate (qdayout(mhru))
   allocate (rch_dakm(mxsubch))
   allocate (rchrg(mhru))
   allocate (rchrg_n(mhru))    !! amount of nitrate getting to the shallow aquifer
   allocate (rchrg_dp(mhru))
   allocate (re(mhru))
   allocate (revapmn(mhru))
   allocate (rhd(mhru))
   allocate (rnd2(mhru))
   allocate (rnd3(mhru))
   allocate (rnd8(mhru))
   allocate (rnd9(mhru))
   allocate (rsdin(mhru))
   allocate (rwt(mhru))
   allocate (sci(mhru))
!    Drainmod tile equations  01/2006
   allocate (sdrain(mhru))
   allocate (sstmaxd(mhru))
!    Drainmod tile equations  01/2006
   allocate (seccip(mhru))
   allocate (secciw(mhru))
   allocate (sed_stl(mhru))
   allocate (sedminpa(mhru))
   allocate (sedminps(mhru))
   allocate (sedorgn(mhru))
   allocate (sedorgp(mhru))
   allocate (sedyld(mhru))

   allocate (sanyld(mhru))
   allocate (silyld(mhru))
   allocate (clayld(mhru))
   allocate (sagyld(mhru))
   allocate (lagyld(mhru))
   allocate (sed_con(mhru))
   allocate (sepbtm(mhru))
   allocate (shallirr(mhru))
   allocate (rchrg_src(mhru))
   allocate (shallst(mhru))
   allocate (shallst_n(mhru))
   allocate (slsoil(mhru))
   allocate (slsubbsn(mhru))
   allocate (smx(mhru))
   allocate (sno_hru(mhru))
   allocate (snotmp(mhru))
   allocate (soln_con(mhru))
   allocate (solp_con(mhru))
   allocate (sol_alb(mhru))
   allocate (sol_avbd(mhru))
   allocate (sol_cnsw(mhru))
   allocate (sol_cov(mhru))
   allocate (sol_crk(mhru))
   allocate (sol_nly(mhru))
   allocate (sol_sumfc(mhru))
   allocate (sol_sumul(mhru))
   allocate (sol_sumwp(mhru))
   allocate (sol_sw(mhru))
   allocate (sol_zmx(mhru))
   allocate (strip_n(20,mhru))
!!    Drainmod tile equations  01/2006
   allocate (strip_cn(20,mhru))
   allocate (strip_p(20,mhru))
   allocate (stmaxd(mhru))
   allocate (itill(mhru))
   allocate (strsa(mhru))
   allocate (strsn(mhru))
   allocate (strstmp(mhru))
   allocate (strsw(mhru))
   allocate (stsol_rd(mhru))
   allocate (subp(mhru))
   allocate (sumix(mhru))
   allocate (surfq(mhru))
   allocate (surqno3(mhru))
   allocate (surqsolp(mhru))
   allocate (swtrg(mhru))
   allocate (t_ov(mhru))
   allocate (tconc(mhru))
   allocate (tdrain(mhru))
   allocate (tc_gwat(mhru))
   allocate (terr_cn(20,mhru))
   allocate (terr_p(20,mhru))
   allocate (terr_sl(20,mhru))
   allocate (tile_ttime(mhru))
   allocate (tileq(mhru))
   allocate (tileno3(mhru))
   allocate (tmn(mhru))
   allocate (tmpav(itempa))
   allocate (tmp_hi(mhru))
   allocate (tmp_lo(mhru))
   allocate (tmx(mhru))
   allocate (trapeff(mhru))
   allocate (twash(mhru))
   allocate (u10(mhru))
   allocate (urblu(mhru))
   allocate (usle_cfac(mhru))
   allocate (usle_eifac(mhru))
   allocate (usle_k(mhru))
   allocate (usle_mult(mhru))
   allocate (usle_ls(mhru))
   allocate (usle_p(mhru))
   allocate (velsetlp(mhru))
   allocate (wtab(mhru))
   allocate (wet_chla(mhru))
   allocate (wet_fr(mhru))
   allocate (iwetgw(mhru))
   allocate (iwetile(mhru))
   allocate (wet_k(mhru))
   allocate (wet_mxsa(mhru))
   allocate (wet_mxvol(mhru))
   allocate (wet_no3(mhru))
   allocate (wet_no3g(mhru))
   allocate (wet_no3s(mhru))
   allocate (wet_nsa(mhru))
   allocate (wet_nsed(mhru))
   allocate (wet_nvol(mhru))
   allocate (wet_orgn(mhru))
   allocate (wet_orgp(mhru))
   allocate (wet_psed(mhru))
   allocate (wet_seci(mhru))
   allocate (wet_sed(mhru))
   allocate (wet_solp(mhru))
   allocate (wet_solpg(mhru))
   allocate (wet_vol(mhru))
   allocate (wfsh(mhru))
   allocate (yldaa(mhru))
   allocate (yldanu(mhru))

   allocate (wet_san(mhru))
   allocate (wet_sil(mhru))
   allocate (wet_cla(mhru))
   allocate (wet_sag(mhru))
   allocate (wet_lag(mhru))

   allocate (frad(mhru,nstep))
!      allocate (hhsubp(mhru,24))

   !     allocate (rhrbsb(24))
   allocate (rstpbsb(nstep))
   allocate (rainsub(mhru,nstep))
   allocate (precipdt(nstep+1))

   allocate (bss(4,mhru))
   allocate (nsetlp(2,mhru))
   allocate (nsetlw(2,mhru))
   allocate (psetlp(2,mhru))
   allocate (psetlw(2,mhru))
   allocate (wrt(2,mhru))
   allocate (wgncur(3,mhru))
   allocate (wgnold(3,mhru))
   allocate (surf_bs(17,mhru))
   allocate (rndseed(10,mhru))
   allocate (pcpband(10,mhru))
   allocate (snoeb(10,mhru))
   allocate (orig_snoeb(10,mhru))
   allocate (snotmpeb(10,mhru))
   allocate (tavband(10,mhru))
   allocate (tmxband(10,mhru))
   allocate (wudeep(12,mhru))
   allocate (wupnd(12,mhru))
   allocate (wushal(12,mhru))
!     allocate (phi(13,msub+1))
   allocate (phi(13,mch))
   allocate (wat_phi1(mhru))
   allocate (wat_phi5(mhru))
   allocate (wat_phi6(mhru))
   allocate (wat_phi9(mhru))

!!    arrays which contain data related to pesticides, HRUs
   allocate (orig_pltpst(mpst,mhru))
   allocate (plt_pst(mpst,mhru))
   allocate (pst_enr(mpst,mhru))
   allocate (pst_sed(mpst,mhru))
   allocate (pst_surq(mpst,mhru))
   allocate (zdb(mpst,mhru))

   allocate (pst_lag(3,mpst,mhru))


!!    arrays which contain data related to HRU output
   allocate (hrupsta(4,mpst,mhru))
   allocate (hrupstd(4,mpst,mhru))
   allocate (hrupstm(4,mpst,mhru))
   allocate (hrupsty(4,mpst,mhru))
   allocate (hrumono(74,mhru))
   allocate (hruyro(74,mhru))
   allocate (hruaao(74,mhru))
   allocate (wtrmon(40,mhru))
   allocate (wtryr(40,mhru))
   allocate (wtraa(40,mhru))

!!    arrays which contain data related to pesticides
   allocate (lat_pst(mpst))
   allocate (npno(mpst))
   allocate (pstsol(mpst))
   allocate (wshd_pstap(mpst))
   allocate (wshd_pstdg(mpst))

!!    arrays which contain data related to years
   allocate (flocnst(mrecc))
   allocate (sedcnst(mrecc))
   allocate (orgncnst(mrecc))
   allocate (orgpcnst(mrecc))
   allocate (no3cnst(mrecc))
   allocate (minpcnst(mrecc))
   allocate (nh3cnst(mrecc))
   allocate (no2cnst(mrecc))
   allocate (bactpcnst(mrecc))
   allocate (bactlpcnst(mrecc))
   allocate (cmtlcnst(3,mrecc))
   allocate (chlacnst(mrecc))
   allocate (disoxcnst(mrecc))
   allocate (cbodcnst(mrecc))
   allocate (solpstcnst(mrecc))
   allocate (srbpstcnst(mrecc))

   allocate (floyr(mrecy,myr))
   allocate (sedyr(mrecy,myr))
   allocate (orgnyr(mrecy,myr))
   allocate (orgpyr(mrecy,myr))
   allocate (no3yr(mrecy,myr))
   allocate (minpyr(mrecy,myr))
   allocate (nh3yr(mrecy,myr))
   allocate (no2yr(mrecy,myr))
   allocate (bactpyr(mrecy,myr))
   allocate (bactlpyr(mrecy,myr))
   allocate (cmtlyr(3,mrecy,myr))
   allocate (chlayr(mrecy,myr))
   allocate (disoxyr(mrecy,myr))
   allocate (cbodyr(mrecy,myr))
   allocate (solpstyr(mrecy,myr))
   allocate (srbpstyr(mrecy,myr))

   allocate (flomon(mrecm,myr,12))
   allocate (sedmon(mrecm,myr,12))
   allocate (orgnmon(mrecm,myr,12))
   allocate (orgpmon(mrecm,myr,12))
   allocate (no3mon(mrecm,myr,12))
   allocate (minpmon(mrecm,myr,12))
   allocate (nh3mon(mrecm,myr,12))
   allocate (no2mon(mrecm,myr,12))
   allocate (bactpmon(mrecm,myr,12))
   allocate (bactlpmon(mrecm,myr,12))
   allocate (cmtlmon(3,mrecm,myr,12))
   allocate (chlamon(mrecm,myr,12))
   allocate (disoxmon(mrecm,myr,12))
   allocate (cbodmon(mrecm,myr,12))
   allocate (solpstmon(mrecm,myr,12))
   allocate (srbpstmon(mrecm,myr,12))

!!    arrays
   allocate (halgae(nstep))
   allocate (hbactlp(nstep))
   allocate (hbactp(nstep))
   allocate (hbod(nstep))
   allocate (hchla(nstep))
   allocate (hdepth(nstep))      ! changed as per nstep  !nstep Mar 19,2008
   allocate (hdisox(nstep))
   allocate (hharea(nstep))  ! changed as per nstep  !nstep Mar 19,2008
   allocate (hhqday(nstep))  ! changed as  nstep  Oct. 18, 2007
   allocate (hhstor(nstep))  ! changed as per nstep   !nstep Mar 19,2008
   allocate (hhtime(nstep))  ! changed as per nstep   !nstep Mar 19,2008
   allocate (hnh4(nstep))
   allocate (hno2(nstep))
   allocate (hno3(nstep))
   allocate (horgn(nstep))
   allocate (horgp(nstep))
   allocate (hrchwtr(nstep))  ! changed as per nstep  !nstep Mar 19,2008
   allocate (hrtwtr(nstep))   ! changed as per nstep  !nstep Mar 19,2008
   allocate (hsdti(nstep))    ! changed as per nstep  !nstep Mar 19,2008
   allocate (hsedst(nstep))
   allocate (hsedyld(nstep))
   allocate (hsolp(nstep))
   allocate (hsolpst(nstep))
   allocate (hsorpst(nstep))

   allocate (wpstaao(4,mpst))
   allocate (wpstmono(4,mpst))
   allocate (wpstyro(4,mpst))
   allocate (wpstdayo(4,mpst))


!!arrays that store initial values
   allocate (wattemp(mch))
!! sj, june 07 modifications to carbon balance routines
   allocate (sol_n(mlyr,mhru))
!! sj, june 07 end

!!Armen Jan 2008
   allocate (tillagef(mlyr,mhru))
!  test rtfr
   allocate (rtfr(mlyr))

!!    added for manure Armen Jan 2009
   allocate (sol_mc(mlyr,mhru))
   allocate (sol_mn(mlyr,mhru))
   allocate (sol_mp(mlyr,mhru))



!!Armen Jan 2008 end
!! sj aug 09 SWAT-C MC stuff
   allocate (cf(mhru))
   allocate (cfh(mhru))
   allocate (cfdec(mhru))
!! sj aug 09 end
   allocate (hhsurf_bs(2,mhru,nstep)) 
   allocate (ubnrunoff(nstep),ubntss(nstep))
   allocate (sub_ubnrunoff(msub,nstep),sub_ubntss(msub,nstep))

!! Arrays for subdaily erosion modeling by Jaehak Jeong
   allocate (hhsedy(mhru,nstep),rhy(nstep))
   allocate (snam(mhru),hydgrp(mhru))
   allocate (dratio(msub))
   allocate (sub_subp_dt(msub,nstep),sub_hhsedy(msub,nstep))
   allocate (sub_atmp(msub,nstep),bmp_recharge(msub))
   allocate (rchhr(mrcho,mch,nstep),hrtevp(nstep),hrttlc(nstep))
   allocate (hhresflwo(nstep),hhressedo(nstep))
!! Arrays for bmp simulation by jaehak jeong
   allocate (bmpdrain(mhru))
   allocate (subdr_km(mhyd),subdr_ickm(mhyd),sub_cn2(msub))
   ! sedimentation-filtration
   allocate (num_sf(msub),sf_fr(10,msub),sf_dim(10,msub),&
      &sf_typ(10,msub),sf_im(10,msub),sf_iy(10,msub),sp_sa(10,msub),&
      &sp_pvol(10,msub),sp_pd(10,msub),sp_sedi(10,msub),&
      &sp_sede(10,msub),ft_sa(10,msub),ft_fsa(10,msub),&
      &ft_dep(10,msub),ft_h(10,msub),ft_pd(10,msub),&
      &ft_k(10,msub),ft_dp(10,msub),ft_dc(10,msub),&
      &ft_por(10,msub),tss_den(10,msub),ft_alp(10,msub),&
      &sp_qi(10,msub),sp_k(10,msub),sp_bpw(10,msub),&
      &ft_bpw(10,msub),sp_dp(10,msub),&
      &ft_qfg(10,msub),sp_qfg(10,msub))
   allocate (sub_ha_imp(msub),ft_qpnd(10,msub),ft_qsw(10,msub),&
      &ft_qin(10,msub),ft_qout(10,msub),ft_sedpnd(10,msub),&
      &sf_ptp(10,msub),ft_fc(10,msub),sub_ha_urb(msub))
!! additional var by Ann
!! Filter Strip variable allocation MJW
   allocate (vfscon(mhru))
   allocate (vfsratio(mhru))
   allocate (vfsch(mhru))
   allocate (vfsi(mhru))
   allocate (filter_i(20,mhru))
   allocate (filter_ratio(20,mhru))
   allocate (filter_con(20,mhru))
   allocate (filter_ch(20,mhru))

   ! detention pond
   allocate(dtp_imo(mhyd),dtp_iyr(mhyd),&
      &dtp_numweir(mhyd),dtp_numstage(mhyd),&
      &dtp_stagdis(mhyd),dtp_reltype(mhyd),dtp_onoff(mhyd))

   allocate(dtp_evrsv(msub),dtp_totwrwid(msub),dtp_lwratio(msub),&
      &dtp_intcept(msub),&
      &dtp_expont(msub),dtp_coef(3,msub),dtp_ivol(msub),dtp_ised(msub))

   allocate(dtp_wdratio(10,msub),dtp_depweir(10,msub),&
      &dtp_diaweir(10,msub),dtp_pcpret(10,msub),&
      &dtp_cdis(10,msub),dtp_flowrate(10,msub),&
      &dtp_wrwid(10,msub),dtp_weirtype(10,msub),dtp_weirdim(10,msub),&
      &dtp_addon(10,msub))
   !! additional var by jeong for nutrient speciation
   allocate (lat_orgn(mhru))
   allocate (lat_orgp(mhru))

!! Variables for soil P and additional operations mjw
   allocate (a_days(mlyr,mhru))
   allocate (psp_store(mlyr,mhru))
   allocate (ssp_store(mlyr,mhru))
   allocate (sol_cal(mlyr,mhru))
   allocate (sol_ph(mlyr,mhru))
   allocate (ro_bmp_flag(20,mhru))
   allocate (ro_bmp_flo(20,mhru))
   allocate (ro_bmp_sed(20,mhru))
   allocate (ro_bmp_pp(20,mhru))
   allocate (ro_bmp_sp(20,mhru))
   allocate (ro_bmp_pn(20,mhru))
   allocate (ro_bmp_sn(20,mhru))
   allocate (ro_bmp_bac(20,mhru))

   allocate (ro_bmp_flos(20,mhru))
   allocate (ro_bmp_seds(20,mhru))
   allocate (ro_bmp_pps(20,mhru))
   allocate (ro_bmp_sps(20,mhru))
   allocate (ro_bmp_pns(20,mhru))
   allocate (ro_bmp_sns(20,mhru))
   allocate (ro_bmp_bacs(20,mhru))

   allocate (ro_bmp_flot(20,mhru))
   allocate (ro_bmp_sedt(20,mhru))
   allocate (ro_bmp_ppt(20,mhru))
   allocate (ro_bmp_spt(20,mhru))
   allocate (ro_bmp_pnt(20,mhru))
   allocate (ro_bmp_snt(20,mhru))
   allocate (ro_bmp_bact(20,mhru))

   allocate (bmp_flo(mhru))
   allocate (bmp_sed(mhru))
   allocate (bmp_pp(mhru))
   allocate (bmp_sp(mhru))
   allocate (bmp_pn(mhru))
   allocate (bmp_sn(mhru))
   allocate (bmp_bac(mhru))

   allocate (bmp_flos(mhru))
   allocate (bmp_seds(mhru))
   allocate (bmp_pps(mhru))
   allocate (bmp_sps(mhru))
   allocate (bmp_pns(mhru))
   allocate (bmp_sns(mhru))
   allocate (bmp_bacs(mhru))

   allocate (bmp_flot(mhru))
   allocate (bmp_sedt(mhru))
   allocate (bmp_ppt(mhru))
   allocate (bmp_spt(mhru))
   allocate (bmp_pnt(mhru))
   allocate (bmp_snt(mhru))
   allocate (bmp_bact(mhru))

   !retention irrigation
   allocate(ri_fr(msub,10),ri_dim(msub,10),&
      &ri_im(msub,10),ri_iy(msub,10),ri_sa(msub,10),ri_vol(msub,10),&
      &ri_qi(msub,10),ri_k(msub,10),ri_dd(msub,10),ri_evrsv(msub,10),&
      &ri_dep(msub,10),ri_nirr(msub,30),&
      &num_noirr(msub),ri_totpvol(nstep),ri_luflg(mhru),&
      &ri_subkm(msub),ri_pumpv(msub,10),ri_sedi(msub,10))
   allocate(num_ri(msub), ri_pmpvol(10,nstep),hrnopcp(msub,0:nstep),&
      &ri_qloss(10,nstep))

   !wet pond
   allocate(wtp_onoff(mhyd),wtp_imo(mhyd),&
      &wtp_iyr(mhyd),wtp_dim(mhyd),wtp_stagdis(mhyd),wtp_sdtype(mhyd),&
      &wtp_pvol(mhyd),wtp_pdepth(mhyd),wtp_sdslope(mhyd),&
      &wtp_lenwdth(mhyd),wtp_extdepth(mhyd),wtp_hydeff(mhyd),&
      &wtp_evrsv(mhyd),wtp_sdintc(mhyd),wtp_sdexp(mhyd),&
      &wtp_pdia(mhyd),wtp_plen(mhyd),&
      &wtp_pmann(mhyd),wtp_ploss(mhyd),wtp_k(mhyd),&
      &wtp_dp(mhyd),wtp_sedi(mhyd),wtp_sede(mhyd),wtp_qi(mhyd))
    allocate(wtp_sdc(3,mhyd))

!!    LID simulations
!!    Common variable
!!    van Genuchten equation's coefficients
   allocate(lid_cuminf_last(4,mhru),lid_sw_last(4,mhru),&
      &lid_f_last(4,mhru),lid_cumr_last(4,mhru),&
      &lid_str_last(4,mhru),lid_farea(4,mhru),lid_qsurf(mhru,4),&
      &lid_sw_add(4,mhru),lid_cumqperc_last(4,mhru),&
      &lid_excum_last(4,mhru))    !!  nbs
   allocate(lid_cumirr_last(mhru))

!!    Green Roof
   allocate(gr_onoff(msub,mudb),&
      &gr_farea(msub,mudb),gr_solop(msub,mudb),gr_etcoef(msub,mudb),&
      &gr_fc(msub,mudb),gr_wp(msub,mudb),gr_ksat(msub,mudb),&
      &gr_por(msub,mudb),gr_hydeff(msub,mudb),gr_soldpt(msub,mudb))

!!    Rain Garden
   allocate(rg_onoff(msub,mudb),&
      &rg_farea(msub,mudb),rg_solop(msub,mudb),rg_etcoef(msub,mudb),&
      &rg_fc(msub,mudb),rg_wp(msub,mudb),rg_ksat(msub,mudb),&
      &rg_por(msub,mudb),rg_hydeff(msub,mudb),rg_soldpt(msub,mudb),&
      &rg_dimop(msub,mudb),rg_sarea(msub,mudb),rg_vol(msub,mudb),&
      &rg_sth(msub,mudb),rg_sdia(msub,mudb),rg_bdia(msub,mudb),&
      &rg_sts(msub,mudb),rg_orifice(msub,mudb),rg_oheight(msub,mudb),&
      &rg_odia(msub,mudb))

!!    CiStern
   allocate(cs_onoff(msub,mudb),&
      &cs_grcon(msub,mudb),cs_farea(msub,mudb),cs_vol(msub,mudb),&
      &cs_rdepth(msub,mudb))

!!    Poropus paVement
   allocate(pv_onoff(msub,mudb),&
      &pv_grvdep(msub,mudb),pv_grvpor(msub,mudb),pv_farea(msub,mudb),&
      &pv_solop(msub,mudb),pv_drcoef(msub,mudb),pv_fc(msub,mudb),&
      &pv_wp(msub,mudb),pv_ksat(msub,mudb),pv_por(msub,mudb),&
      &pv_hydeff(msub,mudb),pv_soldpt(msub,mudb))

!!    LID general
   allocate(lid_onoff(msub,mudb))

   !! By Zhang for C/N cycling
   !! ============================
   allocate(sol_BMC(mlyr,mhru))
   allocate(sol_BMN(mlyr,mhru))
   allocate(sol_HSC(mlyr,mhru))
   allocate(sol_HSN(mlyr,mhru))
   allocate(sol_HPC(mlyr,mhru))
   allocate(sol_HPN(mlyr,mhru))
   allocate(sol_LM(mlyr,mhru))
   allocate(sol_LMC(mlyr,mhru))
   allocate(sol_LMN(mlyr,mhru))
   allocate(sol_LS(mlyr,mhru))
   allocate(sol_LSC(mlyr,mhru))
   allocate(sol_LSN(mlyr,mhru))
   allocate(sol_LSL(mlyr,mhru))
   allocate(sol_LSLC(mlyr,mhru))
   allocate(sol_LSLNC(mlyr,mhru))
   allocate(sol_WOC(mlyr,mhru))
   allocate(sol_WON(mlyr,mhru))

   !!for print out at daily, monthly, and annual scale
   allocate(sedc_d(mhru))
   allocate(surfqc_d(mhru))
   allocate(latc_d(mhru))
   allocate(percc_d(mhru))
   allocate(NPPC_d(mhru))
   allocate(rsdc_d(mhru))
   allocate(grainc_d(mhru))
   allocate(stoverc_d(mhru))
   allocate(emitc_d(mhru))
   allocate(rspc_d(mhru))

   !Tillage factor on SOM decomposition
   allocate(tillage_switch(mhru))
   allocate(tillage_depth(mhru))
   allocate(tillage_days(mhru))
   allocate(tillage_factor(mhru))
   tillage_switch = 0
   tillage_depth = 0.
   tillage_days = 0
   tillage_factor = 0.
   !! By Zhang for C/N cycling
   !! ============================

   !FLOOD ROUTING
   allocate(QHY(nstep+1,mhyd,4), NHY(4*msub))
   allocate(RCHX(msub),RCSS(msub),QCAP(msub),CHXA(msub),CHXP(msub))

   call zero0
   call zero1
   call zero2
   call zeroini
   call zero_urbn

   return
end
