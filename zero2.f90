subroutine zero2

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine zeros all array values

   use parm
   implicit none

   hrumono = 0.
   wtrmon = 0.
   submono = 0.
   rchmono = 0.
   resoutm = 0.
   yldaa = 0.
   bio_aams = 0.
   lai_aamx = 0.

   aairr = 0.
   algae = 0.
   anion_excl = 0.
   alpha_bnk = 0.
   alpha_bnke = 0.
   fr_curb = 0.
   bactlpcnst = 0.
   bactlpmon = 0.
   bactlpq = 0.
   bactlps = 0.
   bactlpyr = 0.
   bactpcnst = 0.
   bactpmon = 0.
   bactpq = 0.
   bactps = 0.
   bactpyr = 0.
   bankst = 0.
   bio_aahv = 0.
   biomix = 0.
   bp1 = 0.
   bp2 = 0.
   br1 = 0.
   br2 = 0.
   bw1 = 0.
   bw2 = 0.
   ch_l1 = 0.
   ch_l2 = 0.
   ch_revap = 0.
   chlap = 0.
   chlar = 0.
   chlaw = 0.
   cht = 0.
   cklsp = 0.
   cmtl1cnst = 0.
   cmtl1mon = 0.
   cmtl1yr = 0.
   cmtl2cnst = 0.
   cmtl2mon = 0.
   cmtl2yr = 0.
   cmtl3cnst = 0.
   cmtl3mon = 0.
   cmtl3yr = 0.
   crdep = 0.
   elevp = 0
   elevt = 0
   erorgn = 0.
   erorgp = 0.
   fcimp = 0.
   flocnst = 0.
   flowfr = 0.
   floyr = 0.
   flwin = 0.
   flwout = 0.
   fsred = 0.
   gwminp = 0.
   gwno3 = 0.
   harveff = 0.
   hhvaroute = 0.
   hru_sub = 0
   hru1 = 0
   hru_seq = 0
   hrutot = 0
   hvstiadj = 0.
   hydgrp = ""
   icanal = 0
   iflod1r = 0
   iflod2r = 0
   ihgage = 0
   ipnd1 = 0
   ipnd2 = 0

   imp_trig = 0   !!Srini pothole
   pot_orgn = 0.
   pot_orgp = 0.
   pot_mps = 0.
   pot_mpa = 0.
   pot_no3 = 0.
   pot_solp = 0.
   pot_evap = 0.
   tile_out = 0.
   tile_sedo = 0.
   tile_no3o = 0.
   tile_solpo = 0.
   tile_orgno = 0.
   tile_orgpo = 0.
   tile_minpso = 0.
   tile_minpao = 0.

   ! irelease = 0 ! not defined
   imp_trig = 1
   ires1 = 0
   ires2 = 0
   iresco = 0
   isgage = 0
   iwgage = 0
   iyres = 0
   laimxfr = 0.
   ldrain = 0
   lkpst_conc = 0.
   lkpst_koc = 0.
   lkpst_mix = 0.
   lkpst_rea = 0.
   lkpst_rsp = 0.
   lkpst_stl = 0.
   lkpst_vol = 0.
   lkspst_act = 0.
   lkspst_bry = 0.
   lkspst_conc = 0.
   lkspst_rea = 0.
   minpcnst = 0.
   minpmon = 0.
   minpyr = 0.
   mores = 0
   ndeat = 0
   ndtargr = 0
   nh3cnst = 0.
   nh3mon = 0.
   nh3yr = 0.
   no2cnst = 0.
   no2mon = 0.
   no2yr = 0.
   no3cnst = 0.
   no3mon = 0.
   no3yr = 0.
   nsetlp = 0.
   nsetlr = 0.
   nsetlw = 0.
   oflowmn = 0.
   oflowmx = 0.
   oflowmn_fps = 0.
   olai = 0.
   orgncnst = 0.
   orgnmon = 0.
   orgnyr = 0.
   orgpcnst = 0.
   orgpmon = 0.
   orgpyr = 0.
   orig_lkpstconc = 0.
   orig_lkspstconc = 0.
   orig_pndno3 = 0.
   orig_pndorgn = 0.
   orig_pndorgp = 0.
   orig_pndsolp = 0.
   orig_resnh3 = 0.
   orig_resno2 = 0.
   orig_resno3 = 0.
   orig_resorgn = 0.
   orig_resorgp = 0.
   orig_ressed = 0.
   orig_ressolp = 0.
   orig_resvol = 0.
   orig_solst = 0.
   orig_solsw = 0.
   orig_soltmp = 0.
   orig_volcr = 0.
   orig_wetno3 = 0.
   orig_wetorgn = 0.
   orig_wetorgp = 0.
   orig_wetsolp = 0.
   ovrlnd = 0.
   pcf = 1.
   phi = 0.
   plt_et = 0.
   plt_pet = 0.
   plantn = 0.
   plantp = 0.
   pltfr_n = 0.
   pltfr_p = 0.
   pnd_chla = 0.
   pnd_no3 = 0.
   pnd_no3g = 0.
   pnd_no3s = 0.
   pnd_orgn = 0.
   pnd_orgp = 0.
   pnd_seci = 0.
   pnd_psed = 0.
   pnd_solp = 0.
   pnd_solpg = 0.
   pot_volx = 0.
   pot_volxmm = 0.
   potflwi = 0.
   potsedi = 0.
   potsani = 0.
   potsili = 0.
   potclai = 0.
   potsagi = 0.
   potlagi = 0.

   psetlp = 0.
   psetlr = 0.
   psetlw = 0.
   rchstor = 0.
   rch_bactp = 0.
   rch_bactlp = 0.
   res_bactlp = 0.
   res_bactp = 0.
   res_chla = 0.
   res_esa = 0.
   res_evol = 0.
   res_k = 0.
   res_nh3 = 0.
   res_no2 = 0.
   res_no3 = 0.
   res_nsed = 0.
   res_orgn = 0.
   res_orgp = 0.
   res_out = 0.
   res_psa = 0.
   res_pvol = 0.
   res_rr = 0.
   res_seci = 0.
   res_sed = 0.
   res_solp = 0.
   res_sub = 0
   res_vol = 0.
   resdata = 0.
   rnd2 = 0.
   rnd3 = 0.
   rnd8 = 0.
   rnd9 = 0.
   rndseed = 0
   rsdco_pl = 0.
   rwt = 0.
   sci = 0.
   seccip = 0.
   seccir = 0.
   secciw = 0.
   sed_stl = 0.
   sed_stlr = 0.
   sedcnst = 0.
   sedmon = 0.
   sedyr = 0.

   sedyld = 0.
   sanyld = 0.
   silyld = 0.
   clayld = 0.
   sagyld = 0.
   lagyld = 0.

   shallirr = 0.
   shyd = 0.
   smx = 0.
   snotmp = 0.
   snotmpeb = 0.
   sol_avbd = 0.
   sol_avpor = 0.
   sol_awc = 0.
   sol_fc = 0.
   sol_hk = 0.
   sol_st = 0.
   sol_sumfc = 0.
   sol_sumul = 0.
   sol_sw = 0.
   sol_tmp = 0.
   sol_ul = 0.
   sol_wpmm = 0.
   starg = 0.
   starg_fps = 0.
   sub_tc = 0.
   surf_bs = 0.
   tmp_hi = 0.
   tmp_lo = 0.
   tmpavp = 0.
   twash = 0.
   urblu = 0
   urbname = ""
   usle_mult = 0.
   values = 0
   varoute = 0.
   vartran = 0.
   volcr = 0.
   welev = 0.
   wet_chla = 0.
   wet_no3 = 0.
   wet_no3g = 0.
   wet_no3s = 0.
   wet_orgn = 0.
   wet_orgp = 0.
   wet_psed = 0.
   wet_seci = 0.
   wet_solp = 0.
   wet_solpg = 0.
   wgnold = 0.
   wlat = 0.
   wpstaao = 0.
   wpstmono = 0.
   wpstyro = 0.
   wrt = 0.
   wshd_pstdg = 0.
   wshdaao = 0.
   wshdmono = 0.
   wshdyro = 0.
   wtraa = 0.
   wtryr = 0.
   wupnd = 0.
   wuresn = 0.
   wurtnf = 0.
   yldn = 0.
   zdb = 0.

   !! MJW
   sol_P_model = 0

   bmp_flag = 0
   !! surface
   bmp_flo = 1.      !! Surface Flow
   bmp_sed = 1.      !! Sediment
   bmp_pp = 1.       !! Particulate P
   bmp_sp = 1.       !! Soluble P
   bmp_pn =  1.      !! Particulate N
   bmp_sn = 1.       !! Soluble N
   bmp_bac = 1.      !! Bacteria
   !! subsurface
   bmp_flos = 1.      !! Subsurface Flow
   bmp_seds = 1.      !! Sediment
   bmp_pps = 1.       !! Particulate P
   bmp_sps = 1.       !! Soluble P
   bmp_pns =  1.      !! Particulate N
   bmp_sns = 1.       !! Soluble N
   bmp_bacs = 1.      !! Bacteria
   !! tile
   bmp_flot = 1.      !! Tile Flow
   bmp_sedt = 1.      !! Sediment
   bmp_ppt = 1.       !! Particulate P
   bmp_spt = 1.       !! Soluble P
   bmp_pnt =  1.      !! Particulate N
   bmp_snt = 1.       !! Soluble N
   bmp_bact = 1.      !! Bacteria

   ro_bmp_flag = 0    !! Flag to turn on or off user BMP

   !! surface runoff removal efficiency
   ro_bmp_flo = 0.    !! Flow
   ro_bmp_sed = 0.    !! Sediment
   ro_bmp_pp = 0.     !! Particulate P
   ro_bmp_sp = 0.     !! Soluble P
   ro_bmp_pn = 0.     !! Particulate N
   ro_bmp_sn = 0.     !! Soluble N
   ro_bmp_bac = 0.    !! Bacteria
   !! subsurface - lateral soil and groundwater
   ro_bmp_flos = 0.   !! Flow
   ro_bmp_seds = 0.   !! Sediment
   ro_bmp_pps = 0.    !! Particulate P
   ro_bmp_sps = 0.    !! Soluble P
   ro_bmp_pns = 0.    !! Particulate N
   ro_bmp_sns = 0.    !! Soluble N
   ro_bmp_bacs = 0.   !! Bacteria
   !! tile flow removal efficiency
   ro_bmp_flot = 0.   !! Flow
   ro_bmp_sedt = 0.   !! Sediment
   ro_bmp_ppt = 0.    !! Particulate P
   ro_bmp_spt = 0.    !! Soluble P
   ro_bmp_pnt = 0.    !! Particulate N
   ro_bmp_snt = 0.    !! Soluble N
   ro_bmp_bact = 0.   !! Bacteria

   ssp_store = 0.
   psp_store = 0.
   a_days = 0
   b_days = 0
   sol_ph = 0
   sol_cal = 0
   bio_init = 0
   lai_init = 0


   return
end
