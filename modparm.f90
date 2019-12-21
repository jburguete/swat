!> @file modparm.f90
!> file containing the module parm
!> @author
!> modified by Javier Burguete Tolosa

!> main module containing the global variables
module parm

!> max number of variables routed through the reach
   integer, parameter :: mvaro = 33
!> maximum number of variables written to HRU output file (output.hru) (none)
   integer, parameter :: mhruo = 79
!> maximum number of variables written to reach output file (.rch) (none)
   integer, parameter :: mrcho = 62
!> maximum number of variables written to subbasin output file (output.sub)
!> (none)
   integer, parameter :: msubo = 24
!> max number of variables summarized in output.std
   integer, parameter :: mstdo = 113
   integer, parameter :: motot = 600 !! (50 years limit)

!> SWAT program header string (name and version)
   character(len=80), parameter :: prog = "SWAT Sep 7    VER 2018/Rev 670"

!> column headers for HRU output file
   character(len=13), dimension(mhruo), parameter ::&
      &heds = (/"  PRECIPmm"," SNOFALLmm"," SNOMELTmm","     IRRmm",&
      &"     PETmm","      ETmm"," SW_INITmm","  SW_ENDmm",&
      &"    PERCmm"," GW_RCHGmm"," DA_RCHGmm","   REVAPmm",&
      &"  SA_IRRmm","  DA_IRRmm","   SA_STmm","   DA_STmm",&
      &"SURQ_GENmm","SURQ_CNTmm","   TLOSSmm"," LATQGENmm",&
      &"    GW_Qmm","    WYLDmm","   DAILYCN"," TMP_AVdgC",&
      &" TMP_MXdgC"," TMP_MNdgC","SOL_TMPdgC","SOLARMJ/m2",&
      &"  SYLDt/ha","  USLEt/ha","N_APPkg/ha","P_APPkg/ha",&
      &"NAUTOkg/ha","PAUTOkg/ha"," NGRZkg/ha"," PGRZkg/ha",&
      &"NCFRTkg/ha","PCFRTkg/ha","NRAINkg/ha"," NFIXkg/ha",&
      &" F-MNkg/ha"," A-MNkg/ha"," A-SNkg/ha"," F-MPkg/ha",&
      &"AO-LPkg/ha"," L-APkg/ha"," A-SPkg/ha"," DNITkg/ha",&
      &"  NUPkg/ha","  PUPkg/ha"," ORGNkg/ha"," ORGPkg/ha",&
      &" SEDPkg/ha","NSURQkg/ha","NLATQkg/ha"," NO3Lkg/ha",&
      &"NO3GWkg/ha"," SOLPkg/ha"," P_GWkg/ha","    W_STRS",&
      &"  TMP_STRS","    N_STRS","    P_STRS","  BIOMt/ha",&
      &"       LAI","   YLDt/ha","  BACTPct ","  BACTLPct",&
      &" WTAB CLIm"," WTAB SOLm","     SNOmm"," CMUPkg/ha",&
      &"CMTOTkg/ha","   QTILEmm"," TNO3kg/ha"," LNO3kg/ha",&
      &"  GW_Q_Dmm"," LATQCNTmm"," TVAPkg/ha"/)

!> space number for beginning of column in HRU output file (none)
   integer, dimension (mhruo), parameter ::&
      &icols = (/43,53,63,73,83,93,103,113,123,133,143,153,&
      &163,173,183,193,203,213,223,233,243,253,263,273,283,&
      &293,303,313,323,333,343,353,363,373,383,393,403,413,&
      &423,433,443,453,463,473,483,493,503,513,523,533,543,&
      &553,563,573,583,593,603,613,623,633,643,653,663,673,&
      &683,693,703,713,723,733,743,753,763,773,783,793,803,&
      &813,823/)

!> column headers for subbasin output file
   character(len=13), dimension(msubo), parameter ::&
      &hedb = (/"  PRECIPmm"," SNOMELTmm","     PETmm","      ETmm",&
      &"      SWmm","    PERCmm","    SURQmm","    GW_Qmm",&
      &"    WYLDmm","  SYLDt/ha"," ORGNkg/ha"," ORGPkg/ha",&
      &"NSURQkg/ha"," SOLPkg/ha"," SEDPkg/ha"," LAT Q(mm)",&
      &"LATNO3kg/h","GWNO3kg/ha","CHOLAmic/L","CBODU mg/L",&
      &" DOXQ mg/L"," TNO3kg/ha","   QTILEmm"," TVAPkg/ha"/)

!> space number for beginning of column in subbasin output file (none)
   integer, dimension (msubo), parameter ::&
      &icolb = (/35,45,55,65,75,85,95,105,115,125,135,145,&
      &155,165,175,185,195,205,215,225,235,245,255,265/)

!> column headers for reach output file
   character(len=13), dimension(mrcho), parameter ::&
      &hedr = (/"  FLOW_INcms"," FLOW_OUTcms","     EVAPcms",&
      &"    TLOSScms","  SED_INtons"," SED_OUTtons",&
      &" SEDCONCmg/L","   ORGN_INkg","  ORGN_OUTkg",&
      &"   ORGP_INkg","  ORGP_OUTkg","    NO3_INkg",&
      &"   NO3_OUTkg","    NH4_INkg","   NH4_OUTkg",&
      &"    NO2_INkg","   NO2_OUTkg","   MINP_INkg",&
      &"  MINP_OUTkg","   CHLA_INkg","  CHLA_OUTkg",&
      &"   CBOD_INkg","  CBOD_OUTkg","  DISOX_INkg",&
      &" DISOX_OUTkg"," SOLPST_INmg","SOLPST_OUTmg",&
      &" SORPST_INmg","SORPST_OUTmg","  REACTPSTmg",&
      &"    VOLPSTmg","  SETTLPSTmg","RESUSP_PSTmg",&
      &"DIFFUSEPSTmg","REACBEDPSTmg","   BURYPSTmg",&
      &"   BED_PSTmg"," BACTP_OUTct","BACTLP_OUTct",&
      &"  CMETAL#1kg","  CMETAL#2kg","  CMETAL#3kg",&
      &"     TOT Nkg","     TOT Pkg"," NO3ConcMg/l",&
      &"    WTMPdegc","            ","            ",&
      &"            ","            ","            ",&
      &"            ","            ","            ",&
      &"            ","            ","            ",&
      &"            ","            ","            ",&
      &"            ","            "/)

!> space number for beginning of column in reach output file (none)
   integer, dimension (mrcho), parameter ::&
      &icolr = (/38,50,62,74,86,98,110,122,134,146,158,170,182,194,206,&
      &218,230,242,254,266,278,290,302,314,326,338,350,362,374,386,398,&
      &410,422,434,446,458,470,482,494,506,518,530,542,554,566,578,590,&
      &602,614,626,638,650,662,674,686,698,710,722,734,746,758,770/)

!> column headers for reservoir output file
   character(len=13), dimension(41), parameter ::&
      &hedrsv = (/"    VOLUMEm3","  FLOW_INcms"," FLOW_OUTcms",&
      &"    PRECIPm3","      EVAPm3","   SEEPAGEm3",&
      &"  SED_INtons"," SED_OUTtons"," SED_CONCppm",&
      &"   ORGN_INkg","  ORGN_OUTkg"," RES_ORGNppm",&
      &"   ORGP_INkg","  ORGP_OUTkg"," RES_ORGPppm",&
      &"    NO3_INkg","   NO3_OUTkg","  RES_NO3ppm",&
      &"    NO2_INkg","   NO2_OUTkg","  RES_NO2ppm",&
      &"    NH3_INkg","   NH3_OUTkg","  RES_NH3ppm",&
      &"   MINP_INkg","  MINP_OUTkg"," RES_MINPppm",&
      &"   CHLA_INkg","  CHLA_OUTkg","SECCHIDEPTHm",&
      &"   PEST_INmg","  REACTPSTmg","    VOLPSTmg",&
      &"  SETTLPSTmg","RESUSP_PSTmg","DIFFUSEPSTmg",&
      &"REACBEDPSTmg","   BURYPSTmg","  PEST_OUTmg",&
      &"PSTCNCWmg/m3","PSTCNCBmg/m3"/)

!> space number for beginning of column in reservoir output file (none)
   integer, dimension (41), parameter ::&
      &icolrsv = (/38,50,62,74,86,98,110,122,134,146,158,170,182,194,&
      &206,218,230,242,254,266,278,290,302,314,326,338,350,362,374,386,&
      &398,410,422,434,446,458,470,482,494,506,518/)

!> column headers for HRU impoundment output file
   character(len=13), dimension(40), parameter ::&
      &hedwtr = (/"  PNDPCPmm","  PND_INmm","PSED_It/ha","  PNDEVPmm",&
      &"  PNDSEPmm"," PND_OUTmm","PSED_Ot/ha"," PNDVOLm^3",&
      &"PNDORGNppm"," PNDNO3ppm","PNDORGPppm","PNDMINPppm",&
      &"PNDCHLAppm","  PNDSECIm","  WETPCPmm","  WET_INmm",&
      &"WSED_It/ha","  WETEVPmm","  WETSEPmm"," WET_OUTmm",&
      &"WSED_Ot/ha"," WETVOLm^3","WETORGNppm"," WETNO3ppm",&
      &"WETORGPppm","WETMINPppm","WETCHLAppm","  WETSECIm",&
      &"  POTPCPmm","  POT_INmm","OSED_It/ha","  POTEVPmm",&
      &"  POTSEPmm"," POT_OUTmm","OSED_Ot/ha"," POTVOLm^3",&
      &"  POT_SAha","HRU_SURQmm","PLANT_ETmm"," SOIL_ETmm"/)

!> forecast region, subbasin, HRU, reach, reservoir or file number (none)
   integer :: i
   integer icalen
!> Basinwide peak rate adjustment factor for sediment routing in the channel.
!> Allows impact of peak flow rate on sediment routing and channel reshaping to
!> be taken into account.
   real*8 :: prf_bsn

!!    srin - co2 (EPA)
   real*8 :: co2_x2, co2_x

   real*8, dimension (:), allocatable :: alph_e
!> denitrification exponential rate coefficient
   real*8, dimension (:), allocatable :: cdn
!> nitrate percolation coefficient (0-1)\n
!> 0:concentration of nitrate in surface runoff is zero\n
!> 1:percolate has same concentration of nitrate as surface runoff
   real*8, dimension (:), allocatable :: nperco
!> Surface runoff lag time. This parameter is needed in subbasins where the time
!> of concentration is greater than 1 day. SURLAG is used to create a "storage"
!> for surface runoff to allow the runoff to take longer than 1 day to reach the
!> subbasin outlet (days)
   real*8, dimension (:), allocatable :: surlag
   real*8, dimension (:), allocatable :: co_p
!> rate factor for humus mineralization on active organic N
   real*8, dimension (:), allocatable :: cmn
!> Phosphorus soil partitioning coefficient. Ratio of soluble phosphorus in
!> surface layer to soluble phosphorus in runoff
   real*8, dimension (:), allocatable :: phoskd
!> Phosphorus availibility index. The fraction of fertilizer P remaining in
!> labile pool after initial rapid phase of P sorption (none)
   real*8, dimension (:), allocatable :: psp
!> denitrification threshold: fraction of field capacity triggering
!> denitrification
   real*8, dimension (:), allocatable :: sdnco

!!   change per JGA 8/31/2011 gsm for output.mgt
!> basinwide retention parameter adjustment factor (greater than 1)
   real*8 :: r2adj_bsn
!> amount of pesticide applied to HRU (kg/ha)
   real*8 :: pst_kg
   real*8 :: yield, burn_frlb
   real*8 :: yieldgrn, yieldbms, yieldtbr, yieldn, yieldp
   real*8 :: hi_bms, hi_rsd, yieldrsd
!!    arrays for Landscape Transport Capacity 5/28/2009 nadia
   real*8, dimension (:), allocatable :: l_k1, l_k2, l_lambda, l_beta
   real*8, dimension (:), allocatable :: l_gama, l_harea, l_vleng
   real*8, dimension (:), allocatable :: l_vslope, l_ktc

!!    arrays for Biofilm variables
   real*8, dimension (:), allocatable :: biofilm_mumax, biofilm_kinv
   real*8, dimension (:), allocatable :: biofilm_klw, biofilm_kla
   real*8, dimension (:), allocatable :: biofilm_cdet, biofilm_bm


!!    new arrays for routing units
   real*8, dimension (:,:), allocatable :: hru_rufr, daru_km, ru_k
   real*8, dimension (:,:), allocatable :: ru_c, ru_eiq, ru_ovsl, ru_a
   real*8, dimension (:,:), allocatable :: ru_ovs, ru_ktc
   real*8, dimension (:), allocatable :: gwq_ru, qdayout
   integer, dimension (:), allocatable :: ils2, ils2flag
!> pesticide identification number from pest.dat (none)
   integer :: ipest
   integer :: iru, mru, irch, isub, mhyd_bsn, ils_nofig
   integer :: mhru1
   integer, dimension (:), allocatable :: mhyd1 , irtun

!! septic variables for output.std
   real*8 :: wshd_sepno3, wshd_sepnh3, wshd_seporgn, wshd_sepfon
   real*8 :: wshd_seporgp, wshd_sepfop, wshd_sepsolp, wshd_sepbod
   real*8 :: wshd_sepmm
   integer, dimension (:), allocatable :: isep_hru
!! septic variables for output.std
   real*8 :: fixco !< nitrogen fixation coefficient
   real*8 :: nfixmx !< maximum daily n-fixation (kg/ha) 
   real*8 :: res_stlr_co !< reservoir sediment settling coefficient
   real*8 :: rsd_covco !< residue cover factor for computing frac of cover
   real*8 :: vcrit !< critical velocity
!> average amount of water stored in snow at the beginning of the simulation for
!> the entire watershed (mm H20)
   real*8 :: wshd_snob
!> average amount of water stored in soil for the entire watershed (mm H2O)
   real*8 :: wshd_sw
!> fraction of watershed area which drains into ponds (none)
   real*8 :: wshd_pndfr
!> total amount of suspended sediment in ponds in the watershed (metric tons)
   real*8 :: wshd_pndsed
!> total volume of water in ponds in the watershed (m^3)
   real*8 :: wshd_pndv
!> pesticide percolation coefficient (0-1)\n
!> 0: concentration of pesticide in surface runoff is zero\n
!> 1: percolate has same concentration of pesticide as surface runoff
   real*8 :: percop
!> fraction of watershed area that drains into reservoirs (none)
   real*8 :: wshd_resfr
!> watershed area in hectares which drains into ponds (ha)
   real*8 :: wshd_pndha
!> watershed area in hectares which drains into reservoirs (ha)
   real*8 :: wshd_resha
!> fraction of watershed area which drains into wetlands (none)
   real*8 :: wshd_wetfr
   real*8 :: wshd_fminp, wshd_ftotn, wshd_fnh3, wshd_fno3, wshd_forgn
   real*8 :: wshd_forgp, wshd_ftotp, wshd_yldn, wshd_yldp, wshd_fixn
   real*8 :: wshd_pup, wshd_wstrs, wshd_nstrs, wshd_pstrs, wshd_tstrs
   real*8 :: wshd_astrs
!> initial soil water content expressed as a fraction of field capacity
   real*8 :: ffcb
   real*8 :: wshd_hmn, wshd_rwn, wshd_hmp, wshd_rmn, wshd_dnit
!> die-off factor for persistent bacteria in soil solution (1/day)
   real*8 :: wdpq
   real*8 :: wshd_rmp, wshd_voln, wshd_nitn, wshd_pas, wshd_pal
!> wash off fraction for persistent bacteria on foliage during a rainfall event
   real*8 :: wof_p
   real*8 :: wshd_plch, wshd_raino3, ressedc, basno3f, basorgnf
   real*8 :: wshd_pinlet, wshd_ptile
   real*8 :: sftmp !< Snowfall temperature (deg C)
!> Minimum melt rate for snow during year (Dec. 21) where deg C refers to the
!> air temperature. (mm/deg C/day)
   real*8 :: smfmn
!> Maximum melt rate for snow during year (June 21) where deg C refers to the
!> air temperature. SMFMX and SMFMN allow the rate of snow melt to vary through
!> the year. These parameters are accounting for the impact of soil temperature
!> on snow melt. (mm/deg C/day)
   real*8 :: smfmx
!> Snow melt base temperature. Mean air temperature at which snow melt will
!> occur. (deg C)
   real*8 :: smtmp
!> growth factor for persistent bacteria in soil solution (1/day)
   real*8 :: wgpq
   real*8 :: basminpf, basorgpf
!> die-off factor for less persistent bacteria in soil solution (1/day)
   real*8 :: wdlpq
!> total amount of suspended sediment in reservoirs in the watershed
!> (metric tons)
   real*8 :: wshd_ressed
!> total volume of water in all reservoirs in the watershed (m^3)
   real*8 :: wshd_resv
!> average amount of phosphorus initially in the mineral P pool in watershed
!> soil (kg P/ha)
   real*8 :: basminpi
!> average amount of nitrogen initially in the nitrate pool in watershed soil
!> (kg N/ha)
   real*8 :: basno3i
!> average amount of nitrogen initially in the organic N pool in watershed soil
!> (kg N/ha)
   real*8 :: basorgni
!> die-off factor for persistent bacteria adsorbed to soil particles (1/day)
   real*8 :: wdps
!> growth factor for less persistent bacteria in soil solution (1/day)
   real*8 :: wglpq
!> average amount of phosphorus initially in the organic P pool in watershed
!> soil (kg P/ha)
   real*8 :: basorgpi
   real*8 :: peakr !< peak runoff rate (m^3/s)
!> albedo, the fraction of the solar radiation reflected at the soil surface
!> back into space (none)
   real*8 :: albday
   real*8 :: pndsedin, sw_excess
!> Snow pack temperature lag factor (0-1)\n
!> 1 = no lag (snow pack temp=current day air temp) as the lag factor goes to
!> zero, the snow pack's temperature will be less influenced by the current
!> day's air temperature
   real*8 :: timp
   real*8 :: wtabelo, tilep, wt_shall
   real*8 :: sq_rto
   real*8 :: qtile !< drainage tile flow in soil layer for the day (mm H2O)
!> amount of precipitation that infiltrates into soil (enters soil) (mm H2O)
   real*8 :: inflpcp
   real*8 :: tloss, snomlt, snofall, fixn, crk, latlyr
   real*8 :: pndloss, wetloss,potloss, lpndloss, lwetloss
   real*8 :: sedrch, fertn, sol_rd, cfertn, cfertp, sepday, bioday
   real*8 :: sepcrk, sepcrktot, fertno3, fertnh3, fertorgn, fertsolp
   real*8 :: fertorgp
!> growth factor for persistent bacteria adsorbed to soil particles (1/day)
   real*8 :: wgps
   real*8 :: qdfr !< fraction of water yield that is surface runoff (none)
   real*8 :: fertp, grazn, grazp, soxy, sdti, rtwtr, ressa
!> die-off factor for less persistent bacteria absorbed to soil particles
!> (1/day)
   real*8 :: wdlps
!> growth factor for less persistent bacteria adsorbed to soil particles
!> (1/day)
   real*8 :: wglps
   real*8 :: da_km !< area of the watershed in square kilometers (km^2)
   real*8 :: rttime, rchdep, rtevp, rttlc, resflwi
!> die-off factor for persistent bacteria in streams (1/day)
   real*8 :: wdprch
   real*8 :: resflwo, respcp, resev, ressep,ressedi,ressedo,dtot
!> phosphorus percolation coefficient. Ratio of soluble phosphorus in surface to
!> soluble phosphorus in percolate
   real*8 :: pperco_bsn
!> basin nitrate percolation coefficient (0-1)\n
!> 0:concentration of nitrate in surface runoff is zero\n
!> 1:percolate has same concentration of nitrate as surface runoff
   real*8 :: nperco_bsn
!> residue decomposition coefficient. The fraction of residue which will
!> decompose in a day assuming optimal moisture, temperature, C:N ratio, and C:P
!> ratio
   real*8 :: rsdco
   real*8 :: phoskd_bsn,voltot
!> weighting factor controling relative importance of inflow rate and outflow
!> rate in determining storage on reach
   real*8 :: msk_x
   real*8 :: volcrmin
!> bacteria soil partitioning coefficient. Ratio of solution bacteria in surface
!> layer to solution bacteria in runoff soluble and sorbed phase in surface
!> runoff.
   real*8 :: bactkdq
!> die-off factor for persistent bacteria on foliage (1/day)
   real*8 :: wdpf
!> amount of water evaporated from canopy storage (mm H2O)
   real*8 :: canev
   real*8 :: precipday !< precipitation for the day in HRU (mm H2O)
   real*8 :: uno3d, usle, rcn, surlag_bsn
   real*8 :: thbact !< temperature adjustment factor for bacteria die-off/growth
!> overall rate change for less persistent bacteria in soil solution (1/day)
   real*8 :: wlpq20
!> overall rate change for less persistent bacteria adsorbed to soil particles
!> (1/day)
   real*8 :: wlps20
!> overall rate change for persistent bacteria in soil solution (1/day)
   real*8 :: wpq20
!> overall rate change for persistent bacteria adsorbed to soil particles
!> (1/day)
   real*8 :: wps20
   real*8 :: bactrop, bactsedp
   real*8 :: wgpf !< growth factor for persistent bacteria on foliage (1/day)
   real*8 :: bactlchp, bactlchlp, enratio, wetpcp, pndpcp, wetsep
   real*8 :: pndsep, wetev, pndev, pndsedo, wetsedo, pndflwi, wetflwi
!> drainage area of watershed in hectares (ha)
   real*8 :: da_ha
   real*8 :: pndflwo, wetflwo, wetsedi, vpd
!> leaf area index at which no evaporation occurs.  This variable is used in
!> ponded HRUs where evaporation from the water surface is restricted by the
!> plant canopy cover. Evaporation from the water surface equals potential ET
!> when LAI = 0 and decreased linearly to O when LAI = EVLAI
   real*8 :: evlai
!> Reach evaporation adjustment factor. Evaporation from the reach is multiplied
!> by EVRCH. This variable was created to limit the evaporation predicted in
!> arid regions.
   real*8 :: evrch
!> die-off factor for less persistent bacteria on foliage (1/day)
   real*8 :: wdlpf
!> actual amount of transpiration that occurs on day in HRU (mm H2O)
   real*8 :: ep_day
!> potential evapotranspiration on current day in HRU (mm H2O)
   real*8 :: pet_day
   real*8 :: bactrolp, bactsedlp
!> peak rate adjustment factor in the subbasin. Used in the MUSLE equation to
!> account for impact of peak flow on erosion (none)
   real*8 :: adj_pkr
!> nitrogen uptake distribution parameter. This parameter controls the amount of
!> nitrogen removed from the different soil layer layers by the plant. In
!> particular, this parameter allows the amount of nitrogen removed from the
!> surface layer via plant uptake to be controlled. While the relationship
!> between UBN and N removed from the surface layer is affected by the depth of
!> the soil profile, in general, as UBN increases the amount of N removed from
!> the surface layer relative to the amount removed from the entire profile
!> increases
   real*8 :: n_updis
!> nitrogen active pool fraction. The fraction of organic nitrogen in the active
!> pool (none)
   real*8 :: nactfr
!> phosphorus uptake distribution parameter This parameter controls the amount
!> of phosphorus removed from the different soil layers by the plant. In
!> particular, this parameter allows the amount of phosphorus removed from the
!> surface layer via plant uptake to be controlled. While the relationship
!> between UBP and P uptake from the surface layer is affected by the depth of
!> the soil profile, in general, as UBP increases the amount of P removed from
!> the surface layer relative to the amount removed from the entire profile
!> increases
   real*8 :: p_updis
   real*8 :: snoev, sno3up, reactw
!> actual amount of evaporation (soil et) that occurs on day in HRU (mm H2O)
   real*8 :: es_day
   real*8 :: sdiegropq, sdiegrolpq, sdiegrops, sdiegrolps
!> wash off fraction for less persistent bacteria on foliage during a rainfall
!> event
   real*8 :: wof_lp
   real*8 :: sbactrop, sbactrolp, sbactsedp, sbactsedlp, ep_max
   real*8 :: sbactlchp, sbactlchlp, psp_bsn, rchwtr, resuspst, setlpst
   real*8 :: bsprev, bssprev, spadyo, spadyev, spadysp, spadyrfv
   real*8 :: spadyosp
!> surface runoff loading to main channel from HRU for day (mm H2O)
   real*8 :: qday
   real*8 :: usle_ei, al5, pndsedc, no3pcp, rcharea, volatpst
!> water uptake distribution parameter. This parameter controls the amount of
!> water removed from the different soil layers by the plant. In particular,
!> this parameter allows the amount of water removed from the surface layer via
!> plant uptake to be controlled. While the relationship between UBW and H2O
!> removed from the surface layer is affected by the depth of the soil profile,
!> in general, as UBW increases the amount of water removed from the surface
!> layer relative to the amount removed from the entire profile increases
   real*8 :: ubw
!> nitrogen uptake normalization parameter. This variable normalizes the
!> nitrogen uptake so that the model can easily verify that upake from the
!> different soil layers sums to 1.0
   real*8 :: uobn
!> phosphorus uptake normalization parameter. This variable normalizes the
!> phosphorus uptake so that the model can easily verify that uptake from the
!> different soil layers sums to 1.0
   real*8 :: uobp
!> water uptake normalization parameter. This variable normalizes the water
!> uptake so that the model can easily verify that uptake from the different
!> soil layers sums to 1.0
   real*8 :: uobw
!> growth factor for less persistent bacteria on foliage (1/day)
   real*8 :: wglpf
   real*8 :: wetsedc, respesti
!> correction coefficient for generated rainfall to ensure that the annual
!> means for generated and observed values are comparable (needed only if
!> IDIST=1)
   real*8 :: rcor
!> value of exponent for mixed exponential rainfall distribution (needed only if
!> IDIST=1)
   real*8 :: rexp
!> 1st shape parameter for snow cover equation. This parameter is determined by
!> solving the equation for 50% snow cover
   real*8 :: snocov1
!> 2nd shape parameter for snow cover equation. This parameter is determined by
!> solving the equation for 95% snow cover
   real*8 :: snocov2
!> Minimum snow water content that corresponds to 100% snow cover. If the snow
!> water content is less than SNOCOVMX, then a certain percentage of the ground
!> will be bare (mm H2O)
   real*8 :: snocovmx, lyrtile, lyrtilex
!> Fraction of SNOCOVMX that corresponds to 50% snow cover. SWAT assumes a
!> nonlinear relationship between snow water and snow cover
   real*8 :: sno50cov
   real*8 :: ai0 !< ratio of chlorophyll-a to algal biomass (ug chla/mg alg)
   real*8 :: ai1 !< fraction of algal biomass that is nitrogen (mg N/mg alg)
   real*8 :: ai2 !< fraction of algal biomass that is phosphorus (mg P/mg alg)
!> the rate of oxygen production per unit of algal photosynthesis (mg O2/mg alg)
   real*8 :: ai3
!> the rate of oxygen uptake per unit of algae respiration (mg O2/mg alg)
   real*8 :: ai4
!> the rate of oxygen uptake per unit of NH3 nitrogen oxidation (mg O2/mg N)
   real*8 :: ai5
!> the rate of oxygen uptake per unit of NO2 nitrogen oxidation (mg O2/mg N)
   real*8 :: ai6
   real*8 :: rhoq !< algal respiration rate (1/day or 1/hr)
!> fraction of solar radiation computed in the temperature heat balance that is
!> photosynthetically active
   real*8 :: tfact
!> half-saturation coefficient for light (MJ/(m2*hr))
   real*8 :: k_l
!> michaelis-menton half-saturation constant for nitrogen (mg N/L)
   real*8 :: k_n
!> michaelis-menton half saturation constant for phosphorus (mg P/L)
   real*8 :: k_p
!> non-algal portion of the light extinction coefficient (1/m)
   real*8 :: lambda0
   real*8 :: lambda1 !< linear algal self-shading coefficient (1/(m*ug chla/L))
!> nonlinear algal self-shading coefficient ((1/m)(ug chla/L)**(-2/3))
   real*8 :: lambda2
   real*8 :: mumax !< maximum specific algal growth rate (1/day or 1/hr)
   real*8 :: p_n !< algal preference factor for ammonia
   real*8 :: rnum1 !< variable to hold value for rnum1s(:) (none)
!> actual evapotranspiration occuring on day in HRU (mm H2O)
   real*8 :: etday
   real*8 :: autop, auton, hmntl, rwntl, hmptl, rmn2tl
   real*8 :: rmptl,wdntl,cmn_bsn,rmp1tl,roctl,gwseep,revapday,reswtr
!> die-off factor for less persistent bacteria in streams (1/day)
   real*8 :: wdlprch
!> die-off factor for persistent bacteria in reservoirs (1/day)
   real*8 :: wdpres
!> potential ET value read in for day (mm H2O)
   real*8 :: petmeas
   real*8 :: bury, difus, reactb, solpesto
!> die-off factor for less persistent bacteria in reservoirs (1/day)
   real*8 :: wdlpres
   real*8 :: sorpesto, spcon_bsn, spexp_bsn, solpesti, sorpesti
!> calibration coefficient to control impact of the storage time constant for
!> the reach at bankfull depth (phi(10,:) upon the storage time constant for the
!> reach used in the Muskingum flow method
   real*8 :: msk_co1
!> calibration coefficient to control impact of the storage time constant for
!> the reach at 0.1 bankfull depth (phi(13,:) upon the storage time constant for
!> the reach used in the Muskingum flow method
   real*8 :: msk_co2
   real*8 :: snoprev, swprev, shallstp, deepstp
   real*8 :: ressolpo, resorgno, resorgpo, resno3o, reschlao, resno2o
   real*8 :: resnh3o, qdbank, potpcpmm, potevmm, potsepmm, potflwo
!> Threshold detection level for less persistent bacteria. When bacteria levels
!> drop to this amount the model considers bacteria in the soil to be
!> insignificant and sets the levels to zero (cfu/m^2)
   real*8 :: bactminlp
!> Threshold detection level for persistent bacteria. When bacteria levels
!> drop to this amount the model considers bacteria in the soil to be
!> insignificant and sets the levels to zero (cfu/m^2)
   real*8 :: bactminp
!> fraction of transmission losses from main channel that enter deep aquifer
   real*8 :: trnsrch
!> overall rate change for persistent bacteria on foliage (1/day)
   real*8 :: wp20p_plt
   real*8 :: potsedo, pest_sol
!> fraction of manure containing active colony forming units (cfu)
   real*8 :: bact_swf
!> bacteria percolation coefficient. Ratio of solution bacteria in surface layer
!> to solution bacteria in percolate
   real*8 :: bactmx
   real*8 :: cncoef !< plant ET curve number coefficient
!> overall rate change for less persistent bacteria on foliage (1/day)
   real*8 :: wp20lp_plt
   real*8 :: cdn_bsn,sdnco_bsn,bactmin
   real*8 :: cn_froz !< drainge coefficient (mm day -1)
   real*8 :: dorm_hr !< time threshold used to define dormant (hours)
   real*8 :: smxco !< adjustment factor for max curve number s factor (0-1)
   real*8 :: tb_adj !< adjustment factor for subdaily unit hydrograph basetime
   real*8 :: chla_subco !< regional adjustment on sub chla_a loading (fraction)
!> depth to impervious layer. Used to model perched water tables in all HRUs in
!> watershed (mm)
   real*8 :: depimp_bsn
   real*8 :: ddrain_bsn !< depth to the sub-surface drain (mm)
   real*8 :: tdrain_bsn !< time to drain soil to field capacity (hours)
   real*8 :: gdrain_bsn
   real*8 :: rch_san, rch_sil, rch_cla, rch_sag, rch_lag, rch_gra


!!    declare mike van liew variables
   real*8 :: hlife_ngw_bsn !< Half-life of nitrogen in groundwater? (days) 
   real*8 :: ch_opco_bsn, ch_onco_bsn
   real*8 :: decr_min !< Minimum daily residue decay
   real*8 :: rcn_sub_bsn !< Concentration of nitrogen in the rainfall (mg/kg)
   real*8 :: bc1_bsn, bc2_bsn, bc3_bsn, bc4_bsn
   real*8 :: anion_excl_bsn
!!    delcare mike van liew variables

!    Drainmod tile equations  01/2006
   real*8, dimension (:), allocatable :: wat_tbl,sol_swpwt
   real*8, dimension (:,:), allocatable :: vwt
   real*8 :: re_bsn !< Effective radius of drains (range 3.0 - 40.0) (mm)
!> Distance bewtween two drain or tile tubes (range 7600.0 - 30000.0) (mm)
   real*8 :: sdrain_bsn
   real*8 :: sstmaxd_bsn
!> Drainage coeffcient (range 10.0 - 51.0) (mm-day-1)
   real*8 :: drain_co_bsn
!> Multiplication factor to determine lateral ksat from SWAT ksat input value
!>for HRU (range 0.01 - 4.0)
   real*8 :: latksatf_bsn
!> Pump capacity (def val = 1.042 mm h-1 or 25 mm day-1) (mm h-1)
   real*8 :: pc_bsn
!    Drainmod tile equations  01/2006
   integer :: i_subhw, imgt, idlast, iwtr, ifrttyp, mo_atmo, mo_atmo1
   integer :: ifirstatmo, iyr_atmo, iyr_atmo1, matmo
   integer :: mch !< maximum number of channels
   integer :: mcr !< maximum number of crops grown per year
!> maximum number of crops/landcover in database file (crop.dat)
   integer :: mcrdb
   integer :: mfcst !< maximum number of forecast stations
   integer :: mfdb !< maximum number of fertilizers in fert.dat
   integer :: mhru !< maximum number of HRUs in watershed
   integer :: mhyd !< maximum number of hydrograph nodes
   integer :: mpdb !< maximum number of pesticides in pest.dat
   integer :: mrg !< maximum number of rainfall/temp gages (none)
   integer :: mcut !< maximum number of cuttings per year
   integer :: mgr !< maximum number of grazings per year
   integer :: mnr !< maximum number of years of rotation
   integer :: myr !< maximum number of years of simulation
!> subbasin water quality code\n
!> 0 do not calculate algae/CBOD
!> 1 calculate algae/CBOD drainmod tile equations
   integer :: isubwq
   integer :: ffcst
!> special project code (none):\n
!> 1 test rewind (run simulation twice)
   integer :: isproj
   integer :: nbyr !< number of calendar years simulated (none)
!> water routing method (none):\n
!> 0 variable storage method\n
!> 1 Muskingum method\n
   integer :: irte
   integer :: nrch !< number of reaches in watershed (none)
   integer :: nres !< number of reservoirs in watershed (none)
!> number of last HRU in previous subbasin or\n
!> number of HRUs in watershed (none)
   integer :: nhru
   integer :: i_mo !< current month being simulated (none)
   integer :: mo, immo
!> wind speed input code (noen)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: wndsim
   integer :: ihru !< HRU number (none)
   integer :: icode !< variable to hold value for icodes(:) (none)
   integer :: ihout !< variable to hold value for ihouts(:) (none)
!> variable to hold value for inum1s(:) (subbasin number) (none)
   integer :: inum1
   integer :: inum2 !< variable to hold value for inum2s(:) (none)
   integer :: inum3 !< variable to hold value for inum3s(:) (none)
   integer :: inum4 !< variable to hold value for inum4s(:) (none)
!> icfac = 0 for C-factor calculation using Cmin (as described in manual)\n
!> = 1 for new C-factor calculation from RUSLE (no minimum needed)
   integer :: icfac
   integer :: inum5, inum6, inum7, inum8
   integer :: mrech !< maximum number of rechour files
   integer :: nrgage !< number of raingage files (none)
   integer :: nrgfil !< number of rain gages per file (none)
   integer :: nrtot !< total number of rain gages (none)
   integer :: ntgage !< number of temperature gage files (none)
   integer :: ntgfil !< number of temperature gages per file (none)
   integer :: nttot !< total number of temperature gages (none)
!> temperature input code (none)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: tmpsim
!> crack flow code\n
!> 1: compute flow in cracks
   integer :: icrk
!> number of pesticide to be routed through the watershed. Redefined to the
!> sequence number of pestcide in NPNO(:) which is to be routed through the
!> watershed (none)
   integer :: irtpest
!> Qual2E option for calculating the local specific growth rate of algae\n
!> 1: multiplicative \f\[u=mumax\,fll\,fnn\,fpp\f\]
!> 2: limiting nutrient \f\[u=mumax\,fll\,\min(fnn,\,fpp)\f\]
!> 3: harmonic mean \f\[u=mumax\,fll\,\frac2{\frac1{fnn}+\frac1{fpp}}\f\]
   integer :: igropt
!> Qual2E light averaging option. Qual2E defines four light averaging options. 
!> The only option currently available in SWAT is #2.
   integer :: lao
!> number of different pesticides used in the simulation (none)
   integer :: npmx
   integer :: curyr !< current year in simulation (sequence) (none)
   integer :: iihru
!    Drainmod tile equations  01/2006
!> tile drainage equations flag/code\n
!> 1 simulate tile flow using subroutine drains(wt_shall)\n
!> 0 simulate tile flow using subroutine origtile(wt_shall,d)
   integer :: itdrn
!> water table depth algorithms flag/code\n
!> 1 simulate wt_shall using subroutine new water table depth routine\n
!> 0 simulate wt_shall using subroutine original water table depth routine
   integer :: iwtdn
!> maximum depressional storage selection flag/code\n
!> 0 = static depressional storage\n
!> 1 = dynamic storage based on tillage and cumulative rainfall
   integer :: ismax
!> not being implemented in this version drainmod tile equations
   integer :: iroutunit
   integer :: ires_nut
!    Drainmod tile equations  01/2006
   integer :: iclb !< auto-calibration flag
   integer :: mrecc !< maximum number of reccnst files
   integer :: mrecd !< maximum number of recday files
   integer :: mrecm !< maximum number of recmon files
   integer :: mtil !< max number of tillage types in till.dat
   integer :: mudb !< maximum number of urban land types in urban.dat
!> rainfall distribution code\n
!>   0 for skewed normal dist\n
!>   1 for mixed exponential distribution
   integer :: idist
   integer :: mrecy !< maximum number of recyear files
   integer :: nyskip !< number of years to not print output
!> solar radiation input code (none)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: slrsim
!> channel degredation code\n
!> 1: compute channel degredation (downcutting and widening)
   integer :: ideg
!> rainfall/runoff code\n
!> 0 daily rainfall/curve number technique
!> 1 sub-daily rainfall/Green&Ampt/hourly routing
!> 3 sub-daily rainfall/Green&Ampt/hourly routing
   integer :: ievent
!> code for potential ET method (none)\n
!> 0 Priestley-Taylor method\n
!> 1 Penman/Monteith method\n
!> 2 Hargreaves method\n
!> 3 read in daily potential ET data
   integer :: ipet
   integer :: iopera
   integer :: idaf !< beginning day of simulation (julian date)
   integer :: idal !< ending day of simulation (julian date)
!> relative humidity input code (none)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: rhsim
!> leap year flag (none)\n
!> 0  leap year\n
!> 1  regular year
   integer :: leapyr
   integer :: id1 !< first day of simulation in current year (julian date)
   integer :: mo_chk
   integer :: nhtot !< total number of relative humidity records in file
   integer :: nstot !< total number of solar radiation records in file (none)
   integer :: nwtot !< total number of wind speed records in file
!> solar radiation data search code (none)\n
!> 0 first day of solar radiation data located in file\n
!> 1 first day of solar radiation data not located in file
   integer :: ifirsts
!> relative humidity data search code (none)\n
!> 0 first day of relative humidity data located in file\n
!> 1 first day of relative humidity data not located in file
   integer :: ifirsth
!> wind speed data search code (none)\n
!> 0 first day of wind speed data located in file\n
!> 1 first day of wind speed data not located in file
   integer :: ifirstw
   integer :: icst
   integer :: ilog !< streamflow print code
   integer :: itotr !< number of output variables printed (output.rch)
   integer :: iyr !< beginning year of simulation (year)
!> stream water quality code\n
!> 0 do not model stream water quality\n
!> 1 model stream water quality (QUAL2E & pesticide transformations)
   integer :: iwq
!> flag for calculations performed only for the first year of simulation (none)
   integer :: iskip
!> potential ET data search code (none)\n
!> 0 first day of potential ET data located in file\n
!> 1 first day of potential ET data not located in file
   integer :: ifirstpet
!> print code for output.pst file\n
!> 0 do not print pesticide output\n
!> 1 print pesticide output
   integer :: iprp
   integer :: itotb !< number of output variables printed (output.sub)
   integer :: itots !< number of output variables printed (output.hru)
   integer :: itoth !< number of HRUs printed (output.hru/output.wtr)
!> rainfall input code (none)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: pcpsim
   integer :: nd_30,iops,iphr,isto,isol
!> number of times forecast period is simulated (using different weather
!> generator seeds each time)
   integer :: fcstcycles
   integer :: fcstday !< beginning date of forecast period (julian date)
   integer :: fcstyr !< beginning year of forecast period
   integer :: iscen !< scenarios counter
   integer :: subtot !< number of subbasins in watershed (none)
   integer :: ogen
   integer :: mapp !< maximum number of applications
   integer :: mlyr !< maximum number of soil layers
   integer :: mpst !< max number of pesticides used in wshed
   integer :: mres !< maximum number of reservoirs
   integer :: msub !< maximum number of subbasins
!> random number generator seed code (none):\n
!> 0: use default numbers\n
!> 1: generate new numbers in every simulation
   integer :: igen
   integer :: iprint !< print code: 0=monthly, 1=daily, 2=annual
   integer :: iida !< day being simulated (current julian date) (julian date)
!> CN method flag (for testing alternative method):\n
!> 0 use traditional SWAT method which bases CN on soil moisture\n
!> 1 use alternative method which bases CN on plant ET
   integer :: icn
!> max half-hour rainfall fraction calc option:\n
!> 0 generate max half-hour rainfall fraction from triangular distribution\n
!> 1 use monthly mean max half-hour rainfall fraction
   integer :: ised_det
   integer :: fcstcnt, mtran, idtill
   integer, dimension(100) :: ida_lup, iyr_lup
   integer :: no_lup, no_up, nostep
!  routing 5/3/2010 gsm per jga
! date
!> date simulation is performed where leftmost eight characters are set to a
!> value of yyyymmdd, where yyyy is the year, mm is the month and dd is the day
   character(len=8) :: date
!> time simulation is performed where leftmost ten characters are set to a value
!> of hhmmss.sss, where hh is the hour, mm is the minutes and ss.sss is the
!> seconds and milliseconds
   character(len=10) :: time
!> time difference with respect to Coordinated Universal Time (ie Greenwich Mean
!> Time)
   character(len=5) :: zone
!> name of file containing calibration parameters
   character(len=13) :: calfile
   character(len=13) :: rhfile !< relative humidity file name (.hmd)
   character(len=13) :: slrfile !< solar radiation file name (.slr)
   character(len=13) :: wndfile !< wind speed file name (.wnd)
   character(len=13) :: petfile !< potential ET file name (.pet)
   character(len=13) :: atmofile, lucfile
!> name of septic tank database file (septwq1.dat)
   character(len=13) :: septdb
   character(len=13) :: dpd_file, wpd_file, rib_file, sfb_file,&
      &lid_file
!> array location of random number seed used for a given process
   integer, dimension (9) :: idg
   integer, dimension (:), allocatable :: ifirstr, ifirsthr
!> values(1): year simulation is performed\n
!> values(2): month simulation is performed\n
!> values(3): day in month simulation is performed\n
!> values(4): time difference with respect to Coordinated Universal Time (ie
!> Greenwich Mean Time)\n
!> values(5): hour simulation is performed\n
!> values(6): minute simulation is performed\n
!> values(7): second simulation is performed\n
!> values(8): millisecond simulation is performed
   integer, dimension (8) :: values
!> julian date for last day of preceding month (where the array location is the
!> number of the month). The dates are for leap years (julian date)
   integer, dimension (13) :: ndays
   integer, dimension (13) :: ndays_noleap, ndays_leap
   integer :: mapex
   real*8, dimension (:), allocatable :: flodaya, seddaya, orgndaya
   real*8, dimension (:), allocatable :: orgpdaya, no3daya, minpdaya
!> harvest index target of cover defined at planting ((kg/ha)/(kg/ha))
   real*8, dimension (:), allocatable :: hi_targ
   real*8, dimension (:), allocatable :: bio_targ !< biomass target (kg/ha)
   real*8, dimension (:), allocatable :: tnyld
   integer, dimension (:), allocatable :: idapa, iypa, ifirsta
   integer, dimension (100) :: mo_transb, mo_transe
   integer, dimension (100) :: ih_tran
!     apex/command output files
!  septic inputs
!! septic change added iseptic 1/28/09 gsm
   integer :: msdb !< maximum number of sept wq data database (none)
   integer :: iseptic
!> flow rate of the septic tank effluent per capita (m3/d)
   real*8, dimension (:), allocatable :: sptqs
   real*8, dimension (:), allocatable :: percp
!> Biological Oxygen Demand of the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptbodconcs
!> concentration of total suspended solid in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: spttssconcs
!> concentration of total nitrogen in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: spttnconcs
!> concentration of total phosphorus of the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptnh4concs
!> concentration of nitrate in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptno3concs
!> concentration of nitrite in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptno2concs
!> concentration of organic nitrogen in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptorgnconcs
!> concentration of total phosphorus in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: spttpconcs
!> concentration of mineral phosphorus in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptminps
!> concentration of organic phosphorus in the septic tank effluent (mg/l)
   real*8, dimension (:), allocatable :: sptorgps
!> concentration of the facel caliform in the septic tank effluent (cfu/100ml)
   real*8, dimension (:), allocatable :: sptfcolis
   real*8, dimension (:), allocatable :: failyr,qstemm
!! septic changes added 1/28/09 gsm
   real*8, dimension (:), allocatable :: bio_amn, bio_bod, biom,rbiom
   real*8, dimension (:), allocatable :: fcoli, bio_ntr, bz_perc
!> number of permanent residents in the hourse (none)
   real*8, dimension (:), allocatable :: sep_cap
   real*8, dimension (:), allocatable :: plqm,bz_area
   real*8, dimension (:), allocatable :: bz_z !< Depth of biozone layer(mm)
   real*8, dimension (:), allocatable :: bz_thk !< thickness of biozone (mm)
   real*8, dimension (:), allocatable :: bio_bd !< density of biomass (kg/m^3)
!! carbon outputs for .hru file
   real*8, dimension (:), allocatable :: cmup_kgh, cmtot_kgh
!! carbon outputs for .hru file
!> denitrification rate coefficient (none)
   real*8, dimension (:), allocatable :: coeff_denitr
!> BOD decay rate coefficient (m^3/day)
   real*8, dimension (:), allocatable :: coeff_bod_dc
!> BOD to live bacteria biomass conversion factor (none)
   real*8, dimension (:), allocatable :: coeff_bod_conv
!> field capacity calibration parameter 1 (none)
   real*8, dimension (:), allocatable :: coeff_fc1
!> field capacity calibration parameter 2 (none)
   real*8, dimension (:), allocatable :: coeff_fc2
!> fecal coliform bacteria decay rate coefficient (m^3/day)
   real*8, dimension (:), allocatable :: coeff_fecal
!> mortality rate coefficient (none)
   real*8, dimension (:), allocatable :: coeff_mrt
!> nitrification rate coefficient (none)
   real*8, dimension (:), allocatable :: coeff_nitr
!> conversion factor for plaque from TDS (none)
   real*8, dimension (:), allocatable :: coeff_plq
!> respiration rate coefficient (none)
   real*8, dimension (:), allocatable :: coeff_rsp
!> slough-off calibration parameter (none)
   real*8, dimension (:), allocatable :: coeff_slg1
!> slough-off calibration parameter (none)
   real*8, dimension (:), allocatable :: coeff_slg2
   real*8, dimension (:), allocatable :: coeff_pdistrb,coeff_solpslp
   real*8, dimension (:), allocatable :: coeff_solpintc,coeff_psorpmax
!! Septic system by Jaehak Jeong
!> septic system type (none)
   integer, dimension (:), allocatable :: isep_typ
   integer, dimension (:), allocatable :: i_sep
!> septic system operation flag (1=active, 2=failing, 3=not operated) (none)
   integer, dimension (:), allocatable :: isep_opt
   integer, dimension (:), allocatable :: sep_tsincefail
   integer, dimension (:), allocatable :: isep_tfail,isep_iyr
   integer, dimension (:), allocatable :: sep_strm_dist,sep_den

   !!   change per JGA 9/8/2011 gsm for output.mgt
   real*8, dimension (:), allocatable :: sol_sumno3, sol_sumsolp
   real*8, dimension (:), allocatable :: strsw_sum, strstmp_sum
   real*8, dimension (:), allocatable :: strsn_sum, strsp_sum
   real*8, dimension (:), allocatable :: strsa_sum


!! New pothole variables
   real*8, dimension (:), allocatable :: spill_hru,tile_out,hru_in
   real*8, dimension (:), allocatable :: spill_precip,pot_seep
   real*8, dimension (:), allocatable :: pot_evap,pot_sedin
!> soluble P loss rate in the pothole (.01 - 0.5) (1/d)
   real*8, dimension (:), allocatable :: pot_solp
   real*8, dimension (:), allocatable :: pot_solpi
   real*8, dimension (:), allocatable :: pot_orgp,pot_orgpi
   real*8, dimension (:), allocatable :: pot_orgn,pot_orgni
   real*8, dimension (:), allocatable :: pot_mps,pot_mpsi
   real*8, dimension (:), allocatable :: pot_mpa,pot_mpai
   real*8, dimension (:), allocatable :: pot_no3i,precip_in
   real*8, dimension (:), allocatable :: tile_sedo,tile_no3o
   real*8, dimension (:), allocatable :: tile_solpo,tile_orgno
   real*8, dimension (:), allocatable :: tile_orgpo,tile_minpso
   real*8, dimension (:), allocatable :: tile_minpao
   integer :: ia_b, ihumus, itemp, isnow
!> output variable codes for output.rch file (none)
   integer, dimension (46) :: ipdvar
!> output varaible codes for output.hru file (none)
   integer, dimension (mhruo) :: ipdvas
!> output variable codes for output.sub file (none)
   integer, dimension (msubo) :: ipdvab
!> HRUs whose output information will be printed to the output.hru and 
!> output.wtr files
   integer, dimension (:), allocatable :: ipdhru
   real*8, dimension (mstdo) :: wshddayo,wshdmono,wshdyro
   real*8, dimension (16) :: fcstaao
   real*8, dimension (mstdo) :: wshdaao
   real*8, dimension (:,:), allocatable :: wpstdayo,wpstmono,wpstyro
   real*8, dimension (:,:), allocatable :: yldkg, bio_hv
!> reach monthly output array (varies)
   real*8, dimension (:,:), allocatable :: rchmono
   real*8, dimension (:,:), allocatable :: wpstaao,rchyro
!> HRU monthly output data array (varies)
   real*8, dimension (:,:), allocatable :: hrumono
   real*8, dimension (:,:), allocatable :: rchaao,rchdy,hruyro
!> subbasin monthly output array (varies)
   real*8, dimension (:,:), allocatable :: submono
   real*8, dimension (:,:), allocatable :: hruaao,subyro,subaao
!> reservoir monthly output array (varies)
   real*8, dimension (:,:), allocatable :: resoutm
   real*8, dimension (:,:), allocatable :: resouty,resouta
   real*8, dimension (12,8) :: wshd_aamon
!> HRU monthly output data array for impoundments (varies)
   real*8, dimension (:,:), allocatable :: wtrmon
   real*8, dimension (:,:), allocatable :: wtryr,wtraa
!> max melt rate for snow during year (June 21) for subbasin(:) where deg C
!> refers to the air temperature. SUB_SMFMX and SMFMN allow the rate of snow
!> melt to vary through the year. These parameters are accounting for the impact
!> of soil temperature on snow melt (range: -5.0/5.0) (mm/deg C/day)
   real*8, dimension (:,:), allocatable :: sub_smfmx
!> min melt rate for snow during year (Dec 21) for subbasin(:) (range: -5.0/5.0)
!> where deg C refers to the air temperature (mm/deg C/day)
   real*8, dimension (:,:), allocatable :: sub_smfmn
   real*8, dimension (:,:,:), allocatable :: hrupstd,hrupsta,hrupstm
   real*8, dimension (:,:,:), allocatable :: hrupsty
!> temperature data search code (none)\n
!> 0 first day of temperature data located in file\n
!> 1 first day of temperature data not located in file
   integer, dimension (:), allocatable :: ifirstt
   integer, dimension (:), allocatable :: ifirstpcp
!> elevation of precipitation gage station (m)
   integer, dimension (:), allocatable :: elevp
!> elevation of temperature gage station (m)
   integer, dimension (:), allocatable :: elevt
!> avg monthly minimum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: ftmpmn
!> avg monthly maximum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: ftmpmx
!> standard deviation for avg monthly minimum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: ftmpstdmn
!> standard deviation for avg monthly maximum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: ftmpstdmx
!> fpcp_stat(:,1,:): average amount of precipitation falling in one day for the
!> month (mm/day)\n
!> fpcp_stat(:,2,:): standard deviation for the average daily precipitation
!> (mm/day)\n
!> fpcp_stat(:,3,:): skew coefficient for the average daily precipitationa
!> (none)
   real*8, dimension (:,:,:), allocatable :: fpcp_stat
!> probability of wet day after dry day in month (none)
   real*8, dimension (:,:), allocatable :: fpr_w1
!> probability of wet day after wet day in month (none)
   real*8, dimension (:,:), allocatable :: fpr_w2
!> proportion of wet days in the month (none)
   real*8, dimension (:,:), allocatable :: fpr_w3
!> average depth of main channel (m)
   real*8, dimension (:), allocatable :: ch_d
   real*8, dimension (:), allocatable :: flwin,flwout,bankst,ch_wi
!> channel organic n concentration (ppm)
   real*8, dimension (:), allocatable :: ch_onco
!> channel organic p concentration (ppm)
   real*8, dimension (:), allocatable :: ch_opco
   real*8, dimension (:), allocatable :: ch_orgn, ch_orgp
   real*8, dimension (:), allocatable :: drift,rch_dox,rch_bactp
!> alpha factor for bank storage recession curve (days)
   real*8, dimension (:), allocatable :: alpha_bnk
!> \f$\exp(-alpha_bnk)\f$ (none)
   real*8, dimension (:), allocatable :: alpha_bnke
   real*8, dimension (:), allocatable :: disolvp,algae,sedst,rchstor
   real*8, dimension (:), allocatable :: organicn,organicp,chlora
!> initial length of main channel (km)
   real*8, dimension (:), allocatable :: ch_li
!> initial slope of main channel (m/m)
   real*8, dimension (:), allocatable :: ch_si
   real*8, dimension (:), allocatable :: nitraten,nitriten

!     Sediment parameters added by Balaji for the new routines

   real*8, dimension (:), allocatable :: ch_bnk_san, ch_bnk_sil
   real*8, dimension (:), allocatable :: ch_bnk_cla, ch_bnk_gra
   real*8, dimension (:), allocatable :: ch_bed_san, ch_bed_sil
   real*8, dimension (:), allocatable :: ch_bed_cla, ch_bed_gra
   real*8, dimension (:), allocatable :: depfp,depsanfp,depsilfp
   real*8, dimension (:), allocatable :: depclafp,depsagfp,deplagfp
   real*8, dimension (:), allocatable :: depch,depsanch,depsilch
   real*8, dimension (:), allocatable :: depclach,depsagch,deplagch
   real*8, dimension (:), allocatable :: depgrach,depgrafp,grast
!> curve number retention parameter adjustment factor to adjust surface runoff
!> for flat slopes (0.5 - 3.0) (dimensionless)
   real*8, dimension (:), allocatable :: r2adj
!> Reach peak rate adjustment factor for sediment routing in the channel. Allows
!> impact of peak flow rate on sediment routing and channel reshaping to be
!> taken into account (none)
   real*8, dimension (:), allocatable :: prf
   real*8, dimension (:), allocatable :: depprch,depprfp
!> linear parameter for calculating sediment reentrained in channel sediment
!> routing
   real*8, dimension (:), allocatable :: spcon
!> exponent parameter for calculating sediment reentrained in channel sediment
!> routing
   real*8, dimension (:), allocatable :: spexp
   real*8, dimension (:), allocatable :: sanst,silst,clast,sagst,lagst
   real*8, dimension (:), allocatable :: pot_san,pot_sil,pot_cla
   real*8, dimension (:), allocatable :: pot_sag,pot_lag
   real*8, dimension (:), allocatable :: potsani,potsili,potclai
   real*8, dimension (:), allocatable :: potsagi,potlagi
   real*8, dimension (:), allocatable :: sanyld,silyld,clayld,sagyld
   real*8, dimension (:), allocatable :: lagyld,grayld
   real*8, dimension (:), allocatable :: res_san,res_sil,res_cla
   real*8, dimension (:), allocatable :: res_sag,res_lag,res_gra
   real*8, dimension (:), allocatable :: pnd_san,pnd_sil,pnd_cla
   real*8, dimension (:), allocatable :: pnd_sag,pnd_lag
   real*8, dimension (:), allocatable :: wet_san,wet_sil,wet_cla
   real*8, dimension (:), allocatable :: wet_lag, wet_sag
   real*8 :: ressano,ressilo,resclao,ressago,reslago, resgrao
   real*8 :: ressani, ressili, resclai, ressagi, reslagi,resgrai
   real*8 :: potsano,potsilo,potclao,potsago,potlago
   real*8 :: pndsanin,pndsilin,pndclain,pndsagin,pndlagin
   real*8 :: pndsano,pndsilo,pndclao,pndsago,pndlago

!> initial depth of main channel (m)
   real*8, dimension (:), allocatable :: ch_di
!> channel erodibility factor (0.0-1.0) (none)\n
!> 0 non-erosive channel\n
!> 1 no resistance to erosion
   real*8, dimension (:), allocatable :: ch_erod
!> length of main channel (km)
   real*8, dimension (:), allocatable :: ch_l2
   real*8, dimension (:), allocatable :: ch_cov
!> bulk density of channel bank sediment (1.1-1.9) (g/cc)
   real*8, dimension (:), allocatable :: ch_bnk_bd
!> bulk density of channel bed sediment (1.1-1.9) (g/cc)
   real*8, dimension (:), allocatable :: ch_bed_bd
!> erodibility of channel bank sediment by jet test (Peter Allen needs to give
!> more info on this)
   real*8, dimension (:), allocatable :: ch_bnk_kd
!> erodibility of channel bed sediment by jet test (Peter Allen needs to give
!> more info on this)
   real*8, dimension (:), allocatable :: ch_bed_kd
!> D50(median) particle size diameter of channel bank sediment (0.001 - 20)
   real*8, dimension (:), allocatable :: ch_bnk_d50
!> D50(median) particle size diameter of channel bed sediment (micrometers)
!> (0.001 - 20)
   real*8, dimension (:), allocatable :: ch_bed_d50
!> channel erodibility factor (0.0-1.0) (none)\n
!> 0 non-erosive channel\n
!> 1 no resistance to erosion
   real*8, dimension (:), allocatable :: ch_cov1
!> channel cover factor (0.0-1.0) (none)\n
!> 0 channel is completely protected from erosion by cover\n
!> 1 no vegetative cover on channel
   real*8, dimension (:), allocatable :: ch_cov2
!> critical shear stress of channel bed (N/m2)
   real*8, dimension (:), allocatable :: tc_bed
!> critical shear stress of channel bank (N/m2)
   real*8, dimension (:), allocatable :: tc_bnk
!> sediment routine methods (DAILY):\n
!>  0 = original SWAT method\n
!>  1 = Bagnold's\n
!>  2 = Kodatie\n
!>  3 = Molinas WU\n
!>  4 = Yang
   integer, dimension (:), allocatable :: ch_eqn
!> pesticide reaction coefficient in reach (1/day)
   real*8, dimension (:), allocatable :: chpst_rea
!> pesticide volatilization coefficient in reach (m/day)
   real*8, dimension (:), allocatable :: chpst_vol
   real*8, dimension (:), allocatable :: chpst_conc
!> pesticide partition coefficient between water and sediment in reach (m^3/g)
   real*8, dimension (:), allocatable :: chpst_koc
!> resuspension velocity in reach for pesticide sorbed to sediment (m/day)
   real*8, dimension (:), allocatable :: chpst_rsp
!> settling velocity in reach for pesticide sorbed to sediment (m/day)
   real*8, dimension (:), allocatable :: chpst_stl
!> channel width to depth ratio (m/m)
   real*8, dimension (:), allocatable :: ch_wdr
!> mixing velocity (diffusion/dispersion) for pesticide in reach (m/day)
   real*8, dimension (:), allocatable :: chpst_mix
!> inital pesticide concentration in river bed sediment (mg/m^3)
   real*8, dimension (:), allocatable :: sedpst_conc
!> pesticide burial velocity in river bed sediment (m/day)
   real*8, dimension (:), allocatable :: sedpst_bry
!> pesticide reaction coefficient in river bed sediment (1/day)
   real*8, dimension (:), allocatable :: sedpst_rea
!> depth of active sediment layer in reach for pesticide (m)
   real*8, dimension (:), allocatable :: sedpst_act
   real*8, dimension (:), allocatable :: rch_cbod,rch_bactlp
!> change in horizontal distance per unit vertical distance (0.0 - 5)\n
!> 0 = for vertical channel bank\n
!> 5 = for channel bank with gentl side slope
   real*8, dimension (:), allocatable :: chside
!> local algal settling rate in reach at 20 deg C (m/day or m/hour)
   real*8, dimension (:), allocatable :: rs1
!> benthos source rate for dissolved phosphorus in reach at 20 deg C
!> ((mg disP-P)/(m^2*day) or (mg disP-P)/(m^2*hour))
   real*8, dimension (:), allocatable :: rs2
!> benthos source rate for ammonia nitrogen in reach at 20 deg C
!> ((mg NH4-N)/(m^2*day) or (mg NH4-N)/(m^2*hour))
   real*8, dimension (:), allocatable :: rs3
!> rate coefficient for organic nitrogen settling in reach at 20 deg C
!> (1/day or 1/hour)
   real*8, dimension (:), allocatable :: rs4
!> organic phosphorus settling rate in reach at 20 deg C (1/day or 1/hour)
   real*8, dimension (:), allocatable :: rs5
!> CBOD deoxygenation rate coefficient in reach at 20 deg C (1/day or 1/hour)
   real*8, dimension (:), allocatable :: rk1
!> reaeration rate in accordance with Fickian diffusion in reach at 20 deg C
!> (1/day or 1/hour)
   real*8, dimension (:), allocatable :: rk2
!> rate of loss of CBOD due to settling in reach at 20 deg C (1/day or 1/hour)
   real*8, dimension (:), allocatable :: rk3
!> sediment oxygen demand rate in reach at 20 deg C
!> (mg O2/(m^2*day) or mg O2/(m^2*hour))
   real*8, dimension (:), allocatable :: rk4
!> coliform die-off rate in reach (1/day)
   real*8, dimension (:), allocatable :: rk5
!> rate coefficient for settling of arbitrary non-conservative constituent in
!> reach (1/day)
   real*8, dimension (:), allocatable :: rs6
!> benthal source rate for arbitrary non-conservative constituent in reach
!> ((mg ANC)/(m^2*day))
   real*8, dimension (:), allocatable :: rs7
!> rate constant for biological oxidation of NH3 to NO2 in reach at 20 deg C
!> (1/day or 1/hour)
   real*8, dimension (:), allocatable :: bc1
!> rate constant for biological oxidation of NO2 to NO3 in reach at 20 deg C
!> (1/day or 1/hour)
   real*8, dimension (:), allocatable :: bc2
!> rate constant for hydrolysis of organic N to ammonia in reach at 20 deg C
!> (1/day or 1/hour)
   real*8, dimension (:), allocatable :: bc3
!> rate constant for the decay of organic P to dissolved P in reach at 20 deg C
!> (1/day or 1/hour)
   real*8, dimension (:), allocatable :: bc4
!> decay rate for arbitrary non-conservative constituent in reach (1/day)
   real*8, dimension (:), allocatable :: rk6
   real*8, dimension (:), allocatable :: ammonian
   real*8, dimension (:), allocatable :: orig_sedpstconc
!> average daily water removal from the reach for the month (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wurch
   integer, dimension (:), allocatable :: icanal
   integer, dimension (:), allocatable :: itb
!> revap coeff: this variable controls the amount of water moving from bank
!> storage to the root zone as a result of soil moisture depletion(none)
   real*8, dimension (:), allocatable :: ch_revap
   real*8, dimension (:), allocatable :: dep_chan
!> coefficient related to radiation used in hargreaves eq
!> (range: 0.0019 - 0.0032)
   real*8, dimension (:), allocatable :: harg_petco
   real*8, dimension (:), allocatable :: subfr_nowtr
!> soil water depletion coefficient used in the new (modified curve number
!> method) same as soil index coeff used in APEX range: 0.5 - 2.0
   real*8, dimension (:), allocatable :: cncoef_sub
   real*8, dimension (:), allocatable :: dr_sub
!> fraction of total watershed area contained in subbasin (km2/km2)
   real*8, dimension (:), allocatable :: sub_fr
   real*8, dimension (:), allocatable :: wcklsp,sub_minp,sub_sw
   real*8, dimension (:), allocatable :: sub_sumfc,sub_gwno3,sub_gwsolp
   real*8, dimension (:), allocatable :: co2 !< CO2 concentration (ppmv)
!> area of subbasin in square kilometers (km^2)
   real*8, dimension (:), allocatable :: sub_km
!> latitude of weather station used to compile data (degrees)
   real*8, dimension (:), allocatable :: wlat
!> time of concentration for subbasin (hour)
   real*8, dimension (:), allocatable :: sub_tc
   real*8, dimension (:), allocatable :: sub_pet
!> elevation of weather station used to compile weather generator data (m)
   real*8, dimension (:), allocatable :: welev
   real*8, dimension (:), allocatable :: sub_orgn,sub_orgp,sub_bd
   real*8, dimension (:), allocatable :: sub_wtmp,sub_sedpa,sub_sedps
!> shortest daylength occurring during the year (hour)
   real*8, dimension (:), allocatable :: daylmn
   real*8, dimension (:), allocatable :: sub_minpa,sub_minps
   real*8, dimension (:), allocatable :: latcos !< \f$\cos(latitude)\f$ (none)
   real*8, dimension (:), allocatable :: latsin !< \f$\sin(latitude)\f$ (none)
!> total potential heat units for year (used when no crop is growing)
!> (heat unit)
   real*8, dimension (:), allocatable :: phutot
!> precipitation lapse rate: precipitation change due to change in elevation
!> (mm H2O/km)
   real*8, dimension (:), allocatable :: plaps
!> temperature lapse rate: temperature change due to change in elevation
!> (deg C/km)
   real*8, dimension (:), allocatable :: tlaps
!> average annual air temperature (deg C)
   real*8, dimension (:), allocatable :: tmp_an
   real*8, dimension (:), allocatable :: sub_precip
!> atmospheric deposition of ammonium values for entire watershed (mg/l)
   real*8, dimension (:), allocatable :: rammo_sub
!> atmospheric deposition of nitrate for entire watershed (mg/l)
   real*8, dimension (:), allocatable :: rcn_sub
   real*8, dimension (:), allocatable :: pcpdays
   real*8, dimension (:), allocatable :: atmo_day
   real*8, dimension (:), allocatable :: sub_snom,sub_qd,sub_sedy
   real*8, dimension (:), allocatable :: sub_tran,sub_no3,sub_latno3
!> snowfall temperature for subbasin(:). Mean air temperature at which precip
!> is equally likely to be rain as snow/freezing rain (range: -5.0/5.0) (deg C)
   real*8, dimension (:,:), allocatable :: sub_sftmp
!> snow melt base temperature for subbasin(:) mean air temperature at which snow
!> melt will occur (range: -5.0/5.0) (deg C)
   real*8, dimension (:,:), allocatable :: sub_smtmp
!> snow pack temperature lag factor (0-1) (none)
   real*8, dimension (:,:), allocatable :: sub_timp
   real*8, dimension (:), allocatable :: sub_tileno3
   real*8, dimension (:), allocatable :: sub_solp,sub_subp,sub_etday
!> average elevation of subbasin (m)
   real*8, dimension (:), allocatable :: sub_elev
   real*8, dimension (:), allocatable :: sub_wyld,sub_surfq
   real*8, dimension (:), allocatable :: qird
   real*8, dimension (:), allocatable :: sub_gwq,sub_sep,sub_chl
   real*8, dimension (:), allocatable :: sub_cbod,sub_dox,sub_solpst
   real*8, dimension (:), allocatable :: sub_sorpst,sub_yorgn,sub_yorgp
!> latitude of HRU/subbasin (degrees)
   real*8, dimension (:), allocatable :: sub_lat
   real*8, dimension (:), allocatable :: sub_bactp,sub_bactlp
   real*8, dimension (:), allocatable :: sub_latq, sub_gwq_d,sub_tileq
   real*8, dimension (:), allocatable :: sub_vaptile
   real*8, dimension (:), allocatable :: sub_dsan, sub_dsil, sub_dcla
   real*8, dimension (:), allocatable :: sub_dsag, sub_dlag

!!!!!! drains
   real*8 :: vap_tile
   real*8, dimension (:), allocatable :: wnan
   real*8, dimension (:,:), allocatable :: sol_stpwt
   real*8, dimension (:,:), allocatable :: sub_pst,sub_hhqd,sub_hhwtmp
!> monthly humidity adjustment. Daily values for relative humidity within the
!> month are rasied or lowered by the specified amount (used in climate change
!> studies) (none)
   real*8, dimension (:,:), allocatable :: huminc
!> monthly solar radiation adjustment. Daily radiation within the month is
!> raised or lowered by the specified amount (used in climate change studies)
!> (MJ/m^2)
   real*8, dimension (:,:), allocatable :: radinc
!> monthly rainfall adjustment. Daily rainfall within the month is adjusted to
!> the specified percentage of the original value (used in climate change
!> studies)(%)
   real*8, dimension (:,:), allocatable :: rfinc
!> monthly temperature adjustment. Daily maximum and minimum temperatures within
!> the month are raised or lowered by the specified amount (used in climate
!> change studies) (deg C)
   real*8, dimension (:,:), allocatable :: tmpinc
!> effective hydraulic conductivity of tributary channel alluvium (mm/hr)
   real*8, dimension (:), allocatable :: ch_k1
!> effective hydraulic conductivity of main channel alluvium (mm/hr)
   real*8, dimension (:), allocatable :: ch_k2
!> elevation at the center of the band (m)
   real*8, dimension (:,:), allocatable :: elevb
!> fraction of subbasin area within elevation band (the same fractions should be
!> listed for all HRUs within the subbasin) (none)
   real*8, dimension (:,:), allocatable :: elevb_fr
!> average wind speed for the month (m/s)
   real*8, dimension (:,:), allocatable :: wndav
!> Manning's "n" value for the tributary channels (none)
   real*8, dimension (:), allocatable :: ch_n1
!> Manning's "n" value for the main channel (none)
   real*8, dimension (:), allocatable :: ch_n2
!> average slope of tributary channels (m/m)
   real*8, dimension (:), allocatable :: ch_s1
!> average slope of main channel (m/m)
   real*8, dimension (:), allocatable :: ch_s2
!> average width of tributary channels (m)
   real*8, dimension (:), allocatable :: ch_w1
!> average width of main channel (m)
   real*8, dimension (:), allocatable :: ch_w2
!> average dew point temperature for the month (deg C)
   real*8, dimension (:,:), allocatable :: dewpt
!> average fraction of total daily rainfall occuring in maximum half-hour period
!> for month (none)
   real*8, dimension (:,:), allocatable :: amp_r
!> average daily solar radiation for the month (MJ/m^2/day)
   real*8, dimension (:,:), allocatable :: solarav
!> standard deviation for avg monthly maximum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: tmpstdmx
!> normalization coefficient for precipitation generated from skewed
!> distribution (none)
   real*8, dimension (:,:), allocatable :: pcf
!> avg monthly minimum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: tmpmn
!> avg monthly maximum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: tmpmx
!> standard deviation for avg monthly minimum air temperature (deg C)
   real*8, dimension (:,:), allocatable :: tmpstdmn
   real*8, dimension (:,:), allocatable :: otmpstdmn,otmpmn,otmpmx
   real*8, dimension (:,:), allocatable :: otmpstdmx, ch_erodmo
   real*8, dimension (:,:), allocatable :: uh, hqdsave, hsdsave
!> probability of wet day after dry day in month (none)
   real*8, dimension (:,:), allocatable :: pr_w1
!> probability of wet day after wet day in month (none)
   real*8, dimension (:,:), allocatable :: pr_w2
!> proportion of wet days in the month (none)
   real*8, dimension (:,:), allocatable :: pr_w3
   real*8, dimension (:,:,:), allocatable :: pcp_stat
   real*8, dimension (:,:), allocatable :: opr_w1
   real*8, dimension (:,:), allocatable :: opr_w2
   real*8, dimension (:,:), allocatable :: opr_w3
   real*8, dimension (:,:,:), allocatable :: opcp_stat
!> precipitation category (none):\n
!> 1 precipitation <= 508 mm/yr\n
!> 2 precipitation > 508 and <= 1016 mm/yr\n
!> 3 precipitation > 1016 mm/yr
   integer, dimension (:), allocatable :: ireg
!> number of HRUs in subbasin (none)
   integer, dimension (:), allocatable :: hrutot
   integer, dimension (:), allocatable :: hru1
!> HRU relative humidity data code (gage # for relative humidity data used in
!> as HRU) (none)
   integer, dimension (:), allocatable :: ihgage
!> HRU solar radiation data code (record # for solar radiation used in HRU)
!> (none)
   integer, dimension (:), allocatable :: isgage
!> HRU wind speed gage data code (gage # for wind speed data used in HRU) (none)
   integer, dimension (:), allocatable :: iwgage
!> GIS code printed to output files (output.sub) (none
   integer, dimension (:), allocatable :: subgis
!> subbasin rain gage data code (gage # for rainfall data used in HRU) (none)
   integer, dimension (:), allocatable :: irgage
!> subbasin temp gage data code (gage # for temperature data used in HRU) (none)
   integer, dimension (:), allocatable :: itgage
!> (none) irelh = 0 (dewpoint)\n
!>        irelh = 1 (relative humidity)\n
!> note:  inputs > 1.0 (dewpoint)\n
!>        inputs < 1.0 (relative hum)
   integer, dimension (:), allocatable :: irelh
   integer, dimension (:), allocatable :: fcst_reg
!> amount of nitrogen stored in the active organic (humic) nitrogen pool
!> (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_aorgn
!> amount of nitrogen stored in the fresh organic (residue) pool (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_fon
   real*8, dimension (:,:), allocatable :: sol_tmp
!> available water capacity of soil layer (mm H20/mm soil)
   real*8, dimension (:,:), allocatable :: sol_awc
!> crack volume for soil layer (mm)
   real*8, dimension (:,:), allocatable :: volcr
   real*8, dimension (:,:), allocatable :: sol_prk
!> subbasin phosphorus percolation coefficient. Ratio of soluble phosphorus in
!> surface to soluble phosphorus in percolate
   real*8, dimension (:,:), allocatable :: pperco_sub
!> amount of phosphorus in the soil layer stored in the stable mineral
!> phosphorus pool(kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_stap
!> factor which converts kg/kg soil to kg/ha (none)
   real*8, dimension (:,:), allocatable :: conv_wt
!> amount of phosphorus stored in the active mineral phosphorus pool (kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_actp
!> soluble P concentration in top soil layer (mg P/kg soil) or\n
!> amount of phosohorus stored in solution. NOTE UNIT CHANGE! (kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_solp
!> maximum or potential crack volume (mm)
   real*8, dimension (:,:), allocatable :: crdep
!> amount of water available to plants in soil layer at field capacity (fc - wp)
!> (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_fc
!> amount of water held in the soil layer at saturation (sat - wp water)
!> (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_ul
!> bulk density of the soil (Mg/m^3)
   real*8, dimension (:,:), allocatable :: sol_bd
!> depth to bottom of soil layer (mm)
   real*8, dimension (:,:), allocatable :: sol_z
!> amount of water stored in the soil layer on any given day (less wp water)
!> (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_st
!> water content of soil at -0.033 MPa (field capacity) (mm H2O/mm soil)
   real*8, dimension (:,:), allocatable :: sol_up
!> percent clay content in soil material (UNIT CHANGE!) (% or none)
   real*8, dimension (:,:), allocatable :: sol_clay
!> beta coefficent to calculate hydraulic conductivity (none)
   real*8, dimension (:,:), allocatable :: sol_hk
   real*8, dimension (:,:), allocatable :: flat,sol_nh3
!> electrical conductivity of soil layer (dS/m)
   real*8, dimension (:,:), allocatable :: sol_ec
!> amount of nitrogen stored in the stable organic N pool. NOTE UNIT CHANGE!
!> (mg N/kg soil or kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_orgn
!> total porosity of soil layer expressed as a fraction of the total volume
!> (none)
   real*8, dimension (:,:), allocatable :: sol_por
!> water content of soil at -1.5 MPa (wilting point) (mm H20/mm soil)
   real*8, dimension (:,:), allocatable :: sol_wp
!> amount of phosphorus stored in the organic P pool. NOTE UNIT CHANGE!
!> (mg P/kg soil or kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_orgp
!> amount of organic matter in the soil layer classified as humic substances
!> (kg humus/ha)
   real*8, dimension (:,:), allocatable :: sol_hum
!> water content of soil at -1.5 MPa (wilting point) (mm H20)
   real*8, dimension (:,:), allocatable :: sol_wpmm
!> amount of nitrogen stored in the nitrate pool. This variable is read in as a
!> concentration and converted to kg/ha (this value is read from the .sol file
!> in units of mg/kg) (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_no3
!> percent organic carbon in soil layer (%)
   real*8, dimension (:,:), allocatable :: sol_cbn
!> saturated hydraulic conductivity of soil layer (mm/hour)
   real*8, dimension (:,:), allocatable :: sol_k
!> amount of organic matter in the soil layer classified as residue (kg/ha)
   real*8, dimension (:,:), allocatable :: sol_rsd
!> amount of phosphorus stored in the fresh organic (residue) pool (kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_fop
!> percent of rock fragments in soil layer (%)
   real*8, dimension (:,:), allocatable :: sol_rock
!> percent silt content in soil material (UNIT CHANGE!) (% or none)
   real*8, dimension (:,:), allocatable :: sol_silt
!> percent sand content of soil material (%)
   real*8, dimension (:,:), allocatable :: sol_sand
   real*8, dimension (:,:), allocatable :: orig_solno3,orig_solorgn
   real*8, dimension (:,:), allocatable :: orig_solsolp,orig_solorgp
   real*8, dimension (:,:), allocatable :: orig_soltmp,orig_solrsd
   real*8, dimension (:,:), allocatable :: orig_solfop,orig_solfon
   real*8, dimension (:,:), allocatable :: orig_solaorgn,orig_solst
   real*8, dimension (:,:), allocatable :: orig_solactp,orig_solstap
   real*8, dimension (:,:), allocatable :: orig_volcr
   real*8, dimension (:,:), allocatable :: conk
!> sol_pst(:,:,1) initial amount of pesticide in first layer read in from .chm
!> file (mg/kg)\n
!> sol_pst(:,:,:) amount of pesticide in layer. NOTE UNIT CHANGE!
!> (kg/ha)
   real*8, dimension (:,:,:), allocatable :: sol_pst
!> pesticide sorption coefficient, Kp; the ratio of the concentration in the
!> solid phase to the concentration in solution ((mg/kg)/(mg/L))
   real*8, dimension (:,:,:), allocatable :: sol_kp
   real*8, dimension (:,:,:), allocatable :: orig_solpst
   real*8, dimension (:), allocatable :: velsetlr, velsetlp
!> 1st shape parameter for reservoir surface area equation (none)
   real*8, dimension (:), allocatable :: br1
!> lake evaporation coefficient (none)
   real*8, dimension (:), allocatable :: evrsv
!> hydraulic conductivity of the reservoir bottom (mm/hr)
   real*8, dimension (:), allocatable :: res_k
!> pesticide concentration in lake water (mg/m^3)
   real*8, dimension (:), allocatable :: lkpst_conc
!> volume of water needed to fill the reservoir to the emergency spillway (read
!> in as 10^4 m^3 and converted to m^3) (m^3)
   real*8, dimension (:), allocatable :: res_evol
!> volume of water needed to fill the reservoir to the principal spillway (read
!> in as 10^4 m^3 and converted to m^3) (m^3)
   real*8, dimension (:), allocatable :: res_pvol
!> reservoir volume (read in as 10^4 m^3 and converted to m^3) (m^3)
   real*8, dimension (:), allocatable :: res_vol
!> reservoir surface area when reservoir is filled to principal spillway (ha)
   real*8, dimension (:), allocatable :: res_psa
!> pesticide reaction coefficient in lake water (1/day)
   real*8, dimension (:), allocatable :: lkpst_rea
!> pesticide volatilization coefficient in lake water (m/day)
   real*8, dimension (:), allocatable :: lkpst_vol
!> 2nd shape parameter for reservoir surface area equation (none)
   real*8, dimension (:), allocatable :: br2
!> average daily principal spillway release volume (read in as a release rate in
!> m^3/s and converted to m^3/day) (m^3/day)
   real*8, dimension (:), allocatable :: res_rr
!> amount of sediment in reservoir (read in as mg/L and converted to kg/L)
!> (kg/L)
   real*8, dimension (:), allocatable :: res_sed
!> pesticide partition coefficient between water and sediment in lake water
!> (m^3/g)
   real*8, dimension (:), allocatable :: lkpst_koc
!> mixing velocity (diffusion/dispersion) in lake water for pesticide (m/day)
   real*8, dimension (:), allocatable :: lkpst_mix
!> resuspension velocity in lake water for pesticide sorbed to sediment (m/day)
   real*8, dimension (:), allocatable :: lkpst_rsp
!> settling velocity in lake water for pesticide sorbed to sediment (m/day)
   real*8, dimension (:), allocatable :: lkpst_stl
!> pesticide concentration in lake bed sediment (mg/m^3)
   real*8, dimension (:), allocatable :: lkspst_conc
!> pesticide reaction coefficient in lake bed sediment (1/day)
   real*8, dimension (:), allocatable :: lkspst_rea
   real*8, dimension (:), allocatable :: theta_n, theta_p, con_nirr
   real*8, dimension (:), allocatable :: con_pirr
!> depth of active sediment layer in lake for for pesticide (m)
   real*8, dimension (:), allocatable :: lkspst_act
!> pesticide burial velocity in lake bed sediment (m/day)
   real*8, dimension (:), allocatable :: lkspst_bry
   real*8, dimension (:), allocatable :: sed_stlr
   real*8, dimension (7) :: resdata
!> normal amount of sediment in reservoir (read in as mg/L and convert to kg/L)
!> (kg/L)
   real*8, dimension (:), allocatable :: res_nsed
!> fraction of water removed from the reservoir via WURESN which is returned and
!> becomes flow from the reservoir outlet (none)
   real*8, dimension (:), allocatable :: wurtnf
!> chlorophyll-a production coefficient for reservoir (none)
   real*8, dimension (:), allocatable :: chlar
!> amount of nitrate in reservoir (kg N)
   real*8, dimension (:), allocatable :: res_no3
!> amount of organic N in reservoir (kg N)
   real*8, dimension (:), allocatable :: res_orgn
!> amount of organic P in reservoir (kg P)
   real*8, dimension (:), allocatable :: res_orgp
!> amount of soluble P in reservoir (kg P)
   real*8, dimension (:), allocatable :: res_solp
   real*8, dimension (:), allocatable :: res_chla,res_seci
!> reservoir surface area when reservoir is filled to emergency spillway (ha)
   real*8, dimension (:), allocatable :: res_esa
!> amount of ammonia in reservoir (kg N)
   real*8, dimension (:), allocatable :: res_nh3
!> amount of nitrite in reservoir (kg N)
   real*8, dimension (:), allocatable :: res_no2
!> water clarity coefficient for reservoir (none)
   real*8, dimension (:), allocatable :: seccir
   real*8, dimension (:), allocatable :: res_bactp, res_bactlp
!> minimum reservoir outflow as a fraction of the principal spillway volume
!> (0-1) (fraction)
   real*8, dimension (:), allocatable :: oflowmn_fps
!> target volume as a fraction of the principal spillway volume (.1-5)
!> (fraction)
   real*8, dimension (:), allocatable :: starg_fps
   real*8, dimension (:), allocatable :: weirc, weirk, weirw
   real*8, dimension (:), allocatable :: acoef, bcoef, ccoef
   real*8, dimension (:), allocatable :: orig_resvol,orig_ressed
   real*8, dimension (:), allocatable :: orig_lkpstconc,orig_lkspstconc
   real*8, dimension (:), allocatable :: orig_ressolp,orig_resorgp
   real*8, dimension (:), allocatable :: orig_resno3,orig_resno2
   real*8, dimension (:), allocatable :: orig_resnh3,orig_resorgn
!> minimum daily ouflow for the month (read in as m^3/s and converted to
!> m^3/day) (m^3/day)
   real*8, dimension (:,:), allocatable :: oflowmn
!> maximum daily ouflow for the month (read in as m^3/s and converted to
!> m^3/day) (m^3/day)
   real*8, dimension (:,:), allocatable :: oflowmx
!> monthly target reservoir storage (needed if IRESCO=2) (read in as 10^4 m^3
!> and converted to m^3) (m^3)
   real*8, dimension (:,:), allocatable :: starg
!> phosphorus settling rate for mid-year period (read in as m/year and converted
!> to m/day) (m/day)
   real*8, dimension (:), allocatable :: psetlr1
!> phosphorus settling rate for remainder of year (read in as m/year and
!> converted to m/day) (m/day)
   real*8, dimension (:), allocatable :: psetlr2
!> nitrogen settling rate for mid-year period (read in as m/year and converted
!> to m/day) (m/day)
   real*8, dimension (:), allocatable :: nsetlr1
!> nitrogen settling rate for remainder of year (read in as m/year and converted
!> to m/day) (m/day)
   real*8, dimension (:), allocatable :: nsetlr2
!> average amount of water withdrawn from reservoir each month for consumptive 
!> water use (read in as 10^4 m^3 and converted to m^3) (m^3)
   real*8, dimension (:,:), allocatable :: wuresn
!> measured average daily outflow from the reservoir for the month (needed if
!> IRESCO=1) (read in as m^3/s and converted to m^3/day) (m^3/day)
   real*8, dimension (:,:,:), allocatable :: res_out
!> number of subbasin reservoir is in (weather for the subbasin is used for the
!> reservoir) (none)
   integer, dimension (:), allocatable :: res_sub
!> beginning of mid-year nutrient settling "season" (none)
   integer, dimension (:), allocatable :: ires1
!> end of mid-year nutrient settling "season" (none)
   integer, dimension (:), allocatable :: ires2
!> outflow simulation code (none):\n
!> 0 compute outflow for uncontrolled reservoir with average annual release
!> rate\n
!> 1 measured monthly outflow\n
!> 2 simulated controlled outflow-target release\n
!> 3 measured daily outflow\n
!> 4 stage/volume/outflow relationship
   integer, dimension (:), allocatable :: iresco
!> year of the simulation that the reservoir becomes operational (none)
   integer, dimension (:), allocatable :: iyres
!> month the reservoir becomes operational (none)
   integer, dimension (:), allocatable :: mores
!> beginning month of non-flood season (needed if IRESCO=2) (none)
   integer, dimension (:), allocatable :: iflod1r
!> ending month of non-flood season (needed if IRESCO=2) (none)
   integer, dimension (:), allocatable :: iflod2r
!> number of days to reach target storage from current reservoir storage (needed
!> if IRESCO=2) (days)
   integer, dimension (:), allocatable :: ndtargr
!> application efficiency (0-1) (none)
   real*8, dimension (:), allocatable :: ap_ef
!> exponential of the rate constant for degradation of the pesticide on foliage
!> (none)
   real*8, dimension (:), allocatable :: decay_f
!> soil adsorption coefficient normalized for soil organic carbon content
!> ((mg/kg)/(mg/L))
   real*8, dimension (:), allocatable :: skoc
!> exponential of the rate constant for degradation of the pesticide in soil
!> (none)
   real*8, dimension (:), allocatable :: decay_s
!> half-life of pesticide on foliage (days)
   real*8, dimension (:), allocatable :: hlife_f
!> half-life of pesticide in soil (days)
   real*8, dimension (:), allocatable :: hlife_s
!> fraction of pesticide on foliage which is washed-off by a rainfall event
!> (none)
   real*8, dimension (:), allocatable :: pst_wof
!> solubility of chemical in water (mg/L (ppm))
   real*8, dimension (:), allocatable :: pst_wsol
   real*8, dimension (:), allocatable :: irramt
   real*8, dimension (:), allocatable :: phusw, phusw_nocrop
!> flag for types of pesticide used in watershed. Array location is pesticide ID
!> number\n
!> 0: pesticide not used\n
!> 1: pesticide used
   integer, dimension (:), allocatable :: pstflg
!> sequence number of pesticide in NPNO(:) (none)
   integer, dimension (:), allocatable :: nope
   integer, dimension (:), allocatable :: nop
   integer, dimension (:), allocatable :: yr_skip, isweep
   integer, dimension (:), allocatable :: icrmx, nopmx
! new management scehduling variables
   integer, dimension (:,:), allocatable :: mgtop, idop
   integer, dimension (:,:), allocatable :: mgt1iop,mgt2iop,mgt3iop
   real*8, dimension (:,:), allocatable ::  mgt4op, mgt5op, mgt6op
   real*8, dimension (:,:), allocatable :: mgt7op, mgt8op, mgt9op
   real*8, dimension (:,:), allocatable :: mgt10iop, phu_op
! mcrdb = maximum number of crops in database
!> fraction of nitrogen in yield (kg N/kg yield)
   real*8, dimension (:), allocatable :: cnyld
!> plant residue decomposition coefficient. The fraction of residue which will
!> decompose in a day assuming optimal moisture, temperature, C:N ratio, and C:P
!> ratio (none)
   real*8, dimension (:), allocatable :: rsdco_pl
!> 1st shape parameter for radiation use efficiency equation (none)
   real*8, dimension (:), allocatable :: wac21
!> 2nd shape parameter for radiation use efficiency equation (none)
   real*8, dimension (:), allocatable :: wac22
!> minimum LAI during winter dormant period \f$(m^2/m^2)\f$
   real*8, dimension (:), allocatable :: alai_min
!> 1st shape parameter for leaf area development equation (none)
   real*8, dimension (:), allocatable :: leaf1
!> 2nd shape parameter for leaf area development equation (none)
   real*8, dimension (:), allocatable :: leaf2
!> Value of harvest index between 0 and HVSTI which represents the lowest value
!> expected due to water stress ((kg/ha)/(kg/ha))
   real*8, dimension (:), allocatable :: wsyf
!> biomass-energy ratio. The potential (unstressed) growth rate per unit of
!> intercepted photosynthetically active radiation.((kg/ha)/(MJ/m**2))
   real*8, dimension (:), allocatable :: bio_e
!> harvest index: crop yield/aboveground biomass ((kg/ha)/(kg/ha))
   real*8, dimension (:), allocatable :: hvsti
!> minimum temperature for plant growth (deg C)
   real*8, dimension (:), allocatable :: t_base
!> optimal temperature for plant growth (deg C)
   real*8, dimension (:), allocatable :: t_opt
   real*8, dimension (:), allocatable :: chtmx !< maximum canopy height (m)
   real*8, dimension (:), allocatable :: cvm !< natural log of USLE_C (none)
!> maximum stomatal conductance (m/s)
   real*8, dimension (:), allocatable :: gsi
!> rate of decline in stomatal conductance per unit increase in vapor pressure
!> deficit ((m/s)*(1/kPa))
   real*8, dimension (:), allocatable :: vpd2
!> rate of decline in radiation use efficiency as a function of vapor pressure
!> deficit (none)
   real*8, dimension (:), allocatable :: wavp
!> fraction of leaf/needle biomass that drops during dormancy (for trees only)
!> (none)
   real*8, dimension (:), allocatable :: bio_leaf
!> maximum (potential) leaf area index (none)
   real*8, dimension (:), allocatable :: blai
!> fraction of phosphorus in yield (kg P/kg yield)
   real*8, dimension (:), allocatable :: cpyld
!> fraction of growing season when leaf area declines (none)
   real*8, dimension (:), allocatable :: dlai
   real*8, dimension (:), allocatable :: rdmx !< maximum root depth of plant (m)
!> 1st shape parameter for plant N uptake equation (none)
   real*8, dimension (:), allocatable :: bio_n1
!> 2nd shape parameter for plant N uptake equation (none)
   real*8, dimension (:), allocatable :: bio_n2
!> 1st shape parameter for plant P uptake equation (none)
   real*8, dimension (:), allocatable :: bio_p1
!> 2st shape parameter for plant P uptake equation (none)
   real*8, dimension (:), allocatable :: bio_p2
!> fraction above ground biomass that dies off at dormancy (fraction)
   real*8, dimension (:), allocatable :: bm_dieoff
   real*8, dimension (:), allocatable :: bmx_trees,ext_coef
!> initial root to shoot ratio at the beg of growing season
   real*8, dimension (:), allocatable :: rsr1
!> root to shoot ratio at the end of the growing season
   real*8, dimension (:), allocatable :: rsr2
!     real*8, dimension (:), allocatable :: air_str
!> nitrogen uptake parameter #1: normal fraction of N in crop biomass at
!> emergence (kg N/kg biomass)
   real*8, dimension (:), allocatable :: pltnfr1
!> nitrogen uptake parameter #2: normal fraction of N in crop biomass at 0.5
!> maturity (kg N/kg biomass)
   real*8, dimension (:), allocatable :: pltnfr2
!> nitrogen uptake parameter #3: normal fraction of N in crop biomass at
!> maturity (kg N/kg biomass)
   real*8, dimension (:), allocatable :: pltnfr3
!> phosphorus uptake parameter #1: normal fraction of P in crop biomass at
!> emergence (kg P/kg biomass)\n
   real*8, dimension (:), allocatable :: pltpfr1
!> phosphorus uptake parameter #2: normal fraction of P in crop biomass at 0.5
!> maturity (kg P/kg biomass)
   real*8, dimension (:), allocatable :: pltpfr2
!> phosphorus uptake parameter #3: normal fraction of P in crop biomass at
!> maturity (kg P/kg biomass)
   real*8, dimension (:), allocatable :: pltpfr3
!> crop/landcover category:\n
!> 1 warm season annual legume\n
!> 2 cold season annual legume\n
!> 3 perennial legume\n
!> 4 warm season annual\n
!> 5 cold season annual\n
!> 6 perennial\n
!> 7 trees
   integer, dimension (:), allocatable :: idc
   integer, dimension (:), allocatable :: mat_yrs
!> concentration of persistent bacteria in manure (fertilizer) (cfu/g manure)
   real*8, dimension (:), allocatable :: bactpdb
!> fraction of mineral N (NO3 + NH3) (kg minN/kg fert)
   real*8, dimension (:), allocatable :: fminn
!> fraction of organic N (kg orgN/kg fert)
   real*8, dimension (:), allocatable :: forgn
!> fraction of organic P (kg orgP/kg fert)
   real*8, dimension (:), allocatable :: forgp
!> bacteria partition coefficient (none):\n
!> 1: all bacteria in solution\n
!> 0: all bacteria sorbed to soil particles
   real*8, dimension (:), allocatable :: bactkddb
!> concentration of less persistent bacteria in manure (fertilizer)
!> (cfu/g manure)
   real*8, dimension (:), allocatable :: bactlpdb
!> fraction of mineral P (kg minP/kg fert)
   real*8, dimension (:), allocatable :: fminp
!> fraction of NH3-N in mineral N (kg NH3-N/kg minN)
   real*8, dimension (:), allocatable :: fnh3n
   character(len=8), dimension (200) :: fertnm !< name of fertilizer
!> curb length density in HRU (km/ha)
   real*8, dimension (:), allocatable :: curbden
!> maximum amount of solids allowed to build up on impervious surfaces
!> (kg/curb km)
   real*8, dimension (:), allocatable :: dirtmx
!> fraction of HRU area that is impervious (both directly and indirectly
!> connected)(fraction)
   real*8, dimension (:), allocatable :: fimp
!> wash-off coefficient for removal of constituents from an impervious surface
!> (1/mm)
   real*8, dimension (:), allocatable :: urbcoef
!> time for the amount of solids on impervious areas to build up to 1/2 the
!> maximum level (days)
   real*8, dimension (:), allocatable :: thalf
!> concentration of total nitrogen in suspended solid load from impervious
!> areas (mg N/kg sed)
   real*8, dimension (:), allocatable :: tnconc
!> concentration of NO3-N in suspended solid load from impervious areas
!> (mg NO3-N/kg sed)
   real*8, dimension (:), allocatable :: tno3conc
!> concentration of total phosphorus in suspended solid load from impervious
!> areas (mg P/kg sed)
   real*8, dimension (:), allocatable :: tpconc
!> fraction of HRU area that is classified as directly connected impervious
!> (fraction)
   real*8, dimension (:), allocatable :: fcimp
!> SCS curve number for moisture condition II in impervious areas (none)
   real*8, dimension (:), allocatable :: urbcn2
!> availability factor, the fraction of the curb length that is sweepable (none)
   real*8 :: fr_curb
   real*8 :: frt_kg !< amount of fertilizer applied to HRU (kg/ha)
!> depth of pesticide in the soil (mm)
   real*8 :: pst_dep
   real*8 :: sweepeff

   real*8, dimension (:), allocatable :: ranrns_hru
   integer, dimension (:), allocatable :: itill
!> depth of mixing caused by operation (mm)
   real*8, dimension (:), allocatable :: deptil
!> mixing efficiency of operation (none)
   real*8, dimension (:), allocatable :: effmix
!> random roughness of a given tillage operation (mm)
   real*8, dimension (:), allocatable :: ranrns
!> 8-character name for the tillage operation
   character(len=8), dimension (550) :: tillnm
!> For ICODES equal to (none)\n
!> 0,1,3,5,9: not used\n
!> 2: Fraction of flow in channel\n
!> 4: amount of water transferred (as defined by INUM4S)\n
!> 7,8,10,11: drainage area in square kilometers associated with the record
!> file\n
!> 12: rearation coefficient
   real*8, dimension (:), allocatable :: rnum1s
!> total drainage area of hydrograph in square kilometers (km^2)
   real*8, dimension (:), allocatable :: hyd_dakm
   real*8, dimension (:,:), allocatable :: varoute,shyd, vartran
   real*8, dimension (:,:,:), allocatable :: hhvaroute
!> routing command code (none):\n
!> 0 = finish\n
!> 1 = subbasin\n
!> 2 = route\n
!> 3 = routres\n
!> 4 = transfer\n
!> 5 = add\n
!> 6 = rechour\n
!> 7 = recmon\n
!> 8 = recyear\n
!> 9 = save\n
!> 10 = recday\n
!> 11 = reccnst\n
!> 12 = structure\n
!> 13 = apex\n
!> 14 = saveconc\n
!> 15 =
   integer, dimension (:), allocatable :: icodes
!> For ICODES equal to (none)\n
!> 0: not used\n
!> 1,2,3,5,7,8,10,11: hydrograph storage location number\n
!> 4: departure type (1=reach, 2=reservoir)\n
!> 9: hydrograph storage location of data to be printed to event file\n
!> 14:hydrograph storage location of data to be printed to saveconc file
   integer, dimension (:), allocatable :: ihouts
!> For ICODES equal to (none)\n
!> 0: not used\n
!> 1: subbasin number\n
!> 2: reach number\n
!> 3: reservoir number\n
!> 4: reach or res # flow is diverted from\n
!> 5: hydrograph storage location of 1st dataset to be added\n
!> 7,8,9,10,11,14: file number
   integer, dimension (:), allocatable :: inum1s
!> For ICODES equal to (none)\n
!> 0,1,7,8,10,11: not used\n
!> 2,3: inflow hydrograph storage location\n
!> 4: destination type (1=reach, 2=reservoir)\n
!> 5: hydrograph storage location of 2nd dataset to be added\n
!> 9,14:print frequency (0=daily, 1=hourly)
   integer, dimension (:), allocatable :: inum2s
!> For ICODES equal to (none)\n
!> 0,1,5,7,8,10,11: not used\n
!> 2,3: subbasin number
!> 4: destination number. Reach or reservoir receiving water\n
!> 9: print format (0=normal, fixed format; 1=txt format for AV interface,
!> recday)
   integer, dimension (:), allocatable :: inum3s
!> For ICODES equal to (none)\n
!> 0,2,3,5,7,8,9,10,11: not used\n
!> 1: GIS code printed to output file (optional)\n
!> 4: rule code governing transfer of water (1=fraction transferred out, 2=min
!> volume or flow left, 3=exact amount transferred)
   integer, dimension (:), allocatable :: inum4s
   integer, dimension (:), allocatable :: inum5s,inum6s,inum7s,inum8s
   integer, dimension (:), allocatable :: subed
   character(len=10), dimension (:), allocatable :: recmonps
   character(len=10), dimension (:), allocatable :: reccnstps
   character(len=5), dimension (:), allocatable :: subnum
   character(len=4), dimension (:), allocatable :: hruno

!> Mannings's n for grassed waterway (none)
   real*8, dimension (:), allocatable :: grwat_n
!> flag for the simulation of grass waterways (none)\n
!> = 0 inactive\n
!> = 1 active
   real*8, dimension (:), allocatable :: grwat_i
!> length of grass waterway (km)
   real*8, dimension (:), allocatable :: grwat_l
!> average width of grassed waterway (m)
   real*8, dimension (:), allocatable :: grwat_w
!> depth of grassed waterway from top of bank to bottom (m)
   real*8, dimension (:), allocatable :: grwat_d
!> average slope of grassed waterway channel (m)
   real*8, dimension (:), allocatable :: grwat_s
!> linear parameter for calculating sediment in grassed waterways (none)
   real*8, dimension (:), allocatable :: grwat_spcon
   real*8, dimension (:), allocatable :: tc_gwat
   real*8, dimension (:), allocatable :: pot_volmm,pot_tilemm,pot_volxmm  !!NUBZ
!> fraction of HRU area that drains into pothole (km^2/km^2)
   real*8, dimension (:), allocatable :: pot_fr
!> average daily outflow to main channel from tile flow if drainage tiles are
!> installed in pothole (needed only if current HRU is IPOT) (m^3/s)
   real*8, dimension (:), allocatable :: pot_tile
!> initial or current volume of water stored in the depression/impounded area
!> (read in as mm and converted to m^3) (needed only if current HRU is IPOT)
!> (mm or m^3 H20)
   real*8, dimension (:), allocatable :: pot_vol
   real*8, dimension (:), allocatable :: potsa
!> maximum volume of water stored in the depression/impounded area (read in as
!> mm and converted to m^3) (needed only if current HRU is IPOT) (mm)
   real*8, dimension (:), allocatable :: pot_volx
!> wetting front matric potential (mm)
   real*8, dimension (:), allocatable :: wfsh
   real*8, dimension (:), allocatable :: potflwi,potsedi
!> nitrate decay rate in impounded area (1/day)
   real*8, dimension (:), allocatable :: pot_no3l
!> normal sediment concentration in impounded water (needed only if current HRU
!> is IPOT)(mg/L)
   real*8, dimension (:), allocatable :: pot_nsed
!> nitrate-N concentration in groundwater loading to reach (mg N/L)
   real*8, dimension (:), allocatable :: gwno3
   real*8, dimension (:), allocatable :: newrti
!> reduction in bacteria loading from filter strip (none)
   real*8, dimension (:), allocatable :: fsred
   real*8, dimension (:), allocatable :: pot_sed,pot_no3,tmpavp
!> average distance to stream (m)
   real*8, dimension (:), allocatable :: dis_stream
!> pothole evaporation coefficient (none)
   real*8, dimension (:), allocatable :: evpot
   real*8, dimension (:), allocatable :: pot_solpl
   real*8, dimension (:), allocatable :: sed_con, orgn_con, orgp_con
!> hydraulic conductivity of soil surface of pothole [defaults to conductivity
!> of upper soil (0.01--10.) layer] (mm/hr)
   real*8, dimension (:), allocatable :: pot_k
   real*8, dimension (:), allocatable :: soln_con, solp_con
!> nitrogen uptake reduction factor (not currently used; defaulted 300.)
   real*8, dimension (:), allocatable :: n_reduc
!> lag coefficient for calculating nitrate concentration in subsurface drains
!> (0.001 - 1.0) (dimensionless)
   real*8, dimension (:), allocatable :: n_lag
!> power function exponent for calculating nitrate concentration in subsurface
!> drains (1.0 - 3.0) (dimensionless)
   real*8, dimension (:), allocatable :: n_ln
!> coefficient for power function for calculating nitrate concentration in
!> subsurface drains (0.5 - 4.0) (dimensionless)
   real*8, dimension (:), allocatable :: n_lnco
   integer, dimension (:), allocatable :: ioper
   integer, dimension (:), allocatable :: ngrwat
!> USLE equation length slope (LS) factor (none)
   real*8, dimension (:), allocatable :: usle_ls
!> filter strip width for bacteria transport (m)
   real*8, dimension (:), allocatable :: filterw
!> fraction of plant heat units accumulated (none)
   real*8, dimension (:), allocatable :: phuacc
!> sum of all tillage mixing efficiencies for HRU operation (none)
   real*8, dimension (:), allocatable :: sumix
!> plant water uptake compensation factor (0-1) (none)
   real*8, dimension (:), allocatable :: epco
!> soil evaporation compensation factor (0-1) (none)
   real*8, dimension (:), allocatable :: esco
!> average slope steepness (m/m)
   real*8, dimension (:), allocatable :: hru_slp
!> average slope length for subbasin (m)
   real*8, dimension (:), allocatable :: slsubbsn
!> organic N enrichment ratio, if left blank the model will calculate for every
!> event (none)
   real*8, dimension (:), allocatable :: erorgn
!> organic P enrichment ratio, if left blank the model will calculate for every
!> event (none)
   real*8, dimension (:), allocatable :: erorgp
!> biological mixing efficiency. Mixing of soil due to activity of earthworms
!> and other soil biota. Mixing is performed at the end of every calendar year
!> (none)
   real*8, dimension (:), allocatable :: biomix
   real*8, dimension (:), allocatable :: pnd_seci
!> maximum canopy storage (mm H2O)
   real*8, dimension (:), allocatable :: canmx
!> maximum daily irrigation diversion from the reach (when IRRSC=1 or IRR=3):
!> when value is positive the units are mm H2O; when the value is negative, the
!> units are (10^4 m^3 H2O) (mm H2O or 10^4 m^3 H2O)
   real*8, dimension (:), allocatable :: divmax
!> minimum instream flow for irrigation diversions when IRRSC=1, irrigation
!> water will be diverted only when streamflow is at or above FLOWMIN (m^3/s)
   real*8, dimension (:), allocatable :: flowmin
!> USLE equation support practice (P) factor (none)
   real*8, dimension (:), allocatable :: usle_p
!> sediment concentration in lateral flow (g/L)
   real*8, dimension (:), allocatable :: lat_sed
!> total drainage area contributing to flow at the outlet (pour point) of the
!> reach in square kilometers (km^2)
   real*8, dimension (:), allocatable :: rch_dakm
   real*8, dimension (:), allocatable :: pnd_no3s,cn1
!> lateral flow travel time or exponential of the lateral flow travel time
!> (days or none)
   real*8, dimension (:), allocatable :: lat_ttime
!> SCS runoff curve number for moisture condition II (none)
   real*8, dimension (:), allocatable :: cn2
!> fraction of available flow in reach that is allowed to be applied to the HRU
!> (none)
   real*8, dimension (:), allocatable :: flowfr
   real*8, dimension (:), allocatable :: sol_zmx !< maximum rooting depth (mm)
!> exponential of the tile flow travel time (none)
   real*8, dimension (:), allocatable :: tile_ttime
!> slope length for lateral subsurface flow (m)
   real*8, dimension (:), allocatable :: slsoil
!> soluble P concentration in groundwater loading to reach (mg P/L)
   real*8, dimension (:), allocatable :: gwminp
!> amount of residue on soil surface (kg/ha)
   real*8, dimension (:), allocatable :: sol_cov
!> fraction of sediment remaining suspended in impoundment after settling for
!> one day (kg/kg)
   real*8, dimension (:), allocatable :: sed_stl
!> Manning's "n" value for overland flow (none)
   real*8, dimension (:), allocatable :: ov_n
!> amount of nitrate in pond (kg N)
   real*8, dimension (:), allocatable :: pnd_no3
!> amount of soluble P in pond (kg P)
   real*8, dimension (:), allocatable :: pnd_solp
!> annual yield (dry weight) in the HRU (metric tons/ha)
   real*8, dimension (:), allocatable :: yldanu
!> coefficient for pesticide drift directly onto stream (none)
   real*8, dimension (:), allocatable :: driftco
!> amount of organic N in pond (kg N)
   real*8, dimension (:), allocatable :: pnd_orgn
!> amount of organic P in pond (kg P)
   real*8, dimension (:), allocatable :: pnd_orgp
   real*8, dimension (:), allocatable :: cn3
   real*8, dimension (:), allocatable :: twlpnd, twlwet               !!srini pond/wet infiltration to shallow gw storage
!> fraction of subbasin area contained in HRU (km^2/km^2)
   real*8, dimension (:), allocatable :: hru_fr
!> amount of water held in soil profile at saturation (mm H2O)
   real*8, dimension (:), allocatable :: sol_sumul
   real*8, dimension (:), allocatable :: pnd_chla
!> area of HRU in square kilometers (km^2)
   real*8, dimension (:), allocatable :: hru_km
   real*8, dimension (:), allocatable :: bio_ms !< cover/crop biomass (kg/ha)
!> albedo when soil is moist (none)
   real*8, dimension (:), allocatable :: sol_alb
   real*8, dimension (:), allocatable :: strsw
!> fraction of HRU/subbasin area that drains into ponds (none)
   real*8, dimension (:), allocatable :: pnd_fr
!> hydraulic conductivity through bottom of ponds (mm/hr)
   real*8, dimension (:), allocatable :: pnd_k
!> surface area of ponds when filled to principal spillway (ha)
   real*8, dimension (:), allocatable :: pnd_psa
!> runoff volume from catchment area needed to fill the ponds to the principal
!> spillway (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: pnd_pvol
!> surface area of ponds when filled to emergency spillway (ha)
   real*8, dimension (:), allocatable :: pnd_esa
!> runoff volume from catchment area needed to fill the ponds to the emergency
!> spillway (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: pnd_evol
!> volume of water in ponds (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: pnd_vol
!> average annual yield in the HRU (metric tons)
   real*8, dimension (:), allocatable :: yldaa
!> normal sediment concentration in pond water (UNIT CHANGE!) (mg/kg or kg/kg)
   real*8, dimension (:), allocatable :: pnd_nsed
!> sediment concentration in pond water (UNIT CHANGE!) (mg/kg or kg/kg)
   real*8, dimension (:), allocatable :: pnd_sed
   real*8, dimension (:), allocatable :: strsa,dep_imp
   real*8, dimension (:), allocatable :: evpnd, evwet
!> fraction of HRU/subbasin area that drains into wetlands (none)
   real*8, dimension (:), allocatable :: wet_fr
!> hydraulic conductivity of bottom of wetlands (mm/hr)
   real*8, dimension (:), allocatable :: wet_k
!> surface area of wetlands in subbasin at normal water level (ha)
   real*8, dimension (:), allocatable :: wet_nsa
!> runoff volume from catchment area needed to fill wetlands to normal water
!> level (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: wet_nvol
   integer, dimension (:), allocatable :: iwetgw, iwetile
!> surface area of wetlands at maximum water level (ha)
   real*8, dimension (:), allocatable :: wet_mxsa
!> runoff volume from catchment area needed to fill wetlands to maximum water
!> level (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: wet_mxvol
!> volume of water in wetlands (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: wet_vol
!> normal sediment concentration in wetland water (UNIT CHANGE!)
!> (mg/kg or kg/kg)
   real*8, dimension (:), allocatable :: wet_nsed
!> sediment concentration in wetland water (UNIT CHANGE!) (mg/L or kg/L)
   real*8, dimension (:), allocatable :: wet_sed
!> 1st shape parameter for pond surface area equation (none)
   real*8, dimension (:), allocatable :: bp1
!> 2nd shape parameter for the pond surface area equation (none)
   real*8, dimension (:), allocatable :: bp2
!> retention coefficient for CN method based on plant ET (none)
   real*8, dimension (:), allocatable :: sci
!> retention coefficient for CN method based on soil moisture (none)
   real*8, dimension (:), allocatable :: smx
!> 1st shape parameter for the wetland surface area equation (none)
   real*8, dimension (:), allocatable :: bw1
!> 2nd shape parameter for the wetland surface area equation (none)
   real*8, dimension (:), allocatable :: bw2
   real*8, dimension (:), allocatable :: bactpq
   real*8, dimension (:), allocatable :: bactp_plt,bactlp_plt,cnday
!> fertilizer application efficiency calculated as the amount of N applied
!> divided by the amount of N removed at harvest (none)
   real*8, dimension (:), allocatable :: auto_eff
!> water clarity coefficient for wetland (none)
   real*8, dimension (:), allocatable :: secciw
!> amount of water stored in soil profile on any given day (mm H2O)
   real*8, dimension (:), allocatable :: sol_sw
   real*8, dimension (:), allocatable :: bactlpq
!> chlorophyll-a production coefficient for wetland (none)
   real*8, dimension (:), allocatable :: chlaw
!> average temperature for the day in HRU (deg C)
   real*8, dimension (:), allocatable :: tmpav
   real*8, dimension (:), allocatable :: bactps,bactlps
!> amount of water stored as snow (mm H2O)
   real*8, dimension (:), allocatable :: sno_hru
!> amount of organic N in wetland (kg N)
   real*8, dimension (:), allocatable :: wet_orgn
!> solar radiation for the day in HRU (MJ/m^2)
   real*8, dimension (:), allocatable :: hru_ra
!> precipitation for the day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: subp
   real*8, dimension (:), allocatable :: rsdin !< initial residue cover (kg/ha)
!> minimum temperature for the day in HRU (deg C)
   real*8, dimension (:), allocatable :: tmn
!> maximum temperature for the day in HRU (deg C)
   real*8, dimension (:), allocatable :: tmx
   real*8, dimension (:), allocatable :: tmp_hi,tmp_lo
!> USLE equation soil erodibility (K) factor (none)
   real*8, dimension (:), allocatable :: usle_k
!> time of concentration for HRU (hour)
   real*8, dimension (:), allocatable :: tconc
!> maximum possible solar radiation for the day in HRU (MJ/m^2)
   real*8, dimension (:), allocatable :: hru_rmx
   real*8, dimension (:), allocatable :: rwt,olai
   real*8, dimension (:), allocatable :: usle_cfac,usle_eifac
!> amount of water held in soil profile at field capacity (mm H2O)
   real*8, dimension (:), allocatable :: sol_sumfc
!> time for flow from farthest point in subbasin to enter a channel (hour)
   real*8, dimension (:), allocatable :: t_ov
!> total amount of NO3 applied during the year in auto-fertilization (kg N/ha)
   real*8, dimension (:), allocatable :: anano3
   real*8, dimension (:), allocatable :: aird
!> amount of organic P in wetland (kg P)
   real*8, dimension (:), allocatable :: wet_orgp
!> average porosity for entire soil profile (none)
   real*8, dimension (:), allocatable :: sol_avpor
!> product of USLE K,P,LS,exp(rock) (none)
   real*8, dimension (:), allocatable :: usle_mult
!> relative humidity for the day in HRU (none)
   real*8, dimension (:), allocatable :: rhd
!> wind speed for the day in HRU (m/s)
   real*8, dimension (:), allocatable :: u10
   real*8, dimension (:), allocatable :: aairr,cht
!> maximum leaf area index for the entire period of simulation in the HRU (none)
   real*8, dimension (:), allocatable :: lai_aamx
   real*8, dimension (:), allocatable :: shallirr,deepirr
!> longest tributary channel length in subbasin (km)
   real*8, dimension (:), allocatable :: ch_l1
!> amount of nitrate in wetland (kg N)
   real*8, dimension (:), allocatable :: wet_no3
   real*8, dimension (:), allocatable :: canstor,ovrlnd
!> maximum irrigation amount per auto application (mm)
   real*8, dimension (:), allocatable :: irr_mx
!> water stress factor which triggers auto irrigation (none or mm)
   real*8, dimension (:), allocatable :: auto_wstr
!> fertilizer/manure id number from database (none)
   real*8, dimension (:), allocatable :: cfrt_id
!> amount of fertilzier applied to HRU on a given day (kg/ha)
   real*8, dimension (:), allocatable :: cfrt_kg
   real*8, dimension (:), allocatable :: cpst_id
   real*8, dimension (:), allocatable :: cpst_kg
   real*8, dimension (:), allocatable :: irr_asq !< surface runoff ratio
   real*8, dimension (:), allocatable :: irr_eff
!> surface runoff ratio (0-1) .1 is 10% surface runoff (frac)
   real*8, dimension (:), allocatable :: irrsq
   real*8, dimension (:), allocatable :: irrefm, irrsalt
!> dry weight of biomass removed by grazing daily ((kg/ha)/day)
   real*8, dimension (:), allocatable :: bio_eat
!> dry weight of biomass removed by trampling daily ((kg/ha)/day)
   real*8, dimension (:), allocatable :: bio_trmp             !!NUBZ
   integer, dimension (:), allocatable :: ifrt_freq,ipst_freq,irr_noa
   integer, dimension (:), allocatable :: irr_sc,irr_no
!> release/impound action code (none):\n
!> 0 begin impounding water\n
!> 1 release impounded water
   integer, dimension (:), allocatable :: imp_trig
   integer, dimension (:), allocatable :: fert_days,irr_sca
!> land cover/crop identification code for first crop grown in HRU (the only
!> crop if there is no rotation) (from crop.dat) (none)
   integer, dimension (:), allocatable :: idplt
   integer, dimension (:), allocatable :: pest_days, wstrs_id
   real*8, dimension (:,:), allocatable :: bio_aahv
   real*8, dimension (:), allocatable :: cumei,cumeira
   real*8, dimension (:), allocatable :: cumrt, cumrai
!> amount of soluble P in wetland (kg P)
   real*8, dimension (:), allocatable :: wet_solp
   real*8, dimension (:), allocatable :: wet_no3s,wet_chla
   real*8, dimension (:), allocatable :: wet_seci,pnd_no3g,pstsol
!> groundwater delay: time required for water leaving the bottom of the root
!> zone to reach the shallow aquifer (days)
   real*8, dimension (:), allocatable :: delay
   real*8, dimension (:), allocatable :: gwht !< groundwater height (m)
!> groundwater contribution to streamflow from HRU on current day (mm H2O)
   real*8, dimension (:), allocatable :: gw_q
   real*8, dimension (:), allocatable :: pnd_solpg
!> alpha factor for groundwater recession curve (1/days)
   real*8, dimension (:), allocatable :: alpha_bf
!> \f$\exp(-alpha_bf)\f$ (none)
   real*8, dimension (:), allocatable :: alpha_bfe
!> specific yield for shallow aquifer (m^3/m^3)
   real*8, dimension (:), allocatable :: gw_spyld
!> alpha factor for groudwater recession curve of the deep aquifer (1/days)
   real*8, dimension (:), allocatable :: alpha_bf_d
!> \f$\exp(-alpha_bf_d)\f$ for deep aquifer (none)
   real*8, dimension (:), allocatable :: alpha_bfe_d
   real*8, dimension (:), allocatable :: gw_qdeep
!> \f$\exp(-1/delay)\f$ (none)
   real*8, dimension (:), allocatable :: gw_delaye
!> revap coeff: this variable controls the amount of water moving from the
!> shallow aquifer to the root zone as a result of soil moisture depletion
!> (none)
   real*8, dimension (:), allocatable :: gw_revap
!> recharge to deep aquifer: the fraction of root zone percolation that reaches
!> the deep aquifer (none)
   real*8, dimension (:), allocatable :: rchrg_dp
!> fraction of porosity from which anions are excluded
   real*8, dimension (:), allocatable :: anion_excl
!> threshold depth of water in shallow aquifer required to allow revap to occur
!> (mm H2O)
   real*8, dimension (:), allocatable :: revapmn
   real*8, dimension (:), allocatable :: rchrg
!> minimum plant biomass for grazing (kg/ha)
   real*8, dimension (:), allocatable :: bio_min
!> initial HRU soil water content expressed as fraction of field capacity (none)
   real*8, dimension (:), allocatable :: ffc
   real*8, dimension (:), allocatable :: surqsolp
!> depth of water in deep aquifer (mm H2O)
   real*8, dimension (:), allocatable :: deepst
!> depth of water in shallow aquifer (mm H2O)
   real*8, dimension (:), allocatable :: shallst
   real*8, dimension (:), allocatable :: cklsp,wet_solpg
   real*8, dimension (:), allocatable :: rchrg_src
!> filter strip trapping efficiency (used for everything but bacteria) (none)
   real*8, dimension (:), allocatable :: trapeff
!> average bulk density for soil profile (Mg/m^3)
   real*8, dimension (:), allocatable :: sol_avbd
   real*8, dimension (:), allocatable :: wet_no3g
!> time to drain soil to field capacity yield used in autofertilization (hours)
   real*8, dimension (:), allocatable :: tdrain
!> threshold depth of water in shallow aquifer required before groundwater flow
!> will occur (mm H2O)
   real*8, dimension (:), allocatable :: gwqmn
   real*8, dimension (:), allocatable :: pplnt,snotmp
!> drain tile lag time: the amount of time between the transfer of water from
!> the soil to the drain tile and the release of the water from the drain tile
!> to the reach (hours)
   real*8, dimension (:), allocatable :: gdrain
!> depth to the sub-surface drain (mm)
   real*8, dimension (:), allocatable :: ddrain
!> crack volume potential of soil (none)
   real*8, dimension (:), allocatable :: sol_crk
!> fraction of surface runoff within the subbasin which takes 1 day or less to
!> reach the subbasin outlet (none)
   real*8, dimension (:), allocatable :: brt
!> day length (hours)
   real*8, dimension (:), allocatable :: dayl
!> static maximum depressional storage; read from .sdr (mm)
   real*8, dimension (:), allocatable :: sstmaxd
   real*8, dimension (:), allocatable :: re !< effective radius of drains (mm)
!> distance between two drain tubes or tiles (mm)
   real*8, dimension (:), allocatable :: sdrain
   real*8, dimension (:), allocatable :: ddrain_hru
!> drainage coefficient (mm/day)
   real*8, dimension (:), allocatable :: drain_co
!> multiplication factor to determine conk(j1,j) from sol_k(j1,j) for HRU (none)
   real*8, dimension (:), allocatable :: latksatf
!> pump capacity (default pump capacity = 1.042mm/hr or 25mm/day) (mm/hr)
   real*8, dimension (:), allocatable :: pc
   real*8, dimension (:), allocatable :: stmaxd
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd3
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd2
   real*8, dimension (:), allocatable :: twash,sol_cnsw,doxq
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd8
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd9
   real*8, dimension (:), allocatable :: percn,sol_sumwp
!> total amount of water entering main channel for day from HRU (mm H2O)
   real*8, dimension (:), allocatable :: qdr
   real*8, dimension (:), allocatable :: tauton,tautop,cbodu,chl_a
   real*8, dimension (:), allocatable :: tfertn,tfertp,tgrazn,tgrazp
!> total lateral flow in soil profile for the day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: latq
   real*8, dimension (:), allocatable :: latno3,minpgw,no3gw,nplnt
   real*8, dimension (:), allocatable :: tileq, tileno3
   real*8, dimension (:), allocatable :: sedminpa,sedminps,sedorgn
!> soil loss for day in HRU (metric tons)
   real*8, dimension (:), allocatable :: sedyld
   real*8, dimension (:), allocatable :: sedorgp,sepbtm,strsn
!> surface runoff generated on day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: surfq
   real*8, dimension (:), allocatable :: strsp,strstmp,surqno3
   real*8, dimension (:), allocatable :: hru_ha !< area of HRU in hectares (ha)
!> fraction of total watershed area contained in HRU (km2/km2)
   real*8, dimension (:), allocatable :: hru_dafr
   real*8, dimension (:), allocatable :: tcfrtn,tcfrtp
!> atmospheric dry deposition of nitrates (kg/ha/yr)
   real*8, dimension (:), allocatable :: drydep_no3
!> atmospheric dry deposition of ammonia (kg/ha/yr)
   real*8, dimension (:), allocatable :: drydep_nh4
!> annual biomass (dry weight) in the HRU (metric tons/ha)
   real*8, dimension (:), allocatable :: bio_yrms
!> base zero total heat units (used when no land cover is growing) (heat units)
   real*8, dimension (:), allocatable :: phubase
   real*8, dimension (:), allocatable :: hvstiadj
   real*8, dimension (:), allocatable :: laiday !< leaf area index (m^2/m^2)
!> chlorophyll-a production coefficient for pond (none)
   real*8, dimension (:), allocatable :: chlap
   real*8, dimension (:), allocatable :: laimxfr,pnd_psed
!> water clarity coefficient for pond (none)
   real*8, dimension (:), allocatable :: seccip
   real*8, dimension (:), allocatable :: wet_psed,plantn,plt_et
!> average annual biomass in the HRU (metric tons)
   real*8, dimension (:), allocatable :: bio_aams
   real*8, dimension (:), allocatable :: plt_pet,plantp
!> time threshold used to define dormant period for plant (when daylength is
!> within the time specified by dl from the minimum daylength for the area, the
!> plant will go dormant) (hour)
   real*8, dimension (:), allocatable :: dormhr
!> maximum leaf area index for the year in the HRU (none)
   real*8, dimension (:), allocatable :: lai_yrmx
   real*8, dimension (:), allocatable :: bio_aamx
   real*8, dimension (:), allocatable :: lat_pst
!> fraction of HRU area that drains into floodplain (km^2/km^2)
   real*8, dimension (:), allocatable :: fld_fr
   real*8, dimension (:), allocatable :: orig_snohru,orig_potvol
   real*8, dimension (:), allocatable :: orig_alai,orig_bioms,pltfr_n
   real*8, dimension (:), allocatable :: orig_phuacc,orig_sumix,pltfr_p
!> total number of heat units to bring plant to maturity (heat units)
   real*8, dimension (:), allocatable :: phu_plt
   real*8, dimension (:), allocatable :: orig_phu
   real*8, dimension (:), allocatable :: orig_shallst,orig_deepst
!> fraction of HRU area that drains into riparian zone (km^2/km^2)
   real*8, dimension (:), allocatable :: rip_fr
   real*8, dimension (:), allocatable :: orig_pndvol,orig_pndsed
   real*8, dimension (:), allocatable :: orig_pndno3,orig_pndsolp
   real*8, dimension (:), allocatable :: orig_pndorgn,orig_pndorgp
   real*8, dimension (:), allocatable :: orig_wetvol,orig_wetsed
   real*8, dimension (:), allocatable :: orig_wetno3,orig_wetsolp
   real*8, dimension (:), allocatable :: orig_wetorgn,orig_wetorgp
   real*8, dimension (:), allocatable :: orig_solcov,orig_solsw
   real*8, dimension (:), allocatable :: orig_potno3,orig_potsed
   real*8, dimension (:), allocatable :: wtab,wtab_mn,wtab_mx
!> nitrate concentration in shallow aquifer converted to kg/ha (ppm NO3-N)
   real*8, dimension (:), allocatable :: shallst_n
   real*8, dimension (:), allocatable :: gw_nloss,rchrg_n
   real*8, dimension (:), allocatable :: det_san, det_sil, det_cla
   real*8, dimension (:), allocatable :: det_sag, det_lag
!> fraction of fertilizer which is applied to top 10 mm of soil (the remaining
!> fraction is applied to first soil layer) (none)
   real*8, dimension (:), allocatable :: afrt_surface
   real*8, dimension (:), allocatable :: tnylda
!> fraction of fertilizer which is applied to the top 10 mm of soil (the
!> remaining fraction is applied to the first soil layer) (none)
   real*8 :: frt_surface
!> maximum NO3-N content allowed to be applied in one year (kg NO3-N/ha)
   real*8, dimension (:), allocatable :: auto_nyr
!> maximum NO3-N content allowed in one fertilizer application (kg NO3-N/ha)
   real*8, dimension (:), allocatable :: auto_napp
!> nitrogen stress factor which triggers auto fertilization (none)
   real*8, dimension (:), allocatable :: auto_nstrs
   real*8, dimension (:), allocatable :: manure_kg
   real*8, dimension (:,:), allocatable :: rcn_mo, rammo_mo
   real*8, dimension (:,:), allocatable :: drydep_no3_mo, drydep_nh4_mo
   real*8, dimension (:), allocatable :: rcn_d, rammo_d
   real*8, dimension (:), allocatable :: drydep_no3_d, drydep_nh4_d
   real*8, dimension (:,:), allocatable :: yldn
   real*8, dimension (:,:), allocatable :: gwati, gwatn, gwatl
   real*8, dimension (:,:), allocatable :: gwatw, gwatd, gwatveg
   real*8, dimension (:,:), allocatable :: gwata, gwats, gwatspcon
   real*8, dimension (:,:), allocatable :: rfqeo_30d,eo_30d
!> phosphorus settling rate for 1st season (m/day)
   real*8, dimension (:), allocatable :: psetlp1
!> phosphorus settling rate for 2nd seaso (m/day)n
   real*8, dimension (:), allocatable :: psetlp2
!> previous value of wgncur(:,:) (none)
   real*8, dimension (:,:), allocatable :: wgnold
!> parameter to predict the impact of precip on other weather attributes
!> (none)\n
!> wgncur(1,:) parameter which predicts impact of precip on daily maximum air
!> temperature\n
!> wgncur(2,:) parameter which predicts impact of precip on daily minimum air
!> temperature\n
!> wgncur(3,:) parameter which predicts impact of precip on daily solar
!> radiation
   real*8, dimension (:,:), allocatable :: wgncur
   real*8, dimension (:,:), allocatable :: wrt
!> pesticide enrichment ratio (none)
   real*8, dimension (:,:), allocatable :: pst_enr
   real*8, dimension (:,:), allocatable :: zdb,pst_surq
!> pesticide on plant foliage (kg/ha)
   real*8, dimension (:,:), allocatable :: plt_pst
!> phosphorus settling rate for 1st season (m/day)
   real*8, dimension (:), allocatable :: psetlw1
!> phosphorus settling rate for 2nd season (m/day)
   real*8, dimension (:), allocatable :: psetlw2
   real*8, dimension (:,:), allocatable :: pst_sed
!> average daily water removal from the pond for the month (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wupnd
!> phi(1,:) cross-sectional area of flow at bankfull depth (m^2)
!> phi(2,:) (none)
!> phi(3,:) (none)
!> phi(4,:) (none)
!> phi(5,:) (none)
!> phi(6,:) bottom width of main channel (m)
!> phi(7,:) depth of water when reach is at bankfull depth (m)
!> phi(8,:) average velocity when reach is at bankfull depth (m/s)
!> phi(9,:) wave celerity when reach is at bankfull depth (m/s)
!> phi(10,:) storage time constant for reach at bankfull depth (ratio of storage
!> to discharge) (hour)
!> phi(11,:) average velocity when reach is at 0.1 bankfull depth (low flow)
!> (m/s)
!> phi(12,:) wave celerity when reach is at 0.1 bankfull depth (low flow) (m/s)
!> phi(13,:) storage time constant for reach at 0.1 bankfull depth (low flow)
!> (ratio of storage to discharge) (hour)
   real*8, dimension (:,:), allocatable :: phi
!> precipitation for the day in band in HRU (mm H2O)
   real*8, dimension (:,:), allocatable :: pcpband
!> average temperature for the day in band in HRU (deg C)
   real*8, dimension (:,:), allocatable :: tavband
   real*8, dimension (:,:), allocatable :: wat_phi
!> initial snow water content in elevation band (mm H2O)
   real*8, dimension (:,:), allocatable :: snoeb
!> average daily water removal from the deep aquifer for the month
!> (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wudeep
!> average daily water removal from the shallow aquifer for the month
!> (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wushal
!> minimum temperature for the day in band in HRU (deg C)
   real*8, dimension (:,:), allocatable :: tmnband
   real*8, dimension (:), allocatable :: bss1
   real*8, dimension (:), allocatable :: bss2
   real*8, dimension (:), allocatable :: bss3
   real*8, dimension (:), allocatable :: bss4
!> nitrogen settling rate for 1st season (m/day)
   real*8, dimension (:), allocatable :: nsetlw1
!> nitrogen settling rate for 2nd season (m/day)
   real*8, dimension (:), allocatable :: nsetlw2
   real*8, dimension (:,:), allocatable :: snotmpeb,surf_bs
!> nitrogen settling rate for 1st season (m/day)
   real*8, dimension (:), allocatable :: nsetlp1
!> nitrogen settling rate for 2nd season (m/day)
   real*8, dimension (:), allocatable :: nsetlp2
!> maximum temperature for the day in band in HRU (deg C)
   real*8, dimension (:,:), allocatable :: tmxband
!> fraction of solar radiation occuring during hour in day in HRU (none)
   real*8, dimension (:,:), allocatable :: frad
!> precipitation for the time step during the day in HRU (mm H2O)
   real*8, dimension (:,:), allocatable :: rainsub
   real*8, dimension (:),   allocatable :: rstpbsb
   real*8, dimension (:,:), allocatable :: orig_snoeb,orig_pltpst
   real*8, dimension (:,:), allocatable :: terr_p, terr_cn, terr_sl
   real*8, dimension (:,:), allocatable :: drain_d, drain_t, drain_g
   real*8, dimension (:,:), allocatable :: drain_idep
   real*8, dimension (:,:), allocatable :: cont_cn, cont_p, filt_w
   real*8, dimension (:,:), allocatable :: strip_n, strip_cn, strip_c
   real*8, dimension (:,:), allocatable :: strip_p, fire_cn
   real*8, dimension (:,:), allocatable :: cropno_upd,hi_upd,laimx_upd
!> fraction of plant heat units at which grazing begins (none)
   real*8, dimension (:,:,:), allocatable :: phug
   real*8, dimension (:,:,:), allocatable :: pst_lag
   !!     integer, dimension (:), allocatable :: ipot
!> pesticide use flag (none)\n
!>  0: no pesticides used in HRU\n
!>  1: pesticides used in HRU
   integer, dimension (:), allocatable :: hrupest
!> sequence number of impound/release operation within the year (none)
   integer, dimension (:), allocatable :: nrelease
   integer, dimension (:), allocatable :: swtrg
!> number of years of rotation (none)
   integer, dimension (:), allocatable :: nrot
!> sequence number of fertilizer application within the year (none)
   integer, dimension (:), allocatable :: nfert
!> sequence number of year in rotation (none)
   integer, dimension (:), allocatable :: nro
!> land cover status code (none). This code informs the model whether or not a
!> land cover is growing at the beginning of the simulation\n
!> 0 no land cover growing\n
!> 1 land cover growing
   integer, dimension (:), allocatable :: igro
!> beginning month of nutrient settling season (none)
   integer, dimension (:), allocatable :: ipnd1
!> ending month of nutrient settling season (none)
   integer, dimension (:), allocatable :: ipnd2
!> sequence number of auto-irrigation application within the year (none)
   integer, dimension (:), allocatable :: nair
!> beginning month of non-flood season (none)
   integer, dimension (:), allocatable :: iflod1
!> ending month of non-flood season (none)
   integer, dimension (:), allocatable :: iflod2
!> number of days required to reach target storage from current pond storage
!> (none)
   integer, dimension (:), allocatable :: ndtarg
!> sequence number of irrigation application within the year (none)
   integer, dimension (:), allocatable :: nirr
   integer, dimension (:), allocatable :: iafrttyp, nstress
   integer, dimension (:), allocatable :: igrotree
   !! burn
   integer, dimension (:), allocatable :: grz_days
!> management code (for GIS output only) (none)
   integer, dimension (:), allocatable :: nmgt
!> sequence number of auto-fert application within the year (none)
   integer, dimension (:), allocatable :: nafert
!> sequence number of street sweeping operation within the year (none)
   integer, dimension (:), allocatable :: nsweep
   integer, dimension (:), allocatable :: icr,ncut
!> irrigation source location (none)\n
!> if IRRSC=1, IRRNO is the number of the reach\n
!> if IRRSC=2, IRRNO is the number of the reservoir\n
!> if IRRSC=3, IRRNO is the number of the subbasin\n
!> if IRRSC=4, IRRNO is the number of the subbasin\n
!> if IRRSC=5, not used
   integer, dimension (:), allocatable :: irrno
!> number of soil in soil profile layers (none)
   integer, dimension (:), allocatable :: sol_nly
!> prior day category (none)\n
!> 1 dry day\n
!> 2 wet day
   integer, dimension (:), allocatable :: npcp
   integer, dimension (:), allocatable :: irn
!> sequence number of continuous fertilization operation within the year (none)
   integer, dimension (:), allocatable :: ncf
!> sequence number of grazing operation within the year (none)
   integer, dimension (:), allocatable :: ngr
   integer, dimension (:), allocatable :: igrz,ndeat
!> subbasin in which HRU is located (none)
   integer, dimension (:), allocatable :: hru_sub
!> urban land type identification number from urban.dat (none)
   integer, dimension (:), allocatable :: urblu
!> soil layer where drainage tile is located (none)
   integer, dimension (:), allocatable :: ldrain
   integer, dimension (:), allocatable :: idorm
   integer, dimension (:), allocatable :: hru_seq
!> urban simulation code (none):\n
!> 0  no urban sections in HRU\n
!> 1  urban sections in HRU, simulate using USGS regression equations\n
!> 2  urban sections in HRU, simulate using build up/wash off algorithm
   integer, dimension (:), allocatable :: iurban
   integer, dimension (:), allocatable :: iday_fert,icfrt
!> number of HRU (in subbasin) that is a floodplain (none)
   integer, dimension (:), allocatable :: ifld
!> number of HRU (in subbasin) that is a riparian zone (none)
   integer, dimension (:), allocatable :: irip
   integer, dimension (:), allocatable :: ndcfrt,hrugis
!> irrigation source code (none):\n
!> 1 divert water from reach\n
!> 2 divert water from reservoir\n
!> 3 divert water from shallow aquifer\n
!> 4 divert water from deep aquifer\n
!> 5 divert water from source outside watershed
   integer, dimension (:), allocatable :: irrsc
   integer, dimension (:), allocatable :: orig_igro,ntil
   integer, dimension (:), allocatable :: iwatable,curyr_mat
   integer, dimension (:), allocatable :: ncpest,icpst,ndcpst
   integer, dimension (:), allocatable :: iday_pest, irr_flag
   integer, dimension (:), allocatable :: irra_flag
!> random number generator seeds array. The seeds in the array are used to
!> generate random numbers for the following purposes (none):\n
!> (1) wet/dry day probability\n
!> (2) solar radiation\n
!> (3) precipitation\n
!> (4) USLE rainfall erosion index\n
!> (5) wind speed\n
!> (6) 0.5 hr rainfall fraction\n
!> (7) relative humidity\n
!> (8) maximum temperature\n
!> (9) minimum temperature\n
!> (10) generate new random numbers
   integer, dimension (:,:), allocatable :: rndseed
   integer, dimension (:,:), allocatable :: iterr, iyterr
   integer, dimension (:,:), allocatable :: itdrain, iydrain, ncrops
!> manure (fertilizer) identification number from fert.dat (none)
   integer, dimension (:), allocatable :: manure_id

!!     gsm added for sdr (drainage) 7/24/08
   integer, dimension (:,:), allocatable :: mgt_sdr,idplrot
   integer, dimension (:,:), allocatable :: icont, iycont
   integer, dimension (:,:), allocatable :: ifilt, iyfilt
   integer, dimension (:,:), allocatable :: istrip, iystrip
   integer, dimension (:,:), allocatable :: iopday, iopyr, mgt_ops
   real*8, dimension (:), allocatable :: wshd_pstap, wshd_pstdg
   integer, dimension (12) :: ndmo
!> array of unique pesticides used in watershed (none)
   integer, dimension (:), allocatable :: npno
   integer, dimension (:), allocatable :: mcrhru
   character(len=13), dimension (18) :: rfile !< rainfall file names (.pcp)
   character(len=13), dimension (18) :: tfile !< temperature file names (.tmp)

   character(len=4), dimension (1000) :: urbname !< name of urban land use

   character(len=1), dimension (:), allocatable :: kirr !< irrigation in HRU
   character(len=1), dimension (:), allocatable :: hydgrp
   character(len=16), dimension (:), allocatable :: snam !< soil series name
   character(len=17), dimension (300) :: pname !< name of pesticide/toxin
   character(len=4) :: title(60) !< description lines in file.cio (1st 3 lines)
   character(len=4) :: cpnm(5000) !< four character code to represent crop name
   character(len=17), dimension(50) :: fname
!> average daily water loading for month (m^3/day)
   real*8, dimension (:,:,:), allocatable :: flomon
!> average daily soluble pesticide loading for month (mg pst/day)
   real*8, dimension (:,:,:), allocatable :: solpstmon
!> average daily sorbed pesticide loading for month (mg pst/day)
   real*8, dimension (:,:,:), allocatable :: srbpstmon
!> average daily organic N loading for month (kg N/day)
   real*8, dimension (:,:,:), allocatable :: orgnmon
!> average daily organic P loading for month (kg P/day)
   real*8, dimension (:,:,:), allocatable :: orgpmon
!> average daily sediment loading for month (metric tons/day)
   real*8, dimension (:,:,:), allocatable :: sedmon
!> average daily mineral P loading for month (kg P/day)
   real*8, dimension (:,:,:), allocatable :: minpmon
!> average amount of NH3-N loaded to stream on a given day in the month
!> (kg N/day)
   real*8, dimension (:,:,:), allocatable :: nh3mon
!> average daily NO3-N loading for month (kg N/day)
   real*8, dimension (:,:,:), allocatable :: no3mon
!> average amount of less persistent bacteria loaded to stream on a given day in
!> the month (# bact/day)
   real*8, dimension (:,:,:), allocatable :: bactlpmon
!> average amount of persistent bacteria loaded to stream on a given day in the
!> month (# bact/day)
   real*8, dimension (:,:,:), allocatable :: bactpmon
!> average amount of NO2-N loaded to stream on a given day in the month
!> (kg N/day)
   real*8, dimension (:,:,:), allocatable :: no2mon
!> average amount of conservative metal #1 loaded to stream on a given day in
!> the month (# bact/day)
   real*8, dimension (:,:,:), allocatable :: cmtl1mon
!> average amount of conservative metal #2 loaded to stream on a given day in
!> the month (# bact/day)
   real*8, dimension (:,:,:), allocatable :: cmtl2mon
!> average amount of conservative metal #3 loaded to stream on a given day in
!> the month (# bact/day)
   real*8, dimension (:,:,:), allocatable :: cmtl3mon
!> average daily loading of CBOD in month (kg/day)
   real*8, dimension (:,:,:), allocatable :: cbodmon
!> average daily loading of chlorophyll-a in month (kg/day)
   real*8, dimension (:,:,:), allocatable :: chlamon
!> average daily loading of dissolved O2 in month (kg/day)
   real*8, dimension (:,:,:), allocatable :: disoxmon
!> average daily water loading for year (m^3/day)
   real*8, dimension (:,:), allocatable :: floyr
!> average daily organic N loading for year (kg N/day)
   real*8, dimension (:,:), allocatable :: orgnyr
!> average daily organic P loading for year (kg P/day)
   real*8, dimension (:,:), allocatable :: orgpyr
!> average daily sediment loading for year (metric tons/day)
   real*8, dimension (:,:), allocatable :: sedyr
!> average daily mineral P loading for year (kg P/day)
   real*8, dimension (:,:), allocatable :: minpyr
!> average daily NH3-N loading for year (kg N/day)
   real*8, dimension (:,:), allocatable :: nh3yr
!> average daily NO2-N loading for year (kg N/day)
   real*8, dimension (:,:), allocatable :: no2yr
!> average daily NO3-N loading for year (kg N/day)
   real*8, dimension (:,:), allocatable :: no3yr
!> average daily loading of less persistent bacteria for year (# bact/day)
   real*8, dimension (:,:), allocatable :: bactlpyr
!> average daily loading of persistent bacteria for year (# bact/day)
   real*8, dimension (:,:), allocatable :: bactpyr
!> average daily loading of conservative metal #1 for year (kg/day)
   real*8, dimension (:,:), allocatable :: cmtl1yr
!> average daily loading of chlorophyll-a in year (kg/day)
   real*8, dimension (:,:), allocatable :: chlayr
!> average daily loading of conservative metal #2 for year (kg/day)
   real*8, dimension (:,:), allocatable :: cmtl2yr
!> average daily loading of conservative metal #3 for year (kg/day)
   real*8, dimension (:,:), allocatable :: cmtl3yr
!> average daily loading of CBOD in year (kg/day)
   real*8, dimension (:,:), allocatable :: cbodyr
!> average daily loading of dissolved O2 in year (kg/day)
   real*8, dimension (:,:), allocatable :: disoxyr
!> average daily soluble pesticide loading for year (mg pst/day)
   real*8, dimension (:,:), allocatable :: solpstyr
!> average daily sorbed pesticide loading for year (mg pst/day)
   real*8, dimension (:,:), allocatable :: srbpstyr
   real*8, dimension (:,:), allocatable :: sol_mc,sol_mn,sol_mp
!average daily water loading to reach (m^3 H2O/day)
   real*8, dimension (:), allocatable :: flocnst
!> average daily organic N loading to reach (kg N/day)
   real*8, dimension (:), allocatable :: orgncnst
!> average daily sediment loading for reach (metric tons/day)
   real*8, dimension (:), allocatable :: sedcnst
!> average daily soluble P loading to reach (kg P/day)
   real*8, dimension (:), allocatable :: minpcnst
!> average daily nitrate loading to reach (kg N/day)
   real*8, dimension (:), allocatable :: no3cnst
!> average daily organic P loading to reach (kg P/day)
   real*8, dimension (:), allocatable :: orgpcnst
!> average daily persistent bacteria loading to reach (# bact/day)
   real*8, dimension (:), allocatable :: bactpcnst
!> average daily ammonia loading to reach (kg N/day)
   real*8, dimension (:), allocatable :: nh3cnst
!> average daily nitrite loading to reach (kg N/day)
   real*8, dimension (:), allocatable :: no2cnst
!> average daily less persistent bacteria loading to reach (# bact/day)
   real*8, dimension (:), allocatable :: bactlpcnst
!> average daily conservative metal #1 loading (kg/day)
   real*8, dimension (:), allocatable :: cmtl1cnst
!> average daily conservative metal #2 loading (kg/day)
   real*8, dimension (:), allocatable :: cmtl2cnst
!> average daily loading of chlorophyll-a (kg/day)
   real*8, dimension (:), allocatable :: chlacnst
!> average daily conservative metal #3 loading (kg/day)
   real*8, dimension (:), allocatable :: cmtl3cnst
!> average daily loading of dissolved O2 (kg/day)
   real*8, dimension (:), allocatable :: disoxcnst
!> average daily loading of CBOD to reach (kg/day)
   real*8, dimension (:), allocatable :: cbodcnst
!> average daily soluble pesticide loading (mg/day)
   real*8, dimension (:), allocatable :: solpstcnst
!> average daily sorbed pesticide loading (mg/day)
   real*8, dimension (:), allocatable :: srbpstcnst

!> max number of time steps per day or number of lines of rainfall data for each
!> day (none)
   integer :: nstep
!> length of time step used to report precipitation data for sub-daily modeling
!> (minutes)
   integer :: idt
   real*8, dimension (:), allocatable :: hrtwtr,hhstor,hdepth,hsdti
   real*8, dimension (:), allocatable :: hrchwtr,halgae,horgn,hnh4
   real*8, dimension (:), allocatable :: hno2,hno3,horgp,hsolp,hbod
   real*8, dimension (:), allocatable :: hdisox,hchla,hsedyld,hsedst
   real*8, dimension (:), allocatable :: hharea,hsolpst,hsorpst
   real*8, dimension (:), allocatable :: hhqday,precipdt
   real*8, dimension (:), allocatable :: hhtime,hbactp,hbactlp
! store initial values
   integer, dimension (10) :: ivar_orig
   real*8, dimension (10) :: rvar_orig
   integer :: nsave !< number of save commands in .fig file
   integer ::  nauto, iatmodep
   real*8, dimension (:), allocatable :: wattemp
   real*8, dimension (:), allocatable :: lkpst_mass, lkspst_mass
   real*8, dimension (:), allocatable :: vel_chan
!> fraction of the total runoff from the entire field entering the most
!> concentrated 10% of the VFS (none)
   real*8, dimension (:), allocatable :: vfscon
!> field area/VFS area ratio (none)
   real*8, dimension (:), allocatable :: vfsratio
!> fraction of flow entering the most concentrated 10% of the VFS which is fully
!> channelized (none)
   real*8, dimension (:), allocatable :: vfsch
   real*8, dimension (:), allocatable :: vfsi
   real*8, dimension (:,:), allocatable :: filter_i,filter_ratio
   real*8, dimension (:,:), allocatable :: filter_con,filter_ch
   real*8, dimension (:,:), allocatable :: sol_n
!> = 0 Static soil carbon (old mineralization routines)\n
!> = 1 C-FARM one carbon pool model\n
!> = 2 Century model
   integer :: cswat
!! sj, june 07 end

!! sj, dec 07 dynamic bulk density
   real*8, dimension (:,:), allocatable :: sol_bdp
!! sj dec 07 end

!! Armen Jan 08
   real*8, dimension (:,:), allocatable :: tillagef
   real*8, dimension (:), allocatable :: rtfr
   real*8, dimension (:), allocatable :: stsol_rd
!! Armen Jan 08 end
   integer:: urban_flag, dorm_flag
   real*8 :: bf_flg, iabstr
   real*8, dimension (:), allocatable :: ubnrunoff,ubntss
   real*8, dimension (:,:), allocatable :: sub_ubnrunoff,sub_ubntss,&
      &ovrlnd_dt
   real*8, dimension (:,:,:), allocatable :: hhsurf_bs

!! subdaily erosion modeling by Jaehak Jeong
!> unit hydrograph method: 1=triangular UH; 2=gamma funtion UH;
   integer:: iuh
!> channel routing for HOURLY; 0=Bagnold; 2=Brownlie; 3=Yang;
   integer:: sed_ch
!> an exponent in the overland flow erosion equation ranges 1.5-3.0
   real*8 :: eros_expo
   real*8 :: eros_spl !< coefficient of splash erosion varing 0.9-3.1
!> Multiplier to USLE_K for soil susceptible to rill erosion, range 0.5-2.0
   real*8 :: rill_mult
   real*8 :: sedprev, c_factor
   real*8 :: ch_d50 !< median particle diameter of channel bed (mm)
!> geometric standard deviation of particle sizes for the main channel. Mean air
!> temperature at which precipitation is equally likely to be rain as
!> snow/freezing rain.
   real*8 :: sig_g
!> alpha coefficient for estimating unit hydrograph using a gamma function
!> (*.bsn)
   real*8 :: uhalpha
   real*8 :: abstinit,abstmax
   real*8, dimension(:,:), allocatable :: hhsedy, sub_subp_dt
   real*8, dimension(:,:), allocatable :: sub_hhsedy,sub_atmp
   real*8, dimension(:), allocatable :: rhy,init_abstrc
   real*8, dimension(:), allocatable :: dratio, hrtevp, hrttlc
   real*8, dimension(:,:,:), allocatable :: rchhr
!! subdaily reservoir modeling by Jaehak Jeong
   real*8, dimension(:), allocatable :: hhresflwi, hhresflwo,&
      &hhressedi, hhressedo

!! bmp modeling by jaehak jeong
   character(len=4), dimension(:), allocatable :: lu_nodrain
   integer, dimension(:), allocatable :: bmpdrain
   real*8, dimension(:), allocatable :: sub_cn2, sub_ha_urb,&
      &bmp_recharge
   !sed-fil
   real*8, dimension(:), allocatable :: sub_ha_imp,subdr_km,subdr_ickm
   real*8, dimension(:,:), allocatable :: sf_im,sf_iy,sp_sa,&
      &sp_pvol,sp_pd,sp_sedi,sp_sede,ft_sa,ft_fsa,&
      &ft_dep,ft_h,ft_pd,ft_k,ft_dp,ft_dc,ft_por,&
      &tss_den,ft_alp,sf_fr,sp_qi,sp_k,ft_qpnd,sp_dp,&
      &ft_qsw,ft_qin,ft_qout,ft_sedpnd,sp_bpw,ft_bpw,&
      &ft_sed_cumul,sp_sed_cumul
   integer, dimension(:), allocatable :: num_sf
   integer, dimension(:,:), allocatable :: sf_typ,sf_dim,ft_qfg,&
      &sp_qfg,sf_ptp,ft_fc
   real*8 :: sfsedmean,sfsedstdev  !Jaehak effluent probability method for urban bmp 2017

   !detention pond
!> month the reservoir becomes operational (none)
   integer, dimension(:), allocatable :: dtp_imo
!> year of the simulation that the reservoir becomes operational (none)
   integer, dimension(:), allocatable :: dtp_iyr
!> total number of stages in the weir (none)
   integer, dimension(:), allocatable :: dtp_numstage
!> total number of weirs in the BMP (none)
   integer, dimension(:), allocatable :: dtp_numweir
!> sub-basin detention pond is associated with (none)
   integer, dimension(:), allocatable :: dtp_onoff
!> equations for stage-discharge relationship (none):\n
!> 1=exponential function,\n
!> 2=linear,\n
!> 3=logarithmic,\n
!> 4=cubic,\n
!> 5=power
   integer, dimension(:), allocatable :: dtp_reltype
!> (none):\n
!> 0=use weir/orifice discharge equation to calculate outflow,\n
!> 1=use stage-dicharge relationship
   integer, dimension(:), allocatable :: dtp_stagdis
   integer, dimension(:), allocatable :: dtp_subnum
!> this parameter controls the response of decomposition to the combined effect
!> of soil temperature and moisture.
   real*8, dimension (:), allocatable :: cf
   real*8, dimension (:), allocatable :: cfh !< maximum humification rate
!> the undisturbed soil turnover rate under optimum soil water and temperature.
!> Increasing it will increase carbon and organic N decomp.
   real*8, dimension (:), allocatable :: cfdec
! additional nutrient variables by jeong for montana bitterroot
   real*8, dimension(:), allocatable :: lat_orgn, lat_orgp

!> weir dimensions (none),\n
!> 1=read user input,\n
!> 0=use model calculation
   integer, dimension(:,:), allocatable :: dtp_weirdim
!> type of weir (none):\n
!> 1=rectangular and\n
!> 2=circular
   integer, dimension(:,:), allocatable :: dtp_weirtype

!> coefficient of 3rd degree in the polynomial equation (none)
   real*8, dimension(:), allocatable :: dtp_coef1
!> coefficient of 2nd degree in the polynomial equation (none)
   real*8, dimension(:), allocatable :: dtp_coef2
!> coefficient of 1st degree in the polynomial equation (none)
   real*8, dimension(:), allocatable :: dtp_coef3
!> detention pond evaporation coefficient (none)
   real*8, dimension(:), allocatable :: dtp_evrsv
!> exponent used in the exponential equation (none)
   real*8, dimension(:), allocatable :: dtp_expont
!> intercept used in regression equations (none)
   real*8, dimension(:), allocatable :: dtp_intcept
!> ratio of length to width of water back up (none)
   real*8, dimension(:), allocatable :: dtp_lwratio
!> total constructed width of the detention wall across the creek (m)
   real*8, dimension(:), allocatable :: dtp_totwrwid
   real*8, dimension(:), allocatable :: &
      &dtp_inflvol,dtp_wdep,dtp_totdep,&
      &dtp_watdepact,dtp_outflow,dtp_totrel,dtp_backoff,dtp_seep_sa,&
      &dtp_evap_sa,dtp_pet_day,dtp_pcpvol,dtp_seepvol,dtp_evapvol,&
      &dtp_flowin,dtp_backup_length,dtp_ivol,dtp_ised

   integer, dimension (:,:),allocatable :: so_res_flag, ro_bmp_flag
   real*8, dimension (:,:),allocatable :: sol_watp, sol_solp_pre
   real*8, dimension (:,:),allocatable :: psp_store, ssp_store, so_res
   real*8, dimension (:,:),allocatable :: sol_cal, sol_ph
   integer:: sol_p_model
   integer, dimension (:,:),allocatable :: a_days, b_days
   real*8, dimension (:), allocatable :: harv_min, fstap, min_res
   real*8, dimension (:,:),allocatable :: ro_bmp_flo, ro_bmp_sed
   real*8, dimension (:,:),allocatable :: ro_bmp_bac
   real*8, dimension (:,:),allocatable :: ro_bmp_pp, ro_bmp_sp
   real*8, dimension (:,:),allocatable :: ro_bmp_pn, ro_bmp_sn

   real*8, dimension (:,:),allocatable :: ro_bmp_flos, ro_bmp_seds
   real*8, dimension (:,:),allocatable :: ro_bmp_bacs
   real*8, dimension (:,:),allocatable :: ro_bmp_pps, ro_bmp_sps
   real*8, dimension (:,:),allocatable :: ro_bmp_pns, ro_bmp_sns

   real*8, dimension (:,:),allocatable :: ro_bmp_flot, ro_bmp_sedt
   real*8, dimension (:,:),allocatable :: ro_bmp_bact
   real*8, dimension (:,:),allocatable :: ro_bmp_ppt, ro_bmp_spt
   real*8, dimension (:,:),allocatable :: ro_bmp_pnt, ro_bmp_snt

   real*8, dimension (:),allocatable :: bmp_flo, bmp_sed, bmp_bac
   real*8, dimension (:),allocatable :: bmp_pp, bmp_sp
   real*8, dimension (:),allocatable :: bmp_pn, bmp_sn, bmp_flag

   real*8, dimension (:),allocatable :: bmp_flos, bmp_seds, bmp_bacs
   real*8, dimension (:),allocatable :: bmp_pps, bmp_sps
   real*8, dimension (:),allocatable :: bmp_pns, bmp_sns

   real*8, dimension (:),allocatable :: bmp_flot, bmp_sedt, bmp_bact
   real*8, dimension (:),allocatable :: bmp_ppt, bmp_spt
   real*8, dimension (:),allocatable :: bmp_pnt, bmp_snt

!> the distance between spillway levels (m)
   real*8, dimension(:,:), allocatable :: dtp_addon
!> discharge coeffieicne for weir/orifice flow (none)
   real*8, dimension(:,:), allocatable :: dtp_cdis
!> depth of rectangular wier at different stages (m)
   real*8, dimension(:,:), allocatable :: dtp_depweir
!> diameter of orifice hole at different stages (m)
   real*8, dimension(:,:), allocatable :: dtp_diaweir
!> maximum discharge from each stage of the weir/hole (m^3/s)
   real*8, dimension(:,:), allocatable :: dtp_flowrate
!> precipitation for different return periods (not used) (mm)
   real*8, dimension(:,:), allocatable :: dtp_pcpret
!> return period at different stages (years)
   real*8, dimension(:,:), allocatable :: dtp_retperd
!> width depth ratio of rectangular weirs (none)
   real*8, dimension(:,:), allocatable :: dtp_wdratio
   real*8, dimension(:,:), allocatable :: dtp_wrwid

   !retention irrigation
   real*8, dimension(:), allocatable :: ri_subkm,ri_totpvol,&
      &irmmdt
   real*8, dimension(:,:), allocatable :: ri_sed,ri_fr,ri_dim,&
      &ri_im,ri_iy,ri_sa,ri_vol,ri_qi,ri_k,ri_dd,ri_evrsv,&
      &ri_dep,ri_ndt,ri_pmpvol,ri_sed_cumul,hrnopcp,ri_qloss,&
      &ri_pumpv,ri_sedi
   character(len=4), dimension(:,:), allocatable :: ri_nirr
   integer, dimension(:), allocatable :: num_ri,ri_luflg,num_noirr

   !wet pond
   integer, dimension(:), allocatable :: wtp_subnum,wtp_onoff,wtp_imo,&
      &wtp_iyr,wtp_dim,wtp_stagdis,wtp_sdtype
   real*8, dimension(:), allocatable :: wtp_pvol,wtp_pdepth,wtp_sdslope,&
      &wtp_lenwdth,wtp_extdepth,wtp_hydeff,wtp_evrsv,wtp_sdintc,&
      &wtp_sdexp,wtp_sdc1,wtp_sdc2,wtp_sdc3,wtp_pdia,wtp_plen,&
      &wtp_pmann,wtp_ploss,wtp_k,wtp_dp,wtp_sedi,wtp_sede,wtp_qi

   real*8 :: lai_init !< initial leaf area index of transplants
   real*8 :: bio_init !< initial biomass of transplants (kg/ha)
!> SCS runoff curve number for moisture condition II (none)
   real*8 :: cnop
!> harvest efficiency: fraction of harvested yield that is removed from HRU; the
!> remainder becomes residue on the soil surface(none)
   real*8 :: harveff
!> harvest index target specified at harvest ((kg/ha)/(kg/ha))
   real*8 :: hi_ovr
   real*8 :: frac_harvk

   ! van Genuchten equation's coefficients
   real*8 :: lid_vgcl,lid_vgcm,lid_qsurf_total,&
      &lid_farea_sum

   ! soil water content and amount of accumulated infiltration
   real*8, dimension(:,:), allocatable :: lid_cuminf_last,&
      &lid_sw_last, interval_last,lid_f_last,lid_cumr_last,lid_str_last,&
      &lid_farea,lid_qsurf,lid_sw_add,lid_cumqperc_last,lid_cumirr_last,&
      &lid_excum_last                                                      !! nbs

   ! Green Roof
   integer, dimension(:,:), allocatable :: gr_onoff,gr_imo,gr_iyr
   real*8, dimension(:,:), allocatable :: gr_farea,gr_solop,gr_etcoef,&
      &gr_fc,gr_wp,gr_ksat,gr_por,gr_hydeff,gr_soldpt

   ! Rain Gerden
   integer, dimension(:,:), allocatable :: rg_onoff,rg_imo,rg_iyr
   real*8, dimension(:,:), allocatable :: rg_farea,rg_solop,rg_etcoef,&
      &rg_fc,rg_wp,rg_ksat,rg_por,rg_hydeff,rg_soldpt,rg_dimop,rg_sarea,&
      &rg_vol,rg_sth,rg_sdia,rg_bdia,rg_sts,rg_orifice,rg_oheight,&
      &rg_odia

   ! CiStern
   integer, dimension(:,:), allocatable :: cs_onoff,cs_imo,cs_iyr,&
      &cs_grcon
   real*8, dimension(:,:), allocatable :: cs_farea,cs_vol,cs_rdepth

   ! Porous paVement
   integer, dimension(:,:), allocatable :: pv_onoff,pv_imo,pv_iyr,&
      &pv_solop
   real*8, dimension(:,:), allocatable :: pv_grvdep,pv_grvpor,pv_farea,&
      &pv_drcoef,pv_fc,pv_wp,pv_ksat,pv_por,pv_hydeff,pv_soldpt

   ! LID general
   integer, dimension(:,:), allocatable :: lid_onoff


!! By Zhang for C/N cycling
   !!SOM-residue C/N state variables -- currently included
   real*8, dimension(:,:), allocatable :: sol_BMC, sol_BMN, sol_HSC,&
      &sol_HSN, sol_HPC, sol_HPN, sol_LM,&
      &sol_LMC, sol_LMN, sol_LS, sol_LSL, sol_LSC, sol_LSN , sol_RNMN,&
      &sol_LSLC, sol_LSLNC, sol_RSPC, sol_WOC, sol_WON, sol_HP, sol_HS,&
      &sol_BM
   ! HSC mass of C present in slow humus (kg ha-1)
   ! HSN mass of N present in slow humus (kg ha-1)
   ! HPC mass of C present in passive humus (kg ha-1)
   ! HPN mass of N present in passive humus (kg ha-1)
   ! LM mass of metabolic litter (kg ha-1)
   ! LMC mass of C in metabolic litter (kg ha-1)
   ! LMN mass of N in metabolic litter (kg ha-1)
   ! LS mass of structural litter (kg ha-1)
   ! LSC mass of C in structural litter (kg ha-1)
   ! LSL mass of lignin in structural litter (kg ha-1)
   ! LSN mass of N in structural litter (kg ha-1)

   !!SOM-residue C/N state variables -- may need to be included
   real*8, dimension(:,:), allocatable :: sol_CAC, sol_CEC

   !!daily updated soil layer associated percolaton and lateral flow Carbon loss
   real*8, dimension(:,:), allocatable :: sol_percc
   real*8, dimension(:,:), allocatable :: sol_latc

   !!Daily carbon change by different means (entire soil profile for each HRU)
   real*8, dimension(:), allocatable :: sedc_d, surfqc_d, latc_d,&
      &percc_d, foc_d, NPPC_d, rsdc_d, grainc_d, stoverc_d, soc_d,&
      &rspc_d, emitc_d
   !!emitc_d include biomass_c eaten by grazing, burnt


   !!Daily carbon change by different means (entire soil profile for each Subbasin)
   !!Only defined the variables, but not used them in the code
   real*8, dimension(:), allocatable :: sub_sedc_d, sub_surfqc_d,&
      &sub_latc_d, sub_percc_d, sub_foc_d, sub_NPPC_d, sub_rsdc_d,&
      &sub_grainc_d, sub_stoverc_d, sub_emitc_d, sub_soc_d, sub_rspc_d


   !!Monthly carbon change by different means (entire soil profile for each HRU)
   real*8, dimension(:), allocatable :: sedc_m, surfqc_m, latc_m, percc_m,&
      &foc_m, NPPC_m, rsdc_m, grainc_m, stoverc_m, emitc_m, soc_m,&
      &rspc_m

   !!Yearly carbon change by different means (entire soil profile for each HRU)
   real*8, dimension(:), allocatable :: sedc_a, surfqc_a, latc_a,&
      &percc_a, foc_a, NPPC_a, rsdc_a, grainc_a, stoverc_a, emitc_a,&
      &soc_a, rspc_a


   !! The following variables are defined and calculated locally
   !! ==================================================================
   ! HSCTP potential transformation of C in slow humus (kg ha-1 day-1)
   ! HSNTP potential transformation of N in slow humus (kg ha.1 day-1)
   ! HPCTP potential transformation of C in passive humus (kg ha-1 day-1)
   ! HPNTP potential transformation of N in passive humus (kg ha-1 day-1)
   ! HPR rate of transformation of passive humus under optimal conditions (subsurface
   ! layers = 0.000012 day-1) (Parton et al.,1993, 1994)
   ! HSR rate of transformation of slow humus under optimal conditions (all layers
   ! = 0.0005 day.1) (Parton et al., 1993, 1994; Vitousek et al., 1993)
   ! KOC liquid C solid partition coefficient for microbial biomass (10^3 m3 Mg-1)
   ! LMF fraction of the litter that is metabolic
   ! LMNF fraction of metabolic litter that is N (kg kg-1)
   ! LMR rate of transformation of metabolic litter under optimal conditions (surface =
   !  0.0405 day-1; all other layers = 0.0507 day-1) (Parton et al., 1994)
   ! LMCTP potential transformation of C in metabolic litter (kg ha-1 day-1)
   ! LMNTP potential transformation of N in metabolic litter (kg ha-1 day-1)
   ! LSCTP potential transformation of C in structural litter (kg ha-1 day-1)
   ! LSF fraction of the litter that is structural
   ! LSLF fraction of structural litter that is lignin (kg kg-1)
   ! LSNF fraction of structural litter that is N (kg kg-1)
   ! LSLCTP potential transformation of C in lignin of structural litter (kg ha-1 day-1)
   ! LSLNCTP potential transformation of C in nonlignin structural litter (kg ha-1 day-1)
   ! LSNTP potential transformation of N in structural litter (kg ha-1 day-1)
   ! LSR rate of potential transformation of structural litter under optimal conditions
   ! (surface = 0.0107 day.1; all other layers = 0.0132 day.1) (Parton et al., 1994)
   ! NCBM N/C ratio of biomass
   ! NCHP N/C ratio passive humus
   ! NCHS N/C ratio of the slow humus
   ! OX oxygen control on biological processes with soil depth
   ! Sf fraction of mineral N sorbed to litter: 0.05 for surface litter, 0.1 for belowground litter

   !!Tillage factor on SOM decomposition
   integer, dimension(:), allocatable :: tillage_switch
   real*8, dimension(:), allocatable :: tillage_depth
   integer, dimension(:), allocatable :: tillage_days
   real*8, dimension(:), allocatable :: tillage_factor
   ! tillage_factor: = 1.6 in 30 days after tillage practices occur; otherwise 1.0;
!! By Zhang for C/N cycling

   !Flood routing variables by Jaehak Jeong 2017
   real*8 :: dthy !< time interval for subdaily routing
   integer, dimension(4) :: IHX
   integer, dimension(:), allocatable :: NHY
   real*8, dimension(:), allocatable :: RCHX,RCSS,QCAP,CHXA,CHXP
   real*8, dimension(:,:,:), allocatable :: QHY

   !!Variables for killop.f90 and harvkillop.f90 files
   real*8 :: ff1, ff2

   ! tdc 2018-03-29 prototypes for some functions that returned REAL implicitly
   ! but which must now return REAL*8 explicitly
   INTERFACE

      function atri(at1,at2,at3,at4i) result (r_atri)
         real*8, intent (in) :: at1, at2, at3
         integer, intent (in out) :: at4i
         real*8 ::r_atri
      end function

      real*8 function aunif (x1) result (unif)
         integer, intent (in out) :: x1
      end function

      real*8 function dstn1(rn1,rn2) result (r_dstn1)
         real*8, intent (in) :: rn1, rn2
      end function

      real*8 function ee(tk) result (r_ee)
         real*8, intent (in) :: tk
      end

      function expo (xx) result (r_expo)
         real*8 :: xx
         real*8 :: r_expo
      end function

      real*8 Function fcgd(xx)
         real*8, intent (in) :: xx
      End function

      real*8 function qman(x1,x2,x3,x4) result (r_qman)
         real*8, intent (in) :: x1, x2, x3, x4
      end function

      real*8 function regres(k) result (r_regres)
         integer, intent (in) :: k
      end

      function tair(hr,jj) result (r_tair)
         integer, intent (in) ::  jj
         real*8, intent(in) :: hr
         real*8 :: r_tair
      end function

      real*8 function theta(r20,thk,tmp) result (r_theta)
         real*8, intent (in) :: r20, thk, tmp
      end

      subroutine ascrv(x1,x2,x3,x4,x5,x6)
         real*8, intent (in) :: x1, x2, x3, x4
         real*8, intent (out) :: x5, x6
      end subroutine

      SUBROUTINE HQDAV(A,CBW,QQ,SSS,ZCH,ZX,CHW,FPW,jrch)
         real*8, intent (in out) :: A, ZX, CHW, FPW
         real*8, intent (in) :: CBW, QQ, SSS, ZCH
         integer, intent (in) :: jrch
      end subroutine

      subroutine layersplit(dep_new)
         real*8, intent(in):: dep_new
      end subroutine

      subroutine ndenit(k,j,cdg,wdn,void)
         integer :: k,j
         real*8 :: cdg, wdn, void
      end subroutine

      subroutine rsedaa(years)
         real*8, intent (in) :: years
      end subroutine

      subroutine vbl(evx,spx,pp,qin,ox,vx1,vy,yi,yo,ysx,vf,vyf,aha)
         real*8, intent (in) :: evx, spx, pp, qin, ox, yi, yo, ysx
         real*8, intent (in) :: vf, vyf, aha
         real*8, intent (in out) :: vx1, vy
      end subroutine

   END INTERFACE

end module parm
