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

!> column headers for subbasin output file
   character(len=13), dimension(msubo), parameter ::&
      &hedb = (/"  PRECIPmm"," SNOMELTmm","     PETmm","      ETmm",&
      &"      SWmm","    PERCmm","    SURQmm","    GW_Qmm",&
      &"    WYLDmm","  SYLDt/ha"," ORGNkg/ha"," ORGPkg/ha",&
      &"NSURQkg/ha"," SOLPkg/ha"," SEDPkg/ha"," LAT Q(mm)",&
      &"LATNO3kg/h","GWNO3kg/ha","CHOLAmic/L","CBODU mg/L",&
      &" DOXQ mg/L"," TNO3kg/ha","   QTILEmm"," TVAPkg/ha"/)

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

!> space number for beginning of column in HRU output file (none)
   integer, dimension (mhruo), parameter ::&
      &icols = (/43,53,63,73,83,93,103,113,123,133,143,153,&
      &163,173,183,193,203,213,223,233,243,253,263,273,283,&
      &293,303,313,323,333,343,353,363,373,383,393,403,413,&
      &423,433,443,453,463,473,483,493,503,513,523,533,543,&
      &553,563,573,583,593,603,613,623,633,643,653,663,673,&
      &683,693,703,713,723,733,743,753,763,773,783,793,803,&
      &813,823/)

!> space number for beginning of column in subbasin output file (none)
   integer, dimension (msubo), parameter ::&
      &icolb = (/35,45,55,65,75,85,95,105,115,125,135,145,&
      &155,165,175,185,195,205,215,225,235,245,255,265/)

!> space number for beginning of column in reach output file (none)
   integer, dimension (mrcho), parameter ::&
      &icolr = (/38,50,62,74,86,98,110,122,134,146,158,170,182,194,206,&
      &218,230,242,254,266,278,290,302,314,326,338,350,362,374,386,398,&
      &410,422,434,446,458,470,482,494,506,518,530,542,554,566,578,590,&
      &602,614,626,638,650,662,674,686,698,710,722,734,746,758,770/)

!> space number for beginning of column in reservoir output file (none)
   integer, dimension (41), parameter ::&
      &icolrsv = (/38,50,62,74,86,98,110,122,134,146,158,170,182,194,&
      &206,218,230,242,254,266,278,290,302,314,326,338,350,362,374,386,&
      &398,410,422,434,446,458,470,482,494,506,518/)

   real*8, parameter :: ab = 0.02083 !< lowest value al5 can have (mm H2O)

   integer, dimension (13), parameter ::&
      &ndays_leap = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
   integer, dimension (13), parameter ::&
      &ndays_noleap = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)

!> code for writing out calendar day or julian day to output.rch, .sub, .hru
!> files;\n
!> icalen = 0 (print julian day), 1 (print month/day/year);\n
!> icalen MUST be == zero if IPRINT == 3 to print subdaily
   integer icalen
!> Basinwide peak rate adjustment factor for sediment routing in the channel.
!> Allows impact of peak flow rate on sediment routing and channel reshaping to
!> be taken into account.
   real*8 :: prf_bsn

!!    srin - co2 (EPA)
   real*8 :: co2_x2, co2_x

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
!> rate factor for humus mineralization on active organic N
   real*8, dimension (:), allocatable :: cmn
!> phosphorus soil partitioning coefficient. Ratio of soluble phosphorus in
!> surface layer attached to sediment to phosphorus dissolved in soil water
   real*8, dimension (:), allocatable :: phoskd
!> phosphorus availibility index. The fraction of fertilizer P remaining in
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
!> yield (dry weight) (kg)
   real*8 :: yield
!> fraction of biomass and residue that burn(input in management file) range
!> (0 - 1.0) (none)
   real*8 :: burn_frlb
   real*8 :: yieldgrn, yieldbms, yieldtbr, yieldn, yieldp
   real*8 :: hi_bms, hi_rsd, yieldrsd

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

!! septic variables for output.std
   real*8 :: wshd_sepno3, wshd_sepnh3, wshd_seporgn, wshd_sepfon
   real*8 :: wshd_seporgp, wshd_sepfop, wshd_sepsolp, wshd_sepbod
   real*8 :: wshd_sepmm
   integer, dimension (:), allocatable :: isep_hru
!! septic variables for output.std
   real*8 :: fixco !< nitrogen fixation coefficient
   real*8 :: nfixmx !< maximum daily n-fixation (kg/ha) 
   real*8 :: res_stlr_co !< reservoir sediment settling coefficient
   real*8 :: rsd_covco !< residue cover factor for computing fraction of cover
   real*8 :: vcrit !< critical velocity
!> average amount of water stored in snow at the beginning of the simulation for
!> the entire watershed (mm H20)
   real*8 :: wshd_snob
!> water in soil at beginning of simulation, or\
!> average amount of water stored in soil for the entire watershed, or\
!> difference between mass balance calculated from watershed averages and actual
!> value for water in soil at end of simulation (goal is to have wshd_sw = 0.)
!> (mm H2O)
   real*8 :: wshd_sw
!> fraction of watershed area which drains into ponds (none)
   real*8 :: wshd_pndfr
!> total amount of suspended sediment in ponds in the watershed (metric tons),\n
!> or mass balance discrepancy for pond sediment expressed as loading per unit
!> hectare of drainage area (metric tons/ha)
   real*8 :: wshd_pndsed
!> total volume of water in ponds in the watershed (m^3),
!> or mass balance discrepancy for pond water volume expressed as depth over
!> drainage area (mm H2O)
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
! fraction of watershed area which drains into wetlands (none)
!   real*8 :: wshd_wetfr !variable not used
!> average annual amount of mineral P applied in watershed (kg P/ha)
   real*8 :: wshd_fminp
!> average annual amount of NH3-N applied in watershed (kg N/ha)
   real*8 :: wshd_fnh3
!> average annual amount of NO3-N applied in watershed (kg N/ha)
   real*8 :: wshd_fno3
!> average annual amount of organic N applied in watershed (kg N/ha)
   real*8 :: wshd_forgn
!> average annual amount of N (mineral & organic) applied in watershed (kg N/ha)
   real*8 :: wshd_ftotn
!> average annual amount of organic P applied in watershed (kg P/ha)
   real*8 :: wshd_forgp
!> average annual amount of P (mineral & organic) applied in watershed (kg P/ha)
   real*8 :: wshd_ftotp
!> amount of nitrogen removed from soil in watershed in the yield (kg N/ha)
   real*8 :: wshd_yldn
!> amount of phosphorus removed from soil in watershed in the yield (kg P/ha)
   real*8 :: wshd_yldp
!> average annual amount of nitrogen added to plant biomass via fixation
!> (kg N/ha)
   real*8 :: wshd_fixn
!> average annual amount of plant uptake of phosphorus (kg P/ha)
   real*8 :: wshd_pup
!> average annual number of nitrogen stress units in watershed (stress units)
   real*8 :: wshd_nstrs
!> average annual number of phosphorus stress units in watershed (stress units)
   real*8 :: wshd_pstrs
!> average annual number of temperature stress units in watershed (stress units)
   real*8 :: wshd_tstrs
!> average annual number of water stress units in watershed (stress units)
   real*8 :: wshd_wstrs
   real*8 :: wshd_astrs
!> initial soil water content expressed as a fraction of field capacity
   real*8 :: ffcb
!> average annual amount of nitrogen lost from nitrate pool due to
!> denitrification in watershed (kg N/ha)
   real*8 :: wshd_dnit
!> average annual amount of nitrogen moving from active organic to nitrate pool
!> in watershed (kg N/ha)
   real*8 :: wshd_hmn
!> average annual amount of phosphorus moving from organic to labile pool in
!> watershed (kg P/ha)
   real*8 :: wshd_hmp
!> average annual amount of nitrogen moving from fresh organic (residue) to
!> nitrate and active organic pools in watershed (kg N/ha)
   real*8 :: wshd_rmn
!> average annual amount of nitrogen moving from active organic to stable
!> organic pool in watershed (kg N/ha)
   real*8 :: wshd_rwn
!> die-off factor for persistent bacteria in soil solution (1/day)
   real*8 :: wdpq
!> average annual amount of phosphorus moving from fresh organic (residue) to
!> labile and organic pools in watershed (kg P/ha)
   real*8 :: wshd_rmp
!> average annual amount of nitrogen moving from the NH3 to the NO3 pool by
!> nitrification in the watershe (kg N/ha)d
   real*8 :: wshd_nitn
!> average annual amount if nitrogen lost by ammonia volatilization in watershed
!> (kg N/ha)
   real*8 :: wshd_voln
!> average annual amount of phosphorus moving from labile mineral to active
!> mineral pool in watershed (kg P/ha)
   real*8 :: wshd_pal
!> average annual amount of phosphorus moving from active mineral to stable
!> mineral pool in watershed (kg P/ha)
   real*8 :: wshd_pas
!> fraction of persistent bacteria on foliage that is washed off by a rainfall
!> event (none)
   real*8 :: wof_p
!> average annual amount of NO3 added to soil by rainfall in watershed (kg N/ha)
   real*8 :: wshd_raino3
!X
!> average annual amount of phosphorus leached into second soil layer (kg P/ha)
   real*8 :: wshd_plch
   real*8 :: ressedc, basno3f, basorgnf
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
   real*8 :: basminpf, basorgpf
!> total amount of suspended sediment in reservoirs in the watershed
!> (metric tons),\n
!> or mass balance discrepancy for reservoir sediment expressed as loading per
!> unit hectare of drainage area (metric tons/ha)
   real*8 :: wshd_ressed
!> total volume of water in all reservoirs in the watershed (m^3),\n
!> or mass balance discrepancy for reservoir water volume expressed as depth
!> over drainage area (mm H2O)
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
!> average amount of phosphorus initially in the organic P pool in watershed
!> soil (kg P/ha)
   real*8 :: basorgpi
   real*8 :: peakr !< peak runoff rate for the day in HRU or channel (m^3/s)
!> albedo of ground for the day in HRU, the fraction of the solar radiation
!> reflected at the soil surface back into space (none)
   real*8 :: albday
!> sediment inflow to the pond from HRU during day (metric tons)
   real*8 :: pndsedin
!> amount of water stored in soil layer on the current day that exceeds field
!> capacity (gravity drained water) (mm H2O)
   real*8 :: sw_excess
!> Snow pack temperature lag factor (0-1)\n
!> 1 = no lag (snow pack temp=current day air temp) as the lag factor goes to
!> zero, the snow pack's temperature will be less influenced by the current
!> day's air temperature
   real*8 :: timp
!> shallow water table depth above the impervious layer (mm H2O)
   real*8 :: wt_shall
   real*8 :: sq_rto
!> amount of water in drainage tile flow in HRU soil layer for the day (mm H2O)
   real*8 :: qtile
!> amount of precipitation that infiltrates into soil (enters soil) (mm H2O)
   real*8 :: inflpcp
!> amount of nitrogen added to the plant biomass via fixation on the day in HRU
!> (kg N/ha)
   real*8 :: fixn
!> amount of water in lateral flow in layer in HRU for the day (mm H2O)
   real*8 :: latlyr
!> amount of precipitation falling as freezing rain/snow on day in HRU (mm H2O)
   real*8 :: snofall
!> amount of water in snow melt for the day in HRU (mm H2O)
   real*8 :: snomlt
!> amount of water removed from surface runoff via transmission losses on day in
!> HRU (mm H2O)
   real*8 :: tloss
   real*8 :: lpndloss, lwetloss
!> biomass generated on current day in HRU (kg)
   real*8 :: bioday
!> total amount of nitrogen applied to soil during continuous fertilizer
!> operation in HRU on day (kg N/ha)
   real*8 :: cfertn
!> amount of phosphorus applied to soil during continuous fertilizer operation
!> in HRU on day (kg P/ha)
   real*8 :: cfertp
!> total amount of nitrogen applied to soil in HRU on day in fertilizer
!> application (kg N/ha)
   real*8 :: fertn
!> micropore percolation from bottom of the soil layer on day in HRU (mm H2O)
   real*8 :: sepday
   real*8 :: sol_rd !< current rooting depth (mm)
!> sediment transported out of channel during time step (metric tons)
   real*8 :: sedrch
   real*8 :: sepcrktot, fertno3, fertnh3, fertorgn, fertsolp
   real*8 :: fertorgp
   real*8 :: qdfr !< fraction of water yield that is surface runoff (none)
!> total amount of phosphorus applied to soil in HRU on day in fertilizer
!> application (kg P/ha)
   real*8 :: fertp
!> amount of nitrogen added to soil in grazing on the day in HRU (kg N/ha)
   real*8 :: grazn
!> amount of phosphorus added to soil in grazing on the day in HRU (kg P/ha)
   real*8 :: grazp
   real*8 :: soxy !< saturation dissolved oxygen concentration (mg/L)
!X
   real*8 :: rtwtr !< water leaving reach on day (m^3 H2O)
   real*8 :: sdti !< average flow rate in reach for day (m^3/s)
   real*8 :: ressa
!> die-off factor for less persistent bacteria absorbed to soil particles
!> (1/day)
   real*8 :: wdlps
!> growth factor for less persistent bacteria adsorbed to soil particles
!> (1/day)
   real*8 :: wglps
   real*8 :: da_km !< area of the watershed in square kilometers (km^2)
   real*8 :: rchdep !< depth of flow on day (m)
   real*8 :: rtevp !< evaporation from reach on day (m^3 H2O)
   real*8 :: rttime !< reach travel time (hour)
   real*8 :: rttlc !< transmission losses from reach on day (m^3 H2O)
   real*8 :: resflwi
   real*8 :: wdprch !< die-off factor for persistent bacteria in streams (1/day)
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
!> total volume of cracks expressed as depth per unit area (mm)
   real*8 :: voltot
   real*8 :: phoskd_bsn
!> weighting factor controling relative importance of inflow rate and outflow
!> rate in determining storage on reach
   real*8 :: msk_x
!> minimum crack volume allowed in any soil layer (mm), or\n
!> minimum soil volume in profile (mm)
   real*8 :: volcrmin
!> bacteria soil partitioning coefficient. Ratio of solution bacteria in surface
!> layer to solution bacteria in runoff soluble and sorbed phase in surface
!> runoff.
   real*8 :: bactkdq
!> die-off factor for persistent bacteria on foliage (1/day)
   real*8 :: wdpf
!> amount of water evaporated from canopy storage (mm H2O)
   real*8 :: canev
!> precipitation, or effective precipitation reaching soil surface, for the
!> current day in HRU (mm H2O)
   real*8 :: precipday
   real*8 :: uno3d !< plant nitrogen deficiency for day in HRU (kg N/ha)
!> daily soil loss predicted with USLE equation (metric tons/ha)
   real*8 :: usle
!> concentration of nitrogen in the rainfall (mg/L)
   real*8 :: rcn
   real*8 :: surlag_bsn
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
!> persistent bacteria transported to main channel with surface runoff
!> (# colonies/ha)
   real*8 :: bactrop
!> persistent bacteria transported with sediment in surface runoff
!> (# colonies/ha)
   real*8 :: bactsedp
   real*8 :: wgpf !< growth factor for persistent bacteria on foliage (1/day)
!> less persistent bacteria removed from soil surface layer by percolation
!> (# colonies/ha)
   real*8 :: bactlchlp
!> persistent bacteria removed from soil surface layer by percolation
!> (# colonies/ha)
   real*8 :: bactlchp
!> enrichment ratio calculated for current day in HRU (none)
   real*8 :: enratio
   real*8 :: pndpcp !< precipitation on pond during day (m^3 H2O)
   real*8 :: wetpcp !< precipitation on wetland for day (m^3 H2O)
   real*8 :: wetsep !< seepage from wetland bottom for day (m^3 H2O)
   real*8 :: pndev !< evaporation from pond on day (m^3 H2O)
   real*8 :: pndflwi !< volume of water flowing into pond on day (m^3 H2O)
   real*8 :: pndsedo !< sediment leaving pond during day (metric tons)
   real*8 :: pndsep !< seepage from pond on day (m^3 H2O)
   real*8 :: wetev !< evaporation from wetland for day (m^3 H2O)
   real*8 :: wetflwi !< volume of water flowing in wetland on day (m^3 H2O)
   real*8 :: wetsedo !< sediment loading from wetland for day (metric tons)
   real*8 :: da_ha !< drainage area of watershed in hectares (ha)
   real*8 :: pndflwo !< volume of water flowing out of pond on day (m^3 H2O)
   real*8 :: vpd !< vapor pressure deficit (kPa)
   real*8 :: wetflwo !< volume of water flowing out wetland on day (m^3 H2O)
   real*8 :: wetsedi !< sediment loading to wetland for day (metric tons)
!> leaf area index at which no evaporation occurs.  This variable is used in
!> ponded HRUs (eg rice) where evaporation from the water surface is restricted
!> by the plant canopy cover. Evaporation from the water surface equals
!> potential ET when LAI = 0 and decreased linearly to O when LAI = EVLAI
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
!> less persistent bacteria transported to main channel with surface runoff
!> (# colonies/ha)
   real*8 :: bactrolp
!> less persistent bacteria transported with sediment in surface runoff
!> (# colonies/ha)
   real*8 :: bactsedlp
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
!> amount of water in snow lost through sublimation on current day in HRU
!> (mm H2O)
   real*8 :: snoev
!> amount of nitrate moving upward in the soil profile in watershed (kg N/ha)
   real*8 :: sno3up
!> amount of pesticide in reach that is lost through reactions (mg pst)
   real*8 :: reactw
!> actual amount of evaporation (soil et) that occurs on day in HRU (mm H2O)
   real*8 :: es_day
!> average annual change in the number of less persistent bacteria colonies in
!> soil solution in watershed (# cfu/m^2)
   real*8 :: sdiegrolpq
!> average annual change in the number of less persistent bacteria colonies on
!> soil particles in watershed (# cfu/m^2)
   real*8 :: sdiegrolps
!> average annual change in the number of persistent bacteria colonies in soil
!> solution in watershed (# cfu/m^2)
   real*8 :: sdiegropq
!> average annual change in the number of persistent bacteria colonies on soil
!> particles in watershed (# cfu/m^2)
   real*8 :: sdiegrops
!> fraction for less persistent bacteria on foliage that is washed off by a
!> rainfall event (none)
   real*8 :: wof_lp
!> maximum amount of transpiration (plant et) that can occur on day in HRU
!> (mm H2O)
   real*8 :: ep_max
!> average annual number of less persistent bacteria transported to main channel
!> with surface runoff in solution (# colonies/ha)
   real*8 :: sbactrolp
!> average annual number of persistent bacteria transported to main channel with
!> surface runoff in solution (# colonies/ha)
   real*8 :: sbactrop
!> average annual number of less persistent bacteria transported with sediment
!> in surface runoff (# colonies/ha)
   real*8 :: sbactsedlp
!> average annual number of persistent bacteria transported with sediment in
!> surface runoff (# colonies/ha)
   real*8 :: sbactsedp
!> average annual number of less persistent bacteria lost from soil surface
!> layer by percolation (# cfu/m^2)
   real*8 :: sbactlchlp
!> average annual number of persistent bacteria lost from soil surface layer by
!> percolation (# cfu/m^2)
   real*8 :: sbactlchp
   real*8 :: rchwtr !< water stored in reach at beginning of day (m^3 H2O)
!> amount of pesticide moving from sediment to reach due to resuspension
!> (mg pst)
   real*8 :: resuspst
!> amount of pesticide moving from water to sediment due to settling (mg pst)
   real*8 :: setlpst
   real*8 :: psp_bsn
!> surface runoff lagged from prior day of simulation (mm H2O)
   real*8 :: bsprev
!> lateral flow lagged from prior day of simulation (mm H2O)
   real*8 :: bssprev
!> average annual amount of water removed from potholes by evaporation in
!> watershed (mm H2O)
   real*8 :: spadyev
!> average annual amount of water released to main channel from potholes in
!> watershed (mm H2O)
   real*8 :: spadyo
!> average annual amount of precipitation on potholes in watershed (mm H2O)
   real*8 :: spadyrfv
!> average annual amount of water removed from potholes by seepage in watershed
!> (mm H2O)
   real*8 :: spadysp
   real*8 :: spadyosp
!> amount of surface runoff loading to main channel from HRU on current day
!> (includes effects of transmission losses) (mm H2O)
   real*8 :: qday
!> fraction of total rainfall that occurs during 0.5h of highest intensity rain
!> (none)
   real*8 :: al5
   real*8 :: no3pcp !< nitrate added to the soil in rainfall (kg N/ha)
   real*8 :: pndsedc !< net change in sediment in pond during day (metric tons)
!> USLE rainfall erosion index on day for HRU (100(ft-tn in)/(acre-hr))
   real*8 :: usle_ei
   real*8 :: rcharea !< cross-sectional area of flow (m^2)
!> amount of pesticide lost from reach by volatilization (mg pst)
   real*8 :: volatpst
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
!> net change in sediment in wetland during day (metric tons)
   real*8 :: wetsedc
   real*8 :: respesti
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
   real*8 :: snocovmx
   real*8 :: lyrtile !< drainage tile flow in soil layer for day in HRU (mm H2O)
   real*8 :: lyrtilex
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
   real*8 :: rhoq !< algal respiration rate at 20 deg C (1/day or 1/hr)
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
!> maximum specific algal growth rate at 20 deg C(1/day or 1/hr)
   real*8 :: mumax
   real*8 :: p_n !< algal preference factor for ammonia
   real*8 :: rnum1 !< variable to hold value for rnum1s(:) (none)
!> actual evapotranspiration occuring on day in HRU (mm H2O)
   real*8 :: etday
!> amount of nitrogen applied in auto-fert application (kg N/ha)
   real*8 :: auton
!> amount of phosphorus applied in auto-fert application (kg P/ha)
   real*8 :: autop
!> amount of nitrogen moving from active organic to nitrate pool in soil profile
!> on current day in HRU (kg N/ha)
   real*8 :: hmntl
!> amount of phosphorus moving from active organic to nitrate pool in soil
!> profile on current day in HRU (kg P/ha)
   real*8 :: hmptl
!> amount of nitrogen moving from the fresh organic (residue) to the
!> nitrate(80%) and active organic(20%) pools in soil profile on current day in
!> HRU (kg N/ha)
   real*8 :: rmn2tl
!> amount of nitrogen moving from active organic to stable organic pool in soil
!> profile on current day in HRU (kg N/ha)
   real*8 :: rwntl
!> amount of water recharging deep aquifer on current day in HRU (mm H2O)
   real*8 :: gwseep
!> amount of water moving from the shallow aquifer into the soil profile or
!> being taken up by plant roots in the shallow aquifer or in the bank storage
!> zone (mm H2O)
   real*8 :: revapday
!> amount of phosphorus moving from the labile mineral pool to the active
!> mineral pool in the soil profile on the current day in the HRU (kg P/ha)
   real*8 :: rmp1tl
!> amount of phosphorus moving from the fresh organic (residue) to the
!> labile(80%) and organic(20%) pools in soil profile on current day in HRU
!> (kg P/ha)
   real*8 :: rmptl
!> amount of phosphorus moving from the active mineral pool to the stable
!> mineral pool in the soil profile on the current day in the HRU (kg P/ha)
   real*8 :: roctl
!> amount of nitrogen lost from nitrate pool by denitrification in soil profile
!> on current day in HRU (kg N/ha)
   real*8 :: wdntl
   real*8 :: cmn_bsn,reswtr
!> die-off factor for less persistent bacteria in streams (1/day)
   real*8 :: wdlprch
!> die-off factor for persistent bacteria in reservoirs (1/day)
   real*8 :: wdpres
!> potential ET value read in for day (mm H2O)
   real*8 :: petmeas
!> loss of pesticide from active sediment layer by burial (mg pst)
   real*8 :: bury
!> diffusion of pesticide from sediment to reach (mg pst)
   real*8 :: difus
!> amount of pesticide in sediment that is lost through reactions (mg pst)
   real*8 :: reactb
!> soluble pesticide concentration in outflow on day (mg pst/m^3)
   real*8 :: solpesto
!> die-off factor for less persistent bacteria in reservoirs (1/day)
   real*8 :: wdlpres
!> sorbed pesticide concentration in outflow on day (mg pst/m^3)
   real*8 :: sorpesto
   real*8 :: spcon_bsn, spexp_bsn, solpesti, sorpesti
!> calibration coefficient to control impact of the storage time constant for
!> the reach at bankfull depth (phi(10,:) upon the storage time constant for the
!> reach used in the Muskingum flow method
   real*8 :: msk_co1
!> calibration coefficient to control impact of the storage time constant for
!> the reach at 0.1 bankfull depth (phi(13,:) upon the storage time constant for
!> the reach used in the Muskingum flow method
   real*8 :: msk_co2
   real*8 :: deepstp !< depth of water in deep aquifer in HRU (mm H2O)
!> depth of water in shallow aquifer in HRU on previous day (mm H2O)
   real*8 :: shallstp
   real*8 :: snoprev !< amount of water stored as snow on previous day (mm H2O)
!> amount of water stored in soil profile in the HRU on the previous day
!> (mm H2O)
   real*8 :: swprev
   real*8 :: ressolpo, resorgno, resorgpo, resno3o, reschlao, resno2o
!> volume of water evaporated from pothole expressed as depth over HRU (mm H2O)
   real*8 :: potevmm
!> volume of water released to main channel from pothole expressed as depth
!> over HRU (mm H2O)
   real*8 :: potflwo
!> precipitation falling on pothole water body expressed as depth over HRU
!> (mm H2O)
   real*8 :: potpcpmm
!> seepage from pothole expressed as depth over HRU (mm H2O)
   real*8 :: potsepmm
   real*8 :: qdbank !< streamflow contribution from bank storage (m^3 H2O)
   real*8 :: resnh3o
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
!> sediment leaving pothole to main channel from HRU on day (metric tons/ha)
   real*8 :: potsedo
   real*8 :: pest_sol
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

   real*8 :: hlife_ngw_bsn !< Half-life of nitrogen in groundwater? (days) 
   real*8 :: ch_opco_bsn, ch_onco_bsn
   real*8 :: decr_min !< Minimum daily residue decay
   real*8 :: rcn_sub_bsn !< Concentration of nitrogen in the rainfall (mg/kg)
   real*8 :: bc1_bsn, bc2_bsn, bc3_bsn, bc4_bsn
   real*8 :: anion_excl_bsn

!> water table based on depth from soil surface (mm)
   real*8, dimension (:), allocatable :: wat_tbl
   real*8, dimension (:), allocatable :: sol_swpwt
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
   integer :: idlast !< number of days simulated in month (none)
   integer :: i_subhw, imgt, iwtr, ifrttyp, mo_atmo, mo_atmo1
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
   integer :: nres !< total number of reservoirs in watershed (none)
!> number of last HRU in previous subbasin or\n
!> number of HRUs in watershed (none)
   integer :: nhru
!> current month being simulated or month of next day of simulation (none)
   integer :: i_mo
   integer :: immo !< current cumulative month of simulation (none)
   integer :: mo
!> wind speed input code (noen)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: wndsim
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
!> 1: simulate crack flow in watershed
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
!> tile drainage equations flag/code\n
!> 1 simulate tile flow using subroutine drains(wt_shall)\n
!> 0 simulate tile flow using subroutine origtile(wt_shall,d)
   integer :: itdrn
!> water table depth algorithms flag/code\n
!> 1 simulate wt_shall using subroutine new water table depth routine\n
!> 0 simulate wt_shall using subroutine original water table depth routine
   integer :: iwtdn
!> maximum depressional storage selection flag/code (none)\n
!> 0 = static depressional storage (stmaxd) read from .bsn for the global value
!> or .sdr for specific HRUs\n
!> 1 = dynamic storage (stmaxd) based on random roughness, tillage and
!> cumulative rainfall intensity by depstor.f90
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
!> number of years to skip output summarization and printing (none)
   integer :: nyskip
!> solar radiation input code (none)\n
!> 1 measured data read for each subbasin\n
!> 2 data simulated for each subbasin
   integer :: slrsim
!> channel degredation code\n
!> 0: do not compute channel degradation\n
!> 1: compute channel degredation (downcutting and widening)
   integer :: ideg
!> rainfall/runoff code (none)\n
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
   integer :: mo_chk !< current month of simulation (none)
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
!> streamflow print code (none)\n
!> 0 print streamflow in reach\n
!> 1 print Log10 streamflow in reach
   integer :: ilog
   integer :: itotr !< number of output variables printed (output.rch)
   integer :: iyr !< current year of simulation (year)
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
   integer :: iprint !< print code (none): 0=monthly, 1=daily, 2=annually
   integer :: iida !< day being simulated (current julian date) (julian date)
!> CN method flag (for testing alternative method):\n
!> 0 use traditional SWAT method which bases CN on soil moisture\n
!> 1 use alternative method which bases CN on plant ET\n
!> 2 use tradtional SWAT method which bases CN on soil moisture but rention is
!> adjusted for mildly-sloped tiled-drained watersheds
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
   integer :: mapex
   real*8, dimension (:), allocatable :: flodaya, seddaya, orgndaya
   real*8, dimension (:), allocatable :: orgpdaya, no3daya, minpdaya
!> harvest index target of cover defined at planting ((kg/ha)/(kg/ha))
   real*8, dimension (:), allocatable :: hi_targ
   real*8, dimension (:), allocatable :: bio_targ !< biomass target (kg/ha)
!> modifier for autofertilization target nitrogen content for plant
!> (kg N/kg yield)
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
!> BOD concentration in biozone (kg/ha)
   real*8, dimension (:), allocatable :: bio_bod
!> biomass of live bacteria in biozone (kg/ha)
   real*8, dimension (:), allocatable :: biom
!> daily change in biomass of live bacteria (kg/ha)
   real*8, dimension (:), allocatable :: rbiom
   real*8, dimension (:), allocatable :: bio_amn
!> concentration of the fecal coliform in the biozone septic tank effluent
!> (cfu/100ml)
   real*8, dimension (:), allocatable :: fcoli
   real*8, dimension (:), allocatable :: bio_ntr, bz_perc
!> number of permanent residents in the hourse (none)
   real*8, dimension (:), allocatable :: sep_cap
   real*8, dimension (:), allocatable :: plqm !< plaque in biozone (kg/ha)
   real*8, dimension (:), allocatable :: bz_area
   real*8, dimension (:), allocatable :: bz_z !< depth of biozone layer (mm)
   real*8, dimension (:), allocatable :: bz_thk !< thickness of biozone (mm)
   real*8, dimension (:), allocatable :: bio_bd !< density of biomass (kg/m^3)
!> current soil carbon for first soil layer (kg/ha)
   real*8, dimension (:), allocatable :: cmup_kgh
!> current soil carbon integrated - aggregating (kg/ha)
   real*8, dimension (:), allocatable :: cmtot_kgh
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
!> septic system type (none)
   integer, dimension (:), allocatable :: isep_typ
!> soil layer where biozone exists (none)
   integer, dimension (:), allocatable :: i_sep
!> septic system operation flag (1=active, 2=failing, 3 or 0=not operated)
!> (none)
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
!> amount of organic P in pothole water body (kg P)
   real*8, dimension (:), allocatable :: pot_orgp
   real*8, dimension (:), allocatable :: pot_orgpi
!> amount of organic N in pothole water body (kg N)
   real*8, dimension (:), allocatable :: pot_orgn
   real*8, dimension (:), allocatable :: pot_orgni
!> amount of stable mineral pool P in pothole water body (kg N)
   real*8, dimension (:), allocatable :: pot_mps
   real*8, dimension (:), allocatable :: pot_mpsi
!> amount of active mineral pool P in pothole water body (kg N)
   real*8, dimension (:), allocatable :: pot_mpa
   real*8, dimension (:), allocatable :: pot_mpai
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
!> wshddayo(1) average amountof precipitation in watershed for the day
!> (mm H20)\n
!> wshddayo(3) surface runoff in watershed for day (mm H20)\n
!> wshddayo(4) lateral flow contribution to streamflow in watershed for day
!> (mm H20)\n
!> wshddayo(5) water percolation past bottom of soil profile in watershed for
!> day (mm H20)\n
!> wshddayo(6) water yield to streamflow from HRUs in watershed for day
!> (mm H20)\n
!> wshddayo(7) actual evapotranspiration in watershed for day (mm H20)\n
!> wshddayo(8) average maximum temperature in watershed for the day (deg C)\n
!> wshddayo(9) average minimum temperature in watershed for the day (deg C)\n
!> wshddayo(12) sediment yield from HRUs in watershed for day
!> (metric tons or metric tons/ha)\n
!> wshddayo(13) sediment loading to ponds in watershed for day (metric tons)\n
!> wshddayo(14) sediment loading from ponds in watershed for day (metric tons)\n
!> wshddayo(15) net change in sediment level in ponds in watershed for day
!> (metric tons)\n
!> wshddayo(16) sediment loading to wetlands for day in watershed
!> (metric tons)\n
!> wshddayo(17) sediment loading to main channels from wetlands for day in
!> watershed (metric tons)\n
!> wshddayo(18) net change in sediment in wetlands for day in watershed
!> (metric tons)\n
!> wshddayo(19) evaporation from ponds in watershed for day (m^3 H2O)\n
!> wshddayo(20) seepage from ponds in watershed for day (m^3 H2O)\n
!> wshddayo(21) precipitation on ponds in watershed for day (m^3 H2O)\n
!> wshddayo(22) volume of water entering ponds in watershed for day (m^3 H2O)\n
!> wshddayo(23) volume of water leaving ponds in watershed for day (m^3 H2O)\n
!> wshddayo(24) evaporation from wetlands for day in watershed (m^3 H2O)\n
!> wshddayo(25) seepage from wetlands for day in watershed (m^3 H2O)\n
!> wshddayo(26) precipitation on wetlands for day in watershed (m^3 H2O)\n
!> wshddayo(27) volume of water entering wetlands on day in watershed
!> (m^3 H2O)\n
!> wshddayo(28) volume of water leaving wetlands on day in watershed (m^3 H2O)\n
!> wshddayo(33) net change in water volume of ponds in watershed for day
!> (m^3 H2O)\n
!> wshddayo(35) amount of water stored in soil profile in watershed at end of
!> day (mm H20)\n
!> wshddayo(36) snow melt in watershed for day (mm H20)\n
!> wshddayo(37) sublimation in watershed for day (mm H20)\n
!> wshddayo(38) average amount of tributary channel transmission losses in
!> watershed on day (mm H20)\n
!> wshddayo(39) freezing rain/snow fall in watershed for day (mm H20)\n
!> wshddayo(40) organic N loading to stream in watershed for day (kg N/ha)\n
!> wshddayo(41) organic P loading to stream in watershed for day (kg P/ha)\n
!> wshddayo(42) nitrate loading to stream in surface runoff in watershed for day
!> (kg N/ha)\n
!> wshddayo(43) soluble P loading to stream in watershed for day (kg P/ha)\n
!> wshddayo(44) plant uptake of N in watershed for day (kg N/ha)\n
!> wshddayo(45) nitrate loading to stream in lateral flow in watershed for day
!> (kg N/ha)\n
!> wshddayo(46) nitrate percolation past bottom of soil profile in watershed for
!> day (kg N/ha)\n
!> wshddayo(104) groundwater contribution to stream in watershed on day
!> (mm H20)\n
!> wshddayo(105) amount of water moving from shallow aquifer to plants/soil
!> profile in watershed on day (mm H2O)\n
!> wshddayo(106) deep aquifer recharge in watershed on day (mm H2O)\n
!> wshddayo(107) total amount of water entering both aquifers in watershed on
!> day (mm H2O)\n
!> wshddayo(108) potential evapotranspiration in watershed on day (mm H20)\n
!> wshddayo(109) drainage tile flow contribution to stream in watershed on day
!> (mm H20)\n
!> wshddayo(110) NO3 yield (gwq) (kg/ha)\n
!> wshddayo(111) NO3 yield (tile) (mm H2O)
   real*8, dimension (mstdo) :: wshddayo
!> watershed monthly output array (see definitions for wshddayo array elements)
!> (varies)\n
!> wshdmono(1) average amount of precipitation in watershed for the month
!> (mm H2O)\n
!> wshdmono(3) surface runoff in watershed for month (mm H2O)\n
!> wshdmono(4) lateral flow contribution to streamflow in watershed for month
!> (mm H2O)\n
!> wshdmono(5) water percolation past bottom of soil profile in watershed for
!> month (mm H2O)\n
!> wshdmono(6) water yield to streamflow from HRUs in watershed for month
!> (mm H2O)\n
!> wshdmono(7) actual evapotranspiration in watershed for month (mm H2O)\n
!> wshdmono(8) average maximum temperature in watershed for the month (deg C)\n
!> wshdmono(9) average minimum temperature in watershed for the month (deg C)\n
!> wshdmono(12) sediment yield from HRUs in watershed for the month
!> (metric tons)\n
!> wshdmono(39) freezing rain/snow fall in watershed for the month (mm H2O)\n
!> wshdmono(40) organic N loading to stream in watershed for the month
!> (kg N/ha)\n
!> wshdmono(41) organic P loading to stream in watershed for the month
!> (kg P/ha)\n
!> wshdmono(42) nitrate loading to stream in surface runoff in watershed for the
!> month (kg N/ha)\n
!> wshdmono(43) soluble P loading to stream in watershed for the month
!> (kg P/ha)\n
!> wshdmono(44) plant uptake of N in watershed for the month (kg N/ha)\n
!> wshdmono(45) nitrate loading to stream in lateral flow in watershed for the
!> month (kg N/ha)\n
!> wshdmono(46) nitrate percolation past bottom of soil profile in watershed
!> for the month (kg N/ha)\n
!> wshdmono(104) groundwater contribution to stream in watershed for the month
!> (mm H2O)\n
!> wshdmono(108) potential evapotranspiration in watershed for the month
!> (mm H2O)\n
!> wshdmono(109) drainage tile flow contribution to stream in watershed for the
!> month (mm H2O)
   real*8, dimension (mstdo) :: wshdmono
!> watershed annual output array (varies)\n
!> wshdyro(1) average amount of precipitation in watershed for the year
!> (mm H2O)\n
!> wshdyro(3) surface runoff in watershed for year (mm H2O)\n
!> wshdyro(4) lateral flow contribution to streamflow in watershed for year
!> (mm H2O)\n
!> wshdyro(5) water percolation past bottom of soil profile in watershed for
!> year (mm H2O)\n
!> wshdyro(6) water yield to streamflow from HRUs in watershed for year
!> (mm H2O)\n
!> wshdyro(7) actual evapotranspiration in watershed for year (mm H2O)\n
!> wshdyro(8) average maximum temperature in watershed for the year (deg C)\n
!> wshdyro(9) average minimum temperature in watershed for the year (deg C)\n
!> wshdyro(12) sediment yield from HRUs in watershed for the year
!> (metric tons)\n
!> wshdyro(40) organic N loading to stream in watershed for the year (kg N/ha)\n
!> wshdyro(41) organic P loading to stream in watershed for the year (kg P/ha)\n
!> wshdyro(42) nitrate loading to stream in surface runoff in watershed for the
!> year (kg N/ha)\n
!> wshdyro(43) soluble P loading to stream in watershed for the year (kg P/ha)\n
!> wshdyro(44) plant uptake of N in watershed for the year\n
!> wshdyro(45) nitrate loading to stream in lateral flow in watershed for the
!> year (kg N/ha)\n
!> wshdyro(46) nitrate percolation past bottom of soil profile in watershed for
!> the year (kg N/ha)\n
!> wshdyro(104) groundwater contribution to stream in watershed for the year
!> (mm H2O)\n
!> wshdyro(108) potential evapotranspiration in watershed for the year
!> (mm H2O)\n
!> wshdyro(109) drainage tile flow contribution to stream in watershed for the
!> year (mm H2O)
   real*8, dimension (mstdo) :: wshdyro
   real*8, dimension (16) :: fcstaao
!> watershed average annual output array (varies)\n
!> wshdaao(1) precipitation in watershed (mm H2O)\n
!> wshdaao(3) surface runoff loading to main channel in watershed (mm H2O)\n
!> wshdaao(4) lateral flow loading to main channel in watershed (mm H2O)\n
!> wshdaao(5) percolation of water out of root zone in watershed (mm H2O)\n
!> wshdaao(7) actual evapotranspiration in watershed (mm H2O)\n
!> wshdaao(13) sediment loading to ponds in watershed (metric tons)\n
!> wshdaao(14) sediment loading from ponds in watershed (metric tons)\n
!> wshdaao(15) net change in sediment level in ponds in watershed (metric tons)\n
!> wshdaao(19) evaporation from ponds in watershed (m^3 H2O)\n
!> wshdaao(20) seepage from ponds in watershed (m^3 H2O)\n
!> wshdaao(21) precipitation on ponds in watershed (m^3 H2O)\n
!> wshdaao(22) volume of water entering ponds in watershed (m^3 H2O)\n
!> wshdaao(23) volume of water leaving ponds in watershed (m^3 H2O)\n
!> wshdaao(38) transmission losses in watershed (mm H2O)
   real*8, dimension (mstdo) :: wshdaao
!> wpstdayo(1,:) amount of pesticide type in surface runoff contribution to
!> stream in watershed on day (in solution) (mg pst/ha)\n
!> wpstdayo(2,:) amount of pesticide type in surface runoff contribution to
!> stream in watershed on day (sorbed to sediment) (mg pst/ha)\n
!> wpstdayo(3,:) amount of pesticide type leached from soil profile in watershed
!> on day (kg pst/ha)\n
!> wpstdayo(4,:) amount of pesticide type in lateral flow contribution to stream
!> in watershed on day (kg pst/ha)
   real*8, dimension (:,:), allocatable :: wpstdayo
   real*8, dimension (:,:), allocatable :: wpstmono,wpstyro
!> harvested biomass (dry weight) (kg/ha)
   real*8, dimension (:,:), allocatable :: bio_hv
!> yield (dry weight) by crop type in the HRU (kg/ha)
   real*8, dimension (:,:), allocatable :: yldkg
!> reach monthly output array (varies)\n
!> rchmono(1,:) flow into reach during month (m^3/s)\n
!> rchmono(2,:) flow out of reach during month (m^3/s)\n
!> rchmono(3,:) sediment transported into reach during month (metric tons)\n
!> rchmono(4,:) sediment transported out of reach during month (metric tons)\n
!> rchmono(5,:) sediment concentration in outflow during month (mg/L)\n
!> rchmono(6,:) organic N transported into reach during month (kg N)\n
!> rchmono(7,:) organic N transported out of reach during month (kg N)\n
!> rchmono(8,:) organic P transported into reach during month (kg P)\n
!> rchmono(9,:) organic P transported out of reach during month (kg P)\n
!> rchmono(10,:) evaporation from reach during month (m^3/s)\n
!> rchmono(11,:) transmission losses from reach during month (m^3/s)\n
!> rchmono(12,:) conservative metal #1 transported out of reach during month
!> (kg)\n
!> rchmono(13,:) conservative metal #2 transported out of reach during month
!> (kg)\n
!> rchmono(14,:) conservative metal #3 transported out of reach during month
!> (kg)\n
!> rchmono(15,:) nitrate transported into reach during month (kg N)\n
!> rchmono(16,:) nitrate transported out of reach during month (kg N)\n
!> rchmono(17,:) soluble P transported into reach during month (kg P)\n
!> rchmono(18,:) soluble P transported out of reach during month (kg P)\n
!> rchmono(19,:) soluble pesticide transported into reach during month
!> (mg pst)\n
!> rchmono(20,:) soluble pesticide transported out of reach during month
!> (mg pst)\n
!> rchmono(21,:) sorbed pesticide transported into reach during month (mg pst)\n
!> rchmono(22,:) sorbed pesticide transported out of reach during month
!> (mg pst)\n
!> rchmono(23,:) amount of pesticide lost through reactions in reach during
!> month (mg pst)\n
!> rchmono(24,:) amount of pesticide lost through volatilization from reach
!> during month (mg pst)\n
!> rchmono(25,:) amount of pesticide settling out of reach to bed sediment
!> during month (mg pst)\n
!> rchmono(26,:) amount of pesticide resuspended from bed sediment to reach
!> during month (mg pst)\n
!> rchmono(27,:) amount of pesticide diffusing from reach to bed sediment during
!> month (mg pst)\n
!> rchmono(28,:) amount of pesticide in sediment layer lost through reactions
!> during month (mg pst)\n
!> rchmono(29,:) amount of pesticide in sediment layer lost through burial
!> during month (mg pst)\n
!> rchmono(30,:) chlorophyll-a transported into reach during month (kg chla)\n
!> rchmono(31,:) chlorophyll-a transported out of reach during month (kg chla)\n
!> rchmono(32,:) ammonia transported into reach during month (kg N)\n
!> rchmono(33,:) ammonia transported out of reach during month (kg N)\n
!> rchmono(34,:) nitrite transported into reach during month (kg N)\n
!> rchmono(35,:) nitrite transported out of reach during month (kg N)\n
!> rchmono(36,:) CBOD transported into reach during month (kg O2)\n
!> rchmono(37,:) CBOD transported out of reach during month (kg O2)\n
!> rchmono(38,:) dissolved oxygen transported into reach during month (kg O2)\n
!> rchmono(39,:) dissolved oxygen transported out of reach during month
!> (kg O2)\n
!> rchmono(40,:) persistent bacteria transported out of reach during month
!> (kg bact)\n
!> rchmono(41,:) less persistent bacteria transported out of reach during month
!> (kg bact)\n
!> rchmono(43,:) total N (org N + no3 + no2 + nh4 outs) (kg)\n
!> rchmono(44,:) total P (org P + sol p outs) (kg)
   real*8, dimension (:,:), allocatable :: rchmono
!> reach annual output array (varies)\n
!> rchyro(1,:) flow into reach during year (m^3/s)\n
!> rchyro(2,:) flow out of reach during year (m^3/s)\n
!> rchyro(3,:) sediment transported into reach during year (metric tons)\n
!> rchyro(4,:) sediment transported out of reach during year (metric tons)\n
!> rchyro(5,:) sediment concentration in outflow during year (mg/L)\n
!> rchyro(6,:) organic N transported into reach during year (kg N)\n
!> rchyro(7,:) organic N transported out of reach during year (kg N)\n
!> rchyro(8,:) organic P transported into reach during year (kg P)\n
!> rchyro(9,:) organic P transported out of reach during year (kg P)\n
!> rchyro(10,:) evaporation from reach during year (m^3/s)\n
!> rchyro(11,:) transmission losses from reach during year (m^3/s)\n
!> rchyro(12,:) conservative metal #1 transported out of reach during year
!> (kg)\n
!> rchyro(13,:) conservative metal #2 transported out of reach during year
!> (kg)\n
!> rchyro(14,:) conservative metal #3 transported out of reach during year
!> (kg)\n
!> rchyro(15,:) nitrate transported into reach during year (kg N)\n
!> rchyro(16,:) nitrate transported out of reach during year (kg N)\n
!> rchyro(17,:) soluble P transported into reach during year (kg P)\n
!> rchyro(18,:) soluble P transported out of reach during year (kg P)\n
!> rchyro(19,:) soluble pesticide transported into reach during year (mg pst)\n
!> rchyro(20,:) soluble pesticide transported out of reach during year
!> (mg pst)\n
!> rchyro(21,:) sorbed pesticide transported into reach during year (mg pst)\n
!> rchyro(22,:) sorbed pesticide transported out of reach during year (mg pst)\n
!> rchyro(23,:) amount of pesticide lost through reactions in reach during
!> year!> (mg pst)\n
!> rchyro(24,:) amount of pesticide lost through volatilization from reach
!> during year (mg pst)\n
!> rchyro(25,:) amount of pesticide settling out of reach to bed sediment during
!> year (mg pst)\n
!> rchyro(26,:) amount of pesticide resuspended from bed sediment to reach
!> during year (mg pst)\n
!> rchyro(27,:) amount of pesticide diffusing from reach to bed sediment during
!> year (mg pst)\n
!> rchyro(28,:) amount of pesticide in sediment layer lost through reactions
!> during year (mg pst)\n
!> rchyro(29,:) amount of pesticide in sediment layer lost through burial during
!> year (mg pst)\n
!> rchyro(30,:) chlorophyll-a transported into reach during year (kg chla)\n
!> rchyro(31,:) chlorophyll-a transported out of reach during year (kg chla)\n
!> rchyro(32,:) ammonia transported into reach during year (kg N)\n
!> rchyro(33,:) ammonia transported out of reach during year (kg N)\n
!> rchyro(34,:) nitrite transported into reach during year (kg N)\n
!> rchyro(35,:) nitrite transported out of reach during year (kg N)\n
!> rchyro(36,:) CBOD transported into reach during year (kg O2)\n
!> rchyro(37,:) CBOD transported out of reach during year (kg O2)\n
!> rchyro(38,:) dissolved oxygen transported into reach during year (kg O2)\n
!> rchyro(39,:) dissolved oxygen transported out of reach during year (kg O2)\n
!> rchyro(40,:) persistent bacteria transported out of reach during year
!> (kg bact)\n
!> rchyro(41,:) less persistent bacteria transported out of reach during year
!> (kg bact)
   real*8, dimension (:,:), allocatable :: rchyro
   real*8, dimension (:,:), allocatable :: wpstaao
!> HRU monthly output data array (varies)\n
!> hrumono(1,:) precipitation in HRU during month (mm H2O)\n
!> hrumono(2,:) amount of precipitation falling as freezing rain/snow in HRU
!> during month (mm H2O)\n
!> hrumono(3,:) amount of snow melt in HRU during month (mm H2O)\n
!> hrumono(4,:) amount of surface runoff to main channel from HRU during month
!> (ignores impact of transmission losses) (mm H2O)\n
!> hrumono(5,:) amount of lateral flow contribution to main channel from HRU
!> during month (mm H2O)\n
!> hrumono(6,:) amount of groundwater flow contribution to main channel from HRU
!> during month (mm H2O)\n
!> hrumono(7,:) amount of water moving from shallow aquifer to plants or soil
!> profile in HRU during mont (mm H2O)h
!> hrumono(8,:) amount of water recharging deep aquifer in HRU during month
!> (mm H2O)\n
!> hrumono(9,:) total amount of water entering both aquifers from HRU during
!> month (mm H2O)\n
!> hrumono(10,:) water yield (total amount of water entering main channel) from
!> HRU during month (mm H2O)\n
!> hrumono(11,:) amount of water percolating out of the soil profile and into
!> the vadose zone in HRU during month (mm H2O)\n
!> hrumono(12,:) actual evapotranspiration in HRU during month (mm H2O)\n
!> hrumono(13,:) amount of transmission losses from tributary channels in HRU
!> for month (mm H2O)\n
!> hrumono(14,:) sediment yield from HRU for month (metric tons/ha)\n
!> hrumono(15,:)  actual amount of transpiration that occurs during month in HRU
!> (mm H2O)\n
!> hrumono(16,:) actual amount of evaporation (from soil) that occurs during
!> month in HRU (mm H2O)\n
!> hrumono(17,:) amount of nitrogen applied in continuous fertilizer operation
!> during month in HRU (kg N/ha)\n
!> hrumono(18,:) amount of phosphorus applied in continuous fertilizer operation
!> during month in HRU (kg P/ha)\n
!> hrumono(19,:) amount of surface runoff generated during month in HRU
!> (mm H2O)\n
!> hrumono(20,:) CN values during month in HRU (none)\n
!> hrumono(21,:) sum of daily soil water values used to calculate the curve
!> number (mm H2O)\n
!> hrumono(22,:) amount of irrigation water applied to HRU during month
!> (mm H2O)\n
!> hrumono(23,:) amount of water removed from shallow aquifer in HRU for
!> irrigation during month (mm H2O)\n
!> hrumono(24,:) amount of water removed from deep aquifer in HRU for irrigation
!> during month (mm H2O)\n
!> hrumono(25,:) potential evapotranspiration in HRU during month (mm H2O)\n
!> hrumono(26,:) monthly amount of N (organic & mineral) applied in HRU during
!> grazing (kg N/ha)\n
!> hrumono(27,:) monthly amount of P (organic & mineral) applied in HRU during
!> grazing (kg P/ha)\n
!> hrumono(28,:) monthly amount of N (organic & mineral) auto-applied in HRU
!> (kg N/ha)\n
!> hrumono(29,:) monthly amount of P (organic & mineral) auto-applied in HRU
!> (kg P/ha)\n
!> hrumono(30,:) sum of daily soil temperature values (deg C)
!> hrumono(31,:) water stress days in HRU during month (stress days)\n
!> hrumono(32,:) temperature stress days in HRU during month (stress days)\n
!> hrumono(33,:) nitrogen stress days in HRU during month (stress days)\n
!> hrumono(34,:) phosphorus stress days in HRU during month (stress days)\n
!> hrumono(35,:) organic nitrogen in surface runoff in HRU during month
!> (kg N/ha)\n
!> hrumono(36,:) organic phosphorus in surface runoff in HRU during month
!> (kg P/ha)\n
!> hrumono(37,:) nitrate in surface runoff in HRU during month (kg N/ha)\n
!> hrumono(38,:) nitrate in lateral flow in HRU during month (kg N/ha)\n
!> hrumono(39,:) soluble phosphorus in surface runoff in HRU during month
!> (kg P/ha)\n
!> hrumono(40,:) amount of nitrogen removed from soil by plant uptake in HRU
!> during month (kg N/ha)\n
!> hrumono(41,:) nitrate percolating past bottom of soil profile in HRU during
!> month (kg N/ha)\n
!> hrumono(42,:) amount of phosphorus removed from soil by plant uptake in HRU
!> during month (kg P/ha)\n
!> hrumono(43,:) amount of phosphorus moving from labile mineral to active
!> mineral pool in HRU during month (kg P/ha)\n
!> hrumono(44,:) amount of phosphorus moving from active mineral to stable
!> mineral pool in HRU during month (kg P/ha)\n
!> hrumono(45,:) amount of nitrogen applied to HRU in fertilizer and grazing
!> operations during month (kg N/ha)\n
!> hrumono(46,:) amount of phosphorus applied to HRU in fertilizer and grazing
!> operations during month (kg P/ha)\n
!> hrumono(47,:) amount of nitrogen added to soil by fixation in HRU during
!> month (kg N/ha)\n
!> hrumono(48,:) amount of nitrogen lost by denitrification in HRU during month
!> (kg N/ha)\n
!> hrumono(49,:) amount of nitrogen moving from active organic to nitrate pool
!> in HRU during month (kg N/ha)\n
!> hrumono(50,:) amount of nitrogen moving from active organic to stable organic
!> pool in HRU during month (kg N/ha)\n
!> hrumono(51,:) amount of phosphorus moving from organic to labile mineral
!> pool in HRU during month (kg P/ha)\n
!> hrumono(52,:) amount of nitrogen moving from fresh organic to nitrate and
!> active organic pools in HRU during month (kg N/ha)\n
!> hrumono(53,:) amount of phosphorus moving from fresh organic to the labile
!> mineral and organic pools in HRU during month (kg P/ha)\n
!> hrumono(54,:) amount of nitrogen added to soil in rain (kg N/ha)\n
!> hrumono(61,:) daily soil loss predicted with USLE equation (metric tons/ha)\n
!> hrumono(62,:) drainage tile flow contribution to main channel from HRU in
!> month (mm H2O)\n
!> hrumono(63,:) less persistent bacteria transported to main channel from HRU
!> during month (#bacteria/ha)\n
!> hrumono(64,:) persistent bacteria transported to main channel from HRU during
!> month (#bacteria/ha)\n
!> hrumono(65,:) nitrate loading from groundwater in HRU to main channel during
!> month (kg N/ha)\n
!> hrumono(66,:) soluble P loading from groundwater in HRU to main channel
!> during month (kg P/ha)\n
!> hrumono(67,:) loading of mineral P attached to sediment in HRU to main
!> channel during month (kg P/ha)
   real*8, dimension (:,:), allocatable :: hrumono
!> rchdy(1,:) flow into reach on day (m^3/s)\n
!> rchdy(2,:) flow out of reach on day (m^3/s)\n
!> rchdy(3,:) evaporation from reach on day (m^3/s)\n
!> rchdy(4,:) transmission losses from reach on day (m^3/s)\n
!> rchdy(5,:) sediment transported into reach on day (metric tons)\n
!> rchdy(6,:) sediment transported out of reach on day (metric tons)\n
!> rchdy(7,:) sediment concentration in outflow (mg/L)\n
!> rchdy(8,:) organic N transported into reach on day (kg N)\n
!> rchdy(9,:) organic N transported out of reach on day (kg N)\n
!> rchdy(10,:) organic P transported into reach on day (kg P)\n
!> rchdy(11,:) organic P transported out of reach on day (kg P)\n
!> rchdy(12,:) nitrate transported into reach on day (kg N)\n
!> rchdy(13,:) nitrate transported out of reach on day (kg N)\n
!> rchdy(14,:) ammonia transported into reach on day (kg N)\n
!> rchdy(15,:) ammonia transported out of reach on day (kg N)\n
!> rchdy(16,:) nitrite transported into reach on day (kg N)\n
!> rchdy(17,:) nitrite transported out of reach on day (kg N)\n
!> rchdy(18,:) soluble P transported into reach on day (kg P)\n
!> rchdy(19,:) soluble P transported out of reach on day (kg P)\n
!> rchdy(20,:) chlorophyll-a transported into reach on day (kg chla)\n
!> rchdy(21,:) chlorophyll-a transported out of reach on day (kg chla)\n
!> rchdy(22,:) CBOD transported into reach on day (kg O2)\n
!> rchdy(23,:) CBOD transported out of reach on day (kg O2)\n
!> rchdy(24,:) dissolved oxygen transported into reach on day (kg O2)\n
!> rchdy(25,:) dissolved oxygen transported out of reach on day (kg O2)\n
!> rchdy(26,:) soluble pesticide transported into reach on day (mg pst)\n
!> rchdy(27,:) soluble pesticide transported out of reach o day (mg pst)\n
!> rchdy(28,:) sorbed pesticide transported into reach on day (mg pst)\n
!> rchdy(29,:) sorbed pesticide transported out of reach on day (mg pst)\n
!> rchdy(30,:) amount of pesticide lost through reactions in reach on day
!> (mg pst)\n
!> rchdy(31,:) amount of pesticide lost through volatilization from reach on
!> day (mg pst)\n
!> rchdy(32,:) amount of pesticide settling out of reach to bed sediment on
!> day (mg pst)\n
!> rchdy(33,:) amount of pesticide resuspended from bed sediment to reach on
!> day (mg pst)\n
!> rchdy(34,:) amount of pesticide diffusing from reach to bed sediment on day
!> (mg pst)\n
!> rchdy(35,:) amount of pesticide in sediment layer lost through reactions on
!> day (mg pst)\n
!> rchdy(36,:) amount of pesticide in sediment layer lost through burial on
!> day (mg pst)\n
!> rchdy(37,:) amount of pesticide stored in river bed sediments (mg pst)\n
!> rchdy(38,:) persistent bacteria transported out of reach on day (kg bact)\n
!> rchdy(39,:) less persistent bacteria transported out of reach on day
!> (kg bact)\n
!> rchdy(40,:) amount of conservative metal #1 transported out of reach on day
!> (kg)\n
!> rchdy(41,:) amount of conservative metal #2 transported out of reach on day
!> (kg)\n
!> rchdy(42,:) amount of conservative metal #3 transported out of reach on day
!> (kg)\n
!> rchdy(43,:) total N (org N + no3 + no2 + nh4 outs) (kg)\n
!> rchdy(44,:) total P (org P + sol p outs) (kg)
   real*8, dimension (:,:), allocatable :: rchdy
!> HRU annual output array (varies)
!> hruyro(1,:) precipitation in HRU during year (mm H2O)\n
!> hruyro(2,:) amount of precipitation falling as freezing rain/snow in HRU
!> during year (mm H2O)\n
!> hruyro(3,:) amount of snow melt in HRU during year (mm H2O)\n
!> hruyro(4,:) amount of surface runoff to main channel from HRU during year
!> (ignores impact of transmission losses) (mm H2O)\n
!> hruyro(5,:) amount of lateral flow contribution to main channel from HRU
!> during year (mm H2O)\n
!> hruyro(6,:) amount of groundwater flow contribution to main channel from HRU
!> during year (mm H2O)\n
!> hruyro(7,:) amount of water moving from shallow aquifer to plants or soil
!> profile in HRU during year (mm H2O)\n
!> hruyro(8,:) amount of water recharging deep aquifer in HRU during year
!> (mm H2O)\n
!> hruyro(9,:) total amount of water entering both aquifers from HRU during year
!> (mm H2O)\n
!> hruyro(10,:) water yield (total amount of water entering main channel) from
!> HRU during year (mm H2O)\n
!> hruyro(11,:) amount of water percolating out of the soil profile and into the
!> vadose zone in HRU during year (mm H2O)\n
!> hruyro(12,:) actual evapotranspiration in HRU during year (mm H2O)\n
!> hruyro(13,:) amount of transmission losses from tributary channels in HRU for
!> year (mm H2O)\n
!> hruyro(14,:) sediment yield from HRU for year (metric tons/ha)\n
!> hruyro(15,:) actual amount of transpiration that occurs during year in HRU
!> (mm H2O)\n
!> hruyro(16,:) actual amount of evaporation (from soil) that occurs during year
!> in HRU (mm H2O)\n
!> hruyro(17,:) amount of nitrogen applied in continuous fertilizer operation
!> during year in HRU (kg N/ha)\n
!> hruyro(18,:) amount of phosphorus applied in continuous fertilizer operation
!> during year in HRU (kg P/ha)\n
!> hruyro(23,:) amount of water removed from shallow aquifer in HRU for
!> irrigation during year (mm H2O)\n
!> hruyro(24,:) amount of water removed from deep aquifer in HRU for irrigation
!> during year (mm H2O)\n
!> hruyro(25,:) potential evapotranspiration in HRU during year (mm H2O)\n
!> hruyro(26,:) annual amount of N (organic & mineral) applied in HRU during
!> grazing (kg N/ha)\n
!> hruyro(27,:) annual amount of P (organic & mineral) applied in HRU during
!> grazing (kg P/ha)\n
!> hruyro(28,:) annual amount of N (organic & mineral) auto-applied in HRU
!> (kg N/ha)\n
!> hruyro(29,:) annual amount of P (organic & mineral) auto-applied in HRU
!> (kg P/ha)\n
!> hruyro(31,:) water stress days in HRU during year (stress days)\n
!> hruyro(32,:) temperature stress days in HRU during year (stress days)\n
!> hruyro(33,:) nitrogen stress days in HRU during year (stress days)\n
!> hruyro(34,:) phosphorus stress days in HRU during year (stress days)\n
!> hruyro(35,:) organic nitrogen in surface runoff in HRU during year
!> (kg N/ha)\n
!> hruyro(36,:) organic phosphorus in surface runoff in HRU during year
!> (kg P/ha)\n
!> hruyro(37,:) nitrate in surface runoff in HRU during year (kg N/ha)\n
!> hruyro(38,:) nitrate in lateral flow in HRU during year (kg N/ha)\n
!> hruyro(39,:) soluble phosphorus in surface runoff in HRU during year
!> (kg P/ha)\n
!> hruyro(40,:) amount of nitrogen removed from soil by plant uptake in HRU
!> during year (kg N/ha)\n
!> hruyro(41,:) nitrate percolating past bottom of soil profile in HRU during
!> year (kg N/ha)\n
!> hruyro(42,:) amount of phosphorus removed from soil by plant uptake in HRU
!> during year (kg P/ha)\n
!> hruyro(43,:) amount of phosphorus moving from labile mineral to active
!> mineral pool in HRU during year (kg P/ha)\n
!> hruyro(44,:) amount of phosphorus moving from active mineral to stable
!> mineral pool in HRU during year (kg P/ha)\n
!> hruyro(45,:) amount of nitrogen applied to HRU in fertilizer and grazing
!> operations during year (kg N/ha)\n
!> hruyro(46,:) amount of phosphorus applied to HRU in fertilizer and grazing
!> operations during year (kg P/ha)\n
!> hruyro(47,:) amount of nitrogen added to soil by fixation in HRU during year
!> (kg N/ha)\n
!> hruyro(48,:) amount of nitrogen lost by denitrification in HRU during year
!> (kg N/ha)\n
!> hruyro(49,:) amount of nitrogen moving from active organic to nitrate pool in
!> HRU during year (kg N/ha)\n
!> hruyro(50,:) amount of nitrogen moving from active organic to stable organic
!> pool in HRU during year (kg N/ha)\n
!> hruyro(51,:) amount of phosphorus moving from organic to labile mineral pool
!> in HRU during year (kg P/ha)\n
!> hruyro(52,:) amount of nitrogen moving from fresh organic to nitrate and
!> active organic pools in HRU during year (kg N/ha)\n
!> hruyro(53,:) amount of phosphorus moving from fresh organic to the labile
!> mineral and organic pools in HRU during year (kg P/ha)\n
!> hruyro(54,:) amount of nitrogen added to soil in rain during year (kg N/ha)\n
!> hruyro(61,:) daily soil loss predicted with USLE equation (metric tons/ha)\n
!> hruyro(63,:) less persistent bacteria transported to main channel from HRU
!> during year (# bacteria/ha)\n
!> hruyro(64,:) persistent bacteria transported to main channel from HRU during
!> year (# bacteria/ha)\n
!> hruyro(65,:) nitrate loading from groundwater in HRU to main channel during
!> year (kg N/ha)\n
!> hruyro(66,:) soluble P loading from groundwater in HRU to main channel during
!> year (kg P/ha)\n
!> hruyro(67,:) loading of mineral P attached to sediment in HRU to main channel
!> during year (kg P/ha)
   real*8, dimension (:,:), allocatable :: hruyro
!> reach average annual output array (varies)\n
!> rchaao(1,:) flow into reach during simulation (m^3/s)\n
!> rchaao(2,:) flow out of reach during simulation (m^3/s)\n
!> rchaao(3,:) sediment transported into reach during simulation
!> (metric tons)\n
!> rchaao(4,:) sediment transported out of reach during simulation
!> (metric tons)\n
!> rchaao(5,:) sediment concentration in outflow during simulation (mg/L)\n
!> rchaao(6,:) organic N transported into reach during simulation (kg N)\n
!> rchaao(7,:) organic N transported out of reach during simulation (kg N)\n
!> rchaao(8,:) organic P transported into reach during simulation (kg P)\n
!> rchaao(9,:) organic P transported out of reach during simulation (kg P)\n
!> rchaao(10,:) evaporation from reach during simulation (m^3/s)\n
!> rchaao(11,:) transmission losses from reach during simulation (m^3/s)\n
!> rchaao(12,:) conservative metal #1 transported out of reach during
!> simulation (kg)\n
!> rchaao(13,:) conservative metal #2 transported out of reach during
!> simulation (kg)\n
!> rchaao(14,:) conservative metal #3 transported out of reach during
!> simulation (kg)\n
!> rchaao(15,:) nitrate transported into reach during simulation (kg N)\n
!> rchaao(16,:) nitrate transported out of reach during simulation (kg N)\n
!> rchaao(17,:) soluble P transported into reach during simulation (kg P)\n
!> rchaao(18,:) soluble P transported out of reach during simulation (kg P)\n
!> rchaao(19,:) soluble pesticide transported into reach during simulation\n
!> rchaao(20,:) soluble pesticide transported out of reach during simulation\n
!> rchaao(21,:) sorbed pesticide transported into reach during simulation\n
!> rchaao(22,:) sorbed pesticide transported out of reach during simulation\n
!> rchaao(23,:) amount of pesticide lost through reactions in reach during
!> simulation\n
!> rchaao(24,:) amount of pesticide lost through volatilization from reach
!> during simulation\n
!> rchaao(25,:) amount of pesticide settling out of reach to bed sediment during
!> simulation\n
!> rchaao(26,:) amount of pesticide resuspended from bed sediment to reach
!> during simulation\n
!> rchaao(27,:) amount of pesticide diffusing from reach to bed sediment during
!> simulation\n
!> rchaao(28,:) amount of pesticide in sediment layer lost through reactions
!> during simulation\n
!> rchaao(29,:) amount of pesticide in sediment layer lost through burial during
!> simulation\n
!> rchaao(30,:) chlorophyll-a transported into reach during simulation
!> (kg chla)\n
!> rchaao(31,:) chlorophyll-a transported out of reach during simulation
!> (kg chla)\n
!> rchaao(32,:) ammonia transported into reach during simuation (kg N)\n
!> rchaao(33,:) ammonia transported out of reach during simuation (kg N)\n
!> rchaao(34,:) nitrite transported into reach during simuation (kg N)\n
!> rchaao(35,:) nitrite transported out of reach during simuation (kg N)\n
!> rchaao(36,:) CBOD transported into reach during simulation (kg O2)\n
!> rchaao(37,:) CBOD transported out of reach during simuation (kg O2)\n
!> rchaao(38,:) dissolved oxygen transported into reach during simuation
!> (kg O2)\n
!> rchaao(39,:) dissolved oxygen transported out of reach during simulation
!> (kg O2)\n
!> rchaao(40,:) persistent bacteria transported out of reach during simulation
!> (kg bact)\n
!> rchaao(41,:) less persistent bacteria transported out of reach during
!> simulation (kg bact)\n
!> rchaao(43,:) Total N (org N + no3 + no2 + nh4 outs) (kg)\n
!> rchaao(44,:) Total P (org P + sol p outs) (kg)
   real*8, dimension (:,:), allocatable :: rchaao
!> subbasin monthly output array (varies)\n
!> submono(1,:) precipitation in subbasin for month (mm H20)\n
!> submono(2,:) snow melt in subbasin for month (mm H20)\n
!> submono(3,:) surface runoff loading in subbasin for month (mm H20)\n
!> submono(4,:) water yield from subbasin for month (mm H20)\n
!> submono(5,:) potential evapotranspiration in subbasin for month (mm H20)\n
!> submono(6,:) actual evapotranspiration in subbasin for month (mm H20)\n
!> submono(7,:) sediment yield from subbasin for month (metric tons/ha)\n
!> submono(8,:) organic N loading from subbasin for month (kg N/ha)\n
!> submono(9,:) organic P loading from subbasin for month (kg P/ha)\n
!> submono(10,:) NO3 loading from surface runoff in subbasin for month
!> (kg N/ha)\n
!> submono(11,:) soluble P loading from subbasin for month (kg P/ha)\n
!> submono(12,:) groundwater loading from subbasin for month (mm H20)\n
!> submono(13,:) percolation out of soil profile in subbasin for month
!> (mm H20)\n
!> submono(14,:) loading to reach of mineral P attached to sediment from
!> subbasin for month (kg P/ha)
   real*8, dimension (:,:), allocatable :: submono
!> subbasin annual output array (varies)\n
!> subyro(1,:) precipitation in subbasin for year (mm H2O)\n
!> subyro(2,:) snow melt in subbasin for year (mm H2O)\n
!> subyro(3,:) surface runoff loading in subbasin for year (mm H2O)\n
!> subyro(4,:) water yield from subbasin for year (mm H2O)\n
!> subyro(5,:) potential evapotranspiration in subbasin for year (mm H2O)\n
!> subyro(6,:) actual evapotranspiration in subbasin for year (mm H2O)\n
!> subyro(7,:) sediment yield from subbasin for year (metric tons/ha)\n
!> subyro(8,:) organic N loading from subbasin for year (kg N/ha)\n
!> subyro(9,:) organic P loading from subbasin for year (kg P/ha)\n
!> subyro(10,:) NO3 loading from surface runoff in subbasin for year (kg N/ha)\n
!> subyro(11,:) soluble P loading from subbasin for year (kg P/ha)\n
!> subyro(12,:) groundwater loading from subbasin for year (mm H2O)\n
!> subyro(13,:) percolation out of soil profile in subbasin for year (mm H2O)\n
!> subyro(14,:) loading to reach of mineral P attached to sediment from subbasin
!> for year (kg P/ha)
   real*8, dimension (:,:), allocatable :: subyro
!> HRU average annual output array (varies)\n
!> hruaao(1,:) precipitation in HRU during simulation (mm H2O)\n
!> hruaao(2,:) amount of precipitation falling as freezing rain/snow in HRU
!> during simulation (mm H2O)\n
!> hruaao(3,:) amount of snow melt in HRU during simulation (mm H2O)\n
!> hruaao(4,:) amount of surface runoff to main channel from HRU during
!> simulation (ignores impact of transmission losses) (mm H2O)\n
!> hruaao(5,:) amount of lateral flow contribution to main channel from HRU
!> during simulation (mm H2O)\n
!> hruaao(6,:) amount of groundwater flow contribution to main channel from HRU
!> during simulation (mm H2O)\n
!> hruaao(7,:) amount of water moving from shallow aquifer to plants or soil
!> profile in HRU during simulation (mm H2O)\n
!> hruaao(8,:) amount of water recharging deep aquifer in HRU during simulation
!> (mm H2O)\n
!> hruaao(9,:) total amount of water entering both aquifers from HRU during
!> simulation (mm H2O)\n
!> hruaao(10,:) water yield (total amount of water entering main channel) from
!> HRU during simulation (mm H2O)\n
!> hruaao(11,:) amount of water percolating out of the soil profile and into the
!> vadose zone in HRU during simulation (mm H2O)\n
!> hruaao(12,:) actual evapotranspiration in HRU during simulation\n
!> hruaao(13,:) amount of transmission losses from tributary channels in HRU for
!> simulation (mm H2O)\n
!> hruaao(14,:) sediment yield from HRU for simulation (metric tons/ha)\n
!> hruaao(17,:) amount of nitrogen applied in continuous fertilizer operation
!> in HRU for simulation (kg N/ha)\n
!> hruaao(18,:) amount of phosphorus applied in continuous fertilizer operation
!> in HRU for simulation (kg P/ha)\n
!> hruaao(23,:) amount of water removed from shallow aquifer in HRU for
!> irrigation during simulation (mm H2O)\n
!> hruaao(24,:) amount of water removed from deep aquifer in HRU for irrigation
!> during simulation (mm H2O)\n
!> hruaao(25,:) potential evapotranspiration in HRU during simulation (mm H2O)\n
!> hruaao(26,:) annual amount of N (organic & mineral) applied in HRU during
!> grazing (kg N/ha)\n
!> hruaao(27,:) annual amount of P (organic & mineral) applied in HRU during
!> grazing (kg P/ha)\n
!> hruaao(28,:) average annual amount of N (organic & mineral) auto-applied in
!> HRU (kg N/ha)\n
!> hruaao(29,:) average annual amount of P (organic & mineral) auto-applied in
!> HRU (kg P/ha)\n
!> hruaao(31,:) water stress days in HRU during simulation (stress days)\n
!> hruaao(32,:) temperature stress days in HRU during simulation (stress days)\n
!> hruaao(33,:) nitrogen stress days in HRU during simulation (stress days)\n
!> hruaao(34,:) phosphorus stress days in HRU during simulation (stress days)\n
!> hruaao(35,:) organic nitrogen in surface runoff in HRU during simulation
!> (kg N/ha)\n
!> hruaao(36,:) organic phosphorus in surface runoff in HRU during simulation
!> (kg P/ha)\n
!> hruaao(37,:) nitrate in surface runoff in HRU during simulation (kg N/ha)\n
!> hruaao(38,:) nitrate in lateral flow in HRU during simulation (kg N/ha)\n
!> hruaao(39,:) soluble phosphorus in surface runoff in HRU during simulation
!> (kg P/ha)\n
!> hruaao(40,:) amount of nitrogen removed from soil by plant uptake in HRU
!> during simulation (kg N/ha)\n
!> hruaao(41,:) nitrate percolating past bottom of soil profile in HRU during
!> simulation (kg N/ha)\n
!> hruaao(42,:) amount of phosphorus removed from soil by plant uptake in HRU
!> during simulation (kg P/ha)\n
!> hruaao(43,:) amount of phosphorus moving from labile mineral to active
!> mineral pool in HRU during simulation (kg P/ha)\n
!> hruaao(44,:) amount of phosphorus moving from active mineral to stable
!> mineral pool in HRU during simulation (kg P/ha)\n
!> hruaao(45,:) amount of nitrogen applied to HRU in fertilizer and grazing
!> operations during simulation (kg N/ha)\n
!> hruaao(46,:) amount of phosphorus applied to HRU in fertilizer and grazing
!> operations during simulation (kg P/ha)\n
!> hruaao(47,:) amount of nitrogen added to soil by fixation in HRU during
!> simulation (kg N/ha)\n
!> hruaao(48,:) amount of nitrogen lost by denitrification in HRU during
!> simulation (kg N/ha)\n
!> hruaao(49,:) amount of nitrogen moving from active organic to nitrate pool
!> in HRU during simulation (kg N/ha)\n
!> hruaao(50,:) amount of nitrogen moving from active organic to stable organic
!> pool in HRU during simulation (kg N/ha)\n
!> hruaao(51,:) amount of phosphorus moving from organic to labile mineral pool
!> in HRU during simulation (kg P/ha)\n
!> hruaao(52,:) amount of nitrogen moving from fresh organic to nitrate and
!> active organic pools in HRU during simulation (kg N/ha)\n
!> hruaao(53,:) amount of phosphorus moving from fresh organic to the labile mineral and organic pools in HRU during simulation (kg P/ha)\n
!> hruaao(54,:) amount of nitrogen added to soil in rain during simulation (kg N/ha)\n
!> hruaao(61,:) daily soil loss predicted with USLE equation (metric tons/ha)\n
!> hruaao(63,:) less persistent bacteria transported to main channel from HRU
!> during simulation (# bacteria/ha)\n
!> hruaao(64,:) persistent bacteria transported to main channel from HRU during
!> simulation (# bacteria/ha)\n
!> hruaao(65,:) nitrate loading from groundwater in HRU to main channel during
!> simulation (kg N/ha)\n
!> hruaao(66,:) soluble P loading from groundwater in HRU to main channel during
!> simulation (kg P/ha)\n
!> hruaao(67,:) loading of mineral P attached to sediment in HRU to main channel
!> during simulation (kg P/ha)
   real*8, dimension (:,:), allocatable :: hruaao
!> subbasin average annual output array (varies)
   real*8, dimension (:,:), allocatable :: subaao
!> reservoir monthly output array (varies)\n
!> resoutm(1,:) flow into reservoir during month (m^3/s)\n
!> resoutm(2,:) flow out of reservoir during month (m^3/s)\n
!> resoutm(3,:) sediment entering reservoir during month (metric tons)\n
!> resoutm(4,:) sediment leaving reservoir during month (metric tons)\n
!> resoutm(5,:) sediment concentration in reservoir during month (mg/L)\n
!> resoutm(6,:) pesticide entering reservoir during month (mg pst)\n
!> resoutm(7,:) pesticide lost from reservoir through reactions during month
!> (mg pst)\n
!> resoutm(8,:) pesticide lost from reservoir through volatilization during
!> month (mg pst)\n
!> resoutm(9,:) pesticide moving from water to sediment through settling during
!> month (mg pst)\n
!> resoutm(10,:) pesticide moving from sediment to water through resuspension
!> during month (mg pst)\n
!> resoutm(11,:) pesticide moving from water to sediment through diffusion
!> during month (mg pst)\n
!> resoutm(12,:) pesticide lost from reservoir sediment layer through reactions
!> during month (mg pst)\n
!> resoutm(13,:) pesticide lost from reservoir sediment layer through burial
!> during month (mg pst)\n
!> resoutm(14,:) pesticide transported out of reservoir during month (mg pst)\n
!> resoutm(15,:) pesticide concentration in reservoir water during month
!> (mg pst/m^3)\n
!> resoutm(16,:) pesticide concentration in reservoir sediment layer during
!> month (mg pst/m^3)\n
!> resoutm(17,:) evaporation from reservoir during month (m^3 H2O)\n
!> resoutm(18,:) seepage from reservoir during month (m^3 H2O)\n
!> resoutm(19,:) precipitation on reservoir during month (m^3 H2O)\n
!> resoutm(22,:) organic N entering reservoir during month (kg N)\n
!> resoutm(23,:) organic N leaving reservoir during month (kg N)\n
!> resoutm(24,:) organic P entering reservoir during month (kg P)\n
!> resoutm(25,:) organic P leaving reservoir during month (kg P)\n
!> resoutm(26,:) nitrate entering reservoir during month (kg N)\n
!> resoutm(27,:) nitrate leaving reservoir during month (kg N)\n
!> resoutm(28,:) nitrite entering reservoir during month (kg N)\n
!> resoutm(29,:) nitrite leaving reservoir during month (kg N)\n
!> resoutm(30,:) ammonia entering reservoir during month (kg N)\n
!> resoutm(31,:) ammonia leaving reservoir during month (kg N)\n
!> resoutm(32,:) mineral P entering reservoir during month (kg P)\n
!> resoutm(33,:) mineral P leaving reservoir during month (kg P)\n
!> resoutm(34,:) chlorophyll-a entering reservoir during month (kg chla)\n
!> resoutm(35,:) chlorophyll-a leaving reservoir during month (kg chla)\n
!> resoutm(36,:) organic P concentration in reservoir water during month
!> (mg P/L)\n
!> resoutm(37,:) mineral P concentration in reservoir water during month
!> (mg P/L)\n
!> resoutm(38,:) organic N concentration in reservoir water during month
!> (mg N/L)\n
!> resoutm(39,:) nitrate concentration in reservoir water during month
!> (mg N/L)\n
!> resoutm(40,:) nitrite concentration in reservoir water during month
!> (mg N/L)\n
!> resoutm(41,:) ammonia concentration in reservoir water during month (mg N/L)
   real*8, dimension (:,:), allocatable :: resoutm
!> reservoir annual output array (varies)\n
!> resouty(1,:) flow into reservoir during year (m^3/s)\n
!> resouty(2,:) flow out of reservoir during year (m^3/s)\n
!> resouty(3,:) sediment entering reservoir during year (metric tons)\n
!> resouty(4,:) sediment leaving reservoir during year (metric tons)\n
!> resouty(5,:) sediment concentration in reservoir during year (mg/L)\n
!> resouty(6,:) pesticide entering reservoir during year (mg pst)\n
!> resouty(7,:) pesticide lost from reservoir through reactions during year
!> (mg pst)\n
!> resouty(8,:) pesticide lost from reservoir through volatilization during year
!> (mg pst)\n
!> resouty(9,:) pesticide moving from water to sediment through settling during
!> year (mg pst)\n
!> resouty(10,:) pesticide moving from sediment to water through resuspension
!> during year (mg pst)\n
!> resouty(11,:) pesticide moving from water to sediment through diffusion
!> during year (mg pst)\n
!> resouty(12,:) pesticide lost from reservoir sediment layer through reactions
!> during year (mg pst)\n
!> resouty(13,:) pesticide lost from reservoir sediment layer through burial
!>during year (mg pst)\n
!> resouty(14,:) pesticide transported out of reservoir during year (mg pst)\n
!> resouty(15,:) pesticide concentration in reservoir water during year
!> (mg pst/m^3)\n
!> resouty(16,:) pesticide concentration in reservoir sediment layer during year
!> (mg pst/m^3)\n
!> resouty(17,:) evaporation from reservoir during year (m^3 H2O)\n
!> resouty(18,:) seepage from reservoir during year (m^3 H2O)\n
!> resouty(19,:) precipitation on reservoir during year (m^3 H2O)\n
!> resouty(22,:) organic N entering reservoir during year (kg N)\n
!> resouty(23,:) organic N leaving reservoir during year (kg N)\n
!> resouty(24,:) organic P entering reservoir during year (kg P)\n
!> resouty(25,:) organic P leaving reservoir during year (kg P)\n
!> resouty(26,:) nitrate entering reservoir during year (kg N)\n
!> resouty(27,:) nitrate leaving reservoir during year (kg N)\n
!> resouty(28,:) nitrite entering reservoir during year (kg N)\n
!> resouty(29,:) nitrite leaving reservoir during year (kg N)\n
!> resouty(30,:) ammonia entering reservoir during year (kg N)\n
!> resouty(31,:) ammonia leaving reservoir during year (kg N)\n
!> resouty(32,:) mineral P entering reservoir during year (kg P)\n
!> resouty(33,:) mineral P leaving reservoir during year (kg P)\n
!> resouty(34,:) chlorophyll-a entering reservoir during year (kg chla)\n
!> resouty(35,:) chlorophyll-a leaving reservoir during year (kg chla)\n
!> resouty(36,:) organic P concentration in reservoir water during year
!> (mg P/L)\n
!> resouty(37,:) mineral P concentration in reservoir water during year
!> (mg P/L)\n
!> resouty(38,:) organic N concentration in reservoir water during year
!> (mg N/L)\n
!> resouty(39,:) nitrate concentration in reservoir water during year (mg N/L)\n
!> resouty(40,:) nitrite concentration in reservoir water during year (mg N/L)\n
!> resouty(41,:) ammonia concentration in reservoir water during year (mg N/L)
   real*8, dimension (:,:), allocatable :: resouty
!> reservoir average annual output array (varies)\n
!> resouta(3,:) sediment entering reservoir during simulation (metric tons)\n
!> resouta(4,:) sediment leaving reservoir during simulation (metric tons)\n
!> resouta(17,:) evaporation from reservoir during simulation (m^3 H2O)\n
!> resouta(18,:) seepage from reservoir during simulation (m^3 H2O)\n
!> resouta(19,:) precipitation on reservoir during simulation (m^3 H2O)\n
!> resouta(20,:) water entering reservoir during simulation (m^3 H2O)\n
!> resouta(21,:) water leaving reservoir during simulation (m^3 H2O)
   real*8, dimension (:,:), allocatable :: resouta
!> wshd_aamon(:,1) average annual precipitation in watershed falling during
!> month (mm H2O)\n
!> wshd_aamon(:,2) average annual freezing rain in watershed falling during
!> month (mm H2O)\n
!> wshd_aamon(:,3) average annual surface runoff in watershed during month
!> (mm H2O)\n
!> wshd_aamon(:,4) average annual lateral flow in watershed during month
!> (mm H2O)\n
!> wshd_aamon(:,5) average annual water yield in watershed during month
!> (mm H2O)\n
!> wshd_aamon(:,6) average annual actual evapotranspiration in watershed during
!> month (mm H2O)\n
!> wshd_aamon(:,7) average annual sediment yield in watershed during month
!> (metric tons)\n
!> wshd_aamon(:,8) average annual potential evapotranspiration in watershed
!> during month (mm H2O)
   real*8, dimension (12,8) :: wshd_aamon
!> HRU monthly output data array for impoundments (varies)\n
!> wtrmon(1,:) evaporation from ponds in HRU for month (mm H2O)\n
!> wtrmon(2,:) seepage from ponds in HRU for month (mm H2O)\n
!> wtrmon(3,:) precipitation on ponds in HRU for month (mm H2O)\n
!> wtrmon(4,:) amount of water entering ponds in HRU for month (mm H2O)\n
!> wtrmon(5,:) sediment entering ponds in HRU for month (metric tons/ha)\n
!> wtrmon(6,:) amount of water leaving ponds in HRU for month (mm H2O)\n
!> wtrmon(7,:) sediment leaving ponds in HRU for month (metric tons/ha)\n
!> wtrmon(8,:) precipitation on wetlands in HRU for month (mm H2O)\n
!> wtrmon(9,:) volume of water entering wetlands from HRU for month (mm H2O)\n
!> wtrmon(10,:) sediment loading to wetlands for month from HRU
!> (metric tons/ha)\n
!> wtrmon(11,:) evaporation from wetlands in HRU for month (mm H2O)\n
!> wtrmon(12,:) seeepage from wetlands in HRU for month (mm H2O)\n
!> wtrmon(13,:) volume of water leaving wetlands in HRU for month (mm H2O)\n
!> wtrmon(14,:) sediment loading from wetlands in HRU to main channel during
!> month (metric tons/ha)\n
!> wtrmon(15,:) precipitation on potholes in HRU for month (mm H2O)\n
!> wtrmon(16,:) evaporation from potholes in HRU for month (mm H2O)\n
!> wtrmon(17,:) seepage from potholes in HRU for month (mm H2O)\n
!> wtrmon(18,:) water leaving potholes in HRU for month (mm H2O)\n
!> wtrmon(19,:) water entering potholes in HRU for month (mm H2O)\n
!> wtrmon(20,:) sediment entering potholes in HRU for month (metric tons/ha)\n
!> wtrmon(21,:) sediment leaving potholes in HRU for month (metric tons/ha)
   real*8, dimension (:,:), allocatable :: wtrmon
!> HRU impoundment annual output array (varies)\n
!> wtryr(1,:) evaporation from ponds in HRU for year (mm H20)\n
!> wtryr(2,:) seepage from ponds in HRU for year (mm H20)\n
!> wtryr(3,:) precipitation on ponds in HRU for year (mm H20)\n
!> wtryr(4,:) amount of water entering ponds in HRU for year (mm H20)\n
!> wtryr(5,:)  sediment entering ponds in HRU for year (metric tons/ha)\n
!> wtryr(6,:) amount of water leaving ponds in HRU for year (mm H20)\n
!> wtryr(7,:)  sediment leaving ponds in HRU for year (metric tons/ha)\n
!> wtryr(8,:) precipitation on wetlands in HRU for year (mm H20)\n
!> wtryr(9,:) volume of water entering wetlands from HRU for year (mm H20)\n
!> wtryr(10,:) sediment loading to wetlands for year from HRU (metric tons/ha)\n
!> wtryr(11,:) evaporation from wetlands in HRU for year (mm H20)\n
!> wtryr(12,:) seeepage from wetlands in HRU for year (mm H20)\n
!> wtryr(13,:) volume of water leaving wetlands in HRU for year (mm H20)\n
!> wtryr(14,:) sediment loading from wetlands in HRU to main channel during year
!> (metric tons/ha)\n
!> wtryr(15,:) precipitation on potholes in HRU during year (mm H20)\n
!> wtryr(16,:) evaporation from potholes in HRU during year (mm H20)\n
!> wtryr(17,:) seepage from potholes in HRU during year (mm H20)\n
!> wtryr(18,:) water leaving potholes in HRU during year (mm H20)\n
!> wtryr(19,:) water entering potholes in HRU during year (mm H20)\n
!> wtryr(20,:) sediment entering potholes in HRU during year (metric tons/ha)\n
!> wtryr(21,:) sediment leaving potholes in HRU during year (metric tons/ha)
   real*8, dimension (:,:), allocatable :: wtryr
!> HRU impoundment average annual output array (varies)
   real*8, dimension (:,:), allocatable :: wtraa
!> max melt rate for snow during year (June 21) for subbasin(:) where deg C
!> refers to the air temperature. SUB_SMFMX and SMFMN allow the rate of snow
!> melt to vary through the year. These parameters are accounting for the impact
!> of soil temperature on snow melt (range: -5.0/5.0) (mm/deg C/day)
   real*8, dimension (:,:), allocatable :: sub_smfmx
!> min melt rate for snow during year (Dec 21) for subbasin(:) (range: -5.0/5.0)
!> where deg C refers to the air temperature (mm/deg C/day)
   real*8, dimension (:,:), allocatable :: sub_smfmn
!> hrupstd(1,:,:) amount of pesticide type in surface runoff contribution to
!> stream from HRU on day (in solution) (mg pst)\n
!> hrupstd(2,:,:) amount of pesticide type in surface runoff contribution to
!> stream from HRU on day (sorbed to sediment) (mg pst)\n
!> hrupstd(3,:,:) total pesticide loading to stream in surface runoff from HRU
!> (mg pst/ha)\n
!> hrupstd(4,:,:) amount of pesticide type in lateral flow contribution to
!> stream from HRU on day (in solution) (mg pst)
   real*8, dimension (:,:,:), allocatable :: hrupstd
!> hrupstm(:,:,:)HRU monthly pesticide output array (varies)\n
!> hrupstm(1,:,:) amount of pesticide type in surface runoff contribution to
!> stream from HRU during month (in solution) (mg pst)\n
!> hrupstm(2,:,:) amount of pesticide type in surface runoff contribution to
!> stream from HRU during month (sorbed to sediment) (mg pst)\n
!> hrupstm(3,:,:) total pesticide loading to stream in surface runoff from HRU
!> during month (mg pst)
   real*8, dimension (:,:,:), allocatable :: hrupstm
!> HRU average annual pesticide output array (varies)
   real*8, dimension (:,:,:), allocatable :: hrupsta
!> hrupsty(:,:,:) HRU annual pesticide output array (varies)\n
!> hrupsty(1,:,:) amount of pesticide type in surface runoff contribution to
!> stream from HRU during year (in solution) (mg pst)\n
!> hrupsty(2,:,:) amount of pesticide type in surface runoff contribution to
!> stream from HRU during year (sorbed to sediment) (mg pst)
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
!> flow into reach on current day (m^3 H2O)
   real*8, dimension (:), allocatable :: flwin
!> flow out of reach on current day (m^3 H2O)
   real*8, dimension (:), allocatable :: flwout
   real*8, dimension (:), allocatable :: bankst,ch_wi
!> channel organic n concentration (ppm)
   real*8, dimension (:), allocatable :: ch_onco
!> channel organic p concentration (ppm)
   real*8, dimension (:), allocatable :: ch_opco
   real*8, dimension (:), allocatable :: ch_orgn, ch_orgp
!> amount of pesticide drifting onto main channel in subbasin (kg)
   real*8, dimension (:), allocatable :: drift
!> dissolved oxygen concentration in reach (mg O2/L)
   real*8, dimension (:), allocatable :: rch_dox
!> persistent bacteria in reach/outflow at end of day (# cfu/100ml)
   real*8, dimension (:), allocatable :: rch_bactp
!> alpha factor for bank storage recession curve (days)
   real*8, dimension (:), allocatable :: alpha_bnk
!> \f$\exp(-alpha_bnk)\f$ (none)
   real*8, dimension (:), allocatable :: alpha_bnke
!> water stored in reach (m^3 H2O)
   real*8, dimension (:), allocatable :: rchstor
!> amount of sediment stored in reach (metric tons)
   real*8, dimension (:), allocatable :: sedst
!> algal biomass concentration in reach (mg alg/L)
   real*8, dimension (:), allocatable :: algae
!> dissolved phosphorus concentration in reach (mg P/L)
   real*8, dimension (:), allocatable :: disolvp
!> chlorophyll-a concentration in reach (mg chl-a/L)
   real*8, dimension (:), allocatable :: chlora
!> organic nitrogen concentration in reach (mg N/L)
   real*8, dimension (:), allocatable :: organicn
!> organic phosphorus concentration in reach (mg P/L)
   real*8, dimension (:), allocatable :: organicp
!> initial length of main channel (km)
   real*8, dimension (:), allocatable :: ch_li
!> initial slope of main channel (m/m)
   real*8, dimension (:), allocatable :: ch_si
!> nitrate concentration in reach (mg N/L)
   real*8, dimension (:), allocatable :: nitraten
!> nitrite concentration in reach (mg N/L)
   real*8, dimension (:), allocatable :: nitriten

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
!> ch_l(1,:) longest tributary channel length in subbasin (km)\n
!> ch_l(2,:) length of main channel (km)
   real*8, dimension (:,:), allocatable :: ch_l
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
!> initial pesticide concentration in reach (mg/(m^3))
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
!> carbonaceous biochemical oxygen demand in reach (mg O2/L)
   real*8, dimension (:), allocatable :: rch_cbod
!> less persistent bacteria in reach/outflow at end of day (# cfu/100ml)
   real*8, dimension (:), allocatable :: rch_bactlp
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
!> ammonia concentration in reach (mg N/L)
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
!> amount of water in soil profile in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_sw
!> amount of phosphorus stored in all mineral pools sorbed to sediment (kg P/ha)
   real*8, dimension (:), allocatable :: sub_minp
   real*8, dimension (:), allocatable :: wcklsp
!> nitrate loading in groundwater from subbasin (kg N/ha)
   real*8, dimension (:), allocatable :: sub_gwno3
!> amount of water in soil at field capacity in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_sumfc
   real*8, dimension (:), allocatable :: sub_gwsolp
   real*8, dimension (:), allocatable :: co2 !< CO2 concentration (ppmv)
!> area of subbasin in square kilometers (km^2)
   real*8, dimension (:), allocatable :: sub_km
!> latitude of weather station used to compile data (degrees)
   real*8, dimension (:), allocatable :: wlat
!> time of concentration for subbasin (hour)
   real*8, dimension (:), allocatable :: sub_tc
!> potential evapotranspiration for day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_pet
!> elevation of weather station used to compile weather generator data (m)
   real*8, dimension (:), allocatable :: welev
!> average bulk density in subbasin for top 10 mm of first soil layer (Mg/m^3)
   real*8, dimension (:), allocatable :: sub_bd
!> amount of nitrogen stored in all organic pools (kg N/ha)
   real*8, dimension (:), allocatable :: sub_orgn
!> amount of phosphorus stored in all organic pools (kg P/ha)
   real*8, dimension (:), allocatable :: sub_orgp
!> amount of active mineral P attached to sediment removed in surface runoff on
!> day in subbasin (kg P/ha)
   real*8, dimension (:), allocatable :: sub_sedpa
!> amount of stable mineral P attached to sediment removed in surface runoff on
!> day in subbasin (kg P/ha)
   real*8, dimension (:), allocatable :: sub_sedps
   real*8, dimension (:), allocatable :: sub_wtmp
!> shortest daylength occurring during the year (hour)
   real*8, dimension (:), allocatable :: daylmn
!> amount of phosphorus stored in active mineral pools sorbed to sediment
!> (kg P/ha)
   real*8, dimension (:), allocatable :: sub_minpa
!> amount of phosphorus stored in stable mineral pools sorbed to sediment
!> (kg P/ha)
   real*8, dimension (:), allocatable :: sub_minps
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
!> effective precipitation (amount of water reaching soil surface) for the day
!> in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_precip
!> atmospheric deposition of ammonium values for entire watershed (mg/l)
   real*8, dimension (:), allocatable :: rammo_sub
!> atmospheric deposition of nitrate for entire watershed (mg/l)
   real*8, dimension (:), allocatable :: rcn_sub
   real*8, dimension (:), allocatable :: pcpdays
   real*8, dimension (:), allocatable :: atmo_day
!> amount of snow melt in subbasin on day (mm H2O)
   real*8, dimension (:), allocatable :: sub_snom
!> surface runoff that reaches main channel during day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_qd
!> sediment yield for the day in subbasin (metric tons)
   real*8, dimension (:), allocatable :: sub_sedy
!> transmission losses on day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_tran
!> NO3-N in surface runoff on day in subbasin (kg N/ha)
   real*8, dimension (:), allocatable :: sub_no3
!> NO3-N in lateral flow on day in subbasin (kg N/ha)
   real*8, dimension (:), allocatable :: sub_latno3
!> snowfall temperature for subbasin(:). Mean air temperature at which precip
!> is equally likely to be rain as snow/freezing rain (range: -5.0/5.0) (deg C)
   real*8, dimension (:,:), allocatable :: sub_sftmp
!> snow melt base temperature for subbasin(:) mean air temperature at which snow
!> melt will occur (range: -5.0/5.0) (deg C)
   real*8, dimension (:,:), allocatable :: sub_smtmp
!> snow pack temperature lag factor (0-1) (none)\n
!> 1 = no lag (snow pack temp=current day air temp) as the lag factor goes to
!> zero, the snow pack's temperature will be less influenced by the current
!> day's air temperature
   real*8, dimension (:,:), allocatable :: sub_timp
   real*8, dimension (:), allocatable :: sub_tileno3
!> actual evapotranspiration on day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_etday
!> soluble P in surface runoff on day in subbasin (kg P/ha)
   real*8, dimension (:), allocatable :: sub_solp
!> precipitation for day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_subp
!> average elevation of HRU (m)
   real*8, dimension (:), allocatable :: sub_elev
!> surface runoff generated on day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_surfq
!> water yield on day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_wyld
   real*8, dimension (:), allocatable :: qird
!> groundwater flow on day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_gwq
!> seepage from bottom of soil profile on day in subbasin (mm H2O)
   real*8, dimension (:), allocatable :: sub_sep
!> chlorophyll-a in water yield on day in subbasin (kg chl-a)
   real*8, dimension (:), allocatable :: sub_chl
!> carbonaceous biological oxygen demand on day for subbasin (kg O2)
   real*8, dimension (:), allocatable :: sub_cbod
!> dissolved oxygen loading on day for subbasin (kg O2)
   real*8, dimension (:), allocatable :: sub_dox
!> pesticide in solution in surface runoff on day in subbasin (mg pst)
   real*8, dimension (:), allocatable :: sub_solpst
!> organic N in surface runoff on day in subbasin (kg P/ha)
   real*8, dimension (:), allocatable :: sub_yorgn
!> organic P in surface runoff on day in subbasin (kg P/ha)
   real*8, dimension (:), allocatable :: sub_yorgp
!> pesticide sorbed to sediment in surface runoff on day in subbasin (mg pst)
   real*8, dimension (:), allocatable :: sub_sorpst
!> latitude of HRU/subbasin (degrees)
   real*8, dimension (:), allocatable :: sub_lat
!> less persistent bacteria in surface runoff for day in subbasin (# cfu/m^2)
   real*8, dimension (:), allocatable :: sub_bactlp
!> persistent bacteria in surface runoff for day in subbasin (# cfu/m^2)
   real*8, dimension (:), allocatable :: sub_bactp
   real*8, dimension (:), allocatable :: sub_latq, sub_gwq_d,sub_tileq
   real*8, dimension (:), allocatable :: sub_vaptile
   real*8, dimension (:), allocatable :: sub_dsan, sub_dsil, sub_dcla
   real*8, dimension (:), allocatable :: sub_dsag, sub_dlag

   real*8 :: vap_tile
   real*8, dimension (:), allocatable :: wnan
   real*8, dimension (:,:), allocatable :: sol_stpwt
!> amount of pesticide in soil layer in subbasin (kg/ha)
   real*8, dimension (:,:), allocatable :: sub_pst
!> water temperature for the time step in subbasin (deg C)
   real*8, dimension (:,:), allocatable :: sub_hhwtmp
   real*8, dimension (:,:), allocatable :: sub_hhqd
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
!> ch_k(1,:) effective hydraulic conductivity of tributary channel alluvium
!> (mm/hr)\n
!> ch_k(2,:) effective hydraulic conductivity of main channel alluvium (mm/hr)
   real*8, dimension (:,:), allocatable :: ch_k
!> elevation at the center of the band in subbasin (m)
   real*8, dimension (:,:), allocatable :: elevb
!> fraction of subbasin area within elevation band (the same fractions should be
!> listed for all HRUs within the subbasin) (none)
   real*8, dimension (:,:), allocatable :: elevb_fr
!> average wind speed for the month (m/s)
   real*8, dimension (:,:), allocatable :: wndav
!> ch_n(1,:) Manning's "n" value for the tributary channels (none)\n
!> ch_n(2,:) Manning's "n" value for the main channel (none)
   real*8, dimension (:,:), allocatable :: ch_n
!> ch_s(1,:) average slope of tributary channels (m/m)\n
!> ch_s(2,:) average slope of main channel (m/m)
   real*8, dimension (:,:), allocatable :: ch_s
!> ch_w(1,:) average width of tributary channels (m)\n
!> ch_w(2,:) average width of main channel (m)
   real*8, dimension (:,:), allocatable :: ch_w
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
!> GIS code printed to output files (output.sub, .rch) (none)
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
!> amount of nitrogen stored in the active organic (humic) nitrogen pool in soil
!> layer (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_aorgn
!> amount of nitrogen stored in the fresh organic (residue) pool in soil layer
!> (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_fon
!> average temperature of soil layer on previous day or\n
!> daily average temperature of soil layer (deg C)
   real*8, dimension (:,:), allocatable :: sol_tmp
!> available water capacity of soil layer (mm H20/mm soil)
   real*8, dimension (:,:), allocatable :: sol_awc
!> crack volume for soil layer (mm)
   real*8, dimension (:,:), allocatable :: volcr
!> percolation storage from soil layer on current day (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_prk
!> subbasin phosphorus percolation coefficient. Ratio of soluble phosphorus in
!> surface to soluble phosphorus in percolate
   real*8, dimension (:,:), allocatable :: pperco_sub
!> amount of phosphorus in the soil layer stored in the stable mineral
!> phosphorus pool (kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_stap
!> factor which converts kg/kg soil to kg/ha (none)
   real*8, dimension (:,:), allocatable :: conv_wt
!> amount of phosphorus stored in the active mineral phosphorus pool (kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_actp
!> soluble P concentration in top soil layer (mg P/kg soil) or\n
!> amount of inorganic phosphorus stored in solution in soil layer. NOTE UNIT
!> CHANGE! (kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_solp
!> maximum or potential crack volume (mm)
   real*8, dimension (:,:), allocatable :: crdep
!> amount of water available to plants in soil layer at field capacity (fc - wp
!> water) (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_fc
!> amount of water held in the soil layer at saturation (sat - wp water)
!> (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_ul
!> bulk density of the soil layer in HRU (Mg/m^3)
   real*8, dimension (:,:), allocatable :: sol_bd
!> depth to bottom of each soil profile layer in a given HRU (mm)
   real*8, dimension (:,:), allocatable :: sol_z
!> amount of water stored in the soil layer on any given day (less wilting point
!> water) (mm H2O)
   real*8, dimension (:,:), allocatable :: sol_st
!> water content of soil at -0.033 MPa (field capacity) (mm H2O/mm soil)
   real*8, dimension (:,:), allocatable :: sol_up
!> percent clay content in soil layer in HRU (UNIT CHANGE!) (% or none)
   real*8, dimension (:,:), allocatable :: sol_clay
!> beta coefficent to calculate hydraulic conductivity (none)
   real*8, dimension (:,:), allocatable :: sol_hk
!> lateral flow storage in soil layer on current day (mm H2O)
   real*8, dimension (:,:), allocatable :: flat
!> amount of nitrogen stored in the ammonium pool in soil layer (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_nh3
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
!> amount of phosphorus stored in the organic P pool in soil layer. NOTE UNIT
!> CHANGE! (mg P/kg soil or kg P/ha)
   real*8, dimension (:,:), allocatable :: sol_orgp
!> amount of organic matter in the soil layer classified as humic substances
!> (kg humus/ha)
   real*8, dimension (:,:), allocatable :: sol_hum
!> water content of soil at -1.5 MPa (wilting point) (mm H20)
   real*8, dimension (:,:), allocatable :: sol_wpmm
!> amount of nitrogen stored in the nitrate pool in the soil layer. This
!> variable is read in as a concentration and converted to kg/ha (this value is
!> read from the .sol file in units of mg/kg) (kg N/ha)
   real*8, dimension (:,:), allocatable :: sol_no3
!> percent organic carbon in soil layer (%)
   real*8, dimension (:,:), allocatable :: sol_cbn
!> saturated hydraulic conductivity of soil layer (mm/hour)
   real*8, dimension (:,:), allocatable :: sol_k
!> amount of organic matter in the soil layer classified as residue (kg/ha)
   real*8, dimension (:,:), allocatable :: sol_rsd
!> amount of phosphorus stored in the fresh organic (residue) pool in soil layer
!> (kg P/ha)
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
!> lateral saturated hydraulic conductivity for each profile layer in a give
!> HRU. For example (conk(2,1) is conductivity of layer from sol_z(1,1)
!> to sol_z(2,1) in HRU1 (mm/hr)
   real*8, dimension (:,:), allocatable :: conk
!> sol_pst(:,:,1) initial amount of pesticide in first layer read in from .chm
!> file (mg/kg)\n
!> sol_pst(:,:,:) amount of pesticide in soil layer. NOTE UNIT CHANGE!
!> (kg/ha)
   real*8, dimension (:,:,:), allocatable :: sol_pst
!> pesticide sorption coefficient, Kp; the ratio of the concentration in the
!> solid phase to the concentration in solution ((mg/kg)/(mg/L) or m^3/ton)
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
!> secchi-disk depth (m)
   real*8, dimension (:), allocatable :: res_seci
   real*8, dimension (:), allocatable :: res_chla
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
!> psetlr(1,:) phosphorus settling rate for mid-year period (read in as m/year
!> and converted to m/day) (m/day)\n
!> psetlr(2,:) phosphorus settling rate for remainder of year (read in as m/year
!> and converted to m/day) (m/day)
   real*8, dimension (:,:), allocatable :: psetlr
!> nsetlr(1,:) nitrogen settling rate for mid-year period (read in as m/year and
!> converted to m/day) (m/day)\n
!> nsetlr(2,:) nitrogen settling rate for remainder of year (read in as m/year
!> and converted to m/day) (m/day)
   real*8, dimension (:,:), allocatable :: nsetlr
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
!> iflodr(1,:) beginning month of non-flood season (needed if IRESCO=2) (none)\n
!> iflodr(2,:) ending month of non-flood season (needed if IRESCO=2) (none)
   integer, dimension (:,:), allocatable :: iflodr
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
!> depth of irrigation water applied to HRU (mm H2O)
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
!> date of street sweeping operation (julian date)
   integer, dimension (:), allocatable :: isweep
   integer, dimension (:), allocatable :: yr_skip
   integer, dimension (:), allocatable :: icrmx, nopmx
   integer, dimension (:,:), allocatable :: mgtop, idop
   integer, dimension (:,:), allocatable :: mgt1iop,mgt2iop,mgt3iop
   real*8, dimension (:,:), allocatable ::  mgt4op, mgt5op, mgt6op
   real*8, dimension (:,:), allocatable :: mgt7op, mgt8op, mgt9op
   real*8, dimension (:,:), allocatable :: mgt10iop, phu_op
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
!> intercepted photosynthetically active radiation ((kg/ha)/(MJ/m**2))
   real*8, dimension (:), allocatable :: bio_e
!> harvest index: crop yield/aboveground biomass ((kg/ha)/(kg/ha))
   real*8, dimension (:), allocatable :: hvsti
!> minimum temperature for plant growth (deg C)
   real*8, dimension (:), allocatable :: t_base
!> optimal temperature for plant growth (deg C)
   real*8, dimension (:), allocatable :: t_opt
   real*8, dimension (:), allocatable :: chtmx !< maximum canopy height (m)
!> natural log of USLE_C (the minimum value of the USLE C factor for the land
!> cover) (none)
   real*8, dimension (:), allocatable :: cvm
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
!> crop/landcover category (none):\n
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
!> fraction of fertilize/manure that is mineral nitrogen (NO3 + NH3)
!> (kg minN/kg fert)
   real*8, dimension (:), allocatable :: fminn
!> fraction of organic nitrogen in fertilizer/manure (kg orgN/kg fert)
   real*8, dimension (:), allocatable :: forgn
!> fraction of fertilizer/manure that is organic phosphorus (kg orgP/kg fert)
   real*8, dimension (:), allocatable :: forgp
!> fraction of bacteria in solution (the remaining fraction is sorbed to soil
!> particles) (none):\n
!> 1: all bacteria in solution\n
!> 0: all bacteria sorbed to soil particles
   real*8, dimension (:), allocatable :: bactkddb
!> concentration of less persistent bacteria in manure (fertilizer)
!> (cfu/g manure)
   real*8, dimension (:), allocatable :: bactlpdb
!> fraction of fertilizer that is mineral phosphorus in fertilizer/manure
!> (kg minP/kg fert)
   real*8, dimension (:), allocatable :: fminp
!> fraction of mineral N content that is NH3-N in fertilizer/manure
!> (kg NH3-N/kg minN)
   real*8, dimension (:), allocatable :: fnh3n
   character(len=8), dimension (200) :: fertnm !< name of fertilizer
!> curb length density in HRU (km/ha)
   real*8, dimension (:), allocatable :: curbden
!> maximum amount of solids allowed to build up on impervious surfaces
!> (kg/curb km)
   real*8, dimension (:), allocatable :: dirtmx
!> fraction of HRU area that is impervious (both directly and indirectly
!> connected) (fraction)
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
   real*8 :: sweepeff !< removal efficiency of sweeping operation (none)

!> random roughness for a given HRU (mm)
   real*8, dimension (:), allocatable :: ranrns_hru
   integer, dimension (:), allocatable :: itill
!> depth of mixing caused by tillage operation (mm)
   real*8, dimension (:), allocatable :: deptil
!> mixing efficiency of tillage operation (none)
   real*8, dimension (:), allocatable :: effmix
!> random roughness of a given tillage operation (mm)
   real*8, dimension (:), allocatable :: ranrns
!> 8-character name for the tillage operation
   character(len=8), dimension (550) :: tillnm
!> For ICODES equal to (none)\n
!> 0,1,3,5,9: not used\n
!> 2: fraction of overland flow in channel\n
!> 4: amount of water transferred (as defined by INUM4S)\n
!> 7,8,10,11: drainage area in square kilometers associated with the record
!> file\n
!> 12: rearation coefficient
   real*8, dimension (:), allocatable :: rnum1s
!> total drainage area of hydrograph in square kilometers (km^2)
   real*8, dimension (:), allocatable :: hyd_dakm
!> shyd(1,:) water (m^3 H2O)\n
!> shyd(2,:) sediment or suspended solid load (metric tons)\n
!> shyd(3,:) organic nitrogen (kg N)\n
!> shyd(4,:) organic phosphorus (kg P)\n
!> shyd(5,:) nitrate (kg N)\n
!> shyd(6,:) soluble phosphorus (kg P)\n
!> shyd(7,:) soluble pesticides (kg P)\n
!> shyd(8,:) sorbed pesticides (kg P)
   real*8, dimension (:,:), allocatable :: shyd
!> varoute(:,:) daily routing storage array (varies):\n
!> varoute(1,:) temperature (deg C)\n
!> varoute(2,:) water (m^3 H2O)\n
!> varoute(3,:) sediment or suspended solid load (metric tons)\n
!> varoute(4,:) organic nitrogen (kg N)\n
!> varoute(5,:) organic phosphorus (kg P)\n
!> varoute(6,:) nitrate (kg N)\n
!> varoute(7,:) mineral phosphorus (kg P)\n
!> varoute(11,:) pesticide in solution (mg pst)\n
!> varoute(12,:) pesticide sorbed to sediment (mg pst)\n
!> varoute(13,:) chlorophyll-a (kg)\n
!> varoute(14,:) ammonium (kg N)\n
!> varoute(15,:) nitrite (kg N)\n
!> varoute(16,:) carbonaceous biological oxygen demand (kg)\n
!> varoute(17,:) dissolved oxygen (kg)\n
!> varoute(18,:) persistent bacteria (# cfu/100ml)\n
!> varoute(19,:) less persistent bacteria (# cfu/100ml)
   real*8, dimension (:,:), allocatable :: varoute
   real*8, dimension (:,:), allocatable :: vartran
!> routing storage array for hourly time step (varies)\n
!> hhvaroute(2,:,:) water (m^3 H2O)\n
!> hhvaroute(4,:,:) organic nitrogen (kg N)\n
!> hhvaroute(5,:,:) organic posphorus (kg P)\n
!> hhvaroute(6,:,:) nitrate (kg N)\n
!> hhvaroute(7,:,:) soluble phosphorus (kg P)\n
!> hhvaroute(13,:,:) chlorophyll-a (kg)\n
!> hhvaroute(14,:,:) ammonium (kg N)\n
!> hhvaroute(15,:,:) nitrite (kg N)\n
!> hhvaroute(16,:,:) carbonaceous biological oxygen demand (kg)\n
!> hhvaroute(17,:,:) dissolved oxygen (kg O2)\n
!> hhvaroute(18,:,:) persistent bacteria (# cfu/100ml)\n
!> hhvaroute(19,:,:) less persistent bacteria (# cfu/100ml)
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
!> 15 =\n
!> 16 = autocal\n
!> 17 = routing unit
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
   integer, dimension (:), allocatable :: grwat_i
!> length of grass waterway (km)
   real*8, dimension (:), allocatable :: grwat_l
!> average width of grassed waterway (m)
   real*8, dimension (:), allocatable :: grwat_w
!> depth of grassed waterway from top of bank to bottom (m)
   real*8, dimension (:), allocatable :: grwat_d
!> average slope of grassed waterway channel (m)
   real*8, dimension (:), allocatable :: grwat_s
!> linear parameter defined by user for calculating sediment transport in
!> grassed waterways (none)
   real*8, dimension (:), allocatable :: grwat_spcon
!> time of concentration for grassed waterway and its drainage area (none)
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
!> surface area of impounded water body (ha)
   real*8, dimension (:), allocatable :: potsa
!> maximum volume of water stored in the depression/impounded area (read in as
!> mm and converted to m^3) (needed only if current HRU is IPOT) (mm)
   real*8, dimension (:), allocatable :: pot_volx
!> wetting front matric potential (average capillary suction at wetting front)
!> (mm)
   real*8, dimension (:), allocatable :: wfsh
!> water entering pothole on day (m^3 H2O)
   real*8, dimension (:), allocatable :: potflwi
!> sediment entering pothole on day (metric tons)
   real*8, dimension (:), allocatable :: potsedi
!> nitrate decay rate in impounded area (1/day)
   real*8, dimension (:), allocatable :: pot_no3l
!> normal sediment concentration in impounded water (needed only if current HRU
!> is IPOT)(mg/L)
   real*8, dimension (:), allocatable :: pot_nsed
!> nitrate-N concentration in groundwater loading to reach (mg N/L)
   real*8, dimension (:), allocatable :: gwno3
!> infiltration rate for last time step from the previous day (mm/hr)
   real*8, dimension (:), allocatable :: newrti
!> reduction in bacteria loading from filter strip (none)
   real*8, dimension (:), allocatable :: fsred
!> amount of nitrate in pothole water body (kg N)
   real*8, dimension (:), allocatable :: pot_no3
!> amount of sediment in pothole water body (metric tons)
   real*8, dimension (:), allocatable :: pot_sed
   real*8, dimension (:), allocatable :: tmpavp
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
!> average slope steepness in HRU (m/m)
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
!> secchi-disk depth of pond (m)
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
!> SCS runoff curve number for moisture condition I (none)
   real*8, dimension (:), allocatable :: cn1
!> amount of nitrate originating from lateral flow in pond at end of day or at
!> beginning of day(kg N)
   real*8, dimension (:), allocatable :: pnd_no3s
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
!> amount of nitrate originating from surface runoff in pond at end of day or at
!> beginning of day (kg N)
   real*8, dimension (:), allocatable :: pnd_no3
!> amount of soluble P originating from surface runoff in pond at end of day or
!> at beginning of day (kg P)
   real*8, dimension (:), allocatable :: pnd_solp
!> annual yield (dry weight) in the HRU (metric tons/ha)
   real*8, dimension (:), allocatable :: yldanu
!> coefficient for pesticide drift directly onto stream (none)
   real*8, dimension (:), allocatable :: driftco
!> amount of organic N originating from surface runoff in pond at end of day or
!> at beginning of day (kg N)
   real*8, dimension (:), allocatable :: pnd_orgn
!> amount of organic P originating from surface runoff in pond at end of day or
!> at beginning of day (kg P)
   real*8, dimension (:), allocatable :: pnd_orgp
!> SCS runoff curve number for moisture condition III (none)
   real*8, dimension (:), allocatable :: cn3
!> water lost through seepage from ponds on day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: twlpnd
!> water lost through seepage from wetlands on day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: twlwet
!> fraction of subbasin area contained in HRU (km^2/km^2)
   real*8, dimension (:), allocatable :: hru_fr
!> amount of water held in soil profile at saturation (mm H2O)
   real*8, dimension (:), allocatable :: sol_sumul
!> amount of chlorophyll-a in pond at end of day (kg chl_a)
   real*8, dimension (:), allocatable :: pnd_chla
!> area of HRU in square kilometers (km^2)
   real*8, dimension (:), allocatable :: hru_km
!> land cover/crop biomass (dry weight) (kg/ha)
   real*8, dimension (:), allocatable :: bio_ms
!> albedo when soil is moist (none)
   real*8, dimension (:), allocatable :: sol_alb
!> fraction of potential plant growth achieved on the day where the reduction is
!> caused by water stress (none)
   real*8, dimension (:), allocatable :: strsw
!> fraction of HRU/subbasin area that drains into ponds (none)
   real*8, dimension (:), allocatable :: pnd_fr
!> hydraulic conductivity through bottom of ponds (mm/hr)
   real*8, dimension (:), allocatable :: pnd_k
!> surface area of ponds when filled to principal spillway (ha)
   real*8, dimension (:), allocatable :: pnd_psa
!> runoff volume of water from catchment area needed to fill the ponds to the
!> principal spillway (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: pnd_pvol
!> surface area of ponds when filled to emergency spillway (ha)
   real*8, dimension (:), allocatable :: pnd_esa
!> runoff volume of water from catchment area needed to fill the ponds to the
!> emergency spillway (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: pnd_evol
!> volume of water in ponds (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: pnd_vol
!> average annual yield (dry weight) in the HRU (metric tons)
   real*8, dimension (:), allocatable :: yldaa
!> normal sediment concentration in pond water (UNIT CHANGE!) (mg/kg or kg/kg)
   real*8, dimension (:), allocatable :: pnd_nsed
!> sediment concentration in pond water (UNIT CHANGE!) (mg/kg or kg/kg)
   real*8, dimension (:), allocatable :: pnd_sed
!> depth to impervious layer (mm)
   real*8, dimension (:), allocatable :: dep_imp
   real*8, dimension (:), allocatable :: strsa
   real*8, dimension (:), allocatable :: evpnd, evwet
!> fraction of HRU/subbasin area that drains into wetlands (none)
   real*8, dimension (:), allocatable :: wet_fr
!> hydraulic conductivity of bottom of wetlands (mm/hr)
   real*8, dimension (:), allocatable :: wet_k
!> surface area of wetlands in subbasin at normal water level (ha)
   real*8, dimension (:), allocatable :: wet_nsa
!> runoff volume of water from catchment area needed to fill wetlands to normal
!> water level (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: wet_nvol
   integer, dimension (:), allocatable :: iwetgw, iwetile
!> surface area of wetlands at maximum water level (ha)
   real*8, dimension (:), allocatable :: wet_mxsa
!> runoff volume of water from catchment area needed to fill wetlands to maximum
!> water level (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: wet_mxvol
!> volume of water in wetlands (UNIT CHANGE!) (10^4 m^3 H2O or m^3 H2O)
   real*8, dimension (:), allocatable :: wet_vol
!> normal sediment concentration in wetland water (UNIT CHANGE!)
!> (mg/kg or kg/kg)
   real*8, dimension (:), allocatable :: wet_nsed
!> sediment concentration in wetland water (UNIT CHANGE!) (mg/L or kg/L)
   real*8, dimension (:), allocatable :: wet_sed
!> bp(1,:) 1st shape parameter for the pond surface area equation (none)\n
!> bp(2,:) 2nd shape parameter for the pond surface area equation (none)
   real*8, dimension (:,:), allocatable :: bp
!> retention coefficient for CN method based on plant ET (none)
   real*8, dimension (:), allocatable :: sci
!> retention coefficient for CN method based on soil moisture (none)
   real*8, dimension (:), allocatable :: smx
!> bw(1,:) 1st shape parameter for the wetland surface area equation (none)\n
!> bw(2,:) 2nd shape parameter for the wetland surface area equation (none)
   real*8, dimension (:,:), allocatable :: bw
!> persistent bacteria in soil solution (# cfu/m^2)
   real*8, dimension (:), allocatable :: bactpq
!> curve number for current day, HRU and at current soil moisture (none)
   real*8, dimension (:), allocatable :: cnday
!> less persistent bacteria on foliage (# cfu/m^2)
   real*8, dimension (:), allocatable :: bactlp_plt
!> persistent bacteria on foliage (# cfu/m^2)
   real*8, dimension (:), allocatable :: bactp_plt
!> fertilizer application efficiency calculated as the amount of N applied
!> divided by the amount of N removed at harvest (none)
   real*8, dimension (:), allocatable :: auto_eff
!> water clarity coefficient for wetland (none)
   real*8, dimension (:), allocatable :: secciw
!> amount of water stored in soil profile at end of any given day (mm H2O)
   real*8, dimension (:), allocatable :: sol_sw
!> less persistent bacteria in soil solution (# cfu/m^2)
   real*8, dimension (:), allocatable :: bactlpq
!> chlorophyll-a production coefficient for wetland (none)
   real*8, dimension (:), allocatable :: chlaw
!> average air temperature on current day in HRU (deg C)
   real*8, dimension (:), allocatable :: tmpav
!> less persistent bacteria attached to soil particles (# cfu/m^2)
   real*8, dimension (:), allocatable :: bactlps
!> persistent bacteria attached to soil particles (# cfu/m^2)
   real*8, dimension (:), allocatable :: bactps
!> amount of water stored as snow in HRU on current day (mm H2O)
   real*8, dimension (:), allocatable :: sno_hru
!> amount of organic N originating from surface runoff in wetland at end of day
!> (kg N)
   real*8, dimension (:), allocatable :: wet_orgn
!> solar radiation for the day in HRU (MJ/m^2)
   real*8, dimension (:), allocatable :: hru_ra
!> precipitation for the day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: subp
   real*8, dimension (:), allocatable :: rsdin !< initial residue cover (kg/ha)
!> minimum air temperature on current day in HRU (deg C)
   real*8, dimension (:), allocatable :: tmn
!> maximum air temperature on current day in HRU (deg C)
   real*8, dimension (:), allocatable :: tmx
!> last maximum temperature in HRU (deg C)
   real*8, dimension (:), allocatable :: tmp_hi
!> last minimum temperature in HRU (deg C)
   real*8, dimension (:), allocatable :: tmp_lo
!> USLE equation soil erodibility (K) factor (none)
   real*8, dimension (:), allocatable :: usle_k
!> time of concentration for HRU (hour)
   real*8, dimension (:), allocatable :: tconc
!> maximum possible solar radiation for the day in HRU (MJ/m^2)
   real*8, dimension (:), allocatable :: hru_rmx
!> fraction of total plant biomass that is in roots (none)
   real*8, dimension (:), allocatable :: rwt
   real*8, dimension (:), allocatable :: olai
   real*8, dimension (:), allocatable :: usle_cfac,usle_eifac
!> amount of water held in soil profile at field capacity (mm H2O)
   real*8, dimension (:), allocatable :: sol_sumfc
!> time for flow from farthest point in subbasin to enter a channel (hour)
   real*8, dimension (:), allocatable :: t_ov
!> total amount of NO3 applied during the year in auto-fertilization (kg N/ha)
   real*8, dimension (:), allocatable :: anano3
!> amount of water applied to HRU on current day (mm H2O)
   real*8, dimension (:), allocatable :: aird
!> amount of organic P originating from surface runoff in wetland at end of day
!> (kg P)
   real*8, dimension (:), allocatable :: wet_orgp
!> average porosity for entire soil profile (none)
   real*8, dimension (:), allocatable :: sol_avpor
!> product of USLE K,P,LS,exp(rock) (none)
   real*8, dimension (:), allocatable :: usle_mult
!> relative humidity for the day in HRU (none)
   real*8, dimension (:), allocatable :: rhd
!> wind speed (measured at 10 meters above surface) for the day in HRU (m/s)
   real*8, dimension (:), allocatable :: u10
   real*8, dimension (:), allocatable :: cht !< canopy height (m)
!> average annual amount of irrigation water applied to HRU (mm H2O)
   real*8, dimension (:), allocatable :: aairr
!> maximum leaf area index for the entire period of simulation in the HRU (none)
   real*8, dimension (:), allocatable :: lai_aamx
!> amount of water removed from deep aquifer for irrigation (mm H2O)
   real*8, dimension (:), allocatable :: deepirr
!> amount of water removed from shallow aquifer for irrigation (mm H2O)
   real*8, dimension (:), allocatable :: shallirr
!> amount of nitrate originating from surface runoff in wetland at end of day
!> (kg N)
   real*8, dimension (:), allocatable :: wet_no3
!> overland flow onto HRU from upstream routing unit (mm H2O)
   real*8, dimension (:), allocatable :: ovrlnd
!> amount of water held in canopy storage (mm H2O)
   real*8, dimension (:), allocatable :: canstor
!> maximum irrigation amount per auto application (mm)
   real*8, dimension (:), allocatable :: irr_mx
!> water stress factor which triggers auto irrigation (none or mm)
   real*8, dimension (:), allocatable :: auto_wstr
!> fertilizer/manure identification number from database (fert.dat) (none)
   real*8, dimension (:), allocatable :: cfrt_id
!> amount of fertilzier/manure applied to HRU on a given day ((kg/ha)/day)
   real*8, dimension (:), allocatable :: cfrt_kg
   real*8, dimension (:), allocatable :: cpst_id
   real*8, dimension (:), allocatable :: cpst_kg
   real*8, dimension (:), allocatable :: irr_asq !< surface runoff ratio
   real*8, dimension (:), allocatable :: irr_eff
!> surface runoff ratio (0-1) .1 is 10% surface runoff (frac)
   real*8, dimension (:), allocatable :: irrsq
!> concentration of salt in irrigation water (mg/kg)
   real*8, dimension (:), allocatable :: irrsalt
   real*8, dimension (:), allocatable :: irrefm
!> dry weight of biomass removed by grazing daily ((kg/ha)/day)
   real*8, dimension (:), allocatable :: bio_eat
!> dry weight of biomass removed by trampling daily ((kg/ha)/day)
   real*8, dimension (:), allocatable :: bio_trmp
!> number of days between applications (days)
   integer, dimension (:), allocatable :: ipst_freq
!> number of days between applications in continuous fertlizer operation (days)
   integer, dimension (:), allocatable :: ifrt_freq
   integer, dimension (:), allocatable :: irr_noa
   integer, dimension (:), allocatable :: irr_sc,irr_no
!> release/impound action code (none):\n
!> 0 begin impounding water\n
!> 1 release impounded water
   integer, dimension (:), allocatable :: imp_trig
!> number of days continuous fertilization will be simulated (none)
   integer, dimension (:), allocatable :: fert_days
   integer, dimension (:), allocatable :: irr_sca
!> land cover/crop identification code for first crop grown in HRU (the only
!> crop if there is no rotation) (from crop.dat) (none)
   integer, dimension (:), allocatable :: idplt
!> water stress identifier (none):\n
!> 1 plant water demand\n
!> 2 soil water deficit
   integer, dimension (:), allocatable :: wstrs_id
   integer, dimension (:), allocatable :: pest_days
   real*8, dimension (:,:), allocatable :: bio_aahv
   real*8, dimension (:), allocatable :: cumei,cumeira
   real*8, dimension (:), allocatable :: cumrt, cumrai
!> amount of soluble P originating from surface runoff in wetland at end of day
!> (kg P)
   real*8, dimension (:), allocatable :: wet_solp
!> amount of chlorophyll-a in wetland at end of day (kg chla)
   real*8, dimension (:), allocatable :: wet_chla
!> amount of nitrate originating from lateral flow in wetland at end of day
!> (kg N)
   real*8, dimension (:), allocatable :: wet_no3s
!> amount of soluble pesticide leached from bottom of soil profile on current
!> day (kg pst/ha)
   real*8, dimension (:), allocatable :: pstsol
!> amount of nitrate originating from groundwater in pond at end of day or at
!> beginning of day (kg N)
   real*8, dimension (:), allocatable :: pnd_no3g
!> secchi-disk depth in wetland at end of day (m)
   real*8, dimension (:), allocatable :: wet_seci
!> groundwater delay: time required for water leaving the bottom of the root
!> zone to reach the shallow aquifer (days)
   real*8, dimension (:), allocatable :: delay
   real*8, dimension (:), allocatable :: gwht !< groundwater height (m)
!> groundwater contribution to streamflow from HRU on current day (mm H2O)
   real*8, dimension (:), allocatable :: gw_q
!> amount of soluble P originating from groundwater in pond at end of day or at
!> beginning of day (kg P)
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
!> groundwater contribution to streamflow from deep aquifer from HRU on current
!> day (mm H2O)
   real*8, dimension (:), allocatable :: gw_qdeep
!> \f$\exp(-1/delay)\f$ where delay(:) is the groundwater delay (time required
!> for water leaving the bottom of the root zone to reach the shallow aquifer;
!> units-days) (none)
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
!> amount of water recharging both aquifers on current day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: rchrg
!> minimum plant biomass for grazing (kg/ha)
   real*8, dimension (:), allocatable :: bio_min
!> initial HRU soil water content expressed as fraction of field capacity (none)
   real*8, dimension (:), allocatable :: ffc
!> amount of soluble phosphorus in surface runoff in HRU for the day (kg P/ha)
   real*8, dimension (:), allocatable :: surqsolp
!> depth of water in deep aquifer (mm H2O)
   real*8, dimension (:), allocatable :: deepst
!> depth of water in shallow aquifer in HRU (mm H2O)
   real*8, dimension (:), allocatable :: shallst
!> amount of soluble P originating from groundwater in wetland at end of day
!> (kg P)
   real*8, dimension (:), allocatable :: wet_solpg
   real*8, dimension (:), allocatable :: cklsp
   real*8, dimension (:), allocatable :: rchrg_src
!> filter strip trapping efficiency (used for everything but bacteria) (none)
   real*8, dimension (:), allocatable :: trapeff
!> average bulk density for soil profile (Mg/m^3)
   real*8, dimension (:), allocatable :: sol_avbd
!> amount of nitrate originating from groundwater in wetland at end of day
!> (kg N)
   real*8, dimension (:), allocatable :: wet_no3g
!> time to drain soil to field capacity yield used in autofertilization (hours)
   real*8, dimension (:), allocatable :: tdrain
!> threshold depth of water in shallow aquifer required before groundwater flow
!> will occur (mm H2O)
   real*8, dimension (:), allocatable :: gwqmn
!> temperature of snow pack in HRU (deg C)
   real*8, dimension (:), allocatable :: snotmp
!> plant uptake of phosphorus in HRU for the day (kg P/ha)
   real*8, dimension (:), allocatable :: pplnt
!> drain tile lag time: the amount of time between the transfer of water from
!> the soil to the drain tile and the release of the water from the drain tile
!> to the reach (hours)
   real*8, dimension (:), allocatable :: gdrain
!> depth of drain tube from the soil surface (mm)
   real*8, dimension (:), allocatable :: ddrain
!> crack volume potential of soil (none)
   real*8, dimension (:), allocatable :: sol_crk
!> fraction of surface runoff within the subbasin which takes 1 day or less to
!> reach the subbasin outlet (none)
   real*8, dimension (:), allocatable :: brt
!> length of the current day (hours)
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
!> maximum surface depressional storage for day in a given HRU (mm)
   real*8, dimension (:), allocatable :: stmaxd
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd3
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd2
!> time that solids have built-up on streets (days)
   real*8, dimension (:), allocatable :: twash
!> dissolved oxygen concentration in the surface runoff on current day in HRU
!> (mg/L)
   real*8, dimension (:), allocatable :: doxq
!> soil water content used to calculate daily CN value (initial soil water
!> content for day) (mm H2O)
   real*8, dimension (:), allocatable :: sol_cnsw
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd8
!> random number between 0.0 and 1.0 (none)
   real*8, dimension (:), allocatable :: rnd9
!> amount of nitrate percolating past bottom of soil profile during the day
!> (kg N/ha)
   real*8, dimension (:), allocatable :: percn
   real*8, dimension (:), allocatable :: sol_sumwp
!> total or net amount of water entering main channel for day from HRU (mm H2O)
   real*8, dimension (:), allocatable :: qdr
!> amount of N applied in autofert operation in year (kg N/ha)
   real*8, dimension (:), allocatable :: tauton
!> amount of P applied in autofert operation in year (kg N/ha)
   real*8, dimension (:), allocatable :: tautop
!> carbonaceous biological oxygen demand of surface runoff on current day in HRU
!> (mg/L)
   real*8, dimension (:), allocatable :: cbodu
!> chlorophyll-a concentration in water yield on current day in HRU
!> (microgram/L)
   real*8, dimension (:), allocatable :: chl_a
   real*8, dimension (:), allocatable :: tfertn,tfertp,tgrazn,tgrazp
!> total amount of water in lateral flow in soil profile for the day in HRU
!> (mm H2O)
   real*8, dimension (:), allocatable :: latq
!> plant uptake of nitrogen in HRU for the day (kg N/ha)
   real*8, dimension (:), allocatable :: nplnt
!> amount of nitrate transported with lateral flow in HRU for the day (kg N/ha)
   real*8, dimension (:), allocatable :: latno3
!> soluble P loading to reach in groundwater (kg P/ha)
   real*8, dimension (:), allocatable :: minpgw
!> nitrate loading to reach in groundwater (kg N/ha)
   real*8, dimension (:), allocatable :: no3gw
   real*8, dimension (:), allocatable :: tileq, tileno3
!> amount of organic nitrogen in surface runoff in HRU for the day (kg N/ha)
   real*8, dimension (:), allocatable :: sedorgn
!> amount of active mineral phosphorus sorbed to sediment in surface runoff in
!> HRU for day (kg P/ha)
   real*8, dimension (:), allocatable :: sedminpa
!> amount of stable mineral phosphorus sorbed to sediment in surface runoff in
!> HRU for day (kg P/ha)
   real*8, dimension (:), allocatable :: sedminps
!> soil loss caused by water erosion for day in HRU (metric tons)
   real*8, dimension (:), allocatable :: sedyld
!> percolation from bottom of soil profile for the day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: sepbtm
!> fraction of potential plant growth achieved on the day where the reduction is
!> caused by nitrogen stress (none)
   real*8, dimension (:), allocatable :: strsn
!> amount of organic phosphorus in surface runoff in HRU for the day (kg P/ha)
   real*8, dimension (:), allocatable :: sedorgp
!> surface runoff generated in HRU on the current day (mm H2O)
   real*8, dimension (:), allocatable :: surfq
!> fraction of potential plant growth achieved on the day in HRU where the
!> reduction is caused by temperature stress (none)
   real*8, dimension (:), allocatable :: strstmp
!> fraction of potential plant growth achieved on the day where the reduction is
!> caused by phosphorus stress (none)
   real*8, dimension (:), allocatable :: strsp
!> amount of nitrate transported in surface runoff in HRU for the day (kg N/ha)
   real*8, dimension (:), allocatable :: surqno3
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
!> optimal harvest index adjusted for water stress for current time during
!> growing season ((kg/ha)/(kg/ha))
   real*8, dimension (:), allocatable :: hvstiadj
!> leaf area index for HRU (m^2/m^2)
   real*8, dimension (:), allocatable :: laiday
!> chlorophyll-a production coefficient for pond (none)
   real*8, dimension (:), allocatable :: chlap
!> amount of mineral P attached to sediment originating from surface runoff in
!> pond at end of day or beginnig of day (kg P)
   real*8, dimension (:), allocatable :: pnd_psed
   real*8, dimension (:), allocatable :: laimxfr
!> water clarity coefficient for pond (none)
   real*8, dimension (:), allocatable :: seccip
!> amount of nitrogen in plant biomass (kg N/ha)
   real*8, dimension (:), allocatable :: plantn
!> actual ET simulated during life of plant (mm H2O)
   real*8, dimension (:), allocatable :: plt_et
!> amount of mineral P attached to sediment originating from surface runoff in
!> wetland at end of day (kg P)
   real*8, dimension (:), allocatable :: wet_psed
!> average annual biomass (dry weight) in the HRU (metric tons)
   real*8, dimension (:), allocatable :: bio_aams
!> amount of phosphorus stored in plant biomass (kg P/ha)
   real*8, dimension (:), allocatable :: plantp
!> potential ET simulated during life of plant (mm H2O)
   real*8, dimension (:), allocatable :: plt_pet
!> time threshold used to define dormant period for plant (when daylength is
!> within the time specified by dl from the minimum daylength for the area, the
!> plant will go dormant) (hour)
   real*8, dimension (:), allocatable :: dormhr
!> maximum leaf area index for the year in the HRU (none)
   real*8, dimension (:), allocatable :: lai_yrmx
   real*8, dimension (:), allocatable :: bio_aamx
!> amount of pesticide in lateral flow in HRU for the day (kg pst/ha)
   real*8, dimension (:), allocatable :: lat_pst
!> fraction of HRU area that drains into floodplain (km^2/km^2)
   real*8, dimension (:), allocatable :: fld_fr
   real*8, dimension (:), allocatable :: orig_snohru,orig_potvol
!> fraction of plant biomass that is nitrogen (none)
   real*8, dimension (:), allocatable :: pltfr_n
   real*8, dimension (:), allocatable :: orig_alai,orig_bioms
!> fraction of plant biomass that is phosphorus (none)
   real*8, dimension (:), allocatable :: pltfr_p
   real*8, dimension (:), allocatable :: orig_phuacc,orig_sumix
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
!> water table based on 30 day antecedent climate (precip,et) (mm)
   real*8, dimension (:), allocatable :: wtab
   real*8, dimension (:), allocatable :: wtab_mn,wtab_mx
!> nitrate concentration in shallow aquifer converted to kg/ha (ppm NO3-N)
   real*8, dimension (:), allocatable :: shallst_n
   real*8, dimension (:), allocatable :: gw_nloss,rchrg_n
   real*8, dimension (:), allocatable :: det_san, det_sil, det_cla
   real*8, dimension (:), allocatable :: det_sag, det_lag
!> fraction of fertilizer which is applied to top 10 mm of soil (the remaining
!> fraction is applied to first soil layer) (none)
   real*8, dimension (:), allocatable :: afrt_surface
!> estimated/target nitrogen content of yield used in autofertilization
!> (kg N/kg yield)
   real*8, dimension (:), allocatable :: tnylda
!> fraction of fertilizer which is applied to the top 10 mm of soil (the
!> remaining fraction is applied to the first soil layer) (none)
   real*8 :: frt_surface
!> maximum NO3-N content allowed to be applied in one year by auto-fertilization
!> (kg NO3-N/ha)
   real*8, dimension (:), allocatable :: auto_nyr
!> maximum NO3-N content allowed in one fertilizer application (kg NO3-N/ha)
   real*8, dimension (:), allocatable :: auto_napp
!> nitrogen stress factor which triggers auto fertilization (none)
   real*8, dimension (:), allocatable :: auto_nstrs
!> dry weight of manure deposited on HRU daily ((kg/ha)/day)
   real*8, dimension (:), allocatable :: manure_kg
   real*8, dimension (:,:), allocatable :: rcn_mo, rammo_mo
   real*8, dimension (:,:), allocatable :: drydep_no3_mo, drydep_nh4_mo
   real*8, dimension (:), allocatable :: rcn_d, rammo_d
   real*8, dimension (:), allocatable :: drydep_no3_d, drydep_nh4_d
   real*8, dimension (:,:), allocatable :: yldn
   integer, dimension (:,:), allocatable :: gwati
   real*8, dimension (:,:), allocatable :: gwatn, gwatl
   real*8, dimension (:,:), allocatable :: gwatw, gwatd, gwatveg
   real*8, dimension (:,:), allocatable :: gwata, gwats, gwatspcon
   real*8, dimension (:,:), allocatable :: rfqeo_30d,eo_30d
!> psetlp(1,:) phosphorus settling rate for 1st season (m/day)\n
!> psetlp(2,:) phosphorus settling rate for 2nd season (m/day)
   real*8, dimension (:,:), allocatable :: psetlp
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
!> 1st shape parameter for calculation of water retention (none)
   real*8, dimension (:), allocatable :: wrt1
!> 2nd shape parameter for calculation of water retention (none)
   real*8, dimension (:), allocatable :: wrt2
!> pesticide enrichment ratio (none)
   real*8, dimension (:,:), allocatable :: pst_enr
!> amount of pesticide type lost in surface runoff on current day in HRU (kg/ha)
   real*8, dimension (:,:), allocatable :: pst_surq
!> division term from net pesticide equation (mm)
   real*8, dimension (:,:), allocatable :: zdb
!> pesticide on plant foliage (kg/ha)
   real*8, dimension (:,:), allocatable :: plt_pst
!> psetlw(1,:) phosphorus settling rate for 1st season (m/day)\n
!> psetlw(2,:) phosphorus settling rate for 2nd season (m/day)
   real*8, dimension (:,:), allocatable :: psetlw
!> pesticide loading from HRU sorbed onto sediment (kg/ha)
   real*8, dimension (:,:), allocatable :: pst_sed
!> average daily water removal from the pond for the month for the HRU within
!> the subbasin (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wupnd
!> phi(1,:) cross-sectional area of flow at bankfull depth (m^2)
!> phi(2,:) (none)
!> phi(3,:) (none)
!> phi(4,:) (none)
!> phi(5,:) flow rate when reach is at bankfull depth (m^3/s)
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
!> cross-sectional area of flow at bankfull depth (m^2)
   real*8, dimension (:), allocatable :: wat_phi1
!> flow rate when reach is at bankfull depth (m^3/s)
   real*8, dimension (:), allocatable :: wat_phi5
!> bottom width of main channel (m)
   real*8, dimension (:), allocatable :: wat_phi6
!> depth of water when reach is at bankfull depth (m)
   real*8, dimension (:), allocatable :: wat_phi7
!> average velocity when reach is at bankfull depth (m/s)
   real*8, dimension (:), allocatable :: wat_phi8
!> wave celerity when reach is at bankfull depth (m/s)
   real*8, dimension (:), allocatable :: wat_phi9
!> storage time constant for reach at bankfull depth (ratio of storage to
!> discharge) (hour)
   real*8, dimension (:), allocatable :: wat_phi10
!> average velocity when reach is at 0.1 bankfull depth (low flow) (m/s)
   real*8, dimension (:), allocatable :: wat_phi11
!> wave celerity when reach is at 0.1 bankfull depth (low flow) (m/s)
   real*8, dimension (:), allocatable :: wat_phi12
!> storage time constant for reach at 0.1 bankfull depth (low flow) (ratio of
!> storage to discharge) (hour)
   real*8, dimension (:), allocatable :: wat_phi13
!> snow water content in elevation band on current day (mm H2O)
   real*8, dimension (:,:), allocatable :: snoeb
!> average daily water removal from the deep aquifer for the month for the HRU
!> within the subbasin (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wudeep
!> average daily water removal from the shallow aquifer for the month for the
!> HRU within the subbasin (10^4 m^3/day)
   real*8, dimension (:,:), allocatable :: wushal
!> minimum temperature for the day in band in HRU (deg C)
   real*8, dimension (:,:), allocatable :: tmnband
!> bss(1,:) amount of lateral flow lagged (mm H2O)\n
!> bss(2,:) amount of nitrate in lateral flow lagged (kg N/ha)\n
!> bss(3,:) amount of tile flow lagged (mm)\n
!> bss(4,:) amount of nitrate in tile flow lagged (kg N/ha)
   real*8, dimension (:,:), allocatable :: bss
!> nsetlw(1,:) nitrogen settling rate for 1st season (m/day)\n
!> nsetlw(2,:) nitrogen settling rate for 2nd season (m/day)
   real*8, dimension (:,:), allocatable :: nsetlw
!> temperature of snow pack in elevation band (deg C)
   real*8, dimension (:,:), allocatable :: snotmpeb
!> surf_bs(1,:) amount of surface runoff lagged over one day (mm H2O)\n
!> surf_bs(2,:) amount of sediment yield lagged over one day (metric tons)\n
!> surf_bs(3,:) amount of organic nitrogen loading lagged over one day
!> (kg N/ha)\n
!> surf_bs(4,:) amount of organic phosphorus loading lagged over one day
!> (kg P/ha)\n
!> surf_bs(5,:) amount of nitrate loading in surface runoff lagged over one day
!> (kg N/ha)\n
!> surf_bs(6,:) amount of soluble phosphorus loading lagged over one day
!> (kg P/ha)\n
!> surf_bs(7,:) amount of active mineral phosphorus loading lagged over one day
!> (kg P/ha)\n
!> surf_bs(8,:) amount of stable mineral phosphorus loading lagged over one day
!> (kg P/ha)\n
!> surf_bs(9,:) amount of less persistent bacteria in solution lagged over one
!> day (# colonies/ha)\n
!> surf_bs(10,:) amount of persistent bacteria in solution lagged over one day
!> (# colonies/ha)\n
!> surf_bs(11,:) amount of less persistent bacteria sorbed lagged over one day
!> (# colonies/ha)\n
!> surf_bs(12,:) amount of persistent bacteria sorbed lagged over one day
!> (# colonies/ha)
   real*8, dimension (:,:), allocatable :: surf_bs
!> nsetlp(1,:) nitrogen settling rate for 1st season (m/day)\n
!> nsetlp(2,:) nitrogen settling rate for 2nd season (m/day)
   real*8, dimension (:,:), allocatable :: nsetlp
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
!> pst_lag(1,:,:) amount of soluble pesticide in surface runoff lagged
!> (kg pst/ha)\n
!> pst_lag(2,:,:) amount of sorbed pesticide in surface runoff lagged
!> (kg pst/ha)\n
!> pst_lag(3,:,:) amount of pesticide lagged (kg pst/ha)
   real*8, dimension (:,:,:), allocatable :: pst_lag
!> pesticide use flag (none)\n
!>  0: no pesticides used in HRU\n
!>  1: pesticides used in HRU
   integer, dimension (:), allocatable :: hrupest
!> sequence number of impound/release operation within the year (none)
   integer, dimension (:), allocatable :: nrelease
!> rainfall event flag (none):\n
!> 0: no rainfall event over midnight\n
!> 1: rainfall event over midnight
   integer, dimension (:), allocatable :: swtrg
!> number of years of rotation (none)
   integer, dimension (:), allocatable :: nrot
!> sequence number of fertilizer application within the year (none)
   integer, dimension (:), allocatable :: nfert
!> sequence number of year in rotation (none)
   integer, dimension (:), allocatable :: nro
!> land cover status code (none). This code informs the model whether or not a
!> land cover is growing at the beginning of the simulation\n
!> 0 no land cover currently growing\n
!> 1 land cover growing
   integer, dimension (:), allocatable :: igro
!> ipnd(1,:) beginning month of 2nd "season" of nutrient settling season
!> (none)\n
!> ipnd(2,:) ending month of 2nd "season" of nutrient settling season (none)
   integer, dimension (:,:), allocatable :: ipnd
!> sequence number of auto-irrigation application within the year (none)
   integer, dimension (:), allocatable :: nair
!> iflod(1,:) beginning month of non-flood season (none)\n
!> iflod(2,:) ending month of non-flood season (none)
   integer, dimension (:,:), allocatable :: iflod
!> number of days required to reach target storage from current pond storage
!> (none)
   integer, dimension (:), allocatable :: ndtarg
!> sequence number of irrigation application within the year (none)
   integer, dimension (:), allocatable :: nirr
!> code for approach used to determine amount of nitrogen to HRU (none):\n
!> 0 nitrogen target approach\n
!> 1 annual max approach
   integer, dimension (:), allocatable :: nstress
   integer, dimension (:), allocatable :: iafrttyp
   integer, dimension (:), allocatable :: igrotree
!> number of days grazing will be simulated (none)
   integer, dimension (:), allocatable :: grz_days
!> management code (for GIS output only) (none)
   integer, dimension (:), allocatable :: nmgt
!> sequence number of auto-fert application within the year (none)
   integer, dimension (:), allocatable :: nafert
!> sequence number of street sweeping operation within the year (none)
   integer, dimension (:), allocatable :: nsweep
!> sequence number of crop grown within the current year (none)
   integer, dimension (:), allocatable :: icr
!> sequence number of harvest operation within a year (none)
   integer, dimension (:), allocatable :: ncut
!> irrigation source location (none)\n
!> if IRRSC=1, IRRNO is the number of the reach\n
!> if IRRSC=2, IRRNO is the number of the reservoir\n
!> if IRRSC=3, IRRNO is the number of the subbasin\n
!> if IRRSC=4, IRRNO is the number of the subbasin\n
!> if IRRSC=5, not used
   integer, dimension (:), allocatable :: irrno
!> number of soil layers in HRU (none)
   integer, dimension (:), allocatable :: sol_nly
!> prior day category (none)\n
!> 1 dry day\n
!> 2 wet day
   integer, dimension (:), allocatable :: npcp
!> average annual number of irrigation applications in HRU (none)
   integer, dimension (:), allocatable :: irn
!> sequence number of continuous fertilization operation within the year (none)
   integer, dimension (:), allocatable :: ncf
!> sequence number of grazing operation within the year (none)
   integer, dimension (:), allocatable :: ngr
!> grazing flag for HRU (none):\n
!> 0 HRU currently not grazed\n
!> 1 HRU currently grazed
   integer, dimension (:), allocatable :: igrz
!> number of days HRU has been grazed (days)
   integer, dimension (:), allocatable :: ndeat
!> subbasin number in which HRU/reach is located (none)
   integer, dimension (:), allocatable :: hru_sub
!> urban land type identification number from urban database (urban.dat) (none)
   integer, dimension (:), allocatable :: urblu
!> soil layer where drainage tile is located (none)
   integer, dimension (:), allocatable :: ldrain
!> dormancy status code (none):\n
!> 0 land cover growing (not dormant)\n
!> 1 land cover dormant
   integer, dimension (:), allocatable :: idorm
   integer, dimension (:), allocatable :: hru_seq
!> urban simulation code (none):\n
!> 0  no urban sections in HRU\n
!> 1  urban sections in HRU, simulate using USGS regression equations\n
!> 2  urban sections in HRU, simulate using build up/wash off algorithm
   integer, dimension (:), allocatable :: iurban
!> continuous fertilizer flag for HRU (none):\n
!> 0 HRU currently not continuously fertilized\n
!> 1 HRU currently continuously fertilized
   integer, dimension (:), allocatable :: icfrt
   integer, dimension (:), allocatable :: iday_fert
!> number of HRU (in subbasin) that is a floodplain (none)
   integer, dimension (:), allocatable :: ifld
!> number of HRU (in subbasin) that is a riparian zone (none)
   integer, dimension (:), allocatable :: irip
!> GIS code printed to output files (output.hru, output.rch) (none)
   integer, dimension (:), allocatable :: hrugis
!> number of days HRU has been continuously fertilized (days)
   integer, dimension (:), allocatable :: ndcfrt
!> irrigation source code (none):\n
!> 1 divert water from reach\n
!> 2 divert water from reservoir\n
!> 3 divert water from shallow aquifer\n
!> 4 divert water from deep aquifer\n
!> 5 divert water from source outside watershed
   integer, dimension (:), allocatable :: irrsc
!> sequence number of tillage operation within current year (none)
   integer, dimension (:), allocatable :: ntil
   integer, dimension (:), allocatable :: orig_igro
!> high water table code (none):\n
!> 0 no high water table\n
!> 1 high water table
   integer, dimension (:), allocatable :: iwatable
   integer, dimension (:), allocatable :: curyr_mat
!> icpst = 0 do not apply\n
!> icpst = 1 application period
   integer, dimension (:), allocatable :: icpst
!> current day within the application period (day)
   integer, dimension (:), allocatable :: ndcpst
   integer, dimension (:), allocatable :: ncpest
!> current day between applications (day)
   integer, dimension (:), allocatable :: iday_pest
   integer, dimension (:), allocatable :: irr_flag
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

   integer, dimension (:,:), allocatable :: mgt_sdr,idplrot
   integer, dimension (:,:), allocatable :: icont, iycont
   integer, dimension (:,:), allocatable :: ifilt, iyfilt
   integer, dimension (:,:), allocatable :: istrip, iystrip
   integer, dimension (:,:), allocatable :: iopday, iopyr, mgt_ops
!> total amount of pesticide type applied in watershed during simulation (kg/ha)
   real*8, dimension (:), allocatable :: wshd_pstap
!> amount of pesticide lost through degradation in watershed (kg pst/ha)
   real*8, dimension (:), allocatable :: wshd_pstdg
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
!> day (depends on model operational time step) (none)
   integer :: nstep
!> length of time step used to report precipitation data for sub-daily modeling
!> (operational time step) (minutes)
   integer :: idt
   real*8, dimension (:), allocatable :: hdepth !< depth of flow during hour (m)
!> water stored in reach at end of hour (m^3 H2O)
   real*8, dimension (:), allocatable :: hhstor
!> water leaving reach in hour (m^3)
   real*8, dimension (:), allocatable :: hrtwtr
!> flow rate in reach for hour (m^3/s)
   real*8, dimension (:), allocatable :: hsdti
!> water stored in reach at beginning of hour (m^3 H2O)
   real*8, dimension (:), allocatable :: hrchwtr
   real*8, dimension (:), allocatable :: halgae,horgn,hnh4
   real*8, dimension (:), allocatable :: hno2,hno3,horgp,hsolp,hbod
   real*8, dimension (:), allocatable :: hdisox,hchla,hsedyld,hsedst
!> cross-sectional area of flow (m^2)
   real*8, dimension (:), allocatable :: hharea
!> soluble pesticide concentration in outflow on day (mg pst/m^3)
   real*8, dimension (:), allocatable :: hsolpst
!> sorbed pesticide concentration in outflow on day (mg pst/m^3)
   real*8, dimension (:), allocatable :: hsorpst
!> surface runoff generated each timestep of day in HRU (mm H2O)
   real*8, dimension (:), allocatable :: hhqday
!> precipitation, or effective precipitation reaching soil surface, in time step
!> for HRU (mm H2O)
   real*8, dimension (:), allocatable :: precipdt
!> travel time of flow in reach for hour (hour)
   real*8, dimension (:), allocatable :: hhtime
!> less persistent bacteria in reach/outflow during hour (# cfu/100mL)
   real*8, dimension (:), allocatable :: hbactlp
!> persistent bacteria in reach/outflow during hour (# cfu/100mL)
   real*8, dimension (:), allocatable :: hbactp
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

   real*8, dimension (:,:), allocatable :: sol_bdp

   real*8, dimension (:,:), allocatable :: tillagef
   real*8, dimension (:), allocatable :: rtfr
!> storing last soil root depth for use in harvestkillop/killop (mm)
   real*8, dimension (:), allocatable :: stsol_rd
   integer:: urban_flag, dorm_flag
   real*8 :: bf_flg, iabstr
!> TSS loading from urban impervious cover (metric tons)
   real*8, dimension (:), allocatable :: ubntss
!> surface runoff from urban impervious cover (mm H2O)
   real*8, dimension (:), allocatable :: ubnrunoff
!> surface runoff from urban impervious cover in subbasin (mm H2O)
   real*8, dimension (:,:), allocatable :: sub_ubnrunoff
!> TSS loading from urban impervious cover in subbasin (metric tons)
   real*8, dimension (:,:), allocatable :: sub_ubntss
   real*8, dimension (:,:), allocatable :: ovrlnd_dt
   real*8, dimension (:,:,:), allocatable :: hhsurf_bs

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
!> sediment yield from HRU drung a time step applied to HRU (tons)
   real*8, dimension(:,:), allocatable :: hhsedy
   real*8, dimension(:,:), allocatable :: sub_subp_dt
!> sediment yield for the time step in subbasin (metric tons)
   real*8, dimension(:,:), allocatable :: sub_hhsedy
   real*8, dimension(:,:), allocatable :: sub_atmp
!> main channel hydraulic radius (m H2O)
   real*8, dimension(:), allocatable :: rhy
   real*8, dimension(:), allocatable :: init_abstrc
!> evaporation losses for hour (m^3 H2O)
   real*8, dimension(:), allocatable :: hrtevp
!> transmission losses for hour (m^3 H2O)
   real*8, dimension(:), allocatable :: hrttlc
   real*8, dimension(:), allocatable :: dratio
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
!> minimum residue allowed due to implementation of residue managment in the OPS
!> file (kg/ha)
   real*8, dimension (:), allocatable :: min_res
   real*8, dimension (:), allocatable :: harv_min, fstap
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
!> depth of rectangular weir at different stages (m)
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

   real*8, dimension(:), allocatable :: ri_subkm,ri_totpvol,&
      &irmmdt
!> total sediment deposited in the pond (tons)
   real*8, dimension(:,:), allocatable :: ri_sed
   real*8, dimension(:,:), allocatable :: ri_fr,ri_dim,&
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

   real*8 :: lid_vgcl !< van Genuchten equation's coefficient, l (none)
   real*8 :: lid_vgcm !< van Genuchten equation's coefficient, m (none)
   real*8 :: lid_qsurf_total,lid_farea_sum

!> cumulative amount of water infiltrated into the amended soil layer at the
!> last time step in a day (mm H2O)
   real*8, dimension(:,:), allocatable :: lid_cuminf_last
!> cumulative amount of rainfall at the last time step in a day (mm H2O)
   real*8, dimension(:,:), allocatable :: lid_cumr_last
!> cumulative amount of excess rainfall at the last time step in a day (mm H2O)
   real*8, dimension(:,:), allocatable :: lid_excum_last
!> potential infiltration rate of the amended soil layer at the last time step
!> in a day (mm/mm H2O)
   real*8, dimension(:,:), allocatable :: lid_f_last
!> soil water content of the amended soil layer at the last time step in a day
!> (mm/mm H2O)
   real*8, dimension(:,:), allocatable :: lid_sw_last
!> depth of runoff generated on a LID in a given time interval (mm H2O)
   real*8, dimension(:,:), allocatable :: lid_qsurf
   real*8, dimension(:,:), allocatable :: interval_last,lid_str_last,&
      &lid_farea,lid_sw_add,lid_cumqperc_last,lid_cumirr_last

   ! Green Roof
   integer, dimension(:,:), allocatable :: gr_onoff,gr_imo,gr_iyr
!> fractional area of a green roof to the HRU (none)
   real*8, dimension(:,:), allocatable :: gr_farea
   real*8, dimension(:,:), allocatable :: gr_solop,gr_etcoef,&
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

!> mass of C present in slow humus (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_HSC
!> mass of N present in slow humus (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_HSN
!> mass of C present in passive humus (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_HPC
!> mass of N present in passive humus (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_HPN
!> mass of metabolic litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LM
!> mass of C in metabolic litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LMC
!> mass of N in metabolic litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LMN
!> mass of structural litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LS
!> mass of C in structural litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LSC
!> mass of lignin in structural litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LSL
!> mass of N in structural litter (kg ha-1)
   real*8, dimension(:,:), allocatable :: sol_LSN

   real*8, dimension(:,:), allocatable :: sol_BMC, sol_BMN, sol_RNMN,&
      &sol_LSLC, sol_LSLNC, sol_RSPC, sol_WOC, sol_WON, sol_HP, sol_HS, sol_BM

   !!SOM-residue C/N state variables -- may need to be included
   real*8, dimension(:,:), allocatable :: sol_CAC, sol_CEC

   !!daily updated soil layer associated percolaton and lateral flow Carbon loss
   real*8, dimension(:,:), allocatable :: sol_percc
   real*8, dimension(:,:), allocatable :: sol_latc

!> amount of C lost with sediment pools (kg C/ha)
   real*8, dimension(:), allocatable :: sedc_d
   real*8, dimension(:), allocatable :: surfqc_d, latc_d,&
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
   real*8 :: dthy !< time interval for subdaily flood routing
   integer, dimension(4) :: IHX
   integer, dimension(:), allocatable :: NHY
   real*8, dimension(:), allocatable :: RCHX,RCSS,QCAP,CHXA,CHXP
   real*8, dimension(:,:,:), allocatable :: QHY

   !!Variables for killop.f90 and harvkillop.f90 files
   real*8 :: ff1, ff2

end module parm
