subroutine readswq

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads parameters from the subbasin instream water
!!    quality file (.swq) and initializes the QUAL2E variables which apply to
!!    the individual subbasins

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |reach number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bc1(:)      |1/day or 1/hr |rate constant for biological oxidation of NH3
!!                               |to NO2 in reach at 20 deg C
!!    bc2(:)      |1/day or 1/hr |rate constant for biological oxidation of NO2
!!                               |to NO3 in reach at 20 deg C
!!    bc3(:)      |1/day or 1/hr |rate constant for hydrolysis of organic N to
!!                               |ammonia in reach at 20 deg C
!!    bc4(:)      |1/day or 1/hr |rate constant for the decay of organic P to
!!                               |dissolved P in reach at 20 deg C
!!    chpst_koc(:)  |m**3/g      |pesticide partition coefficient between
!!                               |water and sediment in reach
!!    chpst_mix(:)  |m/day       |mixing velocity (diffusion/dispersion) for
!!                               |pesticide in reach
!!    chpst_rea(:)  |1/day       |pesticide reaction coefficient in reach
!!    chpst_rsp(:)  |m/day       |resuspension velocity in reach for pesticide
!!                               |sorbed to sediment
!!    chpst_stl(:)  |m/day       |settling velocity in reach for pesticide
!!                               |sorbed to sediment
!!    chpst_vol(:)  |m/day       |pesticide volatilization coefficient in reach
!!    rk1(:)      |1/day or 1/hr |CBOD deoxygenation rate coefficient in reach
!!                               |at 20 deg C
!!    rk2(:)      |1/day or 1/hr |reaeration rate in accordance with Fickian
!!                               |diffusion in reach at 20 deg C
!!    rk3(:)      |1/day or 1/hr |rate of loss of CBOD due to settling in reach
!!                               |at 20 deg C
!!    rk4(:)      |mg O2/        |sediment oxygen demand rate in reach
!!                |  ((m**2)*day)|at 20 deg C
!!                |or mg O2/((m**2)*hr)
!!    rk5(:)      |1/day         |coliform die-off rate in reach
!!    rk6(:)      |1/day         |decay rate for arbitrary non-conservative
!!                               |constituent in reach
!!    rs1(:)      |m/day or m/hr |local algal settling rate in reach at 20 deg C
!!    rs2(:)      |(mg disP-P)/  |benthos source rate for dissolved phosphorus
!!                |  ((m**2)*day)|in reach at 20 deg C
!!                |or (mg disP-P)/((m**2)*hr)|
!!    rs3(:)      |(mg NH4-N)/   |benthos source rate for ammonia nitrogen in
!!                |  ((m**2)*day)|reach at 20 deg C
!!                |or (mg NH4-N)/((m**2)*hr)|
!!    rs4(:)      |1/day or 1/hr |rate coefficient for organic nitrogen
!!                               |settling in reach at 20 deg C
!!    rs5(:)      |1/day or 1/hr |organic phosphorus settling rate in reach at
!!                               |20 deg C
!!    rs6(:)      |1/day         |rate coefficient for settling of arbitrary
!!                               |non-conservative constituent in reach
!!    rs7(:)      |(mg ANC)/     |benthal source rate for arbitrary
!!                   ((m**2)*day)|non-conservative constituent in reach
!!    sedpst_act(:) |m           |depth of active sediment layer in reach for
!!                               |pesticide
!!    sedpst_bry(:) |m/day       |pesticide burial velocity in river bed
!!                               |sediment
!!    sedpst_conc(:)|mg/(m**3)   |inital pesticide concentration in river bed
!!                               |sediment
!!    sedpst_rea(:) |1/day       |pesticide reaction coefficient in river bed
!!                               |sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag
!!    j           |none          |counter
!!    titldum     |NA            |title line in .wq file (not used)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

   character (len=80) :: titldum
   integer :: eof, j

   eof = 0
   j = irch

   do
      read (104,5100,iostat=eof) titldum
      if (eof < 0) exit
      read (104,5100,iostat=eof) titldum
      if (eof < 0) exit
      read (104,*,iostat=eof) rs1(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs2(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs3(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs4(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs5(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs6(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs7(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk1(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk2(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk3(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk4(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk5(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk6(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc1(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc2(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc3(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc4(j)
      if (eof < 0) exit
      read (104,5100,iostat=eof) titldum
      if (eof < 0) exit
      read (104,*,iostat=eof) chpst_rea(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) chpst_vol(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) chpst_koc(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) chpst_stl(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) chpst_rsp(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) chpst_mix(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) sedpst_conc(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) sedpst_rea(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) sedpst_bry(j)
      if (eof < 0) exit
      read (104,*,iostat=eof) sedpst_act(j)
      if (eof < 0) exit
!      read (104,*,iostat=eof) biofilm_mumax(j)
!      if (eof < 0) exit
!      read (104,*,iostat=eof) biofilm_kinv(j)
!      if (eof < 0) exit
!      read (104,*,iostat=eof) biofilm_klw(j)
!      if (eof < 0) exit
!      read (104,*,iostat=eof) biofilm_kla(j)
!      if (eof < 0) exit
!      read (104,*,iostat=eof) biofilm_cdet(j)
!      if (eof < 0) exit
!      read (104,*,iostat=eof) biofilm_bm(j)
      exit
   end do

!!    set default values for undefined parameters
   if (rs1(j) <= 0.) rs1(j) = 1.0
   if (rs2(j) <= 0.) rs2(j) = 0.05
   if (rs3(j) <= 0.) rs3(j) = 0.5
   if (rs4(j) <= 0.) rs4(j) = 0.05
   if (rs5(j) <= 0.) rs5(j) = 0.05
   if (rs6(j) <= 0.) rs6(j) = 2.5
   if (rs7(j) <= 0.) rs7(j) = 2.5
   if (rk1(j) <= 0.) rk1(j) = 1.71
   if (rk2(j) <= 0.) rk2(j) = 1.0    ! previous 50.0
   if (rk4(j) <= 0.) rk4(j) = 2.0
   if (rk5(j) <= 0.) rk5(j) = 2.0
   if (rk6(j) <= 0.) rk6(j) = 1.71
   if (bc1(j) <= 0.) bc1(j) = 0.55
   if (bc2(j) <= 0.) bc2(j) = 1.1
   if (bc3(j) <= 0.) bc3(j) = 0.21
   if (bc4(j) <= 0.) bc4(j) = 0.35
   if (chpst_rea(j) <= 1.e-6) chpst_rea(j) = 0.007
   if (chpst_vol(j) <= 1.e-6) chpst_vol(j) = 0.01
   if (chpst_koc(j) <= 1.e-6) chpst_koc(j) = 0.
   if (chpst_stl(j) <= 1.e-6) chpst_stl(j) = 1.
   if (chpst_rsp(j) <= 1.e-6) chpst_rsp(j) = 0.002
   if (chpst_mix(j) <= 1.e-6) chpst_mix(j) = 0.001
   if (sedpst_conc(j) <= 1.e-6) sedpst_conc(j) = 0.
   if (sedpst_rea(j) <= 1.e-6) sedpst_rea(j) = 0.05
   if (sedpst_bry(j) <= 1.e-6) sedpst_bry(j) = 0.002
   if (sedpst_act(j) <= 1.e-6) sedpst_act(j) = 0.030

!!  set default values for mike van liew
   if (bc1(j) <= 0.) bc1(j) = bc1_bsn
   if (bc2(j) <= 0.) bc2(j) = bc2_bsn
   if (bc3(j) <= 0.) bc3(j) = bc3_bsn
   if (bc4(j) <= 0.) bc4(j) = bc4_bsn
!!  set default values for mike van liew

!! change units from day to hour if hourly routing is performed
   if (ievent > 0) then
      rs1(j) = rs1(j) / 24.
      rs2(j) = rs2(j) / 24.
      rs3(j) = rs3(j) / 24.
      rs4(j) = rs4(j) / 24.
      rs5(j) = rs5(j) / 24.
      rk1(j) = rk1(j) / 24.
      rk2(j) = rk2(j) / 24.
      rk3(j) = rk3(j) / 24.
      rk4(j) = rk4(j) / 24.
      bc1(j) = bc1(j) / 24.
      bc2(j) = bc2(j) / 24.
      bc3(j) = bc3(j) / 24.
      bc4(j) = bc4(j) / 24.
   end if

   close (104)
   return
5100 format (a)
end
