!> @file readswq.f90
!> file containing the subroutine readswq
!> @author
!> modified by Javier Burguete

!> this subroutine reads parameters from the subbasin instream water
!> quality file (.swq) and initializes the QUAL2E variables which apply to
!> the individual subbasins
subroutine readswq

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |reach number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bc(1,:)      |1/day or 1/hr |rate constant for biological oxidation of NH3
!!                               |to NO2 in reach at 20 deg C
!!    bc(2,:)      |1/day or 1/hr |rate constant for biological oxidation of NO2
!!                               |to NO3 in reach at 20 deg C
!!    bc(3,:)      |1/day or 1/hr |rate constant for hydrolysis of organic N to
!!                               |ammonia in reach at 20 deg C
!!    bc(4,:)      |1/day or 1/hr |rate constant for the decay of organic P to
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
!!    rk(1,:)      |1/day or 1/hr |CBOD deoxygenation rate coefficient in reach
!!                               |at 20 deg C
!!    rk(2,:)      |1/day or 1/hr |reaeration rate in accordance with Fickian
!!                               |diffusion in reach at 20 deg C
!!    rk(3,:)      |1/day or 1/hr |rate of loss of CBOD due to settling in reach
!!                               |at 20 deg C
!!    rk(4,:)      |mg O2/        |sediment oxygen demand rate in reach
!!                |  ((m**2)*day)|at 20 deg C
!!                |or mg O2/((m**2)*hr)
!!    rk(5,:)      |1/day         |coliform die-off rate in reach
!!    rk(6,:)      |1/day         |decay rate for arbitrary non-conservative
!!                               |constituent in reach
!!    rs(1,:)      |m/day or m/hr |local algal settling rate in reach at 20 deg C
!!    rs(2,:)      |(mg disP-P)/  |benthos source rate for dissolved phosphorus
!!                |  ((m**2)*day)|in reach at 20 deg C
!!                |or (mg disP-P)/((m**2)*hr)|
!!    rs(3,:)      |(mg NH4-N)/   |benthos source rate for ammonia nitrogen in
!!                |  ((m**2)*day)|reach at 20 deg C
!!                |or (mg NH4-N)/((m**2)*hr)|
!!    rs(4,:)      |1/day or 1/hr |rate coefficient for organic nitrogen
!!                               |settling in reach at 20 deg C
!!    rs(5,:)      |1/day or 1/hr |organic phosphorus settling rate in reach at
!!                               |20 deg C
!!    rs(6,:)      |1/day         |rate coefficient for settling of arbitrary
!!                               |non-conservative constituent in reach
!!    rs(7,:)      |(mg ANC)/     |benthal source rate for arbitrary
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
      read (104,*,iostat=eof) rs(1,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs(2,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs(3,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs(4,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs(5,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs(6,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rs(7,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk(1,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk(2,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk(3,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk(4,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk(5,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) rk(6,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc(1,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc(2,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc(3,j)
      if (eof < 0) exit
      read (104,*,iostat=eof) bc(4,j)
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
   if (rs(1,j) <= 0.) rs(1,j) = 1.0
   if (rs(2,j) <= 0.) rs(2,j) = 0.05
   if (rs(3,j) <= 0.) rs(3,j) = 0.5
   if (rs(4,j) <= 0.) rs(4,j) = 0.05
   if (rs(5,j) <= 0.) rs(5,j) = 0.05
   if (rs(6,j) <= 0.) rs(6,j) = 2.5
   if (rs(7,j) <= 0.) rs(7,j) = 2.5
   if (rk(1,j) <= 0.) rk(1,j) = 1.71
   if (rk(2,j) <= 0.) rk(2,j) = 1.0    ! previous 50.0
   if (rk(4,j) <= 0.) rk(4,j) = 2.0
   if (rk(5,j) <= 0.) rk(5,j) = 2.0
   if (rk(6,j) <= 0.) rk(6,j) = 1.71
   if (bc(1,j) <= 0.) bc(1,j) = 0.55
   if (bc(2,j) <= 0.) bc(2,j) = 1.1
   if (bc(3,j) <= 0.) bc(3,j) = 0.21
   if (bc(4,j) <= 0.) bc(4,j) = 0.35
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
   if (bc(1,j) <= 0.) bc(1,j) = bc1_bsn
   if (bc(2,j) <= 0.) bc(2,j) = bc2_bsn
   if (bc(3,j) <= 0.) bc(3,j) = bc3_bsn
   if (bc(4,j) <= 0.) bc(4,j) = bc4_bsn
!!  set default values for mike van liew

!! change units from day to hour if hourly routing is performed
   if (ievent > 0) then
      rs(1,j) = rs(1,j) / 24.
      rs(2,j) = rs(2,j) / 24.
      rs(3,j) = rs(3,j) / 24.
      rs(4,j) = rs(4,j) / 24.
      rs(5,j) = rs(5,j) / 24.
      rk(1,j) = rk(1,j) / 24.
      rk(2,j) = rk(2,j) / 24.
      rk(3,j) = rk(3,j) / 24.
      rk(4,j) = rk(4,j) / 24.
      bc(1,j) = bc(1,j) / 24.
      bc(2,j) = bc(2,j) / 24.
      bc(3,j) = bc(3,j) / 24.
      bc(4,j) = bc(4,j) / 24.
   end if

   close (104)
   return
5100 format (a)
end
