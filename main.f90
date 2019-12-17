!> @file main.f90
!> file containing the main program that reads input, calls the main simulation
!> model, and writes output.
!> @author
!> modified by Javier Burguete Tolosa
include 'modparm.f90'

!> this is the main program that reads input, calls the main simulation
!> model, and writes output.
program main
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!         ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    fcstcycles  |none        |number of times forecast period is simulated
!!                             |(using different weather generator seeds each
!!                             |time)
!!    fcstday     |julian date |beginning date of forecast period
!!    fcstyr      |year        |beginning year of forecast period
!!    hi_targ(:,:,:)|(kg/ha)/(kg/ha)|harvest index target of cover defined
!!                             |at planting
!!    hruno       |            |
!!    iclb        |none        |auto-calibration flag
!!    iscen       |none        |scenarios counter
!!    isproj      |none        |special project code:
!!                             |1 test rewind (run simulation twice)
!!    prog        |NA          |SWAT program header string
!!    subnum      |NA          |
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    prog        |NA           |program name and version
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    scenario    |none          |scenarios number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: date_and_time
!!    SWAT: getallo, allocate_parms, readfile, readfig
!!    SWAT: readbsn, std1, readwwq, readinpt, std2, storeinitial
!!    SWAT: openwth, headout, simulate, finalbal, writeaa, pestw
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none
   integer :: ii, scenario

   prog = "SWAT Sep 7    VER 2018/Rev 670"
   write (*,1000)
1000 format(1x,"               SWAT2018               ",/,&
   &"               Rev. 670               ",/,&
   &"      Soil & Water Assessment Tool    ",/,&
   &"               PC Version             ",/,&
   &" Program reading from file.cio . . . executing",/)

!! process input

   call getallo
   call allocate_parms
   call readfile
   call readbsn
   call readwwq
   if (fcstyr > 0 .and. fcstday > 0) call readfcst
   call readplant             !! read in the landuse/landcover database
   call readtill              !! read in the tillage database
   call readpest              !! read in the pesticide database
   call readfert              !! read in the fertilizer/nutrient database
   call readurban             !! read in the urban land types database
   call readseptwq            !! read in the septic types database
   call readlup
   call readfig
   call readatmodep
   call readinpt
   call std1
   call std2
   call openwth
   call headout

   !! convert integer to string for output.mgt file
   subnum = ""
   hruno = ""
   do ii = 1, mhru
      write (subnum(ii),fmt=' (i5.5)') hru_sub(ii)
      write (hruno(ii),fmt=' (i4.4)') hru_seq(ii)
   end do

   if (isproj == 2) then
      hi_targ = 0.0
   end if

!! save initial values
   if (isproj == 1) then
      scenario = 2
      call storeinitial
   else if (fcstcycles > 1) then
      scenario =  fcstcycles
      call storeinitial
   else
      scenario = 1
   endif
   if (iclb /= 4) then
      do iscen = 1, scenario

         !! simulate watershed processes
         call simulate

         !! perform summary calculations
         call finalbal
         call writeaa
         call pestw

         !!reinitialize for new scenario
         if (scenario > iscen) call rewind_init
      end do
   end if
   do ii = 101, 109
      close (ii)
   end do
   close(124)
   write (*,1001)
1001 format (/," Execution successfully completed ")

   iscen=1
!! file for Mike White to review to ensure simulation executed normally
   open (9999,file='fin.fin')
   write (9999,*) 'Execution successful'
   close (9999)

   stop
end
