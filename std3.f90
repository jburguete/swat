!> @file std3.f90
!> file containing the subroutine std3
!> @author
!> modified by Javier Burguete

!> this subroutine writes the annual table header to the standard output
!> file
subroutine std3

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    curyr         |none        |current year in simulation
!!    isproj        |none        |special project code:
!!                               |1 test rewind (run simulation twice)
!!    prog          |NA          |program name and version
!!    title         |NA          |title lines from file.cio
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


   use parm
   implicit none

!!    standard output file
   if (iscen == 1) then
      write (26,1000) prog
      write (26,1100) title
      write (26,1200) curyr
      write (26,1300)
   else if (isproj == 1) then
      write (19,1000) prog
      write (19,1100) title
      write (19,1200) curyr
      write (19,1300)
   endif

   return
1000 format ('1',/t5,a80,t105,2(a2,'/'),a2,5x,2(i2,':'),i2)
1100 format (/(t5,20a4))
1200 format ('Annual Summary for Watershed in year ',i4,&
   &' of simulation',/)
1300 format ('UNIT',t41,'PERCO',t50,'TILE',t81,'WATER',t91,'SED',t99,&
   &'NO3',t107,'NO3',t115,'NO3',t123,'NO3',t134,'N',t142,'P',t150,'P',&
   &/,'TIME',t10,'PREC',t18,'SURQ',t26,'LATQ',t35,'GWQ',t42,'LATE',&
   &t53,'Q',t59,'SW',T68,'ET',t75,'PET',t81,'YIELD',t89,'YIELD',t98,&
   &'SURQ',t106,'LATQ',t114,'PERC',t122,'CROP',t128,'ORGANIC',t136,&
   &'SOLUBLE',t144,'ORGANIC',t152,'TILENO3'/,t10,'(mm)',t18,'(mm)',&
   &t26,'(mm)',t34,'(mm)',t42,'(mm)',t50,'(mm)',t58,'(mm)',t66,'(mm)',&
   &t74,'(mm)',t82,'(mm)',t89,'(t/ha)',t97'------------------(kg nutri&
   &ent/ha)--------------------',t152,'(kg/ha)')
end
