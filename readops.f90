!> @file readops.f90
!> file containing the subroutine readops
!> @author
!> modified by Javier Burguete

!> this subroutine reads data from the HRU/subbasin management input file
!> (.mgt). This file contains data related to management practices used in
!> the HRU/subbasin.
subroutine readops

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name            |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    vfsch(:)    |none          |Fraction of flow entering the most concentrated 10% of the VFS.
!!                               |which is fully channelized
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    day         |none          |day operation occurs
!!    eof         |none          |end of file flag (=-1 if eof, else = 0)
!!    iyear
!!    j
!!    k
!!    l
!!    mgt_op      |none          |operation code number
!!                               |0 end of rotation year
!!                               |1 plant/beginning of growing season
!!    mgt         |none          |array of management parameters out of .mgt
!!                               |file (definition changes depending on
!!                               |mgt_op)
!!    mon
!!    titldum     |NA            |title line from input dataset
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: Jdt

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   real*8, dimension(23) :: mgt
   character (len=80) :: titldum
   integer :: day, eof, iyear, j, k, l, mgt_op, mon
   integer :: Jdt

   read (111,5000,end=999) titldum

!!      read scheduled operations

   j = 0
   k = ihru

   do
      mon = 0
      day = 0
      iyear = 0
      mgt_op = 0
      mgt = 0.

      read (111,5200,iostat=eof) mon, day, iyear, mgt_op, (mgt(l), l = 1, 23)
      if (eof < 0) exit
      j = j + 1
      iopday(j,k) = Jdt (ndays,day,mon)
      iopyr(j,k) = iyear
      mgt_ops(j,k) = mgt_op

      select case (mgt_op)

       case (0)  !! end of operations file

       case (1)  !! terracing
         terr_p(j,k) = mgt(4)
         terr_cn(j,k) = mgt(5)
         terr_sl(j,k) = mgt(6)

       case (2)  !! tile drainage
         drain_d(j,k) = mgt(4)
         drain_t(j,k) = mgt(5)
         drain_g(j,k) = mgt(6)
         if (mgt(7) < 1.e-6) mgt(7) = 6000.0
         drain_idep(j,k) = mgt(7)

       case (3)  !! contouring
         cont_cn(j,k) = mgt(4)
         cont_p(j,k) =  mgt(5)

       case (4)  !! filter
         filter_i(j,k) = mgt(1)  !! on off flag
         filter_ratio(j,k) = mgt(3)
         filter_con(j,k) = mgt(4)
         filter_ch(j,k) = mgt(5)

       case (5)  !! strip cropping
         strip_n(j,k) = mgt(4)
         strip_cn(j,k) = mgt(5)
         strip_c(j,k) = mgt(6)
         strip_p(j,k) = mgt(7)

       case (6)  !! fire
         fire_cn(j,k) = mgt(4)

       case (7)  !! grass waterways
         gwati(j,k) = mgt(1)
         gwatn(j,k) = mgt(3)
         gwatspcon(j,k) = mgt(4)
         gwatd(j,k) = mgt(5)
         gwatw(j,k) = mgt(6)
         gwatl(j,k) = mgt(7)
         gwats(j,k) = mgt(8)
!! Set defaults
!! Mannings via Fiener, 2006
         if (gwatn(j,k) <=0.) gwatn(j,k) = 0.35
!! length based on one side of a square HRU
         if (gwatl(j,k) <=0.) gwatl(j,k) = Sqrt(hru_km(k))
!! default to a workable depth based on with and 8:1 sideslopes
         if (gwatd(j,k) <= 0.) then
            gwatd(j,k) = 3. / 64. * gwatw(j,k)
         end if
!! Default to 3/4 HRU slope
         if (gwats(j,k) <=0.) gwats(j,k) = hru_slp(k)*.75
!! default sed coeff to 0.005
         if (gwatspcon(j,k) <= 0.) gwatspcon(j,k) = 0.005


       case (8) !! plant parameter update
         cropno_upd(j,k) = mgt(1)
         hi_upd(j,k) = mgt(4)
         laimx_upd(j,k) = mgt(5)


         !! case (9) !! Residue Managment  mjw
         !! Force residue to a minimum value regardless of tillage.  mjw
         !! so_res_flag(j,k) = mgt(1) !!mjw
         !! so_res(j,k) = mgt(4) !!mjw

       case (10) !! Generic Conservation Practice  mjw
         !! Get user defined removal eff and use these  mjw
         ro_bmp_flag(j,k) = Int(mgt(1)) !! Flag to turn on or off user BMP

         !! surface runoff removal efficiency
         ro_bmp_flo(j,k) = mgt(3)    !! Flow
         ro_bmp_sed(j,k) = mgt(4)    !! Sediment
         ro_bmp_pp(j,k) = mgt(5)     !! Particulate P
         ro_bmp_sp(j,k) = mgt(6)     !! Soluble P
         ro_bmp_pn(j,k) = mgt(7)     !! Particulate N
         ro_bmp_sn(j,k) = mgt(8)     !! Soluble N
         ro_bmp_bac(j,k) = mgt(9)    !! Bacteria
         !! subsurface - lateral soil and groundwater
         ro_bmp_flos(j,k) = mgt(10)  !! Flow
         ro_bmp_seds(j,k) = mgt(11)  !! Sediment
         ro_bmp_pps(j,k) = mgt(12)   !! Particulate P
         ro_bmp_sps(j,k) = mgt(13)   !! Soluble P
         ro_bmp_pns(j,k) = mgt(14)   !! Particulate N
         ro_bmp_sns(j,k) = mgt(15)   !! Soluble N
         ro_bmp_bacs(j,k) = mgt(16)  !! Bacteria
         !! tile flow removal efficiency
         ro_bmp_flot(j,k) = mgt(17)  !! Flow
         ro_bmp_sedt(j,k) = mgt(18)  !! Sediment
         ro_bmp_ppt(j,k) = mgt(19)   !! Particulate P
         ro_bmp_spt(j,k) = mgt(20)   !! Soluble P
         ro_bmp_pnt(j,k) = mgt(21)   !! Particulate N
         ro_bmp_snt(j,k) = mgt(22)   !! Soluble N
         ro_bmp_bact(j,k) = mgt(23)  !! Bacteria

      end select
   end do

   close (111)

   iops = j

999 return
5000 format (a)
5200 format (1x,i2,1x,i2,5x,i4,1x,i2,1x,i4,1x,i3,1x,f6.2,1x,f12.5,1x,&
   &f6.2,1x,f11.5,1x,f8.2,1x,f6.2,1x,16f5.2)
end
