!> @file schedule_ops.f90
!> file containing the subroutine schedule_ops
!> @author
!> modified by Javier Burguete

!> this subroutine controls the simulation of the land phase of the
!> hydrologic cycle
!> @param[in] j HRU number
subroutine schedule_ops(j)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j              |none          |HRU number
!!    iida           |julian date   |day being simulated (current julian date)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    b
!!    iops
!!    jj          |none          |counter
!!    sin_sl
!!    tch
!!    tover
!!    xm
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    INTRINSIC: Exp, Sin, Atan
!!    SWAT: curno, ttcoef_wway

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: j
   real*8 :: b, sin_sl, tch, tover, xm
   integer :: iops, jj
 
   iops = ioper(j)

   do while (iida == iopday(iops,j) .and. iyr == iopyr(iops,j))

      select case (mgt_ops(iops,j))

       case (0)

       case (1)
         xm = 0.6 * (1. - Exp(-35.835 * hru_slp(j)))
         sin_sl = Sin(Atan(hru_slp(j)))
         usle_ls(j) = (terr_sl(iops,j)/22.128) ** xm * (65.41 *&
            &sin_sl * sin_sl + 4.56 * sin_sl + .065)
         usle_mult(j) = sol_rock(1,j) * usle_k(j) *&
            &terr_p(iops,j) * usle_ls(j) * 11.8
         if (terr_cn(iops,j) > 1.e-6) then
            call curno(terr_cn(iops,j),j)
         end if

       case (2)
         ddrain(j) = drain_d(iops,j)
         gdrain(j) = drain_g(iops,j)
         tdrain(j) = drain_t(iops,j)
         dep_imp(j) = drain_idep(iops,j)

         !!    define soil layer that the drainage tile is in
         if (ddrain(j) > 0) then
            do jj = 1, sol_nly(j)
               if (ddrain(j) < sol_z(jj,j)) ldrain(j) = jj
               if (ddrain(j) < sol_z(jj,j)) exit
            end do
         else
            ldrain(j) = 0
         endif
!!    setting tile lage time
         if (ldrain(j) > 0 .and. gdrain(j) > 0.01) then
            tile_ttime(j) = 1. - Exp(-24. / gdrain(j))
         else
            tile_ttime(j) = 0.
         end if

       case (3)
         usle_mult(j) = usle_mult(j) * cont_p(iops,j) /&
            &usle_p(j)
         call curno(cont_cn(iops,j),j)

       case (4) !! filter strip
         vfsi(j) = filter_i(iops,j) !! on off flag
         vfsratio(j) = filter_ratio(iops,j)
         vfscon(j) = filter_con(iops,j)
         vfsch(j) = filter_ch(iops,j)

!! Set some defaults if needed
         if (vfsratio(j) <= 0.) vfsratio(j) = 0.
!! minimum value for vfsratio is 0 max is 300
         if (vfsratio(j) <= 0.) vfsratio(j) = 0.
         if (vfsratio(j) > 300) vfsratio(j) = 300
!! minimum value for vfscon is 0.1 default is 0.5 max is 0.95
         if (vfscon(j) <= 0) vfscon(j) = 0.5
         if (vfscon(j) <= 0.1) vfscon(j) = 0.1
         if (vfscon(j) > 0.95) vfscon(j) = 0.95
!! minimum value for vfsch is 0 max is .95
         if (vfsch(j) <= 0.) vfsch(j) = 0.
         if (vfsch(j) > .95) vfsch(j) = .95


       case (5)
         call curno(strip_cn(iops,j),j)
         usle_mult(j) = usle_mult(j) * strip_p(iops,j) / usle_p(j)
         tover = .0556 * (slsubbsn(j) * strip_n(iops,j)) ** .6 /&
            &hru_slp(j) ** .3

         tconc(j) = tconc(j) + tover - t_ov(j)

       case (6)
         call curno (fire_cn(iops,j),j)

       case (7)
         !if (ngrwat(j) < 0)  ngrwat(j) = 0 ! not used
         !ngrwat(j) = ngrwat(j) + 1 ! not used
         grwat_n(j) = gwatn(iops,j)
         grwat_i(j) = gwati(iops,j)
         grwat_l(j) = gwatl(iops,j)
         grwat_w(j) = gwatw(iops,j)
         grwat_d(j) = gwatd(iops,j)
         grwat_s(j) = gwats(iops,j)
         grwat_spcon(j) = gwatspcon(iops,j)
!! Calculate time of concentration for waterway similar to hydroinit.f
         tch = .62 * grwat_l(j) * grwat_n(j) ** .6 /&
            &(hru_km(j) ** .125 * grwat_s(j) ** .375)
         tc_gwat(j) = tch + t_ov(j)
!! Set counter
         !k = mhru + ngrwat(j) ! not used
!!Check the channel to make sure the enter width and depth will work with 8:1 trap channel, assume width is correct
         b = (grwat_w(j) - 2.) * grwat_d(j) * 8
!! Depth and Width not possible with 8:1 sideslope and trapazoidal channel assume b =.25*width
         if (b <= 0.) grwat_d(j) = 3. / 64. * grwat_w(j)
         call ttcoef_wway(j)

       case (8)
         jj = cropno_upd(iops,j)
         blai(jj) = laimx_upd(iops,j)
         hvsti(jj) = hi_upd(iops,j)

       case (9)
         !! Implement Residue Management MJW
         !! so_res and so_res_flag are not initialized, assuming 0 values
         !if (so_res_flag(iops,j) == 1)  then
         !   min_res(j) = so_res(iops,j)
         !else
         !   min_res(j) = 0.
         !end if
         !min_res(j) = 0. ! not used (it seems to be always 0)

       case (10) !! User defined Upland CP removal MJW
         if (ro_bmp_flag (iops,j) == 1) then
            !bmp_flag(j) = 1 ! not used
            !! surface
            bmp_flo(j) = (1. - ro_bmp_flo(iops,j) / 100.)    !! Surface Flow
            bmp_sed(j) = (1. - ro_bmp_sed(iops,j) / 100.)    !! Sediment
            bmp_pp(j) = (1. - ro_bmp_pp(iops,j) / 100.)      !! Particulate P
            bmp_sp(j) = (1. - ro_bmp_sp(iops,j) / 100.)      !! Soluble P
            bmp_pn(j) =  (1. - ro_bmp_pn(iops,j) / 100.)     !! Particulate N
            bmp_sn(j) = (1. - ro_bmp_sn(iops,j) / 100.)      !! Soluble N
            bmp_bac(j) = (1. - ro_bmp_bac(iops,j) / 100.)    !! Bacteria
            !! subsurface
            bmp_flos(j) = (1. - ro_bmp_flos(iops,j) / 100.)    !! Subsurface Flow
            bmp_seds(j) = (1. - ro_bmp_seds(iops,j) / 100.)    !! Sediment
            bmp_pps(j) = (1. - ro_bmp_pps(iops,j) / 100.)      !! Particulate P
            bmp_sps(j) = (1. - ro_bmp_sps(iops,j) / 100.)      !! Soluble P
            bmp_pns(j) =  (1. - ro_bmp_pns(iops,j) / 100.)     !! Particulate N
            bmp_sns(j) = (1. - ro_bmp_sns(iops,j) / 100.)      !! Soluble N
            bmp_bacs(j) = (1. - ro_bmp_bacs(iops,j) / 100.)    !! Bacteria
            !! tile
            bmp_flot(j) = (1. - ro_bmp_flot(iops,j) / 100.)    !! Tile Flow
            bmp_sedt(j) = (1. - ro_bmp_sedt(iops,j) / 100.)    !! Sediment
            bmp_ppt(j) = (1. - ro_bmp_ppt(iops,j) / 100.)      !! Particulate P
            bmp_spt(j) = (1. - ro_bmp_spt(iops,j) / 100.)      !! Soluble P
            bmp_pnt(j) =  (1. - ro_bmp_pnt(iops,j) / 100.)     !! Particulate N
            bmp_snt(j) = (1. - ro_bmp_snt(iops,j) / 100.)      !! Soluble N
            bmp_bact(j) = (1. - ro_bmp_bact(iops,j) / 100.)    !! Bacteria
         else
            !bmp_flag(j) = 0 ! not used
            !! surface
            bmp_flo(j) = 1.      !! Surface Flow
            bmp_sed(j) = 1.      !! Sediment
            bmp_pp(j) = 1.       !! Particulate P
            bmp_sp(j) = 1.       !! Soluble P
            bmp_pn(j) =  1.      !! Particulate N
            bmp_sn(j) = 1.       !! Soluble N
            bmp_bac(j) = 1.      !! Bacteria
            !! subsurface
            bmp_flos(j) = 1.      !! Subsurface Flow
            bmp_seds(j) = 1.      !! Sediment
            bmp_pps(j) = 1.       !! Particulate P
            bmp_sps(j) = 1.       !! Soluble P
            bmp_pns(j) =  1.      !! Particulate N
            bmp_sns(j) = 1.       !! Soluble N
            bmp_bacs(j) = 1.      !! Bacteria
            !! tile
            bmp_flot(j) = 1.      !! Tile Flow
            bmp_sedt(j) = 1.      !! Sediment
            bmp_ppt(j) = 1.       !! Particulate P
            bmp_spt(j) = 1.       !! Soluble P
            bmp_pnt(j) =  1.      !! Particulate N
            bmp_snt(j) = 1.       !! Soluble N
            bmp_bact(j) = 1.      !! Bacteria

         end if

      end select

      ioper(j) = ioper(j) + 1
      iops = ioper(j)

   end do

   return
end
