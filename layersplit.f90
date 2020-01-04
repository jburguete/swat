!> @file layersplit.f90
!> file containing the subroutine layersplit
!> @author
!> modified by Javier Burguete

subroutine layersplit(dep_new, k)

   use parm
   implicit none

   real*8, intent(in) :: dep_new
   integer, intent(in) :: k
   real*8 :: xx
   integer :: flag, j, jj, n, nly

   nly = sol_nly(k)

!!    create a septic layer
!! changed all sol_zmx(k) in subroutine to dep_new 1/27/09 gsm
   flag = 0

   do j = 2, nly
      xx = abs(dep_new - sol_z(j,k))
      !! if values are within 10 mm of one another, reset boundary
      if (xx < 10.) then
         sol_z(j,k) = dep_new
         exit
      end if

      !! set a soil layer at dep_new and adjust all lower layers
      if (sol_z(j,k) > dep_new) then
         flag = 1
         sol_nly(k) = sol_nly(k) + 1
         nly = nly + 1
         jj = j + 1
         do n = nly, jj, -1
            sol_z(n,k) = sol_z(n-1,k)
            sol_bd(n,k) = sol_bd(n-1,k)
            sol_awc(n,k) = sol_awc(n-1,k)
            sol_k(n,k) = sol_k(n-1,k)
            sol_cbn(n,k) = sol_cbn(n-1,k)
            sol_rock(n,k) = sol_rock(n-1,k) !!! Armen 13 Jan 2008 MJW rev 490
            sol_clay(n,k) = sol_clay(n-1,k)
            sol_sand(n,k) = sol_sand(n-1,k) !!! Claire 2 Dec 2009 MJW rev 490
            sol_silt(n,k) = sol_silt(n-1,k) !!! Claire 2 Dec 2009 MJW rev 490
            sol_ec(n,k) = sol_ec(n-1,k)
            sol_no3(n,k) = sol_no3(n-1,k)
            sol_orgn(n,k) = sol_orgn(n-1,k)
            sol_orgp(n,k) = sol_orgp(n-1,k)
            sol_solp(n,k) = sol_solp(n-1,k)
            sol_mc(n,k) = sol_mc(n-1,k)
            sol_mn(n,k) = sol_mn(n-1,k)
            sol_mp(n,k) = sol_mp(n-1,k)
            sol_n(n,k) = sol_n(n-1,k)
            sol_ph(n,k) = sol_ph(n-1,k) !! MJW rev 490
            sol_cal(n,k) = sol_cal(n-1,k) !! MJW rev 490

         end do
         sol_z(j,k) = dep_new
      end if
      if (flag == 1) exit
   end do

   iseptic = j
end
