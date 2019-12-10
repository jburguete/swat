subroutine chkcst(nopt,xi,bl,bu,ibound)
!
!     This subroutine check if the trial point satisfies all
!     constraints.
!
!     ibound - violation indicator
!            = -1 initial value
!            = 0  no violation
!            = 1  violation
!     nopt = number of optimizing variables
!     ii = the ii'th variable of the arrays x, bl, and bu
!
   implicit none
   dimension xi(nopt),bl(nopt),bu(nopt)
!
   ibound = -1
!
!     Check if explicit constraints are violated
!
   do ii=1, nopt
      if (xi(ii) .lt. bl(ii) .or. xi(ii) .gt. bu(ii)) go to 10
   end do
   if (nopt .eq. 1) go to 9
!
!     Check if implicit constraints are violated
!     (no implicit constraints for this function)
!
!     No constraints are violated
!
9  ibound = 0
   return
!
!     At least one of the constraints are violated
!
10 ibound = 1
   return
end
