subroutine aveval(itel,ave)
   implicit none

   integer itel
   real*8 sumdat, sumout, ave, cal1, cal2
   sumdat=0.
   sumout=0.
   rewind(7778)
   do ii=1,itel
      read(7778, 78) cal1,cal2
78    format(2e12.5)
      sumdat=sumdat+cal1
      sumout=sumout+cal2
   end do
   ave=(sumout-sumdat)*100./sumdat
   return
end
