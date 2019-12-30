!> @file fcgd.f90
!> file containing the function fcgd
!> @author
!> modified by Javier Burguete

real*8 function fcgd(xx)
   implicit none
   real*8, intent (in) :: xx
   real*8, parameter :: tn = -5., top = 35., tx = 50.
   real*8 :: qq
   qq = (tn - top)/(top - tx)
   fcgd = ((xx - tn)**qq) * (tx - xx) / (((top - tn)**qq) * (tx -top ))
   if (fcgd < 0.) fcgd = 0.
end function


