!> @file oxygen_saturation.f90
!> file containing the function oxygen_saturation
!> @author
!> modified by Javier Burguete

!> this function calculates saturation concentration for dissolved oxygen
!> QUAL2E section 3.6.1 equation III-29
!> @param[in] t temperature (deg C)
!> @return saturation concentration for dissolved oxygen
real*8 function oxygen_saturation(t) result(o)
   real*8, intent(in) :: t
   real*8 :: ki ! 1 / (kelvin temperature)
   ki = 1. / (t + 273.15)
   o = Max(0., Exp(-139.34410 + ki * (1.575701e05 + ki * (-6.642308e07&
      + ki * (1.243800e10 + ki * (-8.621949e11))))))
end function
