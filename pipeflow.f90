!> @file pipeflow.f90
!> file containing the function pipeflow
!> @author
!> modified by Javier Burguete

!> this function calculates orifice pipe flow and returns flow rate (m^3/s)
!> @param[in] d diameter (mm)
!> @param[in] h depth (mm)
!> @return flow rate (m^3/s)
real*8 function pipeflow(d, h)
   real*8, intent(in):: d, h
   real*8:: dia, hdepth, area

   dia = d / 1000. !m
   hdepth = h / 1000. !m

   area = dia * dia / 4 * 3.14159 !m2
   pipeflow = 0.6 * area * Sqrt(19.6 * hdepth) !m3/s
end function

