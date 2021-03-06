!> @file resinit.f90
!> file containing the subroutine resinit
!> @author
!> modified by Javier Burguete

!> this subroutine initializes variables for the daily simulation of the
!> channel routing command loop
!> @param[in] jres reservoir number
!> @param[in] k inflow hydrograph storage location number (none)
subroutine resinit(jres, k)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    jres         |none         |reservoir number
!!    k           |none          |inflow hydrograph storage location number
!!    ihout        |none         |outflow hydrograph storage location number
!!    mvaro        |none         |max number of variables routed through the
!!                               |reach
!!    res_sub(:)   |none         |number of subbasin reservoir is in
!!    sub_pet(:)   |mm H2O       |potential evapotranspiration for day in
!!                               |subbasin
!!    varoute(2,:) |m^3 H2O      |water
!!    varoute(3,:) |metric tons  |sediment or suspended solid load
!!    varoute(11,:)|mg pst       |pesticide in solution
!!    varoute(12,:)|mg pst       |pesticide sorbed to sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bury        |mg pst        |loss of pesticide from active sediment layer
!!                               |by burial
!!    difus       |mg pst        |diffusion of pesticide from sediment to lake
!!                               |water
!!    pet_day     |mm H2O        |potential evapotranspiration on day
!!    reactb      |mg pst        |amount of pesticide in sediment that is lost
!!                               |through reactions
!!    reactw      |mg pst        |amount of pesticide in lake water lost through
!!                               |reactions
!!    reschlao    |kg chl-a      |chlorophyll-a leaving reservoir on day
!!    resev       |m^3 H2O       |evaporation from reservoir on day
!!    resflwi     |m^3 H2O       |water entering reservoir on day
!!    resflwo     |m^3 H2O       |water leaving reservoir on day
!!    resnh3o     |kg N          |ammonia leaving reservoir on day
!!    resno2o     |kg N          |nitrite leaving reservoir on day
!!    resno3o     |kg N          |nitrate leaving reservoir on day
!!    resorgno    |kg N          |organic N leaving reservoir on day
!!    resorgpo    |kg P          |organic P leaving reservoir on day
!!    respcp      |m^3 H2O       |precipitation on reservoir for day
!!    respesti    |mg pst        |pesticide entering reservoir on day
!!    ressa       |ha            |surface area of reservoir on day
!!    ressedc     |metric tons   |net change in sediment in reservoir during day
!!    ressedi     |metric tons   |sediment entering reservoir during time step
!!    ressedo     |metric tons   |sediment leaving reservoir during time step
!!    ressep      |m^3 H2O       |seepage from reservoir on day
!!    ressolpo    |kg P          |soluble P leaving reservoir on day
!!    resuspst    |mg pst        |amount of pesticide moving from sediment to
!!                               |lake water due to resuspension
!!    setlpst     |mg pst        |amount of pesticide moving from water to
!!                               |sediment due to settling
!!    solpesti    |mg pst        |soluble pesticide entering reservoir
!!    solpesto    |mg pst        |soluble pesticide in outflow on day
!!    sorpesti    |mg pst        |sorbed pesticide entering reservoir
!!    sorpesto    |mg pst        |sorbed pesticide in outflow on day
!!    varoute(:,:)|varies        |routing storage array
!!    volatpst    |mg pst        |amount of pesticide lost from lake water
!!                               |by volatilization
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

   use parm
   implicit none

   integer, intent(in) :: jres, k
   integer :: ii

!! add transfer amounts
   do ii = 2, mvaro
      varoute(ii,k) = varoute(ii,k) + vartran(ii,jres)
      vartran(ii,jres) = 0.
   end do


!! zero flow out variables
   do ii = 1, mvaro
      varoute(ii,ihout) = 0.
   end do

!! initialize daily variables
   bury = 0.
   difus = 0.
   pet_day = 0.
   pet_day = sub_pet(res_sub(jres))
   reactb = 0.
   reactw = 0.
   reschlao = 0.
   resev = 0.
   if(ievent == 0) then  !!urban modeling by J.Jeong
      resflwi = varoute(2,k)
   else
      resflwi = hhvaroute(2,k,1)
   endif
   !!     resflwi = varoute(2,k)
   resflwo = 0.
   respcp = 0.
   resnh3o = 0.
   resno2o = 0.
   resno3o = 0.
   resorgno = 0.
   resorgpo = 0.
   respesti = 0.
   ressa = 0.
   ressedc = 0.
   if (varoute(3,k) < 1.e-6) varoute(3,k) = 0.0
   ressedi = varoute(3,k)
   ressani = varoute(23,k)
   ressili = varoute(24,k)
   resclai = varoute(25,k)
   ressagi = varoute(26,k)
   reslagi = varoute(27,k)
   resgrai = varoute(28,k)

   if (varoute(3,k) < 1.e-6) varoute(3,k) = 0.0
   if(ievent == 0) then  !!urban modeling by J.Jeong
      ressedi = varoute(3,k)
   else
      ressedi = hhvaroute(3,k,1)
   endif
   ressedo = 0.

   ressep = 0.
   ressolpo = 0.
   resuspst = 0.
   setlpst = 0.
   solpesti = varoute(11,k)
   solpesto = 0.
   sorpesti = varoute(12,k)
   sorpesto = 0.
   volatpst = 0.

   return
end
