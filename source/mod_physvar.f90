module mod_physvar
    use mod_atparam

    implicit none

    private
    public precnv, precls, snowcv, snowls, cbmf, tsr, ssrd, ssr, slrd, slr,&
        & olr, slru, ustr, vstr, shf, evap, hfluxn

    ! precnv = convective precipitation  [g/(m^2 s)], total
    ! precls = large-scale precipitation [g/(m^2 s)], total
    ! snowcv = convective precipitation  [g/(m^2 s)], snow only
    ! snowls = large-scale precipitation [g/(m^2 s)], snow only
    ! cbmf   = cloud-base mass flux
    ! tsr    = top-of-atm. shortwave radiation (downward)
    ! ssrd   = surface shortwave radiation (downward-only)
    ! ssr    = surface shortwave radiation (net downward)
    ! slrd   = surface longwave radiation  (downward-only)
    ! slr    = surface longwave radiation  (net upward)
    ! olr    = outgoing longwave radiation (upward)
    ! slru   = surface longwave emission   (upward)
    !                                   (1:land, 2:sea, 3: wgt. average)
    ! ustr   = u-stress                 (1:land, 2:sea, 3: wgt. average)
    ! vstr   = v-stress                 (1:land, 2:sea, 3: wgt. average)
    ! shf    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
    ! evap   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
    ! hfluxn = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)
    real, dimension(ix,il) :: precnv, precls, snowcv, snowls, cbmf, tsr, ssrd,&
        & ssr, slrd, slr, olr
    real, dimension(ix,il,3) :: slru, ustr, vstr, shf, evap, hfluxn
end module
