module mod_flx_land
    use mod_atparam

    implicit none

    private
    public prec_l, snowf_l, evap_l, ustr_l, vstr_l, ssr_l, slr_l, shf_l, ehf_l,&
        & hflux_l

    ! Fluxes at land surface (all downward, except evaporation)
    ! Precipitation (land)
    real :: prec_l(ix*il)

    ! Snowfall (land)
    real :: snowf_l(ix*il)

    ! Evaporation (land)
    real :: evap_l(ix*il)

    ! u-wind stress (land)
    real :: ustr_l(ix*il)

    ! v-wind stress (land)
    real :: vstr_l(ix*il)

    ! Sfc short-wave radiation (land)
    real :: ssr_l(ix*il)

    ! Sfc long-wave radiation (land)
    real :: slr_l(ix*il)

    ! Sensible heat flux (land)
    real :: shf_l(ix*il)

    ! Latent heat flux (land)
    real :: ehf_l(ix*il)

    ! Net heat flux into land sfc.end module
    real :: hflux_l(ix*il)
end module
