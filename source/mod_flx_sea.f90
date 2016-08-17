module mod_flx_sea
    use mod_atparam

    implicit none

    private
    public prec_s, snowf_s, evap_s, ustr_s, vstr_s, ssr_s, slr_s, shf_s, ehf_s,&
        & hflux_s, hflux_i

    ! Fluxes at sea surface (all downward, except evaporation)
    ! Precipitation (sea)
    real :: prec_s(ix*il)
    
    ! Snowfall (sea)
    real :: snowf_s(ix*il)

    ! Evaporation (sea)
    real :: evap_s(ix*il)

    ! u-wind stress (sea)
    real :: ustr_s(ix*il)

    ! v-wind stress (sea)
    real :: vstr_s(ix*il)

    ! Sfc short-wave radiation (sea)
    real :: ssr_s(ix*il)

    ! Sfc long-wave radiation (sea)
    real :: slr_s(ix*il)

    ! Sensible heat flux (sea)
    real :: shf_s(ix*il)

    ! Latent heat flux (sea)
    real :: ehf_s(ix*il)

    ! Net heat flux into sea sfc.
    real :: hflux_s(ix*il)

    ! Net heat flux into sea-ice sfc.
    real :: hflux_i(ix*il)
end module
