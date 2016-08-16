module mod_cplcon_sea
    use mod_atparam

    implicit none

    private
    public rhcaps, rhcapi, cdsea, cdice, beta

    ! Constant parameters and fields in sea/ice model
    ! 1./heat_capacity (sea)
    real :: rhcaps(ix,il)

    ! 1./heat_capacity (ice)
    real :: rhcapi(ix,il)

    ! 1./dissip_time (sea)
    real :: cdsea(ix,il)

    ! 1./dissip_time (ice)
    real :: cdice(ix,il)

    ! Heat flux coef. at sea/ice int.
    real :: beta = 1.0
end module
