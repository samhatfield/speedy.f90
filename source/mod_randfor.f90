module mod_randfor
    use mod_atparam

    implicit none

    private
    public randfh, randfv

    ! Random diabatic forcing (initial. in INIRDF, modified by XS_RDF))
    real :: randfh(ix,il,2), randfv(il,kx,2)
end module
