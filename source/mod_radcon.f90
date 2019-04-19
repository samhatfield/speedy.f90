module mod_radcon
    use mod_atparam

    implicit none

    private
    public albsea, albice, albsn, epslw, emisfc, ablco2_ref
    public fband
    public alb_l, alb_s, albsfc, snowc
    public tau2, st4a, stratc, flux

    ! Radiation and cloud constants

    ! albsea = Albedo over sea
    ! albice = Albedo over sea ice (for ice fraction = 1)
    ! albsn  = Albedo over snow (for snow cover = 1)

    ! epslw  = fraction of blackbody spectrum absorbed/emitted by PBL only
    ! emisfc = longwave surface emissivity

    real :: albsea = 0.07
    real :: albice = 0.60
    real :: albsn  = 0.60

    real :: epslw  =  0.05
    real :: emisfc =  0.98

    real :: ablco2_ref

    ! Time-invariant fields (initial. in radset)
    ! fband  = energy fraction emitted in each LW band = f(T)
    real :: fband(100:400,4)

    ! Radiative properties of the surface (updated in fordate)
    ! alb_l  = daily-mean albedo over land (bare-land + snow)
    ! alb_s  = daily-mean albedo over sea  (open sea + sea ice)
    ! albsfc = combined surface albedo (land + sea)
    ! snowc  = effective snow cover (fraction)
    real, dimension(ix,il) :: alb_l, alb_s, albsfc, snowc

    ! Transmissivity and blackbody rad. (updated in radsw/radlw)
    ! tau2   = transmissivity of atmospheric layers
    ! st4a   = blackbody emission from full and half atmospheric levels
    ! stratc = stratospheric correction term
    ! flux   = radiative flux in different spectral bands
    real :: tau2(ix,il,kx,4), st4a(ix,il,kx,2), stratc(ix,il,2), flux(ix,il,4)
end module
