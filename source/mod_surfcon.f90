module mod_surfcon
    use mod_atparam

    implicit none

    private
    public fmask, fmask1, phi0, phis0, alb0, swcap, swwil, sd2sc

    ! Land-sea masks (initial. in INBCON)
    ! Original (fractional) land-sea mask
    real :: fmask(ix,il)
    ! Model-defined land fraction
    real :: fmask1(ix,il)
									
    ! Time invariant surface fields 
    ! (initial. in INBCON, phis0 initial. in INVARS)
    ! Unfiltered surface geopotential
    real :: phi0(ix,il)

    ! Spectrally-filtered sfc. geopotential
    real :: phis0(ix,il)

    ! Bare-land annual-mean albedo
    real :: alb0(ix,il)

    ! Soil moisture parameters
    ! Soil wetness at field capacity (volume fraction)
    real :: swcap = 0.30 

    ! Soil wetness at wilting point  (volume fraction)
    real :: swwil = 0.17

    ! Snow depth (mm water) corresponding to snow cover = 1
    real :: sd2sc = 60.0
end module
