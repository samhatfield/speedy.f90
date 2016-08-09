module mod_dyncon1
    use mod_atparam

    implicit none

    private
    public rearth, omega, grav, akap, rgas, pi, a, g
    public hsg, dhs, fsg, dhsr, fsgr
    public radang, gsin, gcos, coriol
    public xgeop1, xgeop2

    ! Physical constants for dynamics
    real, parameter :: rearth = 6.371e+6
    real, parameter :: omega  = 7.292e-05
    real, parameter :: grav   = 9.81
    real, parameter :: akap   = 2./7.
    real, parameter :: rgas   = akap*1004.
    real, parameter :: pi = 4.*atan(1.)
    real, parameter :: a  = rearth
    real, parameter :: g  = grav

    ! Vertical level parameters (initial. in indyns)
    real :: hsg(kxp), dhs(kx), fsg(kx), dhsr(kx), fsgr(kx)

    ! Functions of lat. and lon. (initial. in indyns)
    real :: radang(il), gsin(il), gcos(il), coriol(il)

    ! Constants for hydrostatic eq. (initial. in indyns)
    real :: xgeop1(kx), xgeop2(kx)
end module
