module mod_dyncon1
    use mod_atparam

    implicit none

    private
    public rearth, omega, grav, akap, rgas
    public xgeop1, xgeop2

    ! Physical constants for dynamics
    real, parameter :: rearth = 6.371e+6
    real, parameter :: omega  = 7.292e-05
    real, parameter :: grav   = 9.81
    real, parameter :: akap   = 2./7.
    real, parameter :: rgas   = akap*1004.

    ! Constants for hydrostatic eq. (initial. in indyns)
    real :: xgeop1(kx), xgeop2(kx)
end module
