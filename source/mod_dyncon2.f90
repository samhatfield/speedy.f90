module mod_dyncon2
    use mod_atparam

    implicit none

    private
    public :: tref, tref1, tref2, tref3
    public :: xa, xb, xc, xd, xe, xf, xg, xh, xj, dhsx, elz

    ! Temp. profile for semi-imp. scheme (initial. in IMPINT)
    real, dimension(kx) :: tref, tref1, tref2, tref3

    real, dimension(kx,kx) :: xa, xb, xc, xd, xe
    real, dimension(kx,kx,lmax) :: xf, xg, xh, xj
    real :: dhsx(kx), elz(mx,nx)                               
end module
