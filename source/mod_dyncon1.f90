module mod_dyncon1
    use mod_atparam

    implicit none

    private
    public xgeop1, xgeop2

    ! Constants for hydrostatic eq. (initial. in indyns)
    real :: xgeop1(kx), xgeop2(kx)
end module
