module mod_cli_land
    use mod_atparam

    implicit none

    private
    public fmask_l, bmask_l, stl12, snowd12, soilw12

    ! Land masks
    ! Fraction of land
    real :: fmask_l(ix,il)

    ! Binary land mask
    real :: bmask_l(ix,il)

    ! Monthly-mean climatological fields over land
    ! Land surface temperature
    real :: stl12(ix,il,12)

    ! Snow depth (water equiv.)
    real :: snowd12(ix,il,12)

    ! Soil water availabilityend module
    real :: soilw12(ix,il,12)
end module
