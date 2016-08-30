module mod_cli_sea
    use mod_atparam

    implicit none

    private
    public fmask_s, bmask_s, deglat_s, sst12, sice12, sstan3, hfseacl, sstom12

    ! Sea masks
    ! Fraction of sea
    real :: fmask_s(ix,il)

    ! Binary sea mask
    real :: bmask_s(ix,il)

    ! Grid latitudes
    real :: deglat_s(il)

    ! Monthly-mean climatological fields over sea
    ! Sea/ice surface temperature
    real :: sst12(ix,il,12)

    ! Sea ice fraction
    real :: sice12(ix,il,12)

    ! SST anomaly fields
    ! SST anomaly in 3 consecutive months
    real :: sstan3(ix,il,3)

    ! Climatological fields from model output
    ! Annual-mean heat flux into sea sfc.
    real :: hfseacl(ix,il)

    ! Ocean model SST climatology
    real :: sstom12(ix,il,12)
end module
