module mod_sflcon
    use mod_atparam

    implicit none

    private
    public fwind0, ftemp0, fhum0, cdl, cds, chl, chs, vgust, ctday, dtheta,&
        & fstab, hdrag, fhdrag, clambda, clambsn, forog

    !  Constants for surface fluxes
    ! Ratio of near-sfc wind to lowest-level wind
    real :: fwind0 = 0.95

    ! Weight for near-sfc temperature extrapolation (0-1) :
    !          1 : linear extrapolation from two lowest levels
    !          0 : constant potential temperature ( = lowest level)
    real :: ftemp0 = 1.0

    ! Weight for near-sfc specific humidity extrapolation (0-1) :
    !            1 : extrap. with constant relative hum. ( = lowest level)
    !            0 : constant specific hum. ( = lowest level)
    real :: fhum0 = 0.0

    ! Drag coefficient for momentum over land
    real :: cdl = 2.4e-3

    ! Drag coefficient for momentum over sea
    real :: cds = 1.0e-3

    ! Heat exchange coefficient over land
    real :: chl = 1.2e-3

    ! Heat exchange coefficient over sea
    real :: chs = 0.9e-3

    ! Wind speed for sub-grid-scale gusts
    real :: vgust = 5.0

    ! Daily-cycle correction (dTskin/dSSRad)
    real :: ctday = 1.0e-2

    ! Potential temp. gradient for stability correction
    real :: dtheta = 3.0

    ! Amplitude of stability correction (fraction)
    real :: fstab = 0.67

    ! Height scale for orographic correction
    real :: hdrag = 2000.0

    ! Amplitude of orographic correction (fraction)
    real :: fhdrag = 0.5

    ! Heat conductivity in skin-to-root soil layer
    real :: clambda = 7.0

    ! Heat conductivity in soil for snow cover = 1
    real :: clambsn = 7.0

    ! Time-invariant fields (initial. in SFLSET)
    real :: forog(ix*il)
end module
