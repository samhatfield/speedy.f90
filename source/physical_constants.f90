module physical_constants
    use params

    implicit none

    private
    public rearth, omega, grav, akap, rgas
    public p0, rd, cp, alhc, alhs, sbc
    public sigl, sigh, grdsig, grdscp, wvi

    ! Physical constants for dynamics
    real, parameter :: rearth = 6.371e+6
    real, parameter :: omega  = 7.292e-05
    real, parameter :: grav   = 9.81
    real, parameter :: akap   = 2.0/7.0
    real, parameter :: rgas   = akap*1004.0

    ! Physical constants for thermodynamics
    real, parameter :: p0   = 1.e+5   ! Reference pressure
    real, parameter :: rd   = 287.0   ! Gas constant for dry air
    real, parameter :: cp   = 1004.0  ! Specific heat at constant pressure
    real, parameter :: alhc = 2501.0  ! Latent heat of condensation, in J/g for consistency with
                                      ! specific humidity in g/Kg
    real, parameter :: alhs = 2801.0  ! Latent heat of sublimation
    real, parameter :: sbc  = 5.67e-8 ! Stefan-Boltzmann constant

    !   Functions of sigma and latitude (initial. in INPHYS)
    real, dimension(kx)   :: sigl   ! Logarithm of full-level sigma
    real, dimension(0:kx) :: sigh   ! Half-level sigma
    real, dimension(kx)   :: grdsig ! g/(d_sigma p0) : to convert fluxes of u,v,q into d(u,v,q)/dt
    real, dimension(kx)   :: grdscp ! g/(d_sigma p0 c_p): to convert energy fluxes into dT/dt
    real, dimension(kx,2) :: wvi    ! Weights for vertical interpolation
end module
