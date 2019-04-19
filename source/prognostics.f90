! Prognostic spectral variables for model dynamics, and geopotential.
module prognostics
    use mod_atparam, only: mx, nx, kx, ntr

    implicit none

    private
    public vor, div, t, ps, tr
    public phi, phis

    ! Prognostic spectral variables (updated in step)
    complex :: vor(mx,nx,kx,2)    ! Vorticity
    complex :: div(mx,nx,kx,2)    ! Divergence
    complex :: t(mx,nx,kx,2)      ! Absolute temperature
    complex :: ps(mx,nx,2)        ! Log of (normalised) surface pressure (p_s/p0)
    complex :: tr(mx,nx,kx,2,ntr) ! Tracers (tr(1): specific humidity in g/kg)

    ! Geopotential (updated in geop)
    complex :: phi(mx,nx,kx) ! Atmospheric geopotential
    complex :: phis(mx,nx)   ! Surface geopotential
end module
