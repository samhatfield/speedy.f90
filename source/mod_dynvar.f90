!> @brief
!> Prognostic spectral variables for model dynamics, and geopotential.
!> Initialised in invars.
module mod_dynvar
    use mod_atparam

    implicit none

    private
    public vor, div, t, ps, tr
    public phi, phis

    ! Prognostic spectral variables (updated in step)
    ! Vorticity
    complex :: vor(MX,NX,KX,2)

    ! Divergence 
    complex :: div(MX,NX,KX,2)

    ! Absolute temperature
    complex :: t(MX,NX,KX,2)

    ! Log of (norm.) sfc pressure (p_s/p0)
    complex :: PS(MX,NX,2)

    ! Tracers (tr.1: spec. humidity in g/kg)
    complex :: TR(MX,NX,KX,2,NTR)

    ! Geopotential (updated in geop)
    ! Atmos. geopotential
    complex :: PHI(MX,NX,KX)

    ! Surface geopotential
    complex :: PHIS(MX,NX)
end module
