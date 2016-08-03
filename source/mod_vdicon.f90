!> @brief
!> Constants for vertical diffusion and shallow convection.
module mod_vdicon
    implicit none

    private
    public trshc, trvdi, trvds, redshc, rhgrad, segrad

    ! Relaxation time (in hours) for shallow convection
    real, parameter :: trshc = 6.0

    ! Relaxation time (in hours) for moisture diffusion
    real, parameter :: trvdi = 24.0

    ! Relaxation time (in hours) for super-adiab. conditions
    real, parameter :: trvds = 6.0

    ! Reduction factor of shallow conv. in areas of deep conv.
    real, parameter :: redshc = 0.5

    ! Maximum gradient of relative humidity (d_RH/d_sigma)
    real, parameter :: rhgrad = 0.5

    ! Minimum gradient of dry static energy (d_DSE/d_phi)
    real, parameter :: segrad = 0.1
end module
