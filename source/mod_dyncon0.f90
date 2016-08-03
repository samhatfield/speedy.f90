!> @brief
!> Constants for initialization of dynamics.
module mod_dyncon0
    implicit none

    private
    public gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    ! Ref. temperature lapse rate (-dT/dz in deg/km)
    real, parameter :: gamma = 6.0

    ! Ref. scale height for pressure (in km)
    real, parameter :: hscale = 7.5

    ! Ref. scale height for spec. humidity (in km)
    real, parameter :: hshum = 2.5

    ! Ref. relative humidity of near-surface air
    real, parameter :: refrh1 = 0.7

    ! Max damping time (in hours) for hor. diffusion (del^6) of temperature and
    ! vorticity
    real, parameter :: thd = 2.4

    ! Max damping time (in hours) for hor. diffusion (del^6)
    ! of divergence
    real, parameter :: thdd = 2.4

    ! Max damping time (in hours) for extra diffusion (del^2)
    ! in the stratosphere 
    real, parameter :: thds = 12.0

    ! Damping time (in hours) for drag on zonal-mean wind
    ! in the stratosphere 
    real, parameter :: tdrs = 24.0 * 30.0
end module
