!> @brief
!> Constants for initialization of dynamics.
module dynamical_constants
    implicit none

    private
    public gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    real, parameter :: gamma  = 6.0       ! Reference temperature lapse rate (-dT/dz in deg/km)
    real, parameter :: hscale = 7.5       ! Reference scale height for pressure (in km)
    real, parameter :: hshum  = 2.5       ! Reference scale height for specific humidity (in km)
    real, parameter :: refrh1 = 0.7       ! Reference relative humidity of near-surface air
    real, parameter :: thd    = 2.4       ! Max damping time (in hours) for horizontal diffusion
                                          ! (del^6) of temperature and vorticity
    real, parameter :: thdd   = 2.4       ! Max damping time (in hours) for horizontal diffusion
                                          ! (del^6) of divergence
    real, parameter :: thds   = 12.0      ! Max damping time (in hours) for extra diffusion (del^2)
                                          ! in the stratosphere
    real, parameter :: tdrs   = 24.0*30.0 ! Damping time (in hours) for drag on zonal-mean wind
                                          ! in the stratosphere
end module
