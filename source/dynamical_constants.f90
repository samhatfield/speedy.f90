!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 01/05/2019
!  Constants for initialization of dynamics.
module dynamical_constants
    use types, only: p

    implicit none

    private
    public gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    real(p), parameter :: gamma  = 6.0       ! Reference temperature lapse rate (-dT/dz in deg/km)
    real(p), parameter :: hscale = 7.5       ! Reference scale height for pressure (in km)
    real(p), parameter :: hshum  = 2.5       ! Reference scale height for specific humidity (in km)
    real(p), parameter :: refrh1 = 0.7       ! Reference relative humidity of near-surface air
    real(p), parameter :: thd    = 2.4       ! Max damping time (in hours) for horizontal diffusion
                                             ! (del^6) of temperature and vorticity
    real(p), parameter :: thdd   = 2.4       ! Max damping time (in hours) for horizontal diffusion
                                             ! (del^6) of divergence
    real(p), parameter :: thds   = 12.0      ! Max damping time (in hours) for extra diffusion
                                             ! (del^2) in the stratosphere
    real(p), parameter :: tdrs   = 24.0*30.0 ! Damping time (in hours) for drag on zonal-mean wind
                                             ! in the stratosphere
end module
