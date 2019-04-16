!> @brief
!> Length of the integration and time stepping constants.
module mod_tsteps
    implicit none

    private
    public nmonts, ndaysl, nsteps, nstdia, nsteps_out
    public iseasc, nstrad, sppt_on, issty0
    public isst0, delt, delt2, rob, wil, alph

    ! Integration length in months
    integer :: nmonts = 3

    ! No. of days in the last month of int. (max=30)
    integer :: ndaysl = 0

    ! No. of time steps in one day
    integer, parameter :: nsteps = 36

    ! Period (no. of steps) for diagnostic print-out
    integer, parameter :: nstdia = 36*5

    ! Number of time steps between outputs
    integer, parameter :: nsteps_out = 1

    ! Seasonal cycle flag (0=no, 1=yes)
    integer, parameter :: iseasc = 1

    ! Period (no. of steps) for shortwave radiation
    integer, parameter :: nstrad = 3

    ! Turn on SPPT?
    logical, parameter :: sppt_on = .false.

    integer, parameter :: issty0 = 1979

    ! Record in SST anomaly file corr. to the initial month
    ! Initialized in agcm_init
    integer :: isst0

    ! Time step in seconds
    real, parameter :: delt = 86400.0 / nsteps

    ! 2 * time step in seconds
    real, parameter :: delt2 = 2 * delt

    ! Damping factor in Robert time filter
    real, parameter :: rob = 0.05

    ! Parameter of Williams filter
    real, parameter :: wil = 0.53

    ! Coefficient for semi-implicit computations
    real :: alph
end module
