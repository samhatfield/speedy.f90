module params
    implicit none

    private
    public trunc, ix, iy, il, kx, nx, mx, ntr
    public nsteps, nstdia, nsteps_out, iseasc, nstrad, sppt_on, issty0, delt, rob, wil, alph

    ! Model geometry parameters
    integer, parameter :: trunc = 30         ! Spectral truncation
    integer, parameter :: ix = 96            ! Number of longitudes
    integer, parameter :: iy = 24, il = 2*iy ! Number of latitudes in hemisphere and full sphere
    integer, parameter :: kx = 8             ! Number of vertical levels
    integer, parameter :: nx = trunc+2       ! Number of total wavenumbers
    integer, parameter :: mx = trunc+1       ! Number of zonal wavenumbers
    integer, parameter :: ntr = 1            ! Number of tracers

    ! Time stepping parameters
    integer, parameter :: nsteps = 36             ! Number of time steps in one day
    real, parameter    :: delt = 86400.0 / nsteps ! Time step in seconds
    real, parameter    :: rob = 0.05              ! Damping factor in Robert time filter
    real, parameter    :: wil = 0.53              ! Parameter of Williams filter
    real, parameter :: alph = 0.5                 ! Coefficient for semi-implicit computations
                                                  ! 0 -> forward step for gravity wave terms,
                                                  ! 1 -> backward implicit, 0.5 -> centered implicit

    ! Output parameters
    integer, parameter :: nstdia = 36*5        ! Period (number of steps) for diagnostic print-out
    integer, parameter :: nsteps_out = 1       ! Number of time steps between outputs

    ! Physics parameters
    integer, parameter :: iseasc = 1           ! Seasonal cycle flag (0=no, 1=yes)
    integer, parameter :: nstrad = 3           ! Period (number of steps) for shortwave radiation
    logical, parameter :: sppt_on = .false.    ! Turn on SPPT?
    integer, parameter :: issty0 = 1979        ! Starting year for SST anomaly file
end module
