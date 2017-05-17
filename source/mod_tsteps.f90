!> @brief
!> Length of the integration and time stepping constants.
module mod_tsteps
    implicit none

    private
    public nmonts, ndaysl, nsteps, nstdia, nstppr, nstout, idout, nmonrs, ihout, sixhrrun
    public iseasc, istart, iyear0, imont0, ipout, nstrad, sppt_on, nstrdf, indrdf, issty0
    public isst0, delt, delt2, rob, wil, alph

    ! Integration length in months
    integer :: nmonts = 3

    ! No. of days in the last month of int. (max=30)
    integer :: ndaysl = 0

    ! No. of time steps in one day
    integer, parameter :: nsteps = 36

    ! Period (no. of steps) for diagnostic print-out
    integer, parameter :: nstdia = 36*5

    ! Period (no. of steps) for post-processing
    integer, parameter :: nstppr = 6

    ! Period (no. of steps) for time-mean output
    integer :: nstout = -1

    ! Daily output flag (0=no, 1=basic (Z500,PREC,MSLP,TEMP0), 2=full)
    integer, parameter :: idout  = 0

    ! Period (no. of months) for restart file update
    integer, parameter :: nmonrs = 3

    ! 6-hourly output flags
    logical, parameter :: ihout = .false.
    logical, parameter :: ipout = .false.
    logical, parameter :: sixhrrun = .false.

    ! Seasonal cycle flag (0=no, 1=yes)
    integer, parameter :: iseasc = 1

    ! Start flag (0: from rest, 1: from restart file)
    integer :: istart

    ! Year of initial date (4-digit, eg 1900)
    integer :: iyear0 = 1978
    
    ! Month of initial date (1 to 12)
    integer :: imont0 = 1

    ! Period (no. of steps) for shortwave radiation 
    integer, parameter :: nstrad = 3

    ! Turn on SPPT?
    logical, parameter :: sppt_on = .false.
    
    ! Duration of random diabatic forcing ( 0 : no forcing, > 0 : no. of
    ! initial steps, < 0 : whole integration)
    integer, parameter :: nstrdf = 0

    ! Initialization index for random diabatic forcing
    integer, parameter :: indrdf = -1

    integer, parameter :: issty0 = 1854

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
