subroutine agcm_init(cexp, inidate, ntimes, irstart, ndays)
    !   subroutine agcm_init (cexp,inidate,ntimes,irstart,
    !  &                      ndays)
    !
    !   purpose: initialization of atmos. model and coupling interface 
    !

    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, ndaytot, iyear, imonth, iday, ihour

    implicit none

    ! input (reset by input/include files if inidate = 0):
    character(len=3), intent(inout) :: cexp        ! experiment identifier
    integer, intent(in) :: inidate     ! initial date yyyymm
    integer, intent(in) :: ntimes      ! integr. length in months (< 0) or days (> 0)
    integer, intent(in) :: irstart     ! restart flag: 0 = no, > 0 = yes

    ! output:
    integer, intent(inout) :: ndays       ! total no. of integration days

    print *, ' hallo from speedy_agcm'

    ! 1. set run initial time, duration, time-stepping and coupling options
    read (2,*) istart

    ! Read date from fort.2 file
    read (2,*) iyear0
    read (2,*) imont0
    read (2,*) iday
    read (2,*) ihour
    iyear = iyear0
    imonth = imont0

    call newdate(0)

    print *, 'start date ', iyear, imonth, iday, ihour

    isst0 = (iyear0 - issty0) * 12 + imont0

    ndays = ndaytot

    ! check consistency of coupling and prescribed SST anomaly flags
    if (icsea >= 4) isstan = 1

    ! 2. initialization of atmospheric model constants and variables 
    call ini_atm(cexp)

    ! 3. initialization of coupled modules (land, sea, ice)
    call ini_coupler(istart)

    ! 4. set up the forcing fields for the first time step
    call fordate(0)

    ! 5. do the initial (2nd-order) time step, initialize the semi-impl. scheme
    call stepone
end subroutine
