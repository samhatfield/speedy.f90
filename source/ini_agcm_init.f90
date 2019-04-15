! Initialization of atmospheric model and coupling interface
subroutine agcm_init()
    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, ndaytot, iyear, imonth, iday, ihour

    implicit none

    ! Read mode from fort.2 file
    read (2,*) istart

    ! Read date from fort.2 file
    read (2,*) iyear0
    read (2,*) imont0
    read (2,*) iday
    read (2,*) ihour
    iyear = iyear0
    imonth = imont0

    call newdate(0)

    print *, 'Start date ', iyear, imonth, iday, ihour

    isst0 = (iyear0 - issty0) * 12 + imont0

    ! Check consistency of coupling and prescribed SST anomaly flags
    if (icsea >= 4) isstan = 1

    ! Initialization of atmospheric model constants and variables
    call ini_atm()

    ! Initialization of coupled modules (land, sea, ice)
    call ini_coupler(istart)

    ! Set up the forcing fields for the first time step
    call fordate(0)

    ! Do the initial (2nd-order) time step, initialize the semi-implicit scheme
    call stepone
end subroutine
