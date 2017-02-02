subroutine agcm_init(cexp, inidate, ntimes, irstart, ndays)
    !   subroutine agcm_init (cexp,inidate,ntimes,irstart,
    !  &                      ndays)
    !
    !   purpose: initialization of atmos. model and coupling interface 
    !

    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, ndaytot

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
    if (inidate <= 0) then
       read (2,*) istart
       read (2,'(a3)') cexp
    else
       istart = irstart
    endif

!   if (istart /= 0) istart = 1
    if (istart == 1110) then ! restart and 6-hourly output
        istart = 1
        ihout = .true.
        ipout = .true.
    else if (istart == 1111) then ! Start with gridded data
        ihout = .true.
        ipout = .true.
    else if (istart == 1112) then ! Start with gridded data
        istart = 1111
        ihout = .true.
        ipout = .false.
    else if (istart /= 0) then
        istart = 1
        ihout = .false.
        ipout = .false.
    else
        istart = 0
        ihout = .false.
        ipout = .false.
    end if

    if (inidate > 0) then
       iyear0 = inidate/100
       imont0 = mod(inidate,100)

       if (ntimes < 0) then
          nmonts = -ntimes
          ndaysl = 0
       else
          nmonts = 0
          ndaysl = ntimes
       endif
    endif

    isst0 = (iyear0 - issty0) * 12 + imont0

    call newdate(0)

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
