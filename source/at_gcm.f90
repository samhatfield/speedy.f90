program agcm_main
    use mod_date, only: ndaytot

    implicit none

    integer :: jday

    ! Initialization
    call agcm_init()

    print *, 'integration length in days: ', ndaytot

    ! Do loop over total number of integration days
    do jday = 1, ndaytot
        ! Run atmospheric model for 1 day
        call agcm_1day(jday)

        ! Exchange data with coupler
        call agcm_to_coupler(jday)
        call coupler_to_agcm(jday)
    enddo
end

! Perform atmospheric model integration for 1 day
subroutine agcm_1day(jday)

    use mod_tsteps, only: nsteps
    use mod_date, only: iyear, imonth, iday

    implicit none

    integer, intent(in) :: jday
    integer :: istep

    if (iday == 1) print *, ' start of year/month = ', iyear, imonth

    istep = 1 + (jday - 1) * nsteps

    ! 1. set forcing terms according to date
    call fordate(1)

    ! 2. set daily-average flux arrays to zero
    call dmflux(0)

    ! 3. integrate the atmospheric model for 1 day
    call stloop(istep)
end
