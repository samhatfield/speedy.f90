program agcm_main
    use mod_tsteps, only: nsteps
    use mod_date, only: ndaytot, iyear, imonth, iday

    implicit none

    integer :: jday, istep

    ! Initialization
    call agcm_init()

    write (*,'(A28, I5)') 'Integration length in days: ', ndaytot

    ! Do loop over total number of integration days
    do jday = 1, ndaytot
        if (iday == 1) write(*,'(A14, I4, I0.2, I0.2)') 'Current date: ', iyear, imonth, iday

		! Increment time step counter
        istep = 1 + (jday - 1) * nsteps

		! Set forcing terms according to date
        call fordate(1)

		! Set daily-average flux arrays to zero
        call dmflux(0)

		! Integrate the atmospheric model for 1 day
        call stloop(istep)

        ! Exchange data with coupler
        call agcm_to_coupler(jday)
        call coupler_to_agcm(jday)
    enddo
end
