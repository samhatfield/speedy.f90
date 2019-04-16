program agcm_main
    use mod_tsteps, only: nsteps, alph, delt2, ihout, nstrad
    use mod_date, only: ihour, newdate, ndaytot, iyear, imonth, iday
	use mod_lflags, only: lradsw

    implicit none

    integer :: jday, istep, jj

    ! Initialization
    call agcm_init()

    write (*,'(A28, I5)') 'Integration length in days: ', ndaytot

	! Initialize time step counter
	istep = 1

    ! Do loop over total number of integration days
    do jday = 1, ndaytot
        if (iday == 1) write(*,'(A14, I4, I0.2, I0.2)') 'Current date: ', iyear, imonth, iday

		! Set forcing terms according to date
        call fordate(1)

		! Set daily-average flux arrays to zero
        call dmflux(0)

		! Integrate the atmospheric model for 1 day
	    do jj = 1, nsteps
            ! Set logical flags
            lradsw = (mod(istep,nstrad) == 1)

            ! Perform one leapfrog time step
            call step(2, 2, delt2, alph)

            ! Do diagnostic, post-processing and I/O tasks
            call diagns(2, istep)

			! Increment time step counter
            istep = istep + 1

	        ! Increment hour timer
			if (floor(istep*24.0/nsteps) > floor((istep-1)*24.0/nsteps)) ihour = mod(ihour + 1, 24)

	        ! If it's a new day, then compute a new date
	        if (mod(istep,nsteps) == 0) call newdate(1)

			! Gridded data output every 6 hours
	        if (ihout .and. mod(istep,nsteps/4) == 0) call iogrid (4)
	    end do

        ! Exchange data with coupler
        call agcm_to_coupler(jday)
        call coupler_to_agcm(jday)
    enddo
end
