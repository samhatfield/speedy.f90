program agcm_main
    use mod_tsteps, only: nsteps, alph, delt2, ihout, nstrad
    use mod_date, only: model_datetime, model_step, newdate, ndaytot
	use mod_lflags, only: lradsw

    implicit none

    integer :: jday, jj

    ! Initialization
    call agcm_init()

    write (*,'(A28, I5)') 'Integration length in days: ', ndaytot

	! Initialize time step counter
	model_step = 1

    ! Do loop over total number of integration days
    do jday = 1, ndaytot
        if (model_datetime%day == 1) write(*,'(A14, I4, I0.2, I0.2)') 'Current date: ', &
			& model_datetime%year, model_datetime%month, model_datetime%day

		! Set forcing terms according to date
        call fordate(1)

		! Set daily-average flux arrays to zero
        call dmflux(0)

		! Integrate the atmospheric model for 1 day
	    do jj = 1, nsteps
            ! Set logical flags
            lradsw = (mod(model_step,nstrad) == 1)

            ! Perform one leapfrog time step
            call step(2, 2, delt2, alph)

            ! Do diagnostic, post-processing and I/O tasks
            call diagns(2, model_step)

			! Increment time step counter
            model_step = model_step + 1

	        ! Increment model datetime
	        call newdate(1)

			! Gridded data output every 6 hours
	        if (ihout .and. mod(model_step-1,nsteps/4) == 0) call iogrid (4)
	    end do

        ! Exchange data with coupler
        call agcm_to_coupler(jday)
        call coupler_to_agcm(jday)
    enddo
end
