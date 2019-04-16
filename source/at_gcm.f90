program agcm_main
    use mod_tsteps, only: nsteps, alph, delt2, nsteps_out, nstrad
    use mod_date, only: model_datetime, end_datetime, newdate, datetime_equal
    use mod_lflags, only: lradsw
    use mod_output, only: output_step

    implicit none

    ! Time step counter
    integer :: model_step = 1

    ! Initialization
    call agcm_init()

    ! Model main loop
    do while (.not. datetime_equal(model_datetime, end_datetime))
        ! Daily tasks
        if (model_datetime%hour == 0 .and. model_datetime%minute == 0) then
            ! Print date once per month
            if (model_datetime%day == 1) then
            	write(*,'(A14, I4, I0.2, I0.2)') 'Current date: ', &
            	& model_datetime%year, model_datetime%month, model_datetime%day
            end if

            ! Set forcing terms according to date
            call fordate(1)

            ! Set daily-average flux arrays to zero
            call dmflux(0)
        end if

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

        ! Output
        if (mod(model_step-1, nsteps_out) == 0) call output_step(model_step)

        ! Exchange data with coupler once per day
        if (model_datetime%hour == 0 .and. model_datetime%minute == 0) then
            call agcm_to_coupler(1+model_step/nsteps)
            call coupler_to_agcm(1+model_step/nsteps)
        end if
    enddo
end
