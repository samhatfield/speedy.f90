program speedy
    use mod_tsteps, only: nsteps, delt2, nsteps_out, nstrad
    use date, only: model_datetime, end_datetime, newdate, datetime_equal
    use shortwave_radiation, only: compute_shortwave
    use mod_output, only: output_step
    use coupler, only: couple_sea_land
    use initialization, only: initialize

    implicit none

    ! Time step counter
    integer :: model_step = 1

    ! Initialization
    call initialize

    ! Model main loop
    do while (.not. datetime_equal(model_datetime, end_datetime))
        ! Daily tasks
        if (mod(model_step-1, nsteps) == 0) then
            ! Set forcing terms according to date
            call fordate(1)

            ! Set daily-average flux arrays to zero
            call dmflux(0)
        end if

        ! Determine whether to compute shortwave radiation on this time step
        compute_shortwave = mod(model_step, nstrad) == 1

        ! Perform one leapfrog time step
        call step(2, 2, delt2)

        ! Do diagnostic, post-processing and I/O tasks
        call diagns(2, model_step)

        ! Increment time step counter
        model_step = model_step + 1

        ! Increment model datetime
        call newdate(1)

        ! Output
        if (mod(model_step-1, nsteps_out) == 0) call output_step(model_step-1)

        ! Exchange data with coupler once per day
        if (mod(model_step-1, nsteps) == 0) then
            call couple_sea_land(1+model_step/nsteps)
        end if
    enddo
end
