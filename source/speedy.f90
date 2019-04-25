program speedy
    use params, only: nsteps, delt, nsteps_out, nstrad
    use date, only: model_datetime, end_datetime, newdate, datetime_equal
    use shortwave_radiation, only: compute_shortwave
    use input_output, only: output
    use coupler, only: couple_sea_land
    use initialization, only: initialize
    use time_stepping, only: step
    use diagnostics, only: check_diagnostics
    use prognostics, only: vor, div, t

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
        call step(2, 2, 2*delt)

        ! Check model diagnostics
        call check_diagnostics(vor(:,:,:,2), div(:,:,:,2), t(:,:,:,2), model_step)

        ! Increment time step counter
        model_step = model_step + 1

        ! Increment model datetime
        call newdate(1)

        ! Output
        if (mod(model_step-1, nsteps_out) == 0) call output(model_step-1)

        ! Exchange data with coupler once per day
        if (mod(model_step-1, nsteps) == 0) then
            call couple_sea_land(1+model_step/nsteps)
        end if
    enddo
end
