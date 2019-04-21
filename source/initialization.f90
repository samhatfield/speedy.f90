module initialization
    implicit none

    private
    public initialize

contains
    ! Initialization of atmospheric model and coupling interface
    subroutine initialize
        use mod_tsteps
        use date, only: newdate, model_datetime, start_datetime, end_datetime
        use coupler, only: initialize_coupler
        use sea_model, only: sea_coupling_flag, sst_anomaly_coupling_flag
        use mod_spectral, only: inifft
        use physics, only: initialize_physics
        use mod_output, only: output_step
        use time_stepping, only: first_step

        ! Read start and end dates from fort.2 file
        read (2,*) start_datetime%year
        read (2,*) start_datetime%month
        read (2,*) start_datetime%day
        read (2,*) start_datetime%hour
        read (2,*) start_datetime%minute
        read (2,*) end_datetime%year
        read (2,*) end_datetime%month
        read (2,*) end_datetime%day
        read (2,*) end_datetime%hour
        read (2,*) end_datetime%minute
        model_datetime = start_datetime

        call newdate(0)

        write(*,'(A12,I4,A,I0.2,A,I0.2,A,I0.2,A,I0.2)') 'Start date: ', &
            & start_datetime%year,'/',start_datetime%month,'/',start_datetime%day,' ', &
            & start_datetime%hour,':',start_datetime%minute
        write(*,'(A12,I4,A,I0.2,A,I0.2,A,I0.2,A,I0.2)') 'End date: ', &
            & end_datetime%year,'/',end_datetime%month,'/',end_datetime%day,' ', &
            & end_datetime%hour,':',end_datetime%minute

        isst0 = (start_datetime%year - issty0) * 12 + start_datetime%month

        ! Check consistency of coupling and prescribed SST anomaly flags
        if (sea_coupling_flag >= 4) sst_anomaly_coupling_flag = 1

        ! =========================================================================
        ! Initialization of atmospheric model constants and variables
        ! =========================================================================

        ! Initialize ffts
        call inifft

        ! Initialize dynamical constants and operators
        call indyns

        ! Initialize constants for physical parametrization
        call initialize_physics

        ! Initialize forcing fields (boundary cond. + random forcing)
        call inbcon

        ! Initialize model variables
        call invars

        ! Initialize time-mean arrays for surface fluxes and output fields
        call dmflux(0)

        ! =========================================================================
        ! Initialization of coupled modules (land, sea, ice)
        ! =========================================================================

        call initialize_coupler

        ! =========================================================================
        ! Initialization of first time step
        ! =========================================================================

        ! Write initial data
        call output_step(0)

        ! Set up the forcing fields for the first time step
        call fordate(0)

        ! Do the initial (2nd-order) time step, initialize the semi-implicit scheme
        call first_step
    end subroutine
end module
