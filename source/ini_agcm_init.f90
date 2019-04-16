! Initialization of atmospheric model and coupling interface
subroutine agcm_init()
    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, ndaytot, model_datetime, start_datetime

    implicit none

    ! Read date from fort.2 file
    read (2,*) start_datetime%year
    read (2,*) start_datetime%month
    read (2,*) start_datetime%day
    read (2,*) start_datetime%hour
    read (2,*) start_datetime%minute
    model_datetime = start_datetime

    call newdate(0)

    write(*,'(A12, I4, I0.2, I0.2)') 'Start date: ', &
        & model_datetime%year, model_datetime%month, model_datetime%day

    isst0 = (start_datetime%year - issty0) * 12 + start_datetime%month

    ! Check consistency of coupling and prescribed SST anomaly flags
    if (icsea >= 4) isstan = 1

    ! Initialization of atmospheric model constants and variables
    call ini_atm()

    ! Initialization of coupled modules (land, sea, ice)
    call ini_coupler()

    ! Set up the forcing fields for the first time step
    call fordate(0)

    ! Do the initial (2nd-order) time step, initialize the semi-implicit scheme
    call stepone
end subroutine
