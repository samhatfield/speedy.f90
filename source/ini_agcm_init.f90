! Initialization of atmospheric model and coupling interface
subroutine agcm_init()
    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, model_datetime, start_datetime, end_datetime
    use coupler, only: initialize_coupler

    implicit none

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
    if (icsea >= 4) isstan = 1

    ! Initialization of atmospheric model constants and variables
    call ini_atm

    ! Initialization of coupled modules (land, sea, ice)
    call initialize_coupler

    ! Set up the forcing fields for the first time step
    call fordate(0)

    ! Do the initial (2nd-order) time step, initialize the semi-implicit scheme
    call stepone
end subroutine
