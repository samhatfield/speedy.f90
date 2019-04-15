! Call initialization routines for all model common blocks
subroutine ini_atm()
    use mod_tsteps, only: ihout

    implicit none

    ! Initialize ffts
    call inifft

    ! Initialize dynamical constants and operators
    call indyns

    ! Initialize constants for physical parametrization
    call inphys

    ! Initialize forcing fields (boundary cond. + random forcing)
    call inbcon

    ! Initialize model variables
    call invars

    ! Initialize time-mean arrays for surface fluxes and output fields
    call dmflux(0)

    ! Create control file for 6-hourly output
    call iogrid(5)

    ! Output files for grid-point fields
    if (ihout .eqv. .false.) call setgrd(0)

    ! Write initial data
    if (ihout) call iogrid(4)
end
