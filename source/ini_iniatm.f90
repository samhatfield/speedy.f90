! Call initialization routines for all model common blocks
subroutine ini_atm()
    use mod_output, only: output_step

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

    ! Write initial data
    call output_step(0)
end
