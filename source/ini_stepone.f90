subroutine stepone
    ! subroutine stepone
    !
    ! purpose : call initialization of semi-implicit scheme
    !           and perform initial time step

    use mod_tsteps, only: delt, delt2

    implicit none

    real :: delth

    delth = 0.5 * delt

    call impint(delth)

    call step(1, 1, delth)

    call impint(delt)

    call step(1, 2, delt)

    call impint(delt2)
end
