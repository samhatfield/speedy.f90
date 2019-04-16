subroutine stepone
    ! subroutine stepone
    !
    ! purpose : call initialization of semi-implicit scheme
    !           and perform initial time step

    use mod_tsteps, only: delt, delt2, alph, rob, wil

    implicit none

    real :: delth

    delth = 0.5 * delt

    call impint(delth, alph)

    call step(1, 1, delth, alph, rob, wil)

    call impint(delt, alph)

    call step(1, 2, delt, alph, rob, wil)

    call impint(delt2, alph)
end
