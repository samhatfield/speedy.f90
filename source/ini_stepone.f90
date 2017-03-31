subroutine stepone
    ! subroutine stepone
    !
    ! purpose : call initialization of semi-implicit scheme 
    !           and perform initial time step

    use mod_tsteps, only: delt, delt2, alph, rob, wil, istart

    implicit none

    integer :: iitest = 1
    real :: delth

    if (iitest == 1) print *, ' instep: initial time step'

    if (istart == 0 .or. istart == 2) then

      delth = 0.5 * delt

      if (iitest == 1) print *, ' semi-impl. initialization'
      call impint(delth, alph)

      if (iitest == 1) print *, ' forward half-step'
      call step(1, 1, delth, alph, rob, wil)

      if (iitest == 1) print *, ' semi-impl. initialization'
      call impint(delt, alph)

      if (iitest == 1) print *, ' leapfrog half-step'
      call step(1, 2, delt, alph, rob, wil)
    end if

    if (iitest == 1) print *, ' semi-impl. initialization'
    call impint(delt2, alph)
end
