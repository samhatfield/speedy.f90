module tendencies
    use mod_atparam

    implicit none

contains
    subroutine get_tendencies(vordt, divdt, tdt, psdt, trdt, j2)
        use implicit, only: implicit_terms
        use mod_tsteps, only: alph

        complex, dimension(mx,nx,kx), intent(inout) ::  vordt, divdt, tdt
        complex, intent(inout) :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
        integer, intent(in) :: j2

        ! =========================================================================
        ! Computation of grid-point tendencies (converted to spectral at the end of
        ! grtend)
        ! =========================================================================

        call grtend(vordt, divdt, tdt, psdt, trdt, 1, j2)

        ! =========================================================================
        ! Computation of spectral tendencies
        ! =========================================================================

        if (alph < 0.5) then
            call sptend(divdt, tdt, psdt, j2)
        else
            call sptend(divdt, tdt, psdt, 1)

            ! Implicit correction
            call implicit_terms(divdt, tdt, psdt)
        end if
    end subroutine
end module