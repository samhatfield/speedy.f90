! Call initialization routines for all model common blocks
subroutine ini_atm()
    use mod_tsteps, only: indrdf, ipout, ihout
    use mod_atparam
    use mod_dyncon1, only: grav, hsg, fsg, radang
    use mod_tmean

    implicit none

    ! Post-processing levels (hpa/1000)
    real :: ppl(kx)
    integer :: is3d = 1, k, nddm, ndm, ndtm, ntm

    ! Initialize ffts
    call inifft

    ! Initialize dynamical constants and operators
    call indyns

    ! Set post-processing levels
    do k = 1, kx
        ppl(k) = prlev(fsg(k))
    end do

    ! Initialize constants for physical parametrization
    call inphys(hsg, ppl, radang)

    ! Initialize forcing fields (boundary cond. + random forcing)
    call inbcon(grav,radang)

    call inirdf(indrdf)

    ! Initialize model variables
    call invars

    ! Initialize time-mean arrays for surface fluxes and output fields
    call dmflux(0)

    call dmout(0)

    ! Create control file for 6-hourly output
    call iogrid(5)
    if (ipout) then
        call iogrid(3)
        call geop(1)
    end if

    ! Output files for grid-point fields
    if (ihout .eqv. .false.) call setgrd(0)

    ! Write initial data
    if (ihout .and. ipout) call iogrid(2)
    if (ihout) call iogrid(4)

    contains
        function prlev(siglev)
            ! function prlev (siglev)
            ! purpose : select the closest standard pressure level for post-proc.
            ! input :   siglev = sigma level
            implicit none

            real, intent(in) :: siglev
            real :: plev(14) = (/ 0.925, 0.850, 0.775, 0.700, 0.600, 0.500, 0.400,&
                & 0.300, 0.250, 0.200, 0.150, 0.100, 0.050, 0.030 /)
            real :: prlev, dif, adif
            integer :: k

            dif = 1.0 - siglev

            prlev = 1.0

            do k = 1, 14
                adif = abs(plev(k) - siglev)
                if (adif <= dif) then
                    dif = adif
                    prlev = plev(k)
                end if
            end do
        end
end
