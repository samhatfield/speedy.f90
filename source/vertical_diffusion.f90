!> Parametrization of vertical diffusion
module vertical_diffusion
    use params

    implicit none

    private
    public get_vertical_diffusion_tend

    ! Constants for vertical diffusion and shallow convection.
    real, parameter :: trshc  = 6.0  !! Relaxation time (in hours) for shallow convection
    real, parameter :: trvdi  = 24.0 !! Relaxation time (in hours) for moisture diffusion
    real, parameter :: trvds  = 6.0  !! Relaxation time (in hours) for super-adiabatic conditions
    real, parameter :: redshc = 0.5  !! Reduction factor of shallow convection in areas of deep
                                     !! convection
    real, parameter :: rhgrad = 0.5  !! Maximum gradient of relative humidity (d_RH/d_sigma)
    real, parameter :: segrad = 0.1  !! Minimum gradient of dry static energy (d_DSE/d_phi)

contains
    !> Compute tendencies of momentum, energy and moisture due to vertical diffusion
    !  and shallow convection
    subroutine get_vertical_diffusion_tend(se, rh, qa, qsat, phi, icnv, utenvd, vtenvd, &
        & ttenvd, qtenvd)
        use physical_constants, only: cp, alhc, sigh
        use geometry, only: fsg, dhs

        real, intent(in)    :: se(ix,il,kx)     !! Dry static energy
        real, intent(in)    :: rh(ix,il,kx)     !! Relative humidity
        real, intent(in)    :: qa(ix,il,kx)     !! Specific humidity [g/kg]
        real, intent(in)    :: qsat(ix,il,kx)   !! Saturated specific humidity [g/kg]
        real, intent(in)    :: phi(ix,il,kx)    !! Geopotential
        integer, intent(in) :: icnv(ix,il)      !! Sigma-level index of deep convection
        real, intent(out)   :: utenvd(ix,il,kx) !! u-wind tendency
        real, intent(out)   :: vtenvd(ix,il,kx) !! v-wind tendency
        real, intent(out)   :: ttenvd(ix,il,kx) !! Temperature tendency
        real, intent(out)   :: qtenvd(ix,il,kx) !! Specific humidity tendency

        integer :: nl1, i, j, k, k1
        real :: cshc, cvdi, fshcq, fshcse, fvdiq, fvdise, drh0, fvdiq2, dmse, drh
        real :: fluxse, fluxq, fcnv, se0
        real, dimension(kx) :: rsig, rsig1

        ! 1. Initalization

        ! N.B. In this routine, fluxes of dry static energy and humidity
        !      are scaled in such a way that:
        !      d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma

        nl1  = kx - 1
        cshc = dhs(kx)/3600.0
        cvdi = (sigh(nl1) - sigh(1))/((nl1 - 1)*3600.0)

        fshcq  = cshc/trshc
        fshcse = cshc/(trshc*cp)

        fvdiq  = cvdi/trvdi
        fvdise = cvdi/(trvds*cp)

        do k = 1, nl1
            rsig(k) = 1.0/dhs(k)
            rsig1(k) = 1.0/(1.0 - sigh(k))
        end do
        rsig(kx)=1.0/dhs(kx)

        utenvd = 0.0
        vtenvd = 0.0
        ttenvd = 0.0
        qtenvd = 0.0

        ! 2. Shallow convection
        drh0   = rhgrad*(fsg(kx) - fsg(nl1))
        fvdiq2 = fvdiq*sigh(nl1)

        do i = 1, ix
            do j = 1, il
                dmse = se(i,j,kx) - se(i,j,nl1) + alhc*(qa(i,j,kx) - qsat(i,j,nl1))
                drh  = rh(i,j,kx) - rh(i,j,nl1)
                fcnv = 1.0

                if (dmse >= 0.0) then
                    if (icnv(i,j) > 0) fcnv = redshc

                    fluxse           =  fcnv*fshcse*dmse
                    ttenvd(i,j,nl1)  =  fluxse*rsig(nl1)
                    ttenvd(i,j,kx)   = -fluxse*rsig(kx)

                    if (drh >= 0.0) then
                        fluxq           =  fcnv*fshcq*qsat(i,j,kx)*drh
                        qtenvd(i,j,nl1) =  fluxq*rsig(nl1)
                        qtenvd(i,j,kx)  = -fluxq*rsig(kx)
                    end if
                else if (drh > drh0) then
                    fluxq           =  fvdiq2*qsat(i,j,nl1)*drh
                    qtenvd(i,j,nl1) =  fluxq*rsig(nl1)
                    qtenvd(i,j,kx)  = -fluxq*rsig(kx)
                end if
            end do
        end do

        ! 3. Vertical diffusion of moisture above the PBL
        do k = 3, kx - 2
            if (sigh(k) > 0.5) then
                drh0   = rhgrad*(fsg(k+1) - fsg(k))
                fvdiq2 = fvdiq*sigh(k)

                do i = 1, ix
                    do j = 1, il
                        drh = rh(i,j,k+1) - rh(i,j,k)
                        if (drh >= drh0) then
                            fluxq           = fvdiq2*qsat(i,j,k)*drh
                            qtenvd(i,j,k)   = qtenvd(i,j,k) + fluxq*rsig(k)
                            qtenvd(i,j,k+1) = qtenvd(i,j,k+1) - fluxq*rsig(k+1)
                        end if
                    end do
                end do
            end if
        end do

        ! 4. Damping of super-adiabatic lapse rate
        do k = 1, nl1
            do i = 1, ix
                do j = 1, il
                    se0 = se(i,j,k+1) + segrad*(phi(i,j,k) - phi(i,j,k+1))

                    if (se(i,j,k) < se0) then
                        fluxse        = fvdise*(se0 - se(i,j,k))
                        ttenvd(i,j,k) = ttenvd(i,j,k) + fluxse*rsig(k)
                        do k1 = k+1, kx
                            ttenvd(i,j,k1) = ttenvd(i,j,k1) - fluxse*rsig1(k)
                        end do
                    end if
                end do
            end do
        end do
    end
end module
