module tendencies
    use mod_atparam

    implicit none

    private
    public get_tendencies

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

        call get_grid_point_tendencies(vordt, divdt, tdt, psdt, trdt, 1, j2)

        ! =========================================================================
        ! Computation of spectral tendencies
        ! =========================================================================

        if (alph < 0.5) then
            call get_spectral_tendencies(divdt, tdt, psdt, j2)
        else
            call get_spectral_tendencies(divdt, tdt, psdt, 1)

            ! Implicit correction
            call implicit_terms(divdt, tdt, psdt)
        end if
    end subroutine

    ! Compute non-linear tendencies in grid point space from dynamics and
    ! physical parametrizations, and convert them to spectral tendencies
    ! dF/dt = T_dyn(F(J2)) + T_phy(F(J1))
    !   Input:  j1 = time level index for physical tendencies
    !           j2 = time level index for dynamical tendencies
    !   Output: vordt = spectral tendency of vorticity
    !           divdt = spectral tendency of divergence
    !           tdt   = spectral tendency of temperature
    !           psdt  = spectral tendency of log(p_s)
    !           trdt  = spectral tendency of tracers
    subroutine get_grid_point_tendencies(vordt, divdt, tdt, psdt, trdt, j1, j2)
        use prognostics
        use physical_constants, only: akap, rgas
        use geometry, only: dhs, fsg, dhsr, fsgr, coriol
        use mod_dyncon2, only: tref, tref3
        use physics, only: get_physical_tendencies
        use spectral, only: grid_to_spec, spec_to_grid, laplacian, grad, uvspec, vdspec

        !** notes ****
        ! -- TG does not have to be computed at both time levels every time step,
        !     I have left it this way to retain parallel structure with subroutine
        !     using latitude loop
        ! -- memory can be reduced considerably eliminating TGG, computing VORG
        !     only when needed, etc -- I have not optimized this subroutine for
        !     routine use on the YMP
        ! -- results from grtend1.F should duplicate results from grtend.F
        !                              -- Isaac
        !************

        complex, dimension(mx,nx,kx), intent(inout) :: vordt, divdt, tdt
        complex, intent(inout) :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
        integer, intent(in) :: j1, j2

        complex :: dumc(mx,nx,2)

        real, dimension(ix,il,kx) :: utend, vtend, ttend
        real :: trtend(ix,il,kx,ntr)

        real, dimension(ix,il,kx) :: ug, vg, tg, vorg, divg, tgg, puv
        real, dimension(ix,il) :: px, py, umean, vmean, dmean
        real :: trg(ix,il,kx,ntr), sigdt(ix,il,kxp)
        real :: temp(ix,il,kxp), sigm(ix,il,kxp)

        integer :: k, i, itr, j

        ! =========================================================================
        ! Convert prognostics to grid point space
        ! =========================================================================

        do k = 1, kx
            vorg(:,:,k) = spec_to_grid(vor(:,:,k,j2), 1)
            divg(:,:,k) = spec_to_grid(div(:,:,k,j2), 1)
            tg(:,:,k)   = spec_to_grid(t(:,:,k,j2), 1)

            do itr = 1, ntr
                trg(:,:,k,itr) = spec_to_grid(tr(:,:,k,j2,itr), 1)
            end do

            call uvspec(vor(:,:,k,j2), div(:,:,k,j2), dumc(:,:,1), dumc(:,:,2))
            vg(:,:,k) = spec_to_grid(dumc(:,:,2), 2)
            ug(:,:,k) = spec_to_grid(dumc(:,:,1), 2)

            do j = 1, il
                do i = 1, ix
                    vorg(i,j,k) = vorg(i,j,k) + coriol(j)
                end do
            end do
        end do

        umean(:,:) = 0.0
        vmean(:,:) = 0.0
        dmean(:,:) = 0.0

        do k = 1, kx
            umean(:,:) = umean(:,:) + ug(:,:,k) * dhs(k)
            vmean(:,:) = vmean(:,:) + vg(:,:,k) * dhs(k)
            dmean(:,:) = dmean(:,:) + divg(:,:,k) * dhs(k)
        end do

        ! Compute tendency of log(surface pressure)
        ! ps(1,1,j2)=zero
        call grad(ps(:,:,j2), dumc(:,:,1), dumc(:,:,2))
        px = spec_to_grid(dumc(:,:,1), 2)
        py = spec_to_grid(dumc(:,:,2), 2)

        psdt = grid_to_spec(-umean*px - vmean*py)
        psdt(1,1) = (0.0, 0.0)

        ! Compute "vertical" velocity
        sigdt(:,:,1) = 0.0
        sigdt(:,:,kxp) = 0.0
        sigm(:,:,1) = 0.0
        sigm(:,:,kxp) = 0.0

        ! (The following combination of terms is utilized later in the
        !     temperature equation)
        do k = 1, kx
            puv(:,:,k) = (ug(:,:,k) - umean) * px + (vg(:,:,k) - vmean) * py
        end do

        do k = 1, kx
            sigdt(:,:,k+1) = sigdt(:,:,k) - dhs(k)*(puv(:,:,k)+divg(:,:,k)-dmean)
            sigm(:,:,k+1) = sigm(:,:,k) - dhs(k)*puv(:,:,k)
        end do

        ! Subtract part of temperature field that is used as reference for
        ! implicit terms
        do k = 1, kx
            tgg(:,:,k) = tg(:,:,k) - tref(k)
        end do

        ! Zonal wind tendency
        temp(:,:,1) = 0.0
        temp(:,:,kxp) = 0.0

        do k = 2, kx
            temp(:,:,k) = sigdt(:,:,k) * (ug(:,:,k) - ug(:,:,k-1))
        end do

        do k = 1, kx
            utend(:,:,k) = vg(:,:,k) * vorg(:,:,k) - tgg(:,:,k)*rgas*px &
                & - (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
        end do

        ! Meridional wind tendency
        do k = 2, kx
            temp(:,:,k) = sigdt(:,:,k) * (vg(:,:,k) - vg(:,:,k-1))
        end do

        do k = 1, kx
            vtend(:,:,k) = -ug(:,:,k)*vorg(:,:,k) - tgg(:,:,k)*rgas*py &
                & - (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
        end do

        ! Temperature tendency
        do k = 2, kx
            temp(:,:,k) = sigdt(:,:,k)*(tgg(:,:,k) - tgg(:,:,k-1)) &
                & + sigm(:,:,k)*(tref(k) - tref(k-1))
        end do

        do k = 1, kx
            ttend(:,:,k) = tgg(:,:,k)*divg(:,:,k) - (temp(:,:,k+1)+temp(:,:,k))*dhsr(k) &
                & + fsgr(k)*tgg(:,:,k)*(sigdt(:,:,k+1) + sigdt(:,:,k)) + tref3(k)*(sigm(:,:,k+1) &
                & + sigm(:,:,k)) + akap*(tg(:,:,k)*puv(:,:,k) - tgg(:,:,k)*dmean(:,:))
        end do

        ! Tracer tendency
        do itr = 1, ntr
            do k = 2, kx
                temp(:,:,k) = sigdt(:,:,k)*(trg(:,:,k,itr) - trg(:,:,k-1,itr))
            end do

            temp(:,:,2:3) = 0.0

            do k = 1, kx
                trtend(:,:,k,itr) = trg(:,:,k,itr)*divg(:,:,k)-(temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
            end do
        end do

        ! =========================================================================
        ! Compute physical tendencies
        ! =========================================================================

        call geop(j1)

        call get_physical_tendencies(vor(:,:,:,j1), div(:,:,:,j1), t(:,:,:,j1), tr(:,:,:,j1,1), &
            & phi, ps(:,:,j1), utend, vtend, ttend, trtend)

        ! =========================================================================
        ! Convert tendencies to spectral space
        ! =========================================================================

        do k = 1, kx
            !  Convert u and v tendencies to vor and div spectral tendencies
            !  vdspec takes a grid u and a grid v and converts them to
            !  spectral vor and div
            call vdspec(utend(:,:,k), vtend(:,:,k), vordt(:,:,k), divdt(:,:,k), 2)

            ! Divergence tendency
            ! add -lapl(0.5*(u**2+v**2)) to div tendency
            divdt(:,:,k) = divdt(:,:,k) &
                & - laplacian(grid_to_spec(0.5*(ug(:,:,k)**2.0 + vg(:,:,k)**2.0)))

            ! Temperature tendency
            ! and add div(vT) to spectral t tendency
            call vdspec(-ug(:,:,k)*tgg(:,:,k), -vg(:,:,k)*tgg(:,:,k), dumc(:,:,1), tdt(:,:,k), 2)
            tdt(:,:,k) = tdt(:,:,k) + grid_to_spec(ttend(:,:,k))

            ! tracer tendency
            do itr = 1, ntr
                call vdspec(-ug(:,:,k)*trg(:,:,k,itr), -vg(:,:,k)*trg(:,:,k,itr), &
                    & dumc(:,:,1), trdt(:,:,k,itr), 2)
                trdt(:,:,k,itr) = trdt(:,:,k,itr) + grid_to_spec(trtend(:,:,k,itr))
            end do
        end do
    end subroutine

    ! Compute spectral tendencies of divergence, temperature  and log(surface pressure)
    ! Input/output : divdt = divergence tendency (spectral)
    !                tdt   = temperature tendency (spectral)
    !                psdt  = tendency of log_surf.pressure (spectral)
    !                j2    = time level index (1 or 2)
    subroutine get_spectral_tendencies(divdt, tdt, psdt, j2)
        use prognostics, only: div, phi, ps
        use physical_constants, only: rgas
        use geometry, only: dhs, dhsr
        use mod_dyncon2, only: tref, tref2, tref3
        use spectral, only: laplacian

        complex, intent(inout) :: psdt(mx,nx), divdt(mx,nx,kx), tdt(mx,nx,kx)
        integer, intent(in) :: j2

        complex :: dumk(mx,nx,kxp), dmeanc(mx,nx), sigdtc(mx,nx,kxp)

        integer :: k

        ! Vertical mean div and pressure tendency
        dmeanc(:,:) = (0.0, 0.0)
        do k = 1, kx
            dmeanc = dmeanc + div(:,:,k,j2) * dhs(k)
        end do

        psdt = psdt - dmeanc
        psdt(1,1) = (0.0, 0.0)

        ! Sigma-dot "velocity" and temperature tendency
        sigdtc(:,:,1) = (0.0, 0.0)
        sigdtc(:,:,kxp) = (0.0, 0.0)

        do k = 1, kxm
            sigdtc(:,:,k+1) = sigdtc(:,:,k) - dhs(k)*(div(:,:,k,j2) - dmeanc)
        end do

        dumk(:,:,1) = (0.0, 0.0)
        dumk(:,:,kxp) = (0.0, 0.0)

        do k = 2, kx
            dumk(:,:,k) = sigdtc(:,:,k) * (tref(k) - tref(k-1))
        end do

        do k = 1, kx
            tdt(:,:,k) = tdt(:,:,k) - (dumk(:,:,k+1) + dumk(:,:,k)) * dhsr(k)&
                & + tref3(k) * (sigdtc(:,:,k+1) + sigdtc(:,:,k))&
                & - tref2(k) * dmeanc
        end do

        ! Geopotential and divergence tendency
        call geop(j2)

        do k = 1, kx
            divdt(:,:,k) = divdt(:,:,k) - laplacian(phi(:,:,k) + rgas*tref(k)*ps(:,:,j2))
        end do
    end subroutine
end module
