module time_stepping
    use params

    implicit none

    private
    public first_step, step

contains
    ! Call initialization of semi-implicit scheme and perform initial time step
    subroutine first_step
        use implicit, only: initialize_implicit

        call initialize_implicit(0.5*delt)

        call step(1, 1, 0.5*delt)

        call initialize_implicit(delt)

        call step(1, 2, delt)

        call initialize_implicit(2*delt)
    end

    ! Perform one time step starting from F(1) and F(2) and using the following scheme:
    ! Fnew = F(1) + DT * [ T_dyn(F(J2)) + T_phy(F(1)) ]
    ! F(1) = (1-2*eps)*F(J1) + eps*[F(1)+Fnew]
    ! F(2) = Fnew
    ! Input:
    ! If j1 == 1, j2 == 1 : forward time step (eps = 0)
    ! If j1 == 1, j2 == 2 : initial leapfrog time step (eps = 0)
    ! If j1 == 2, j2 == 2 : leapfrog time step with time filter (eps = ROB)
    ! dt = time step
    subroutine step(j1, j2, dt)
        use dynamical_constants, only: tdrs
        use prognostics
        use horizontal_diffusion, only: do_horizontal_diffusion, &
            & dmp, dmpd, dmps, dmp1, dmp1d, dmp1s, tcorv, qcorv, tcorh, qcorh
        use tendencies, only: get_tendencies

        integer, intent(in) :: j1, j2
        real, intent(in) :: dt
        complex, dimension(mx,nx,kx) ::  vordt, divdt, tdt
        complex :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
        real :: eps, sdrag

        complex :: ctmp(mx,nx,kx)

        integer :: n, itr, k, m

        ! =========================================================================
        ! Compute tendencies of prognostic variables
        ! =========================================================================

        call get_tendencies(vordt, divdt, tdt, psdt, trdt, j2)

        ! =========================================================================
        ! Horizontal diffusion
        ! =========================================================================

        ! Diffusion of wind and temperature
        vordt = do_horizontal_diffusion(vor(:,:,:,1), vordt, dmp,  dmp1)
        divdt = do_horizontal_diffusion(div(:,:,:,1), divdt, dmpd, dmp1d)

        do k = 1, kx
            do m = 1, mx
                do n = 1, nx
                    ctmp(m,n,k) = t(m,n,k,1) + tcorh(m,n)*tcorv(k)
                end do
            end do
        end do

        tdt = do_horizontal_diffusion(ctmp, tdt, dmp, dmp1)

        ! Stratospheric diffusion and zonal wind damping
        sdrag = 1.0/(tdrs*3600.0)
        do n = 1, nx
            vordt(1,n,1) = vordt(1,n,1) - sdrag*vor(1,n,1,1)
            divdt(1,n,1) = divdt(1,n,1) - sdrag*div(1,n,1,1)
        end do

        vordt = do_horizontal_diffusion(vor(:,:,:,1),  vordt, dmps, dmp1s)
        divdt = do_horizontal_diffusion(div(:,:,:,1),  divdt, dmps, dmp1s)
        tdt   = do_horizontal_diffusion(ctmp, tdt,   dmps, dmp1s)

        ! Diffusion of tracers
        do k = 1, kx
            do m = 1, mx
                do n = 1, nx
                    ctmp(m,n,k) = tr(m,n,k,1,1) + qcorh(m,n)*qcorv(k)
                end do
            end do
        end do

        trdt(:,:,:,1) = do_horizontal_diffusion(ctmp, trdt(:,:,:,1), dmpd, dmp1d)

        if (ntr > 1) then
            do itr = 2, ntr
                trdt(:,:,:,1) = do_horizontal_diffusion(tr(:,:,:,1,itr), trdt(:,:,:,itr), dmp, dmp1)
            enddo
        endif

        ! =========================================================================
        ! Time integration with Robert filter
        ! =========================================================================

        if (j1 == 1) then
            eps = 0.0
        else
            eps = rob
        endif

        ps  = step_field_2d(j1, dt, eps, ps, psdt)
        vor = step_field_3d(j1, dt, eps, vor, vordt)
        div = step_field_3d(j1, dt, eps, div, divdt)
        t   = step_field_3d(j1, dt, eps, t, tdt)

        do itr = 1, ntr
            tr(:,:,:,:,itr) = step_field_3d(j1, dt, eps, tr(:,:,:,:,itr), trdt(:,:,:,itr))
        end do
    end



    ! Perform time integration of field across all model levels using tendency fdt
    function step_field_3d(j1, dt, eps, input, fdt) result(output)
        use spectral, only: trunct

        integer, intent(in) :: j1
        real, intent(in) :: dt, eps
        complex, intent(inout) :: fdt(mx,nx,kx)
        complex, intent(in) :: input(mx,nx,kx,2)
        complex :: output(mx,nx,kx,2)
        integer :: k

        do k = 1, kx
            output(:,:,k,:) = step_field_2d(j1, dt, eps, input(:,:,k,:), fdt(:,:,k))
        end do
    end

    function step_field_2d(j1, dt, eps, input, fdt) result(output)
        use spectral, only: trunct

        integer, intent(in) :: j1
        real, intent(in) :: dt, eps
        complex, intent(inout) :: fdt(mx,nx)
        complex, intent(in) :: input(mx,nx,2)
        complex :: output(mx,nx,2)
        real :: eps2
        complex :: fnew(mx,nx)

        output = input

        eps2 = 1.0 - 2.0*eps

        if (ix == iy*4) then
            call trunct(fdt)
        end if

        ! The actual leap frog with the Robert filter
        fnew = output(:,:,1) + dt*fdt
        output(:,:,1) = output(:,:,j1) + wil*eps*(output(:,:,1) - 2*output(:,:,j1) + fnew)

        ! Williams' innovation to the filter
        output(:,:,2) = fnew - (1.0 - wil)*eps*(output(:,:,1) - 2.0*output(:,:,j1) + fnew)
    end
end module
