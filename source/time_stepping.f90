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
        use mod_dyncon0, only: tdrs
        use prognostics
        use mod_hdifcon
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
        call hordif(kx, vor, vordt, dmp,  dmp1)
        call hordif(kx, div, divdt, dmpd, dmp1d)

        do k = 1, kx
            do m = 1, mx
                do n = 1, nx
                    ctmp(m,n,k) = t(m,n,k,1) + tcorh(m,n)*tcorv(k)
                end do
            end do
        end do

        call hordif(kx,ctmp,tdt,dmp,dmp1)

        ! Stratospheric diffusion and zonal wind damping
        sdrag = 1.0/(tdrs*3600.0)
        do n = 1, nx
            vordt(1,n,1) = vordt(1,n,1) - sdrag*vor(1,n,1,1)
            divdt(1,n,1) = divdt(1,n,1) - sdrag*div(1,n,1,1)
        end do

        call hordif(1, vor,  vordt, dmps, dmp1s)
        call hordif(1, div,  divdt, dmps, dmp1s)
        call hordif(1, ctmp, tdt,   dmps, dmp1s)

        ! Diffusion of tracers
        do k = 1, kx
            do m = 1, mx
                do n = 1, nx
                    ctmp(m,n,k) = tr(m,n,k,1,1) + qcorh(m,n)*qcorv(k)
                end do
            end do
        end do

        call hordif(kx, ctmp, trdt, dmpd, dmp1d)

        if (ntr > 1) then
            do itr = 2, ntr
                call hordif(kx, tr(:,:,:,1,itr), trdt(:,:,:,itr), dmp, dmp1)
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

    ! Add horizontal diffusion tendency of FIELD to spectral tendency FDT at NLEV
    ! levels using damping coefficients DMP and DMP1
    subroutine hordif(nlev, field, fdt, dmp, dmp1)
        integer, intent(in) :: nlev
        complex, intent(in) :: field(mx*nx,kx)
        complex, intent(inout) :: fdt(mx*nx,kx)
        real, intent(in) :: dmp(mx*nx), dmp1(mx*nx)
        integer :: k, m

        do k = 1, nlev
            do m = 1, mx*nx
                fdt(m,k) = (fdt(m,k) - dmp(m)*field(m,k))*dmp1(m)
            end do
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
