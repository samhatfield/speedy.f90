!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 04/07/2019
!  For storing and initializing prognostic spectral variables for model dynamics, and geopotential.
module prognostics
    use types, only: p
    use params

    implicit none

    private
    public vor, div, t, ps, tr
    public phi, phis
    public initialize_prognostics

    ! Prognostic spectral variables
    complex(p) :: vor(mx,nx,kx,2)    !! Vorticity
    complex(p) :: div(mx,nx,kx,2)    !! Divergence
    complex(p) :: t(mx,nx,kx,2)      !! Absolute temperature
    complex(p) :: ps(mx,nx,2)        !! Log of (normalised) surface pressure (p_s/p0)
    complex(p) :: tr(mx,nx,kx,2,ntr) !! Tracers (tr(1): specific humidity in g/kg)

    ! Geopotential
    complex(p) :: phi(mx,nx,kx) !! Atmospheric geopotential
    complex(p) :: phis(mx,nx)   !! Surface geopotential

    character(len=*), parameter :: initfile = "init.nc" !! Restart file

contains
    !> Initializes all spectral variables starting from either a reference
    !  atmosphere or a restart file.
    subroutine initialize_prognostics
        logical :: init_file_exists
        inquire(file=initfile, exist=init_file_exists)
        if (init_file_exists) then
            call initialize_from_file
        else
            call initialize_from_rest_state
        end if
    end subroutine

    !> Initializes all spectral variables starting from a reference atmosphere.
    subroutine initialize_from_rest_state
        use dynamical_constants, only: gamma, hscale, hshum, refrh1
        use physical_constants, only: grav, rgas
        use geometry, only: fsg
        use boundaries, only: phis0
        use diagnostics, only: check_diagnostics
        use spectral, only: grid_to_spec, trunct
        use input_output, only: output
        use geopotential, only: get_geopotential

        complex(p) :: surfs(mx,nx)
        real(p) :: surfg(ix,il)
        real(p) :: gam1, esref, gam2, qexp, qref, rgam, rgamr, rlog0, tref, ttop
        integer :: i, j, k

        gam1 = gamma/(1000.0*grav)

        ! 1. Compute spectral surface geopotential
        phis = grid_to_spec(phis0)

        ! 2. Start from reference atmosphere (at rest)
        write (*,'(A)') 'Starting from rest'

        ! 2.1 Set vorticity, divergence and tracers to zero
        vor(:,:,:,1) = (0.0, 0.0)
        div(:,:,:,1) = (0.0, 0.0)
        tr(:,:,:,1,:) = (0.0, 0.0)

        ! 2.2 Set reference temperature :
        !     tropos:  T = 288 degK at z = 0, constant lapse rate
        !     stratos: T = 216 degK, lapse rate = 0
        tref  = 288.0
        ttop  = 216.0
        gam2  = gam1/tref
        rgam  = rgas*gam1
        rgamr = 1.0/rgam

        ! Surface and stratospheric air temperature
        t(:,:,1,1) = (0.0, 0.0)
        t(:,:,2,1) = (0.0, 0.0)
        surfs = -gam1 * phis

        t(1,1,1,1) = sqrt(2.0)*(1.0, 0.0)*ttop
        t(1,1,2,1) = sqrt(2.0)*(1.0, 0.0)*ttop
        surfs(1,1) = sqrt(2.0)*(1.0, 0.0)*tref - gam1*phis(1,1)

        ! Temperature at tropospheric levels
        do k = 3, kx
            t(:,:,k,1) = surfs*fsg(k)**rgam
        end do

        ! 2.3 Set log(ps) consistent with temperature profile
        !     p_ref = 1013 hPa at z = 0
        rlog0 = log(1.013)

        do j=1,il
            do i=1,ix
                surfg(i,j) = rlog0 + rgamr*log(1.0 - gam2*phis0(i,j))
            end do
        end do

        ps(:,:,1) = grid_to_spec(surfg)
        if (ix == iy*4) call trunct(ps)

        ! 2.4 Set tropospheric specific humidity in g/kg
        !     Qref = RHref * Qsat(288K, 1013hPa)
        esref = 17.0
        qref = refrh1*0.622*esref
        qexp = hscale/hshum

        ! Specific humidity at the surface
        do j = 1, il
            do i = 1, ix
                surfg(i,j)=qref*exp(qexp*surfg(i,j))
            end do
        end do

        surfs = grid_to_spec(surfg)
        if (ix == iy*4) call trunct (surfs)

        ! Specific humidity at tropospheric levels
        do k = 3, kx
            tr(:,:,k,1,1) = surfs*fsg(k)**qexp
        end do

        ! Print diagnostics from initial conditions
        call check_diagnostics(vor(:,:,:,1), div(:,:,:,1), t(:,:,:,1), 0)

        ! Geopotential is just a diagnostic variable so we need to compute it here
        phi = get_geopotential(t, phis)

        ! Write initial data
        call output(0, vor, div, t, ps, tr, phi)
    end subroutine

    !> Initializes all spectral variables from a file.
    subroutine initialize_from_file
        use input_output, only: load_init_file
        use physical_constants, only: p0
        use boundaries, only: phis0
        use diagnostics, only: check_diagnostics
        use spectral, only: grid_to_spec, trunct, vds
        use geometry, only: coa
        use input_output, only: output
        use geopotential, only: get_geopotential

        complex(p), dimension(mx,nx) :: ucosm, vcosm
        real(p), dimension(ix,il,kx) :: u_grid, v_grid, s_grid
        integer :: j, k

        ! 1. Compute spectral surface geopotential
        phis = grid_to_spec(phis0)

        ! 2. Start from init.nc
        write (*,'(A)') 'Starting from file'

        ! 2.1 Set vorticity, divergence
        vor(:,:,:,1) = (0.0, 0.0)
        div(:,:,:,1) = (0.0, 0.0)
        tr(:,:,:,1,:) = (0.0, 0.0)

        u_grid = load_init_file(initfile, "u")
        v_grid = load_init_file(initfile, "v")

        do j = 1, il
            u_grid(:,j,:) = u_grid(:,j,:) * coa(j)
            v_grid(:,j,:) = v_grid(:,j,:) * coa(j)
        end do
        do k = 1, kx
            ucosm(:,:) = grid_to_spec(u_grid(:,:,k))
            vcosm(:,:) = grid_to_spec(v_grid(:,:,k))
            call vds(ucosm(:,:),vcosm(:,:),vor(:,:,k,1),div(:,:,k,1))
        end do

        ! 2.2 Set temperature :
        t(:,:,:,1) = (0.0, 0.0)
        s_grid = load_init_file(initfile, "t")
        do k = 1, kx
            t(:,:,k,1) = grid_to_spec(s_grid(:,:,k))
        end do

        ! 2.3 Set log(ps)
        s_grid = load_init_file(initfile, "ps")
        ps(:,:,1) = grid_to_spec(log(s_grid(:,:,1)/p0))
        if (ix == iy*4) call trunct(ps)

        ! 2.4 Set specific humidity in g/kg
        ! Specific humidity
        s_grid = load_init_file(initfile, "q")
        do k = 1, kx
            tr(:,:,k,1,1) = grid_to_spec(s_grid(:,:,k) * 1.0e3)
        end do

        ! Print diagnostics from initial conditions
        call check_diagnostics(vor(:,:,:,1), div(:,:,:,1), t(:,:,:,1), 0)

        ! Geopotential is just a diagnostic variable so we need to compute it here
        phi = get_geopotential(t, phis)

        ! Write initial data
        call output(0, vor, div, t, ps, tr, phi)
    end subroutine
end module
