subroutine invars
    !   Purpose : initialize all spectral variables starting from
    !             either a reference atmosphere or a restart file

    use mod_dyncon0, only: gamma, hscale, hshum, refrh1
    use params
    use prognostics
    use physical_constants, only: grav, rgas
    use geometry, only: fsg
    use boundaries, only: phis0
    use diagnostics, only: check_diagnostics
    use spectral, only: grid_to_spec, trunct

    implicit none

    complex :: zero, ccon, surfs(mx,nx)
    real :: surfg(ix,il)
    real :: gam1, esref, factk, gam2, qexp, qref, rgam, rgamr, rlog0, tref, ttop
    integer :: i, j, k

    gam1 = gamma/(1000.*grav)
    zero = (0.,0.)
    ccon = (1.,0.)*sqrt(2.)

    ! 1. Compute spectral surface geopotential
    phis = grid_to_spec(phis0, 1)

    ! 2. Start from reference atmosphere (at rest)
    print*, 'Starting from rest'

    ! 2.1 Set vorticity, divergence and tracers to zero
    vor(:,:,:,1) = zero
    div(:,:,:,1) = zero
    tr(:,:,:,1,:) = zero

    ! 2.2 Set reference temperature :
    !     tropos:  T = 288 degK at z = 0, constant lapse rate
    !     stratos: T = 216 degK, lapse rate = 0
    tref  = 288.
    ttop  = 216.
    gam2  = gam1/tref
    rgam  = rgas*gam1
    rgamr = 1./rgam

    ! Surface and stratospheric air temperature
    t(:,:,1,1) = zero
    t(:,:,2,1) = zero
    surfs = -gam1 * phis

    t(1,1,1,1) = ccon*ttop
    t(1,1,2,1) = ccon*ttop
    surfs(1,1) = ccon*tref - gam1*phis(1,1)

    ! Temperature at tropospheric levels
    do k=3,kx
        factk=fsg(k)**rgam
        t(:,:,k,1) = surfs * factk
    end do

    ! 2.3 Set log(ps) consistent with temperature profile
    !     p_ref = 1013 hPa at z = 0
    rlog0 = log(1.013)

    do j=1,il
        do i=1,ix
            surfg(i,j) = rlog0 + rgamr*log(1.-gam2*phis0(i,j))
        end do
    end do

    ps(:,:,1) = grid_to_spec(surfg, 1)
    if (ix.eq.iy*4) call trunct(ps)

    ! 2.4 Set tropospheric spec. humidity in g/kg
    !     Qref = RHref * Qsat(288K, 1013hPa)
    esref = 17.
    qref = refrh1*0.622*esref
    qexp = hscale/hshum

    ! Spec. humidity at the surface
    do j=1,il
        do i=1,ix
            surfg(i,j)=qref*exp(qexp*surfg(i,j))
        end do
    end do

    surfs = grid_to_spec(surfg, 1)
    if (ix.eq.iy*4) call trunct (surfs)

    ! Spec. humidity at tropospheric levels
    do k=3,kx
        factk=fsg(k)**qexp
        tr(:,:,k,1,1) = surfs * factk
    end do

    ! Print diagnostics from initial conditions
    call check_diagnostics(vor(:,:,:,1), div(:,:,:,1), t(:,:,:,1), 0)
end
