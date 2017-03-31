subroutine invars 
    !   subroutine invars (istart)
    !
    !   Purpose : initialize all spectral variables starting from
    !             either a reference atmosphere or a restart file
    !   Input :   istart = 0 : reference atmosphere (at rest)
    !                    = 1 : restart file
    !   Initialized common blocks : date1, dynsp1, dynsp2 (phis only),
    !                               sfcanom, sfcflux

    use mod_tsteps, only: iyear0, imont0, istart
    use mod_dyncon0, only: gamma, hscale, hshum, refrh1
    use mod_atparam
    use mod_dynvar
    use mod_dyncon1, only: grav, rgas, fsg
    use mod_surfcon, only: phi0, phis0
    use mod_date, only: iyear, imonth, iday, ihour

    implicit none

    complex :: zero, ccon, surfs(mx,nx)
    real :: surfg(ix,il)
    real :: gam1, esref, factk, gam2, qexp, qref, rgam, rgamr, rlog0, tref, ttop
    integer :: i, j, k

    gam1 = gamma/(1000.*grav)
    zero = (0.,0.)
    ccon = (1.,0.)*sqrt(2.)

    ! 1. Compute spectral surface geopotential
    call spec(phi0,phis)
    if (ix.eq.iy*4) call trunct(phis)

    call grid(phis,phis0,1)

    if (istart.eq.0) then
        ! 2. Start from reference atmosphere (at rest) 
        print*, ' starting from rest'

        iyear  = iyear0
        imonth = imont0
        iday = 1
        ihour = 0

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

        call spec(surfg,ps)
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
            print *, ' q0 jlat = ', j, surfg(1,j)
        end do

        call spec(surfg,surfs)
        if (ix.eq.iy*4) call trunct (surfs)

        ! Spec. humidity at tropospheric levels      
        do k=3,kx
            factk=fsg(k)**qexp
            print *, 'vertical scale factor at level ', k, factk
            tr(:,:,k,1,1) = surfs * factk
        end do

        ! Print diagnostics from initial conditions
        call diagns (1,0)
    else if (istart .eq. 1) then
        ! 3. Start from restart file 
        print*,' reading a restart file'

        call restart(0)

        ! Print diagnostics from initial conditions
        call diagns (2,0)
    else if (istart .eq. 2) then
        ! 4. Start from grid initial condition

        call iogrid(1)

        call diagns(1,0)
    else
        print *, 'IMPOSSIBLE!! check the fort.2 file!'
        stop
    endif
end
