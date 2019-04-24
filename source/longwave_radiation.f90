module longwave_radiation
    implicit none

    private
    public get_longwave_rad_fluxes, radset

contains
    ! Compute the absorption of longwave radiation
    ! Input:   imode  = index for operation mode
    !                   -1 : downward flux only
    !                    0 : downward + upward flux
    !                   +1 : upward flux only
    !          ta     = absolute temperature (3-dim)
    !          ts     = surface temperature                    [if imode=0]
    !          fsfcd  = downward flux of lw rad. at the sfc.   [if imode=1]
    !          fsfcu  = surface blackbody emission (upward)    [if imode=1]
    !          dfabs  = DFABS output from RADLW(-1,... )       [if imode=1]
    ! Output:  fsfcd  = downward flux of lw rad. at the sfc.[if imode=-1,0]
    !          fsfcu  = surface blackbody emission (upward)  [if imode=  0]
    !          fsfc   = net upw. flux of lw rad. at the sfc. [if imode=0,1]
    !          ftop   = outgoing flux of lw rad. at the top  [if imode=0,1]
    !          dfabs  = flux of lw rad. absorbed by each atm. layer (3-dim)
    subroutine get_longwave_rad_fluxes(imode, ta, ts, fsfcd, fsfcu, fsfc, ftop, dfabs)
        use mod_atparam
        use physical_constants, only: sbc, wvi
        use geometry, only: dhs
        use mod_radcon, only: epslw, emisfc, fband, tau2, st4a, stratc, flux

        integer, intent(in) :: imode

        ! Number of radiation bands with tau < 1
        integer, parameter :: nband = 4

        real, intent(in) :: ta(ix,il,kx), ts(ix,il)
        real, intent(inout) :: fsfcd(ix,il), fsfcu(ix,il), ftop(ix,il), fsfc(ix,il)
        real, intent(inout) :: dfabs(ix,il,kx)

        integer :: i, j, jb, k, nl1
        real :: anis, brad, corlw(ix,il), corlw1(ix,il), corlw2(ix,il), emis, refsfc
        real :: st3a(ix,il), tsq

        nl1 = kx - 1

        refsfc = 1.0 - emisfc

        if (imode == 1) go to 410
        ! 1. Blackbody emission from atmospheric levels.
        ! The linearized gradient of the blakbody emission is computed
        ! from temperatures at layer boundaries, which are interpolated
        ! assuming a linear dependence of T on log_sigma.
        ! Above the first (top) level, the atmosphere is assumed isothermal.

        ! Temperature at level boundaries
        do k = 1, nl1
            st4a(:,:,k,1) = ta(:,:,k) + wvi(k,2)*(ta(:,:,k+1) - ta(:,:,k))
        end do

        ! Mean temperature in stratospheric layers
        st4a(:,:,1,2) = 0.75*ta(:,:,1) + 0.25*st4a(:,:,1,1)
        st4a(:,:,2,2) = 0.50*ta(:,:,2) + 0.25*(st4a(:,:,1,1) + st4a(:,:,2,1))

        ! Temperature gradient in tropospheric layers
        anis  = 1.0

        do k = 3, nl1
            st4a(:,:,k,2) = 0.5*anis*max(st4a(:,:,k,1) - st4a(:,:,k-1,1), 0.0)
        end do

        st4a(:,:,kx,2) = anis*max(ta(:,:,kx) - st4a(:,:,nl1,1), 0.0)

        ! Blackbody emission in the stratosphere
        do k = 1, 2
            st4a(:,:,k,1) = sbc*st4a(:,:,k,2)**4.0
            st4a(:,:,k,2) = 0.0
        end do

        ! Blackbody emission in the troposphere
        do k = 3, kx
            st3a = sbc*ta(:,:,k)**3.0
            st4a(:,:,k,1) = st3a*ta(:,:,k)
            st4a(:,:,k,2) = 4.0*st3a*st4a(:,:,k,2)
        end do

        ! 2. Initialization of fluxes
        fsfcd = 0.0
        dfabs = 0.0

        ! 3. Emission ad absorption of longwave downward flux.
        !    For downward emission, a correction term depending on the
        !    local temperature gradient and on the layer transmissivity is
        !    added to the average (full-level) emission of each layer.

        ! 3.1  Stratosphere
        k=1
        do jb = 1, 2
            do i = 1, ix
                do j = 1, il
                    emis = 1.0 - tau2(i,j,k,jb)
                    brad = fband(nint(ta(i,j,k)),jb)*(st4a(i,j,k,1) + emis*st4a(i,j,k,2))
                    flux(i,j,jb) = emis*brad
                    dfabs(i,j,k) = dfabs(i,j,k) - flux(i,j,jb)
                end do
            end do
        end do

        flux(:,:,3:nband) = 0.0

        ! 3.2  Troposphere
        do jb = 1, nband
            do k = 2, kx
                do i = 1, ix
                    do j = 1, il
                        emis = 1.0 - tau2(i,j,k,jb)
                        brad = fband(nint(ta(i,j,k)),jb)*(st4a(i,j,k,1) + emis*st4a(i,j,k,2))
                        dfabs(i,j,k) = dfabs(i,j,k) + flux(i,j,jb)
                        flux(i,j,jb) = tau2(i,j,k,jb)*flux(i,j,jb) + emis*brad
                        dfabs(i,j,k) = dfabs(i,j,k) - flux(i,j,jb)
                    end do
                end do
            end do
        end do

        ! 3.3 Surface downward flux
        do jb = 1, nband
            fsfcd = fsfcd + emisfc*flux(:,:,jb)
        end do

        ! 3.4 Correction for "black" band (incl. surface reflection)
        corlw = epslw*emisfc*st4a(:,:,kx,1)
        dfabs(:,:,kx) = dfabs(:,:,kx) - corlw
        fsfcd = fsfcd + corlw

        if (imode == -1) return

        ! 4. Emission ad absorption of longwave upward flux.
        !    For upward emission, a correction term depending on the
        !    local temperature gradient and on the layer transmissivity is
        !    subtracted from the average (full-level) emission of each layer.

        ! 4.1  Surface

        ! Black-body (or grey-body) emission
        fsfcu = emisfc*sbc*ts**4.0

        ! Entry point for upward-only mode (imode = 1)
     410  continue

        fsfc = fsfcu - fsfcd

        do jb = 1, nband
            do i = 1, ix
                do j = 1, il
                    flux(i,j,jb) = fband(nint(ts(i,j)),jb)*fsfcu(i,j) + refsfc*flux(i,j,jb)
                end do
            end do
        end do

        ! 4.2  Troposphere

        ! Correction for "black" band
        dfabs(:,:,kx) = dfabs(:,:,kx) + epslw*fsfcu

        do jb = 1, nband
            do k = kx, 2, -1
                do i = 1, ix
                    do j = 1, il
                        emis = 1.0 - tau2(i,j,k,jb)
                        brad = fband(nint(ta(i,j,k)),jb)*(st4a(i,j,k,1) - emis*st4a(i,j,k,2))
                        dfabs(i,j,k) = dfabs(i,j,k) + flux(i,j,jb)
                        flux(i,j,jb) = tau2(i,j,k,jb)*flux(i,j,jb) + emis*brad
                        dfabs(i,j,k) = dfabs(i,j,k) - flux(i,j,jb)
                    end do
                end do
            end do
        end do

        ! 4.3  Stratosphere
        k = 1
        do jb = 1, 2
            do i = 1, ix
                do j = 1, il
                    emis = 1.0 - tau2(i,j,k,jb)
                    brad = fband(nint(ta(i,j,k)),jb)*(st4a(i,j,k,1) - emis*st4a(i,j,k,2))
                    dfabs(i,j,k) = dfabs(i,j,k) + flux(i,j,jb)
                    flux(i,j,jb) = tau2(i,j,k,jb)*flux(i,j,jb) + emis*brad
                    dfabs(i,j,k) = dfabs(i,j,k) - flux(i,j,jb)
                end do
            end do
        end do

        ! Correction for "black" band and polar night cooling
        corlw1 = dhs(1)*stratc(:,:,2)*st4a(:,:,1,1) + stratc(:,:,1)
        corlw2 = dhs(2)*stratc(:,:,2)*st4a(:,:,2,1)
        dfabs(:,:,1) = dfabs(:,:,1) - corlw1
        dfabs(:,:,2) = dfabs(:,:,2) - corlw2
        ftop = corlw1 + corlw2

        ! 4.4  Outgoing longwave radiation
        do jb = 1, nband
            ftop = ftop + flux(:,:,jb)
        end do
    end

    ! Compute energy fractions in longwave bands as a function of temperature
    subroutine radset
        use mod_radcon, only: epslw, fband

        integer :: jb, jtemp
        real :: eps1

        eps1 = 1.0 - epslw

        do jtemp = 200, 320
            fband(jtemp,2) = (0.148 - 3.0e-6*(jtemp - 247)**2)*eps1
            fband(jtemp,3) = (0.356 - 5.2e-6*(jtemp - 282)**2)*eps1
            fband(jtemp,4) = (0.314 + 1.0e-5*(jtemp - 315)**2)*eps1
            fband(jtemp,1) = eps1 - (fband(jtemp,2) + fband(jtemp,3) + fband(jtemp,4))
        end do

        do jb = 1, 4
            do jtemp = 100, 199
                fband(jtemp,jb) = fband(200,jb)
            end do
            do jtemp = 321, 400
                fband(jtemp,jb) = fband(320,jb)
            end do
        end do
    end
end module
