module spectral
    use mod_atparam

    implicit none

    private
    public el2
    public initialize_spectral
    public laplacian, inverse_laplacian, spec_to_grid, grid_to_spec
    public grad, vds, uvspec, vdspec, trunct

    real, dimension(mx,nx) :: el2, elm2, el4, trfilt
    real :: gradx(mx), gradym(mx,nx), gradyp(mx,nx)
    real, dimension(mx,nx) :: uvdx, uvdym, uvdyp
    real, dimension(mx,nx) :: vddym, vddyp

contains
    ! Initialize spectral transforms
    subroutine initialize_spectral
        use physical_constants, only: rearth
        use fourier, only: initialize_fourier
        use legendre, only: initialize_legendre, epsi

        real :: el1
        integer :: m, m1, m2, n, l2(mx,nx), wavenum_tot(mx,nx), mm(mx)

        ! Initialize Fourier transforms
        call initialize_fourier

        ! Initialize Legendre transforms
        call initialize_legendre

        !  MM = zonal wavenumber = m
        !     ISC=3 implies that only wavenumber 0,3,6,9,etc are included in model
        !  wavenum_tot = total wavenumber of spherical harmonic = l
        !  L2 = l*(l+1)
        !  EL2 = l*(l+1)/(a**2)
        !  EL4 = EL2*EL2 ; for biharmonic diffusion
        !  ELM2 = 1./EL2
        !  TRFILT used to filter out "non-triangular" part of rhomboidal truncation
        do n = 1, nx
            do m = 1, mx
                mm(m) = isc*(m - 1)
                wavenum_tot(m,n) = mm(m) + n - 1
                el2(m,n) = float(wavenum_tot(m,n)*(wavenum_tot(m,n) + 1))/(rearth**2.0)
                el4(m,n) = el2(m,n)**2.0
                if (wavenum_tot(m,n) <= ntrun) then
                    trfilt(m,n) = 1.0
                else
                    trfilt(m,n) = 0.0
                end if
            end do
        end do

        elm2(1,1) = 0.0
        elm2(2:mx,:) = 1.0/el2(2:mx,:)
        elm2(1,2:nx) = 1.0/el2(1,2:nx)

        ! quantities required by subroutines GRAD, UVSPEC, and VDS
        ! saved in spectral
        do m = 1, mx
            do n = 1, nx
                m1 = mm(m)
                m2 = m1 + 1
                el1 = float(wavenum_tot(m,n))
                if (n == 1) then
                    gradx(m)   = float(m1)/rearth
                    uvdx(m,1)  = -rearth/float(m1+1)
                    uvdym(m,1) = 0.0
                    vddym(m,1) = 0.0
                else
                    uvdx(m,n)   = -rearth*float(m1)/(el1*(el1 + 1))
                    gradym(m,n) = (el1 - 1.0)*epsi(m2,n)/rearth
                    uvdym(m,n)  = -rearth*epsi(m2,n)/el1
                    vddym(m,n)  = (el1 + 1)*epsi(m2,n)/rearth
                end if
                gradyp(m,n) = (el1 + 2.0)*epsi(m2,n+1)/rearth
                uvdyp(m,n)  = -rearth*epsi(m2,n+1)/(el1 + 1.0)
                vddyp(m,n)  = el1*epsi(m2,n+1)/rearth
            end do
        end do
    end

    function laplacian(input) result(output)
        complex, intent(in) :: input(mx,nx)
        complex :: output(mx,nx)

        output = -input*el2
    end function

    function inverse_laplacian(input) result(output)
        complex, intent(in) :: input(mx,nx)
        complex :: output(mx,nx)

        output = -input*elm2
    end function

    function spec_to_grid(vorm, kcos) result(vorg)
        use legendre, only: legendre_inv
        use fourier, only: gridx

        complex, intent(in) :: vorm(mx,nx)
        integer, intent(in) :: kcos

        real :: vorg(ix,il)
        real :: vorm_r(mx2,nx), varm(mx2,il)

        vorm_r = reshape(transfer(vorm, vorm_r), (/ mx2, nx /))
        call legendre_inv(vorm_r,varm)
        call gridx(varm,vorg,kcos)
    end function

    function grid_to_spec(vorg) result(vorm)
        use legendre, only: legendre_dir
        use fourier, only: specx

        real, intent(in) :: vorg(ix,il)
        complex :: vorm(mx,nx)
        real :: vorm_r(mx2,nx), varm(mx2,il)

        call specx(vorg,varm)
        call legendre_dir(varm,vorm_r)
        vorm = reshape(transfer(vorm_r, vorm), (/ mx, nx /))
    end function

    subroutine grad(psi,psdx,psdy)
        complex, dimension(mx,nx), intent(inout) :: psi
        complex, dimension(mx,nx), intent(inout) :: psdx, psdy

        integer :: k, n, m

        do n = 1, nx
            psdx(:,n) = gradx*psi(:,n)*(0.0, 1.0)
        end do

        do m=1,mx
            psdy(m,1)  =  gradyp(m,1)*psi(m,2)
            psdy(m,nx) = -gradym(m,nx)*psi(m,ntrun1)
        end do

        do n=2,ntrun1
            do m=1,mx
                psdy(m,n) = -gradym(m,n)*psi(m,n-1) + gradyp(m,n)*psi(m,n+1)
            end do
        end do
    end

    subroutine vds(ucosm,vcosm,vorm,divm)
        complex, dimension(mx,nx) :: ucosm, vcosm
        complex, dimension(mx,nx), intent(inout) :: vorm, divm
        complex, dimension(mx,nx) :: zc, zp

        integer :: n, m, k

        do n=1,nx
            zp(:,n) = gradx*ucosm(:,n)*(0.0, 1.0)
            zc(:,n) = gradx*vcosm(:,n)*(0.0, 1.0)
        end do

        do m=1,mx
            vorm(m,1)  = zc(m,1) - vddyp(m,1)*ucosm(m,2)
            vorm(m,nx) = vddym(m,nx)*ucosm(m,ntrun1)
            divm(m,1)  = zp(m,1) + vddyp(m,1)*vcosm(m,2)
            divm(m,nx) = -vddym(m,nx)*vcosm(m,ntrun1)
        end do

        do n=2,ntrun1
            do m=1,mx
                vorm(m,n) =  vddym(m,n)*ucosm(m,n-1) - vddyp(m,n)*ucosm(m,n+1) + zc(m,n)
                divm(m,n) = -vddym(m,n)*vcosm(m,n-1) + vddyp(m,n)*vcosm(m,n+1) + zp(m,n)
            end do
        end do
    end

    subroutine uvspec(vorm,divm,ucosm,vcosm)
        complex, dimension(mx,nx), intent(in) :: vorm,divm
        complex, dimension(mx,nx), intent(inout) :: ucosm,vcosm
        complex, dimension(mx,nx) :: zc,zp

        integer :: k, n, m

        zp = uvdx*vorm*(0.0, 1.0)
        zc = uvdx*divm*(0.0, 1.0)

        do m=1,mx
            ucosm(m,1)  =  zc(m,1) - uvdyp(m,1)*vorm(m,2)
            ucosm(m,nx) =  uvdym(m,nx)*vorm(m,ntrun1)
            vcosm(m,1)  =  zp(m,1) + uvdyp(m,1)*divm(m,2)
            vcosm(m,nx) = -uvdym(m,nx)*divm(m,ntrun1)
        end do

        do n=2,ntrun1
            do m=1,mx
              vcosm(m,n) = -uvdym(m,n)*divm(m,n-1) + uvdyp(m,n)*divm(m,n+1) + zp(m,n)
              ucosm(m,n) =  uvdym(m,n)*vorm(m,n-1) - uvdyp(m,n)*vorm(m,n+1) + zc(m,n)
            end do
        end do
    end

    subroutine vdspec(ug,vg,vorm,divm,kcos)
        use geometry, only: cosgr, cosgr2

        real, intent(in) :: ug(ix,il), vg(ix,il)
        complex, intent(out) :: vorm(mx,nx), divm(mx,nx)
        integer, intent(in) :: kcos
        integer :: i, j
        real :: ug1(ix,il), vg1(ix,il)
        complex :: specu(mx,nx), specv(mx,nx)

        if (kcos.eq.2) then
            do j=1,il
                do i=1,ix
                    ug1(i,j)=ug(i,j)*cosgr(j)
                    vg1(i,j)=vg(i,j)*cosgr(j)
                end do
            end do
        else
            do j=1,il
                do i=1,ix
                    ug1(i,j)=ug(i,j)*cosgr2(j)
                    vg1(i,j)=vg(i,j)*cosgr2(j)
                end do
            end do
        end if

        specu = grid_to_spec(ug1)
        specv = grid_to_spec(vg1)
        call vds(specu,specv,vorm,divm)
    end

    subroutine trunct(vor)
        complex, intent(inout) :: vor(mx,nx)

        vor = vor * trfilt
    end
end module
