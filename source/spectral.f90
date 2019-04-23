module spectral
    use mod_atparam

    implicit none

    private
    public el2, elm2, el4, trfilt, nsh2, sia, coa, wt, wght, cosg,&
        & cosgr, cosgr2, gradx, gradym, gradyp, consq, epsi, repsi,&
        & cpol, uvdx, uvdym, uvdyp, vddym, vddyp
    public wsave
    public initialize_spectral
    public laplacian, inverse_laplacian, spec_to_grid, grid_to_spec

    real, dimension(mx,nx) :: el2, elm2, el4, trfilt
    integer :: nsh2(nx)
    real, dimension(iy) :: sia, coa, wt, wght
    real, dimension(il) :: cosg, cosgr, cosgr2
    real :: gradx(mx), gradym(mx,nx), gradyp(mx,nx)
    real :: consq(mxp), epsi(mxp,nxp), repsi(mxp,nxp)
    real :: cpol(mx2,nx,iy)
    real, dimension(mx,nx) :: uvdx, uvdym, uvdyp
    real, dimension(mx,nx) :: vddym, vddyp

    real, dimension(2*ix+15) :: wsave(2*ix+15)

contains
    ! Initialize spectral transforms
    subroutine initialize_spectral
        use mod_dyncon1, only: rearth

        real :: el1, ell2, emm2, poly(mx,nx)

        integer :: j, jj, m, m1, m2, n, l2(mx,nx), wavenum_tot(mx,nx), mm(mx)

        call rffti(ix,wsave)

        ! Initializes Legendre transforms and constants used for other
        ! subroutines that manipulate spherical harmonics
        !
        ! first compute Gaussian latitudes and weights at the IY points from
        !     pole to equator
        ! SIA(IY) is sin of latitude, WT(IY) are Gaussian weights for quadratures,
        call gaussl(sia, wt, iy)

        ! COA(IY) = cos(lat); WGHT needed for transforms,
        coa = sqrt(1.0 - sia**2.0)
        wght = wt/(rearth*coa**2.0)

        ! expand cosine and its reciprocal to cover both hemispheres,
        do j = 1, iy
            jj = il+1-j
            cosg(j)   = coa(j)
            cosg(jj)  = coa(j)
            cosgr(j)  = 1.0/coa(j)
            cosgr(jj) = 1.0/coa(j)
            cosgr2(j) = 1.0/(coa(j)*coa(j))
            cosgr2(jj)= 1.0/(coa(j)*coa(j))
        end do

        !  MM = zonal wavenumber = m
        !     ISC=3 implies that only wavenumber 0,3,6,9,etc are included in model
        !  wavenum_tot = total wavenumber of spherical harmonic = l
        !  L2 = l*(l+1)
        !  EL2 = l*(l+1)/(a**2)
        !  EL4 = EL2*EL2 ; for biharmonic diffusion
        !  ELM2 = 1./EL2
        !  TRFILT used to filter out "non-triangular" part of rhomboidal truncation
        do n = 1, nx
            nsh2(n) = 0
            do m = 1, mx
                mm(m) = isc*(m - 1)
                wavenum_tot(m,n) = mm(m) + n - 1
                el2(m,n) = float(wavenum_tot(m,n)*(wavenum_tot(m,n) + 1))/(rearth**2.0)
                el4(m,n) = el2(m,n)**2.0
                if (wavenum_tot(m,n) <= ntrun1 .or. ix /= 4*iy) nsh2(n) = nsh2(n) + 2
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

        ! quantities needed to generate and differentiate Legendre polynomials
        ! all m values up to MXP = ISC*MTRUN+1 are needed by recursion relation
        ! saved in spectral
        do m = 1, mxp
            do n = 1, nxp
                emm2 = float(m - 1)**2
                ell2 = float(n + m - 2)**2
                if (n == nxp) then
                    epsi(m,n) = 0.0
                else if(n == 1 .and. m == 1) then
                    epsi(m,n) = 0.0
                else
                    epsi(m,n)=sqrt((ell2 - emm2)/(4.0*ell2 - 1.0))
                end if
                repsi(m,n) = 0.0
                if (epsi(m,n) > 0.) repsi(m,n) = 1.0/epsi(m,n)
            end do
        end do

        do m = 2, mxp
            consq(m) = sqrt(0.5*(2.0*float(m - 1) + 1.0)/float(m - 1))
        end do

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

        !  generate associated Legendre polynomial
        !  LGNDRE computes the polynomials at a particular latitiude, POLY(MX,NX)
        !  polynomials and 'clones' stored in spectral
        do j = 1, iy
            call lgndre(j, poly)
            do n = 1, nx
                do m = 1, mx
                    m1 = 2*m - 1
                    m2 = 2*m
                    cpol(m1,n,j) = poly(m,n)
                    cpol(m2,n,j) = poly(m,n)
                end do
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
        complex, intent(in) :: vorm(mx,nx)
        integer, intent(in) :: kcos

        real :: vorg(ix,il)
        real :: vorm_r(mx2,nx), varm(mx2,il)

        vorm_r = reshape(transfer(vorm, vorm_r), (/ mx2, nx /))
        call legendre_inv(vorm_r,varm)
        call gridx(varm,vorg,kcos)
    end function

    function grid_to_spec(vorg) result(vorm)
        real, intent(in) :: vorg(ix,il)
        complex :: vorm(mx,nx)
        real :: vorm_r(mx2,nx), varm(mx2,il)

        call specx(vorg,varm)
        call legendre_dir(varm,vorm_r)
        vorm = reshape(transfer(vorm_r, vorm), (/ mx, nx /))
    end function
end module
