module legendre
    use mod_atparam

    implicit none

    private
    public initialize_legendre, legendre_dir, legendre_inv
    public epsi

    real :: cpol(mx2,nx,iy)
    real :: consq(mxp)
    real :: epsi(mxp,nxp), repsi(mxp,nxp)
    integer :: nsh2(nx)
    real, dimension(iy) :: wt, wght

contains
    subroutine initialize_legendre
        use physical_constants, only: rearth
        use geometry, only: coa_half

        ! Initializes Legendre transforms and constants used for other
        ! subroutines that manipulate spherical harmonics

        real :: emm2, ell2, poly(mx,nx)
        integer :: j, n, m, m1, m2, mm(mx), wavenum_tot(mx,nx)

        ! first compute Gaussian latitudes and weights at the IY points from
        !     pole to equator
        ! WT(IY) are Gaussian weights for quadratures,
        !   saved in spectral
        call gaussl(wt,iy)

        ! WGHT needed for transforms saved in spectral
        do j=1,iy
            wght(j)=wt(j)/(rearth*coa_half(j))
        end do

        do n = 1, nx
            nsh2(n) = 0
            do m = 1, mx
                mm(m) = isc*(m - 1)
                wavenum_tot(m,n) = mm(m) + n - 1
                if (wavenum_tot(m,n) <= ntrun1 .or. ix /= 4*iy) nsh2(n) = nsh2(n) + 2

            end do
        end do

        do m = 2, mxp
            consq(m) = sqrt(0.5*(2.0*float(m - 1) + 1.0)/float(m - 1))
        end do

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

    end subroutine

    ! Computes inverse Legendre transformation
    subroutine legendre_inv(v_in,v_out)
        ! mx2 = 2*mx because these arrays actually represent complex variables
        real, intent(in) :: v_in(mx2,nx)
        real, intent(inout) :: v_out(mx2,il)
        real :: even(mx2),odd(mx2)

        integer :: j, j1, m, n

        ! Loop over Northern Hemisphere, computing odd and even decomposition of
        ! incoming field
        do j=1,iy
            j1=il+1-j

            ! Initialise arrays
            even = 0.0
            odd = 0.0

            ! Compute even decomposition
            do n=1,nx,2
                do m=1,nsh2(n)
                    even(m) = even(m) + v_in(m,n)*cpol(m,n,j)
                end do
            end do

            ! Compute odd decomposition
            do n=2,nx,2
                do m=1,nsh2(n)
                    odd(m) = odd(m) + v_in(m,n)*cpol(m,n,j)
                end do
            end do

            ! Compute Southern Hemisphere
            v_out(:,j1) = even + odd

            ! Compute Northern Hemisphere
            v_out(:,j)  = even - odd
        end do
    end
    !******************************************************************
    ! Computes direct Legendre transformation
    subroutine legendre_dir(v_in,v_out)
        ! mx2 = 2*mx because these arrays actually represent complex variables
        real, intent(in) :: v_in(mx2,il)
        real, intent(inout) :: v_out(mx2,nx)
        real :: even(mx2,iy), odd(mx2,iy)

        integer :: j, j1, m, n

        ! Initialise output array
        v_out = 0.0

        ! Loop over Northern Hemisphere, computing odd and even decomposition of
        ! incoming field. The Legendre weights (wt) are applied here
        do j=1,iy
            ! Corresponding Southern Hemisphere latitude
            j1=il+1-j

            even(:,j) = (v_in(:,j1) + v_in(:,j)) * wt(j)
            odd(:,j)  = (v_in(:,j1) - v_in(:,j)) * wt(j)
        end do

        ! The parity of an associated Legendre polynomial is the same
        ! as the parity of n' - m'. n', m' are the actual total wavenumber and zonal
        ! wavenumber, n and m are the indices used for SPEEDY's spectral packing.
        ! m' = m - 1 and n' = m + n - 2, therefore n' - m' = n - 1

        ! Loop over coefficients corresponding to even associated Legendre
        ! polynomials
        do n=1,ntrun1,2
            do m=1,nsh2(n)
                v_out(m,n) = dot_product(cpol(m,n,:iy), even(m,:iy))
            end do
        end do

        ! Loop over coefficients corresponding to odd associated Legendre
        ! polynomials
        do n=2,ntrun1,2
            do m=1,nsh2(n)
                v_out(m,n) = dot_product(cpol(m,n,:iy), odd(m,:iy))
            end do
        end do
    end

    subroutine gaussl(w,m)
        !   a slightly modified version of a program in Numerical Recipes
        !       (Cambridge Univ. Press, 1989)
        !   input:
        !      m    = number of gaussian latitudes between pole and equator
        !   output:
        !      w(m) = weights in gaussian quadrature (sum should equal 1.0)
        real, intent(inout) :: w(m)
        integer, intent(in) :: m
        double precision :: z,z1,p1,p2,p3,pp
        double precision, parameter :: eps=3.d-14
        integer :: n, j, i

        n = 2*m

        z1 = 2.0

        do i=1,m
            z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
            do while (abs(z-z1).gt.eps)
                p1=1.d0
                p2=0.d0

                do j=1,n
                  p3=p2
                  p2=p1
                  p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
                end do

                pp=n*(z*p1-p2)/(z*z-1.d0)
                z1=z
                z=z1-p1/pp
            end do

            w(i)=2.d0/((1.d0-z*z)*pp*pp)
        end do
    end

    subroutine lgndre(j, poly)
        use geometry, only: sia_half, coa_half

        integer, intent(in) :: j
        real, intent(inout) :: poly(mx,nx)
        real, parameter :: small = 1.e-30

        integer :: m, n, mm2
        real :: alp(mxp,nx), x, y
        y = coa_half(j)
        x = sia_half(j)

        ! start recursion with N=1 (M=L) diagonal
        alp(1,1) = sqrt(0.5)
        do m=2,mxp
            alp(m,1) = consq(m)*y*alp(m-1,1)
        end do

        ! continue with other elements
        do m=1,mxp
            alp(m,2)=(x*alp(m,1))*repsi(m,2)
        end do

        do n=3,nx
            do m=1,mxp
              alp(m,n)=(x*alp(m,n-1)-epsi(m,n-1)*alp(m,n-2))*repsi(m,n)
            end do
        end do

        ! zero polynomials with absolute values smaller than 10**(-30)
        do n=1,nx
            do m=1,mxp
                if(abs(alp(m,n)) .le. small) alp(m,n)=0.0
            end do
        end do

        ! pick off the required polynomials
        do n=1,nx
            do m=1,mx
                mm2=isc*(m-1)+1
                poly(m,n)=alp(mm2,n)
            end do
        end do
    end
end module
