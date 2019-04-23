!******************************************************************
subroutine gaussl(x,w,m)
    !   a slightly modified version of a program in Numerical Recipes
    !       (Cambridge Univ. Press, 1989)
    !   input:
    !      m    = number of gaussian latitudes between pole and equator
    !   output:
    !      x(m) = sin(gaussian latitude)
    !      w(m) = weights in gaussian quadrature (sum should equal 1.0)

    implicit none

    real, intent(inout) :: x(m),w(m)
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

        x(i)=z
        w(i)=2.d0/((1.d0-z*z)*pp*pp)
    end do
end
!****************************************************************
subroutine lgndre(j, poly)
    ! follows Leith Holloways code

    use mod_atparam
    use spectral, only: sia, coa, consq, repsi, epsi

    implicit none

    integer, intent(in) :: j
    real, intent(inout) :: poly(mx,nx)
    real, parameter :: small = 1.e-30

    integer :: m, n, mm2
    real :: alp(mxp,nx), x, y
    y = coa(j)
    x = sia(j)

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
!*********************************************************************
subroutine grad(psi,psdx,psdy)
    use mod_atparam
    use spectral, only: gradx, gradyp, gradym

    implicit none

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
!******************************************************************
subroutine vds(ucosm,vcosm,vorm,divm)
    use mod_atparam
    use spectral, only: gradx, vddyp, vddym

    implicit none

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
!******************************************************************
subroutine uvspec(vorm,divm,ucosm,vcosm)
    use mod_atparam
    use spectral, only: uvdx, uvdyp, uvdym

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
!*********************************************************************
subroutine vdspec(ug,vg,vorm,divm,kcos)
    use mod_atparam
    use spectral, only: cosgr, cosgr2, grid_to_spec

    implicit none

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
!*********************************************************************
! Computes inverse Legendre transformation
subroutine legendre_inv(v_in,v_out)
    use mod_atparam
    use spectral, only: cpol, nsh2

    implicit none

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
    use mod_atparam
    use spectral, only: wt, cpol, nsh2

    implicit none

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
!******************************************************************
subroutine trunct(vor)
    use mod_atparam
    use spectral, only: trfilt

    implicit none

    complex, intent(inout) :: vor(mx,nx)

    vor = vor * trfilt
end
