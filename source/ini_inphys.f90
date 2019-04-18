! Initialize physical parametrization routines
subroutine inphys
    use mod_atparam
    use mod_physcon
    use mod_dyncon1, only: grav, hsg, radang

    implicit none

    integer :: j, k

    ! 1.2 Functions of sigma and latitude
    sigh(0) = hsg(1)

    do k = 1, kx
        sig(k)  = 0.5*(hsg(k+1)+hsg(k))
        sigl(k) = log(sig(k))
        sigh(k) = hsg(k+1)
        dsig(k) = hsg(k+1)-hsg(k)
        grdsig(k) = grav/(dsig(k)*p0)
        grdscp(k) = grdsig(k)/cp
    end do

    ! Weights for vertical interpolation at half-levels(1,kx) and surface
    ! Note that for phys.par. half-lev(k) is between full-lev k and k+1
    ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
    ! Fsurf = Ffull(kx)+WVI(kx,2)*(Ffull(kx)-Ffull(kx-1))
    do k = 1, kx-1
        wvi(k,1) = 1./(sigl(k+1)-sigl(k))
        wvi(k,2) = (log(sigh(k))-sigl(k))*wvi(k,1)
    end do

    wvi(kx,1) = 0.
    wvi(kx,2) = (log(0.99)-sigl(kx))*wvi(kx-1,1)

    do j = 1, il
        slat(j) = sin(radang(j))
        clat(j) = cos(radang(j))
    end do
end
