subroutine inphys(hsg,ppl,rlat)
    !
    ! subroutine inphys (hsg,ppl,rlat)
    !
    ! Purpose: Initialize common blocks for physical parametrization routines 
    ! Input :  hsg  : sigma at half levels
    !          ppl  : pressure levels for post-processing
    !          rlat : gaussian-grid latitudes
    ! Initialized common blocks: phycon, fsiglt, forcon, 
    !                            cnvcon, lsccon, radcon, sflcon, vdicon

    use mod_atparam
    use mod_physcon

    implicit none

    integer, parameter :: nlon = ix, nlat = il, nlev = kx, ngp = nlon*nlat

    real :: hsg(0:nlev), ppl(nlev), rlat(nlat)  
    integer :: j, k
    
    ! 1.2 Functions of sigma and latitude
    sigh(0) = hsg(0)

    do k = 1, nlev
        sig(k)  = 0.5*(hsg(k)+hsg(k-1))
        sigl(k) = log(sig(k))
        sigh(k) = hsg(k)
        dsig(k) = hsg(k)-hsg(k-1)
        pout(k) = ppl(k)
        grdsig(k) = gg/(dsig(k)*p0)
        grdscp(k) = grdsig(k)/cp
    end do

    ! Weights for vertical interpolation at half-levels(1,nlev) and surface
    ! Note that for phys.par. half-lev(k) is between full-lev k and k+1 
    ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
    ! Fsurf = Ffull(nlev)+WVI(nlev,2)*(Ffull(nlev)-Ffull(nlev-1))
    do k = 1, nlev-1
        wvi(k,1) = 1./(sigl(k+1)-sigl(k))
        wvi(k,2) = (log(sigh(k))-sigl(k))*wvi(k,1)
    end do

    wvi(nlev,1) = 0.
    wvi(nlev,2) = (log(0.99)-sigl(nlev))*wvi(nlev-1,1)

    do j = 1, nlat
        slat(j) = sin(rlat(j))
        clat(j) = cos(rlat(j))
    end do
end
