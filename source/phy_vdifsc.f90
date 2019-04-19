! Compute tendencies of momentum, energy and moisture due to vertical diffusion
! and shallow convection
! Input:   ua     = u-wind                           (3-dim)
!          va     = v-wind                           (3-dim)
!          se     = dry static energy                (3-dim)
!          rh     = relative humidity [0-1]          (3-dim)
!          qa     = specific humidity [g/kg]         (3-dim)
!          qsat   = saturation sp. humidity [g/kg]   (3-dim)
!          phi    = geopotential                     (3-dim)
!          icnv   = index of deep convection         (2-dim)
! Output:  utenvd = u-wind tendency                  (3-dim)
!          vtenvd = v-wind tendency                  (3-dim)
!          ttenvd = temperature tendency             (3-dim)
!          qtenvd = sp. humidity tendency [g/(kg s)] (3-dim)
subroutine vdifsc(ua, va, se, rh, qa, qsat, phi, icnv, utenvd, vtenvd, ttenvd, qtenvd)
    use mod_atparam
    use mod_vdicon
    use mod_physcon, only: cp, alhc, sig, sigh, dsig

    implicit none

    real, dimension(ix,il,kx), intent(in) :: ua, va, se, rh, qa, qsat, phi
    integer, intent(in) :: icnv(ix,il)
    real, dimension(ix,il,kx), intent(inout) :: utenvd, vtenvd, ttenvd, qtenvd

    integer :: nl1, i, j, k, k1
    real :: cshc, cvdi, fshcq, fshcse, fvdiq, fvdise, drh0, fvdiq2, dmse, drh
    real :: fluxse, fluxq, fcnv, se0
    real, dimension(kx) :: rsig, rsig1

    ! 1. Initalization

    ! N.B. In this routine, fluxes of dry static energy and humidity
    !      are scaled in such a way that:
    !      d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma

    nl1  = kx - 1
    cshc = dsig(kx)/3600.0
    cvdi = (sigh(nl1) - sigh(1))/((nl1 - 1)*3600.0)

    fshcq  = cshc/trshc
    fshcse = cshc/(trshc*cp)

    fvdiq  = cvdi/trvdi
    fvdise = cvdi/(trvds*cp)

    do k = 1, nl1
        rsig(k) = 1.0/dsig(k)
	    rsig1(k) = 1.0/(1.0 - sigh(k))
    end do
    rsig(kx)=1.0/dsig(kx)

    utenvd = 0.0
    vtenvd = 0.0
    ttenvd = 0.0
    qtenvd = 0.0

    ! 2. Shallow convection
    drh0   = rhgrad*(sig(kx) - sig(nl1))
    fvdiq2 = fvdiq*sigh(nl1)

    do i = 1, ix
        do j = 1, il
            dmse = se(i,j,kx) - se(i,j,nl1) + alhc*(qa(i,j,kx) - qsat(i,j,nl1))
            drh  = rh(i,j,kx) - rh(i,j,nl1)
            fcnv = 1.0

            if (dmse >= 0.0) then
                if (icnv(i,j) > 0) fcnv = redshc

                fluxse           =  fcnv*fshcse*dmse
                ttenvd(i,j,nl1)  =  fluxse*rsig(nl1)
                ttenvd(i,j,kx)   = -fluxse*rsig(kx)

                if (drh >= 0.0) then
                    fluxq           =  fcnv*fshcq*qsat(i,j,kx)*drh
                    qtenvd(i,j,nl1) =  fluxq*rsig(nl1)
                    qtenvd(i,j,kx)  = -fluxq*rsig(kx)
                end if
            else if (drh > drh0) then
                fluxq           =  fvdiq2*qsat(i,j,nl1)*drh
                qtenvd(i,j,nl1) =  fluxq*rsig(nl1)
                qtenvd(i,j,kx)  = -fluxq*rsig(kx)
            end if
        end do
    end do

    ! 3. Vertical diffusion of moisture above the PBL
    do k = 3, kx - 2
        if (sigh(k) > 0.5) then
            drh0   = rhgrad*(sig(k+1) - sig(k))
            fvdiq2 = fvdiq*sigh(k)

            do i = 1, ix
                do j = 1, il
                    drh = rh(i,j,k+1) - rh(i,j,k)
                    if (drh >= drh0) then
                        fluxq           = fvdiq2*qsat(i,j,k)*drh
                        qtenvd(i,j,k)   = qtenvd(i,j,k) + fluxq*rsig(k)
                        qtenvd(i,j,k+1) = qtenvd(i,j,k+1) - fluxq*rsig(k+1)
                    end if
                end do
            end do
        end if
    end do

    ! 4. Damping of super-adiabatic lapse rate
    do k = 1, nl1
        do i = 1, ix
            do j = 1, il
                se0 = se(i,j,k+1) + segrad*(phi(i,j,k) - phi(i,j,k+1))

                if (se(i,j,k) < se0) then
                    fluxse        = fvdise*(se0 - se(i,j,k))
                    ttenvd(i,j,k) = ttenvd(i,j,k) + fluxse*rsig(k)
                    do k1 = k+1, kx
                        ttenvd(i,j,k1) = ttenvd(i,j,k1) - fluxse*rsig1(k)
                    end do
                end if
            end do
        end do
    end do
end
