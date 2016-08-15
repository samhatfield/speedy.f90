subroutine vdifsc(ua,va,se,rh,qa,qsat,phi,icnv,utenvd,vtenvd,ttenvd,qtenvd)
    !   subroutine vdifsc (ua,va,se,rh,qa,qsat,phi,icnv,
    !  &                   utenvd,vtenvd,ttenvd,qtenvd)
    !
    !   Purpose: Compute tendencies of momentum, energy and moisture
    !            due to vertical diffusion and shallow convection
    !   Input:   ua     = u-wind                           (3-dim)
    !            va     = v-wind                           (3-dim)
    !            se     = dry static energy                (3-dim)
    !            rh     = relative humidity [0-1]          (3-dim)
    !            qa     = specific humidity [g/kg]         (3-dim)
    !            qsat   = saturation sp. humidity [g/kg]   (3-dim)
    !            phi    = geopotential                     (3-dim)
    !            icnv   = index of deep convection         (2-dim)
    !   Output:  utenvd = u-wind tendency                  (3-dim)
    !            vtenvd = v-wind tendency                  (3-dim)
    !            ttenvd = temperature tendency             (3-dim)
    !            qtenvd = sp. humidity tendency [g/(kg s)] (3-dim)
    !

    use mod_atparam
    use mod_vdicon
    use mod_physcon, only: cp, alhc, sig, sigh, dsig

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real, dimension(ngp,nlev), intent(in) :: ua, va, se, rh, qa, qsat, phi
    integer, intent(in) :: icnv(ngp)
    real, dimension(ngp,nlev), intent(inout) :: utenvd, vtenvd, ttenvd, qtenvd

    integer :: nl1, j, k, k1
    real :: cshc, cvdi, fshcq, fshcse, fvdiq, fvdise, drh0, fvdiq2, dmse, drh
    real :: fluxse, fluxq, fcnv, se0
    real, dimension(nlev) :: rsig, rsig1

    ! 1. Initalization

    ! N.B. In this routine, fluxes of dry static energy and humidity
    !      are scaled in such a way that:
    !      d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma

    nl1  = nlev-1
    cshc = dsig(nlev)/3600.
    cvdi = (sigh(nl1)-sigh(1))/((nl1-1)*3600.)

    fshcq  = cshc/trshc
    fshcse = cshc/(trshc*cp)

    fvdiq  = cvdi/trvdi
    fvdise = cvdi/(trvds*cp)

    do k=1,nl1
        rsig(k)=1./dsig(k)
	    rsig1(k)=1./(1.-sigh(k))
    end do
    rsig(nlev)=1./dsig(nlev)

    utenvd = 0.0
    vtenvd = 0.0
    ttenvd = 0.0
    qtenvd = 0.0
   
    ! 2. Shallow convection
    drh0   = rhgrad*(sig(nlev)-sig(nl1))
    fvdiq2 = fvdiq*sigh(nl1)

    do j=1,ngp
        dmse = (se(j,nlev)-se(j,nl1))+alhc*(qa(j,nlev)-qsat(j,nl1))
        drh  = rh(j,nlev)-rh(j,nl1)
        fcnv = 1.

        if (dmse.ge.0.0) then
            if (icnv(j).gt.0) fcnv = redshc
  
            fluxse         = fcnv*fshcse*dmse
            ttenvd(j,nl1)  = fluxse*rsig(nl1)
            ttenvd(j,nlev) =-fluxse*rsig(nlev)
  
            if (drh.ge.0.0) then
                fluxq          = fcnv*fshcq*qsat(j,nlev)*drh
                qtenvd(j,nl1)  = fluxq*rsig(nl1) 
                qtenvd(j,nlev) =-fluxq*rsig(nlev)
            end if
        else if (drh.ge.drh0) then
          fluxq          = fvdiq2*qsat(j,nl1)*drh
          qtenvd(j,nl1)  = fluxq*rsig(nl1) 
          qtenvd(j,nlev) =-fluxq*rsig(nlev)
        end if
    end do

    ! 3. Vertical diffusion of moisture above the PBL
    do k=3,nlev-2
        if (sigh(k).gt.0.5) then
            drh0   = rhgrad*(sig(k+1)-sig(k))
            fvdiq2 = fvdiq*sigh(k)

            do j=1,ngp
                drh=rh(j,k+1)-rh(j,k)
                if (drh.ge.drh0) then
                    fluxq        = fvdiq2*qsat(j,k)*drh
                    qtenvd(j,k)  = qtenvd(j,k)  +fluxq*rsig(k)
                    qtenvd(j,k+1)= qtenvd(j,k+1)-fluxq*rsig(k+1)
                end if
            end do
        end if
    end do

    ! 4. Damping of super-adiabatic lapse rate
    do k=1,nl1
        do j=1,ngp
            se0 = se(j,k+1)+segrad*(phi(j,k)-phi(j,k+1))
  
            if (se(j,k).lt.se0) then
                fluxse      = fvdise*(se0-se(j,k))
                ttenvd(j,k) = ttenvd(j,k)+fluxse*rsig(k)
                do k1=k+1,nlev
                    ttenvd(j,k1) = ttenvd(j,k1)-fluxse*rsig1(k)
                end do
            end if
        end do
    end do
end
