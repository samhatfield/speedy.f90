! Compute convective fluxes of dry static energy and moisture using a
! simplified mass-flux scheme
!  Input:   psa    = norm. surface pressure [p/p0]            (2-dim)
!           se     = dry static energy                        (3-dim)
!           qa     = specific humidity [g/kg]                 (3-dim)
!           qsat   = saturation spec. hum. [g/kg]             (3-dim)
!  Output:  itop   = top of convection (layer index)          (2-dim)
!           cbmf   = cloud-base mass flux                     (2-dim)
!           precnv = convective precipitation [g/(m^2 s)]     (2-dim)
!           dfse   = net flux of d.s.en. into each atm. layer (3-dim)
!           dfqa   = net flux of sp.hum. into each atm. layer (3-dim)
subroutine convmf (psa, se, qa, qsat, itop, cbmf, precnv, dfse, dfqa)
    use mod_cnvcon
    use mod_atparam
    use mod_physcon, only: p0, alhc, alhs, sig, dsig, wvi
    use mod_dyncon1, only: grav

    implicit none

    real, intent(in) :: psa(ix,il), se(ix,il,kx), qa(ix,il,kx), qsat(ix,il,kx)
    integer, intent(inout) :: itop(ix,il)
    real, intent(inout) :: cbmf(ix,il), precnv(ix,il), dfse(ix,il,kx), dfqa(ix,il,kx)

    integer :: i, j, k, k1, ktop1, ktop2, nl1, nlp
    real :: mss(ix,il,2:kx), mse0, mse1, mss0, mss2, msthr, qdif(ix,il)
    real :: entr(2:kx-1), delq, enmass, fdq, fds, fm0, fmass, fpsa, fqmax
    real :: fsq, fuq, fus, qb, qmax, qsatb, qthr0, qthr1, rdps, rlhc, sb, sentr
    logical :: lqthr

    ! 1. Initialization of output and workspace arrays
    nl1 = kx - 1
    nlp = kx + 1
    fqmax = 5.0

    fm0 = p0*dsig(kx)/(grav*trcnv*3600.0)
    rdps=2.0/(1.0 - psmin)

    dfse = 0.0
    dfqa = 0.0

    cbmf = 0.0
    precnv = 0.0

    ! Saturation moist static energy
    do k = 2, kx
        mss(:,:,k) = se(:,:,k) + alhc*qsat(:,:,k)
    end do

    ! Entrainment profile (up to sigma = 0.5)
    sentr = 0.0
    do k = 2, nl1
        entr(k) = (max(0.0, sig(k) - 0.5))**2.0
        sentr = sentr + entr(k)
    end do

    sentr = entmax/sentr
    entr(2:nl1) = entr(2:nl1) * sentr

    ! 2. Check of conditions for convection
    rlhc = 1.0/alhc

    do i = 1, ix
        do j = 1, il
            itop(i,j) = nlp

            if (psa(i,j) > psmin) then
                ! Minimum of moist static energy in the lowest two levels
                mse0 = se(i,j,kx) + alhc*qa(i,j,kx)
                mse1 = se(i,j,nl1) + alhc*qa(i,j,nl1)
                mse1 = min(mse0, mse1)

                ! Saturation (or super-saturated) moist static energy in PBL
                mss0 = max(mse0, mss(i,j,kx))

                ktop1 = kx
                ktop2 = kx

                do k = kx-3, 3, -1
                    mss2 = mss(i,j,k) + wvi(k,2)*(mss(i,j,k+1) - mss(i,j,k))

                    ! Check 1: conditional instability
                    !          (MSS in PBL > MSS at top level)
                    if (mss0 > mss2) then
                       ktop1 = k
                    end if

                    ! Check 2: gradient of actual moist static energy
                    !          between lower and upper troposphere
                    if (mse1 > mss2) then
                       ktop2 = k
                       msthr = mss2
                    end if
                end do

                if (ktop1 < kx) then
                    ! Check 3: RH > RH_c at both k=NLEV and k=NL1
                    qthr0 = rhbl*qsat(i,j,kx)
                    qthr1 = rhbl*qsat(i,j,nl1)
                    lqthr = (qa(i,j,kx) > qthr0 .and. qa(i,j,nl1) > qthr1)

                    if (ktop2 < kx) then
                       itop(i,j) = ktop1
                       qdif(i,j) = max(qa(i,j,kx) - qthr0, (mse0 - msthr)*rlhc)
                    else if (lqthr) then
                       itop(i,j) = ktop1
                       qdif(i,j) = qa(i,j,kx) - qthr0
                    end if
                end if
            end if
        end do
    end do

    ! 3. Convection over selected grid-points
    do i = 1, ix
        do j = 1, il
            if (itop(i,j) == nlp) cycle

            ! 3.1 Boundary layer (cloud base)
            k = kx
            k1 = k - 1

            ! Maximum specific humidity in the PBL
            qmax = max(1.01*qa(i,j,k), qsat(i,j,k))

            ! Dry static energy and moisture at upper boundary
            sb = se(i,j,k1) + wvi(k1,2)*(se(i,j,k) - se(i,j,k1))
            qb = qa(i,j,k1) + wvi(k1,2)*(qa(i,j,k) - qa(i,j,k1))
            qb = min(qb, qa(i,j,k))

            ! Cloud-base mass flux, computed to satisfy:
            ! fmass*(qmax-qb)*(g/dp)=qdif/trcnv
            fpsa = psa(i,j)*min(1.0, (psa(i,j) - psmin)*rdps)
            fmass = fm0*fpsa*min(fqmax, qdif(i,j)/(qmax - qb))
            cbmf(i,j) = fmass

            ! Upward fluxes at upper boundary
            fus = fmass*se(i,j,k)
            fuq = fmass*qmax

            ! Downward fluxes at upper boundary
            fds = fmass*sb
            fdq = fmass*qb

            ! Net flux of dry static energy and moisture
            dfse(i,j,k) = fds - fus
            dfqa(i,j,k) = fdq - fuq

            ! 3.2 Intermediate layers (entrainment)
            do k = kx - 1, itop(i, j) + 1, -1
                k1 = k - 1

                ! Fluxes at lower boundary
                dfse(i,j,k) = fus - fds
                dfqa(i,j,k) = fuq - fdq

                ! Mass entrainment
                enmass = entr(k)*psa(i,j)*cbmf(i,j)
                fmass = fmass + enmass

                ! Upward fluxes at upper boundary
                fus = fus + enmass*se(i,j,k)
                fuq = fuq + enmass*qa(i,j,k)

                ! Downward fluxes at upper boundary
                sb = se(i,j,k1) + wvi(k1,2)*(se(i,j,k) - se(i,j,k1))
                qb = qa(i,j,k1) + wvi(k1,2)*(qa(i,j,k) - qa(i,j,k1))
                fds = fmass*sb
                fdq = fmass*qb

                ! Net flux of dry static energy and moisture
                dfse(i,j,k) = dfse(i,j,k) + fds - fus
                dfqa(i,j,k) = dfqa(i,j,k) + fdq - fuq

                ! Secondary moisture flux
                delq = rhil*qsat(i,j,k) - qa(i,j,k)
                if (delq > 0.0) then
                    fsq = smf*cbmf(i,j)*delq
                    dfqa(i,j,k)   = dfqa(i,j,k) + fsq
                    dfqa(i,j,kx) = dfqa(i,j,kx) - fsq
                end if
            end do

            ! 3.3 Top layer (condensation and detrainment)
            k = itop(i,j)

            ! Flux of convective precipitation
            qsatb = qsat(i,j,k) + wvi(k,2)*(qsat(i,j,k+1) - qsat(i,j,k))

            precnv(i,j) = max(fuq - fmass*qsatb, 0.0)

            ! Net flux of dry static energy and moisture
            dfse(i,j,k) = fus - fds + alhc*precnv(i,j)
            dfqa(i,j,k) = fuq - fdq - precnv(i,j)
        end do
    end do
end
