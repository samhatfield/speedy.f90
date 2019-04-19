subroutine sol_oz(tyear)
    !  subroutine sol_oz (tyear)
    !
    !  Purpose: Compute zonally-averaged fields to be used
    !           in the computation of SW absorption:
    !           fsol   = flux of incoming solar radiation
    !           ozone  = flux absorbed by ozone (lower stratos.)
    !           ozupp  = flux absorbed by ozone (upper stratos.)
    !           zenit  = function of solar zenith angle
    !  Input:   tyear  = time as fraction of year (0-1, 0 = 1jan.h00)

    use mod_atparam
    use mod_physcon, only: slat, clat
    use mod_radcon, only: solc, epssw, fsol, ozone, ozupp, zenit, stratz

    implicit none

    real, intent(in) :: tyear
    real :: topsr(ix), alpha, azen, coz1, coz2, dalpha, flat2, fs0
    real :: nzen, rzen
    integer :: j

    ! alpha = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
    alpha = 4.0*asin(1.0)*(tyear + 10.0/365.0)
    dalpha = 0.0

    coz1 = 1.0*max(0.0, cos(alpha - dalpha))
    coz2 = 1.8

    azen = 1.0
    nzen = 2

    rzen = -cos(alpha)*23.45*asin(1.0)/90.0

    fs0 = 6.0

    ! Solar radiation at the top
    call solar(tyear, 4.0*solc, ix, clat, slat, topsr)

    do j = 1, il
        flat2 = 1.5*slat(j)**2 - 0.5

        ! Solar radiation at the top
        fsol(:,j) = topsr(j)

        ! Ozone depth in upper and lower stratosphere
        ozupp(:,j) = 0.5*epssw
        ozone(:,j) = 0.4*epssw*(1.0 + coz1*slat(j) + coz2*flat2)

        ! Zenith angle correction to (downward) absorptivity
        zenit(:,j) = 1.0 + azen*(1.0 - (clat(j)*cos(rzen) + slat(j)*sin(rzen)))**nzen

        ! Ozone absorption in upper and lower stratosphere
        ozupp(:,j) = fsol(:,j)*ozupp(:,j)*zenit(:,j)
        ozone(:,j) = fsol(:,j)*ozone(:,j)*zenit(:,j)

        ! Polar night cooling in the stratosphere
        stratz(:,j) = max(fs0 - fsol(:,j), 0.0)
    end do
end

subroutine solar(tyear,csol,nlat,clat,slat,topsr)
    ! Average daily flux of solar radiation, from Hartmann (1994)

    implicit none

    real, intent(in) :: tyear, csol
    integer, intent(in) :: nlat
    real, dimension(nlat), intent(in) :: clat, slat
    real, intent(inout) :: topsr(nlat)

    integer :: j
    real :: ca1, ca2, ca3, cdecl, ch0, csolp, decl, fdis, h0, alpha, pigr, sa1
    real :: sa2, sa3, sdecl, sh0, tdecl

    ! 1. Compute declination angle and Earth-Sun distance factor
    pigr  = 2.*asin(1.)
    alpha = 2.*pigr*tyear

    ca1 = cos(alpha)
    sa1 = sin(alpha)
    ca2 = ca1*ca1-sa1*sa1
    sa2 = 2.*sa1*ca1
    ca3 = ca1*ca2-sa1*sa2
    sa3 = sa1*ca2+sa2*ca1

    decl = 0.006918-0.399912*ca1+0.070257*sa1-0.006758*ca2+0.000907*sa2&
        & -0.002697*ca3+0.001480*sa3

    fdis = 1.000110+0.034221*ca1+0.001280*sa1+0.000719*ca2+0.000077*sa2

    cdecl = cos(decl)
    sdecl = sin(decl)
    tdecl = sdecl/cdecl

    ! 2. Compute daily-average insolation at the atm. top
    csolp=csol/pigr

    do j=1,nlat
        ch0 = min(1.,max(-1.,-tdecl*slat(j)/clat(j)))
        h0  = acos(ch0)
        sh0 = sin(h0)

        topsr(j) = csolp*fdis*(h0*slat(j)*sdecl+sh0*clat(j)*cdecl)
    end do
end

subroutine cloud(qa,rh,precnv,precls,iptop,gse,fmask,icltop,cloudc,clstr)
    !  subroutine cloud (qa,rh,precnv,precls,iptop,gse,fmask,
    ! &                  icltop,cloudc,clstr)
    !
    !  Purpose: Compute cloud-top level and cloud cover
    !  Input:   qa     = specific humidity [g/kg]                (3-dim)
    !           rh     = relative humidity                       (3-dim)
    !           precnv = convective precipitation                (2-dim)
    !           precls = large-scale precipitation               (2-dim)
    !           iptop  = top level of precipitating cloud        (2-dim)
    !           gse    = gradient of dry st. energy (dSE/dPHI)   (2-dim)
    !           fmask  = fractional land-sea mask                (2-dim)
    !  Output:  icltop = cloud top level (all clouds)            (2-dim)
    !           cloudc = total cloud cover                       (2-dim)
    !           clstr  = stratiform cloud cover                  (2-dim)

    use mod_atparam
    use mod_radcon, only: rhcl1, rhcl2, qacl, wpcl, pmaxcl, clsmax, clsminl,&
        & gse_s0, gse_s1, albcl, qcloud

    implicit none

    integer :: iptop(ix,il)
    real, intent(in) :: qa(ix,il,kx), rh(ix,il,kx), precnv(ix,il), precls(ix,il), gse(ix,il),&
        & fmask(ix,il)
    real, intent(inout) :: cloudc(ix,il), clstr(ix,il)
    integer, intent(inout) :: icltop(ix,il)

    integer :: inew, i, j, k, nl1, nlp
    real :: clfact, clstrl, drh, fstab, pr1, rgse, rrcl

    nl1  = kx-1
    nlp  = kx+1
    rrcl = 1./(rhcl2-rhcl1)

    ! 1.  Cloud cover, defined as the sum of:
    !     - a term proportional to the square-root of precip. rate
    !     - a quadratic function of the max. relative humidity
    !       in tropospheric layers above PBL where Q > QACL :
    !       ( = 0 for RHmax < RHCL1, = 1 for RHmax > RHCL2 )
    !     Cloud-top level: defined as the highest (i.e. least sigma)
    !       between the top of convection/condensation and
    !       the level of maximum relative humidity.

    do i = 1, ix
        do j = 1, il
            if (rh(i,j,nl1) > rhcl1) then
                cloudc(i,j) = rh(i,j,nl1) - rhcl1
                icltop(i,j) = nl1
            else
                cloudc(i,j) = 0.0
                icltop(i,j) = nlp
            end if
        end do
    end do

    do k = 3, kx - 2
        do i = 1, ix
            do j = 1, il
                drh = rh(i,j,k) - rhcl1
                if (drh > cloudc(i,j) .and. qa(i,j,k) > qacl) then
                    cloudc(i,j) = drh
                    icltop(i,j) = k
                end if
            end do
        end do
    end do

    do i = 1, ix
        do j = 1, il
            pr1 = min(pmaxcl, 86.4*(precnv(i,j) + precls(i,j)))
            cloudc(i,j) = min(1.0, wpcl*sqrt(pr1) + min(1.0, cloudc(i,j)*rrcl)**2.0)
            icltop(i,j) = min(iptop(i,j), icltop(i,j))
        end do
    end do

    ! 2.  Equivalent specific humidity of clouds
    qcloud = qa(:,:,nl1)

    ! 3. Stratiform clouds at the top of PBL
    inew = 1

    if (inew > 0) then
        clfact = 1.2
        rgse   = 1.0/(gse_s1 - gse_s0)

        do i = 1, ix
            do j = 1, il
                ! Stratocumulus clouds over sea
                fstab    = max(0.0, min(1.0, rgse*(gse(i,j) - gse_s0)))
                clstr(i,j) = fstab*max(clsmax - clfact*cloudc(i,j), 0.0)

                ! Stratocumulus clouds over land
                clstrl     = max(clstr(i,j), clsminl)*rh(i,j,kx)
                clstr(i,j) = clstr(i,j) + fmask(i,j)*(clstrl - clstr(i,j))
            end do
        end do
    else
        clsmax  = 0.3
        clsminl = 0.1

        ! Stratocumulus clouds over sea
        clstr = max(clsmax - cloudc, 0.0)

        ! Rescale for consistency with previous albedo values
        clstr = clstr*albcl/0.5

        ! Correction for aerosols over land
        clstr = clstr + fmask*(clsminl - clstr)
    end if
end

subroutine radsw(psa,qa,icltop,cloudc,clstr,fsfcd,fsfc,ftop,dfabs)
    !  subroutine radsw (psa,qa,icltop,cloudc,clstr,
    ! &                  fsfcd,fsfc,ftop,dfabs)
    !
    !  purpose: compute the absorption of shortwave radiation and
    !           initialize arrays for longwave-radiation routines
    !  input:   psa    = norm. surface pressure [p/p0]           (2-dim)
    !           qa     = specific humidity [g/kg]                (3-dim)
    !           icltop = cloud top level                         (2-dim)
    !           cloudc = total cloud cover                       (2-dim)
    !           clstr  = stratiform cloud cover                  (2-dim)
    !  output:  fsfcd  = downward-only flux of sw rad. at the surface (2-dim)
    !           fsfc   = net (downw.) flux of sw rad. at the surface  (2-dim)
    !           ftop   = net (downw.) flux of sw rad. at the atm. top (2-dim)
    !           dfabs  = flux of sw rad. absorbed by each atm. layer  (3-dim)

    use mod_atparam
    use mod_physcon, only: sig, dsig
    use mod_radcon

    implicit none

    integer, intent(in) :: icltop(ix,il)
    real, intent(in) :: psa(ix,il), qa(ix,il,kx), cloudc(ix,il), clstr(ix,il)
    real, intent(inout) :: ftop(ix,il), fsfc(ix,il), fsfcd(ix,il), dfabs(ix,il,kx)

    integer :: i, j, k, nl1
    real :: acloud(ix,il), psaz(ix,il), abs1, acloud1, deltap, eps1
    real :: fband1, fband2

    nl1 = kx - 1

    fband2 = 0.05
    fband1 = 1.0 - fband2

    ! 1.  Initialization
    tau2 = 0.0

    do i = 1, ix
        do j = 1, il
            if (icltop(i,j) <= kx) then
                tau2(i,j,icltop(i,j),3) = albcl*cloudc(i,j)
            end if
            tau2(i,j,kx,3) = albcls*clstr(i,j)
        end do
    end do

    ! 2. Shortwave transmissivity:
    ! function of layer mass, ozone (in the statosphere),
    ! abs. humidity and cloud cover (in the troposphere)
    psaz = psa*zenit
    acloud = cloudc*min(abscl1*qcloud, abscl2)
    tau2(:,:,1,1) = exp(-psaz*dsig(1)*absdry)

    do k = 2, nl1
        abs1 = absdry + absaer*sig(k)**2

        if (k >= icltop(i,j)) then
            tau2(:,:,k,1) = exp(-psaz*dsig(k)*(abs1 + abswv1*qa(:,:,k) + acloud))
        else
            tau2(:,:,k,1) = exp(-psaz*dsig(k)*(abs1 + abswv1*qa(:,:,k)))
        endif
    end do

    abs1 = absdry + absaer*sig(kx)**2
    tau2(:,:,kx,1) = exp(-psaz*dsig(kx)*(abs1 + abswv1*qa(:,:,kx)))

    do k = 2, kx
        tau2(:,:,k,2) = exp(-psaz*dsig(k)*abswv2*qa(:,:,k))
    end do

    ! 3. Shortwave downward flux
    ! 3.1 Initialization of fluxes
    ftop = fsol
    flux(:,:,1) = fsol*fband1
    flux(:,:,2) = fsol*fband2

    ! 3.2 Ozone and dry-air absorption in the stratosphere
    k = 1
    dfabs(:,:,k) = flux(:,:,1)
    flux(:,:,1)  = tau2(:,:,k,1)*(flux(:,:,1) - ozupp*psa)
    dfabs(:,:,k) = dfabs(:,:,k) - flux(:,:,1)

    k = 2
    dfabs(:,:,k) = flux(:,:,1)
    flux(:,:,1)  = tau2(:,:,k,1)*(flux(:,:,1) - ozone*psa)
    dfabs(:,:,k) = dfabs(:,:,k) - flux(:,:,1)

    ! 3.3  Absorption and reflection in the troposphere
    do k = 3, kx
        tau2(:,:,k,3) = flux(:,:,1)*tau2(:,:,k,3)
        flux (:,:,1)  = flux(:,:,1) - tau2(:,:,k,3)
        dfabs(:,:,k)  = flux(:,:,1)
        flux (:,:,1)  = tau2(:,:,k,1)*flux(:,:,1)
        dfabs(:,:,k)  = dfabs(:,:,k) - flux(:,:,1)
    end do

    do k = 2, kx
        dfabs(:,:,k) = dfabs(:,:,k) + flux(:,:,2)
        flux(:,:,2)  = tau2(:,:,k,2)*flux(:,:,2)
        dfabs(:,:,k) = dfabs(:,:,k) - flux(:,:,2)
    end do

    ! 4. Shortwave upward flux
    ! 4.1  Absorption and reflection at the surface
    fsfcd       = flux(:,:,1) + flux(:,:,2)
    flux(:,:,1) = flux(:,:,1)*albsfc
    fsfc        = fsfcd - flux(:,:,1)

    ! 4.2  Absorption of upward flux
    do k=kx,1,-1
        dfabs(:,:,k) = dfabs(:,:,k) + flux(:,:,1)
        flux(:,:,1)  = tau2(:,:,k,1)*flux(:,:,1)
        dfabs(:,:,k) = dfabs(:,:,k) - flux(:,:,1)
        flux(:,:,1)  = flux(:,:,1) + tau2(:,:,k,3)
    end do

    ! 4.3  Net solar radiation = incoming - outgoing
    ftop = ftop - flux(:,:,1)

    ! 5.  Initialization of longwave radiation model
    ! 5.1  Longwave transmissivity:
    ! function of layer mass, abs. humidity and cloud cover.

    ! Cloud-free levels (stratosphere + PBL)
    k = 1
    tau2(:,:,k,1) = exp(-psa*dsig(k)*ablwin)
    tau2(:,:,k,2) = exp(-psa*dsig(k)*ablco2)
    tau2(:,:,k,3) = 1.0
    tau2(:,:,k,4) = 1.0

    do k = 2, kx, kx - 2
        tau2(:,:,k,1) = exp(-psa*dsig(k)*ablwin)
        tau2(:,:,k,2) = exp(-psa*dsig(k)*ablco2)
        tau2(:,:,k,3) = exp(-psa*dsig(k)*ablwv1*qa(:,:,k))
        tau2(:,:,k,4) = exp(-psa*dsig(k)*ablwv2*qa(:,:j,k))
    end do

    ! Cloudy layers (free troposphere)
    acloud = cloudc * ablcl2

    do k = 3, nl1
        do i = 1, ix
            do j = 1, il
                 deltap = psa(i,j)*dsig(k)

                 if (k < icltop(i,j)) then
                   acloud1 = acloud(i,j)
                 else
                   acloud1 = ablcl1*cloudc(i,j)
                 endif

                 tau2(i,j,k,1) = exp(-deltap*(ablwin+acloud1))
                 tau2(i,j,k,2) = exp(-deltap*ablco2)
                 tau2(i,j,k,3) = exp(-deltap*max(ablwv1*qa(i,j,k), acloud(i,j)))
                 tau2(i,j,k,4) = exp(-deltap*max(ablwv2*qa(i,j,k), acloud(i,j)))
            end do
        end do
    end do

    ! 5.2  Stratospheric correction terms
    eps1 = epslw/(dsig(1) + dsig(2))
    stratc(:,:,1) = stratz*psa
    stratc(:,:,2) = eps1*psa
end

subroutine radlw(imode,ta,ts,fsfcd,fsfcu,fsfc,ftop,dfabs)
    !  subroutine radlw(imode,ta,ts,
    ! &                  fsfcd,fsfcu,
    ! &                  fsfc,ftop,dfabs)
    !
    !  Purpose: Compute the absorption of longwave radiation
    !  Input:   imode  = index for operation mode
    !                    -1 : downward flux only
    !                     0 : downward + upward flux
    !                    +1 : upward flux only
    !           ta     = absolute temperature (3-dim)
    !           ts     = surface temperature                    [if imode=0]
    !           fsfcd  = downward flux of lw rad. at the sfc.   [if imode=1]
    !           fsfcu  = surface blackbody emission (upward)    [if imode=1]
    !           dfabs  = DFABS output from RADLW(-1,... )       [if imode=1]
    !  Output:  fsfcd  = downward flux of lw rad. at the sfc.[if imode=-1,0]
    !           fsfcu  = surface blackbody emission (upward)  [if imode=  0]
    !           fsfc   = net upw. flux of lw rad. at the sfc. [if imode=0,1]
    !           ftop   = outgoing flux of lw rad. at the top  [if imode=0,1]
    !           dfabs  = flux of lw rad. absorbed by each atm. layer (3-dim)
    !

    use mod_atparam
    use mod_physcon, only: sbc, dsig, wvi
    use mod_radcon, only: epslw, emisfc, fband, tau2, st4a, stratc, flux

    implicit none

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
    corlw1 = dsig(1)*stratc(:,:,2)*st4a(:,:,1,1) + stratc(:,:,1)
    corlw2 = dsig(2)*stratc(:,:,2)*st4a(:,:,2,1)
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

    implicit none

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
