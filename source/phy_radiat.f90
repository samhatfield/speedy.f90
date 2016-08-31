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
    !  Updated common blocks: radzon

    use mod_atparam
    use mod_physcon, only: slat, clat
    use mod_radcon, only: solc, epssw, fsol, ozone, ozupp, zenit, stratz

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real, intent(in) :: tyear
    real :: topsr(nlat), alpha, azen, coz1, coz2, czen, dalpha, flat2, fs0
    real :: nzen, rzen, szen
    integer :: i, j, j0

    ! alpha = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
    alpha=4.*asin(1.)*(tyear+10./365.)
    dalpha=0.
    !DALPHA=ASIN(0.5)

    coz1= 1.0*max(0.,cos(alpha-dalpha))
    coz2= 1.8

    azen=1.0
    nzen=2

    rzen=-cos(alpha)*23.45*asin(1.)/90.
    czen=cos(rzen)
    szen=sin(rzen)

    fs0=6.

    ! Solar radiation at the top
    call solar(tyear,4.*solc,nlat,clat,slat,topsr)

    do j=1,nlat
        j0=1+nlon*(j-1)
        flat2=1.5*slat(j)**2-0.5

        ! Solar radiation at the top
        fsol(j0)=topsr(j)

        ! Ozone depth in upper and lower stratosphere 
        ozupp(j0)=0.5*epssw
        ozone(j0)=0.4*epssw*(1.0+coz1*slat(j)+coz2*flat2)

        ! Zenith angle correction to (downward) absorptivity 
        zenit(j0)=1.+azen*(1.-(clat(j)*czen+slat(j)*szen))**nzen

        ! Ozone absorption in upper and lower stratosphere 
        ozupp(j0)=fsol(j0)*ozupp(j0)*zenit(j0)
        ozone(j0)=fsol(j0)*ozone(j0)*zenit(j0)

        ! Polar night cooling in the stratosphere
        stratz(j0)=max(fs0-fsol(j0),0.)

        do i=1,nlon-1
            fsol  (i+j0) = fsol  (j0)
            ozone (i+j0) = ozone (j0)
            ozupp (i+j0) = ozupp (j0)
            zenit (i+j0) = zenit (j0)
            stratz(i+j0) = stratz(j0)
        end do
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

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer :: iptop(ngp)
    real, intent(in) :: qa(ngp,nlev), rh(ngp,nlev), precnv(ngp), precls(ngp), gse(ngp),&
        & fmask(ngp)
    real, intent(inout) :: cloudc(ngp), clstr(ngp)
    integer, intent(inout) :: icltop(ngp)

    integer :: inew, j, k, nl1, nlp
    real :: albcor, cl1, clfact, clstrl, drh, fstab, pr1, rgse, rrcl
      
    nl1  = nlev-1
    nlp  = nlev+1
    rrcl = 1./(rhcl2-rhcl1)

    ! 1.  Cloud cover, defined as the sum of:
    !     - a term proportional to the square-root of precip. rate 
    !     - a quadratic function of the max. relative humidity
    !       in tropospheric layers above PBL where Q > QACL :
    !       ( = 0 for RHmax < RHCL1, = 1 for RHmax > RHCL2 )
    !     Cloud-top level: defined as the highest (i.e. least sigma)
    !       between the top of convection/condensation and
    !       the level of maximum relative humidity. 

    do j=1,ngp
        if (rh(j,nl1).gt.rhcl1) then
            cloudc(j) = rh(j,nl1)-rhcl1
            icltop(j) = nl1
        else
            cloudc(j) = 0.
            icltop(j) = nlp
        end if
    end do

    do k=3,nlev-2
        do j=1,ngp
            drh = rh(j,k)-rhcl1
            if (drh.gt.cloudc(j).and.qa(j,k).gt.qacl) then
                cloudc(j) = drh
                icltop(j) = k
            end if
        end do
    end do

    do j=1,ngp
        cl1 = min(1.,cloudc(j)*rrcl)
        pr1 = min(pmaxcl,86.4*(precnv(j)+precls(j)))
        cloudc(j) = min(1.,wpcl*sqrt(pr1)+cl1*cl1)
        icltop(j) = min(iptop(j),icltop(j))
    end do

    ! 2.  Equivalent specific humidity of clouds 
    qcloud = qa(:,nl1)

    ! 3. Stratiform clouds at the top of PBL
    inew = 1

    if (inew.gt.0) then
        !        CLSMAX  = 0.6
        !        CLSMINL = 0.15
        !        GSE_S0  = 0.25
        !        GSE_S1  = 0.40

        clfact = 1.2
        rgse   = 1./(gse_s1-gse_s0)

        do j=1,ngp
            ! Stratocumulus clouds over sea
            fstab    = max(0.,min(1.,rgse*(gse(j)-gse_s0)))
            clstr(j) = fstab*max(clsmax-clfact*cloudc(j),0.)
            ! Stratocumulus clouds over land
            clstrl   = max(clstr(j),clsminl)*rh(j,nlev)
            clstr(j) = clstr(j)+fmask(j)*(clstrl-clstr(j))
        end do
    else
        clsmax  = 0.3
        clsminl = 0.1
        albcor  = albcl/0.5
 
        do j=1,ngp
            ! stratocumulus clouds over sea
            clstr(j) = max(clsmax-cloudc(j),0.)
            ! rescale for consistency with previous albedo values
            clstr(j) = clstr(j)*albcor
            ! correction for aerosols over land
            clstr(j) = clstr(j)+fmask(j)*(clsminl-clstr(j))
        end do
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

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer, intent(in) :: icltop(ngp)
    real, intent(in) :: psa(ngp), qa(ngp,nlev), cloudc(ngp), clstr(ngp)
    real, intent(inout) :: ftop(ngp), fsfc(ngp), fsfcd(ngp), dfabs(ngp,nlev)

    integer :: j, k, nl1
    real :: acloud(ngp), psaz(ngp), abs1, acloud1, deltap, eps1
    real :: fband1, fband2

    nl1 = nlev-1

    fband2 = 0.05
    fband1 = 1.-fband2

    ! ALBMINL=0.05
    ! ALBCLS = 0.5
    
    ! 1.  Initialization
    tau2 = 0.0

    do j=1,ngp
        !fk-- change to ensure only icltop <= nlev used
        if(icltop(j) .le. nlev) then
          tau2(j,icltop(j),3)= albcl*cloudc(j)
        endif
        !fk-- end change
        tau2(j,nlev,3)     = albcls*clstr(j)
    end do

    ! 2. Shortwave transmissivity:
    ! function of layer mass, ozone (in the statosphere),
    ! abs. humidity and cloud cover (in the troposphere)

    do j=1,ngp
        psaz(j)=psa(j)*zenit(j)
        acloud(j)=cloudc(j)*min(abscl1*qcloud(j),abscl2)
    end do

    do j=1,ngp
        deltap=psaz(j)*dsig(1)
        tau2(j,1,1)=exp(-deltap*absdry)
    end do

    do k=2,nl1
        abs1=absdry+absaer*sig(k)**2
        do j=1,ngp
            deltap=psaz(j)*dsig(k)
            if (k.ge.icltop(j)) then
                tau2(j,k,1)=exp(-deltap*(abs1+abswv1*qa(j,k)+acloud(j)))
            else
              tau2(j,k,1)=exp(-deltap*(abs1+abswv1*qa(j,k)))
            endif
        end do
    end do

    abs1=absdry+absaer*sig(nlev)**2
    do j=1,ngp
        deltap=psaz(j)*dsig(nlev)
        tau2(j,nlev,1)=exp(-deltap*(abs1+abswv1*qa(j,nlev)))
    end do

    do k=2,nlev
        do j=1,ngp
          deltap=psaz(j)*dsig(k)
          tau2(j,k,2)=exp(-deltap*abswv2*qa(j,k))
        end do
    end do

    ! 3. Shortwave downward flux 
    ! 3.1 Initialization of fluxes 
    ftop = fsol
    flux(:,1) = fsol * fband1
    flux(:,2) = fsol * fband2

    ! 3.2 Ozone and dry-air absorption in the stratosphere
    K=1
    do j=1,ngp
        dfabs(j,k)=flux(j,1)
        flux (j,1)=tau2(j,k,1)*(flux(j,1)-ozupp(j)*psa(j))
        dfabs(j,k)=dfabs(j,k)-flux(j,1)
    end do

    k=2
    do j=1,ngp
        dfabs(j,k)=flux(j,1)
        flux (j,1)=tau2(j,k,1)*(flux(j,1)-ozone(j)*psa(j))
        dfabs(j,k)=dfabs(j,k)-flux(j,1)
    end do

    ! 3.3  Absorption and reflection in the troposphere
    do k=3,nlev
        do j=1,ngp
            tau2(j,k,3)=flux(j,1)*tau2(j,k,3)
            flux (j,1)=flux(j,1)-tau2(j,k,3)
            dfabs(j,k)=flux(j,1)
            flux (j,1)=tau2(j,k,1)*flux(j,1)
            dfabs(j,k)=dfabs(j,k)-flux(j,1)
        end do
    end do

    do k=2,nlev
        do j=1,ngp
          dfabs(j,k)=dfabs(j,k)+flux(j,2)
          flux (j,2)=tau2(j,k,2)*flux(j,2)
          dfabs(j,k)=dfabs(j,k)-flux(j,2)
        end do
    end do

    ! 4. Shortwave upward flux 
    ! 4.1  Absorption and reflection at the surface
    do j=1,ngp
        fsfcd(j)  = flux(j,1)+flux(j,2)
        flux(j,1) = flux(j,1)*albsfc(j)
        fsfc(j)   = fsfcd(j)-flux(j,1)
    end do

    ! 4.2  Absorption of upward flux
    do k=nlev,1,-1
        do j=1,ngp
            dfabs(j,k)=dfabs(j,k)+flux(j,1)
            flux (j,1)=tau2(j,k,1)*flux(j,1)
            dfabs(j,k)=dfabs(j,k)-flux(j,1)
            flux (j,1)=flux(j,1)+tau2(j,k,3)
        end do
    end do

    ! 4.3  Net solar radiation = incoming - outgoing
    ftop = ftop - flux(:,1)

    ! 5.  Initialization of longwave radiation model
    ! 5.1  Longwave transmissivity:
    ! function of layer mass, abs. humidity and cloud cover.

    ! Cloud-free levels (stratosphere + PBL)
    k=1
    do j=1,ngp
        deltap=psa(j)*dsig(k)
        tau2(j,k,1)=exp(-deltap*ablwin)
        tau2(j,k,2)=exp(-deltap*ablco2)
        tau2(j,k,3)=1.
        tau2(j,k,4)=1.
    end do

    do k=2,nlev,nlev-2
        do j=1,ngp
            deltap=psa(j)*dsig(k)
            tau2(j,k,1)=exp(-deltap*ablwin)
            tau2(j,k,2)=exp(-deltap*ablco2)
            tau2(j,k,3)=exp(-deltap*ablwv1*qa(j,k))
            tau2(j,k,4)=exp(-deltap*ablwv2*qa(j,k))
        end do
    end do

    ! Cloudy layers (free troposphere)
    acloud = cloudc * ablcl2

    do k=3,nl1
       do j=1,ngp
         deltap=psa(j)*dsig(k)
         if (k.lt.icltop(j)) then
           acloud1=acloud(j)
         else
           acloud1=ablcl1*cloudc(j)
         endif
         tau2(j,k,1)=exp(-deltap*(ablwin+acloud1))
         tau2(j,k,2)=exp(-deltap*ablco2)
         tau2(j,k,3)=exp(-deltap*max(ablwv1*qa(j,k),acloud(j)))
         tau2(j,k,4)=exp(-deltap*max(ablwv2*qa(j,k),acloud(j)))
       end do
    end do

    ! 5.2  Stratospheric correction terms
    eps1=epslw/(dsig(1)+dsig(2))
    do j=1,ngp
        stratc(j,1)=stratz(j)*psa(j)
        stratc(j,2)=eps1*psa(j)
    end do
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
    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    ! Number of radiation bands with tau < 1
    integer, parameter :: nband=4

    real, intent(in) :: ta(ngp,nlev), ts(ngp)
    real, intent(inout) :: fsfcd(ngp), fsfcu(ngp), ftop(ngp), fsfc(ngp)
    real, intent(inout) :: dfabs(ngp,nlev)

    integer :: j, jb, k, nl1
    real :: anis, anish, brad, corlw, corlw1, corlw2, emis, eps1, esbc, refsfc
    real :: st3a, tsq

    nl1=nlev-1

    refsfc=1.-emisfc

    if (imode.eq.1) go to 410
    ! 1. Blackbody emission from atmospheric levels.
    ! The linearized gradient of the blakbody emission is computed
    ! from temperatures at layer boundaries, which are interpolated 
    ! assuming a linear dependence of T on log_sigma.
    ! Above the first (top) level, the atmosphere is assumed isothermal.

    ! Temperature at level boundaries
    do k=1,nl1
        do j=1,ngp
            st4a(j,k,1)=ta(j,k)+wvi(k,2)*(ta(j,k+1)-ta(j,k))
        end do
    end do 

    ! Mean temperature in stratospheric layers 
    do j=1,ngp
        st4a(j,1,2)=0.75*ta(j,1)+0.25* st4a(j,1,1)
        st4a(j,2,2)=0.50*ta(j,2)+0.25*(st4a(j,1,1)+st4a(j,2,1))
    end do

    ! Temperature gradient in tropospheric layers 
    anis =1.0
    anish=0.5*anis

    do k=3,nl1
        do j=1,ngp
            st4a(j,k,2)=anish*max(st4a(j,k,1)-st4a(j,k-1,1),0.)
        end do
    end do

    do j=1,ngp
        st4a(j,nlev,2)=anis*max(ta(j,nlev)-st4a(j,nl1,1),0.)
    end do

    ! Blackbody emission in the stratosphere
    do k=1,2
        do j=1,ngp
            st4a(j,k,1)=sbc*st4a(j,k,2)**4
            st4a(j,k,2)=0.
        end do
    end do

    ! Blackbody emission in the troposphere
    do k=3,nlev
        do j=1,ngp
            st3a=sbc*ta(j,k)**3
            st4a(j,k,1)=st3a*ta(j,k)
            st4a(j,k,2)=4.*st3a*st4a(j,k,2)
        end do
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
    do jb=1,2
        do j=1,ngp
            emis=1.-tau2(j,k,jb)
            brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1)+emis*st4a(j,k,2))
            flux(j,jb)=emis*brad
            dfabs(j,k)=dfabs(j,k)-flux(j,jb)
        end do
    end do

    flux(:,3:nband) = 0.0
	
    ! 3.2  Troposphere
    do jb=1,nband
        do k=2,nlev
            do j=1,ngp
                emis=1.-tau2(j,k,jb)
                brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1)+emis*st4a(j,k,2))
                dfabs(j,k)=dfabs(j,k)+flux(j,jb)
                flux(j,jb)=tau2(j,k,jb)*flux(j,jb)+emis*brad
                dfabs(j,k)=dfabs(j,k)-flux(j,jb)
            end do
        end do
    end do

    ! 3.3 Surface downward flux
    do jb=1,nband
        do j=1,ngp
            fsfcd(j)=fsfcd(j)+emisfc*flux(j,jb)
        end do
    end do

    ! 3.4 Correction for "black" band (incl. surface reflection)
    eps1=epslw*emisfc
    do j=1,ngp
        corlw=eps1*st4a(j,nlev,1)
        dfabs(j,nlev)=dfabs(j,nlev)-corlw
        fsfcd(j)     =fsfcd(j)     +corlw
    end do

    if (imode.eq.-1) return

    ! 4. Emission ad absorption of longwave upward flux. 
    !    For upward emission, a correction term depending on the      
    !    local temperature gradient and on the layer transmissivity is  
    !    subtracted from the average (full-level) emission of each layer. 
	
    ! 4.1  Surface

    ! Black-body (or grey-body) emission 
    esbc=emisfc*sbc
    do j=1,ngp
        tsq=ts(j)*ts(j)
        fsfcu(j)=esbc*tsq*tsq
    end do

    ! Entry point for upward-only mode (IMODE=1)
 410  continue

    fsfc = fsfcu - fsfcd

    do jb=1,nband
        do j=1,ngp
            flux(j,jb)=fband(nint(ts(j)),jb)*fsfcu(j)+refsfc*flux(j,jb)
        end do
    end do
	
    ! 4.2  Troposphere

    ! Correction for "black" band
    do j=1,ngp
        dfabs(j,nlev)=dfabs(j,nlev)+epslw*fsfcu(j)
    end do

    do jb=1,nband
        do k=nlev,2,-1
            do j=1,ngp
                emis=1.-tau2(j,k,jb)
                brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1)-emis*st4a(j,k,2))
                dfabs(j,k)=dfabs(j,k)+flux(j,jb)
                flux(j,jb)=tau2(j,k,jb)*flux(j,jb)+emis*brad
                dfabs(j,k)=dfabs(j,k)-flux(j,jb)
            end do
        end do
    end do
	
    ! 4.3  Stratosphere
    k=1
    do jb=1,2
        do j=1,ngp
            emis=1.-tau2(j,k,jb)
            brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1)-emis*st4a(j,k,2))
            dfabs(j,k)=dfabs(j,k)+flux(j,jb)
            flux(j,jb)=tau2(j,k,jb)*flux(j,jb)+emis*brad
            dfabs(j,k)=dfabs(j,k)-flux(j,jb)
        end do
    end do

    ! Correction for "black" band and polar night cooling
    do j=1,ngp
        corlw1=dsig(1)*stratc(j,2)*st4a(j,1,1)+stratc(j,1)
        corlw2=dsig(2)*stratc(j,2)*st4a(j,2,1)
        dfabs(j,1)=dfabs(j,1)-corlw1
        dfabs(j,2)=dfabs(j,2)-corlw2
        ftop(j)   =corlw1+corlw2
    end do

    ! 4.4  Outgoing longwave radiation 
    do jb=1,nband
        do j=1,ngp
            ftop(j)=ftop(j)+flux(j,jb)
        end do
    end do
end

subroutine radset
    ! subroutine radset
    !
    ! Purpose: compute energy fractions in LW bands
    !          as a function of temperature

    use mod_atparam
    use mod_radcon, only: epslw, fband

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer :: jb, jtemp
    real :: eps1

    eps1=1.-epslw

    do jtemp=200,320
        fband(jtemp,2)=(0.148-3.0e-6*(jtemp-247)**2)*eps1
        fband(jtemp,3)=(0.356-5.2e-6*(jtemp-282)**2)*eps1
        fband(jtemp,4)=(0.314+1.0e-5*(jtemp-315)**2)*eps1
        fband(jtemp,1)=eps1-(fband(jtemp,2)+fband(jtemp,3)+fband(jtemp,4))
    end do

    do jb=1,4
        do jtemp=100,199
            fband(jtemp,jb)=fband(200,jb)
        end do
        do jtemp=321,400
            fband(jtemp,jb)=fband(320,jb)
        end do
    end do
end
