!> Parametrization of short-wave radiation
module shortwave_radiation
    use params

    implicit none

    private
    public ablco2
    public get_shortwave_rad_fluxes, get_zonal_average_fields, clouds
    public increase_co2, compute_shortwave

    ! Shortwave radiation and cloud constants
    real, parameter :: solc    = 342.0 !! Solar constant (area averaged) in W/m^2
    real, parameter :: rhcl1   = 0.30  !! Relative humidity threshold corresponding to
                                       !! cloud cover = 0
    real, parameter :: rhcl2   = 1.00  !! Relative humidity correponding to cloud cover = 1
    real, parameter :: qacl    = 0.20  !! Specific humidity threshold for cloud cover
    real, parameter :: wpcl    = 0.2   !! Cloud cover weight for the square-root of precipitation
                                       !! (for p = 1 mm/day)
    real, parameter :: pmaxcl  = 10.0  !! Maximum value of precipitation (mm/day) contributing to
                                       !! cloud cover
    real, parameter :: clsmax  = 0.60  !! Maximum stratiform cloud cover
    real, parameter :: clsminl = 0.15  !! Minimum stratiform cloud cover over land (for RH = 1)
    real, parameter :: gse_s0  = 0.25  !! Gradient of dry static energy corresponding to stratiform
                                       !! cloud cover = 0
    real, parameter :: gse_s1  = 0.40  !! Gradient of dry static energy corresponding to stratiform
                                       !! cloud cover = 1
    real, parameter :: albcl   = 0.43  !! Cloud albedo (for cloud cover = 1)
    real, parameter :: albcls  = 0.50  !! Stratiform cloud albedo (for st. cloud cover = 1)
    real, parameter :: epssw   = 0.020 !! Fraction of incoming solar radiation absorbed by ozone

    ! Shortwave absorptivities (for dp = 10^5 Pa)
    real, parameter :: absdry =  0.033 !! Absorptivity of dry air (visible band)
    real, parameter :: absaer =  0.033 !! Absorptivity of aerosols (visible band)
    real, parameter :: abswv1 =  0.022 !! Absorptivity of water vapour
                                       !! (visible band, for dq = 1 g/kg)
    real, parameter :: abswv2 = 15.000 !! Absorptivity of water vapour
                                       !! (near IR band, for dq = 1 g/kg)
    real, parameter :: abscl1 =  0.015 !! Absorptivity of clouds (visible band, maximum value)
    real, parameter :: abscl2 =  0.15  !! Absorptivity of clouds
                                       !! (visible band, for dq_base = 1 g/kg)

    ! Longwave absorptivities (per dp = 10^5 Pa)
    real, parameter :: ablwin =  0.3 !! Absorptivity of air in "window" band
    real            :: ablco2 =  6.0 !! Absorptivity of air in CO2 band
    real, parameter :: ablwv1 =  0.7 !! Absorptivity of water vapour in H2O band 1 (weak),
                                    !! (for dq = 1 g/kg)
    real, parameter :: ablwv2 = 50.0 !! Absorptivity of water vapour in H2O band 2 (strong),
                                    !! (for dq = 1 g/kg)
    real, parameter :: ablcl1 = 12.0 !! Absorptivity of "thick" clouds in window band
                                    !! (below cloud top)
    real, parameter :: ablcl2 =  0.6 !! Absorptivity of "thin" upper clouds in window and H2O bands

    ! Zonally-averaged fields for SW/LW scheme (updated in sol_oz)
    real, dimension(ix,il) :: fsol   !! Flux of incoming solar radiation
    real, dimension(ix,il) :: ozone  !! Flux absorbed by ozone (lower stratosphere)
    real, dimension(ix,il) :: ozupp  !! Flux absorbed by ozone (upper stratosphere)
    real, dimension(ix,il) :: zenit  !! Optical depth ratio (function of solar zenith angle)
    real, dimension(ix,il) :: stratz !! Stratospheric correction for polar night

    real, dimension(ix,il) :: qcloud !! Equivalent specific humidity of clouds

    ! Logical flags to control shortwave radiation behaviour
    logical, parameter :: increase_co2 = .false. !! Flag for CO2 optical thickness increase
    logical :: compute_shortwave       = .true.  !! Flag for shortwave radiation routine (turned on
                                                 !! and off in main loop depending on the value of
                                                 !! nstrad)

contains
    !> Compute the absorption of shortwave radiation and initialize arrays
    !  for longwave-radiation routines
    subroutine get_shortwave_rad_fluxes(psa, qa, icltop, cloudc, clstr, fsfcd, fsfc, ftop, dfabs)
        use geometry, only: fsg, dhs
        use mod_radcon

        real, intent(in) :: psa(ix,il)        !! Normalised surface pressure [p/p0]
        real, intent(in) :: qa(ix,il,kx)      !! Specific humidity [g/kg]
        integer, intent(in) :: icltop(ix,il)  !! Cloud top level
        real, intent(in) :: cloudc(ix,il)     !! Total cloud cover
        real, intent(in) :: clstr(ix,il)      !! Stratiform cloud cover
        real, intent(out) :: fsfcd(ix,il)     !! Total downward flux of short-wave radiation at the
                                              !! surface
        real, intent(out) :: fsfc(ix,il)      !! Net downward flux of short-wave radiation at the
                                              !! surface
        real, intent(out) :: ftop(ix,il)      !! Net downward flux of short-wave radiation at the
                                              !! top of the atmosphere
        real, intent(out) :: dfabs(ix,il,kx)  !! Flux of short-wave radiation absorbed in each
                                              !! atmospheric layer

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
        tau2(:,:,1,1) = exp(-psaz*dhs(1)*absdry)

        do k = 2, nl1
            abs1 = absdry + absaer*fsg(k)**2

            do i = 1, ix
                do j = 1, il
                    if (k >= icltop(i,j)) then
                        tau2(i,j,k,1) = exp(-psaz(i,j)*dhs(k)*(abs1 + abswv1*qa(i,j,k) + acloud(i,j)))
                    else
                        tau2(i,j,k,1) = exp(-psaz(i,j)*dhs(k)*(abs1 + abswv1*qa(i,j,k)))
                    end if
                end do
            end do
        end do

        abs1 = absdry + absaer*fsg(kx)**2
        tau2(:,:,kx,1) = exp(-psaz*dhs(kx)*(abs1 + abswv1*qa(:,:,kx)))

        do k = 2, kx
            tau2(:,:,k,2) = exp(-psaz*dhs(k)*abswv2*qa(:,:,k))
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
        tau2(:,:,k,1) = exp(-psa*dhs(k)*ablwin)
        tau2(:,:,k,2) = exp(-psa*dhs(k)*ablco2)
        tau2(:,:,k,3) = 1.0
        tau2(:,:,k,4) = 1.0

        do k = 2, kx, kx - 2
            tau2(:,:,k,1) = exp(-psa*dhs(k)*ablwin)
            tau2(:,:,k,2) = exp(-psa*dhs(k)*ablco2)
            tau2(:,:,k,3) = exp(-psa*dhs(k)*ablwv1*qa(:,:,k))
            tau2(:,:,k,4) = exp(-psa*dhs(k)*ablwv2*qa(:,:,k))
        end do

        ! Cloudy layers (free troposphere)
        acloud = cloudc * ablcl2

        do k = 3, nl1
            do i = 1, ix
                do j = 1, il
                     deltap = psa(i,j)*dhs(k)

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
        eps1 = epslw/(dhs(1) + dhs(2))
        stratc(:,:,1) = stratz*psa
        stratc(:,:,2) = eps1*psa
    end

    !> Compute zonally-averaged fields to be used in the computation of
    !  short-wave absorption
    subroutine get_zonal_average_fields(tyear)
        use geometry, only: sia, coa

        real, intent(in) :: tyear !! time as fraction of year (0-1, 0 = 1jan.h00)

        real :: topsr(il), alpha, azen, coz1, coz2, dalpha, flat2, fs0
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
        call solar(tyear, 4.0*solc, topsr)

        do j = 1, il
            flat2 = 1.5*sia(j)**2 - 0.5

            ! Solar radiation at the top
            fsol(:,j) = topsr(j)

            ! Ozone depth in upper and lower stratosphere
            ozupp(:,j) = 0.5*epssw
            ozone(:,j) = 0.4*epssw*(1.0 + coz1*sia(j) + coz2*flat2)

            ! Zenith angle correction to (downward) absorptivity
            zenit(:,j) = 1.0 + azen*(1.0 - (coa(j)*cos(rzen) + sia(j)*sin(rzen)))**nzen

            ! Ozone absorption in upper and lower stratosphere
            ozupp(:,j) = fsol(:,j)*ozupp(:,j)*zenit(:,j)
            ozone(:,j) = fsol(:,j)*ozone(:,j)*zenit(:,j)

            ! Polar night cooling in the stratosphere
            stratz(:,j) = max(fs0 - fsol(:,j), 0.0)
        end do
    end

    ! Average daily flux of solar radiation, from Hartmann (1994)
    subroutine solar(tyear, csol, topsr)
        use geometry, only: coa, sia

        real, intent(in) :: tyear      !! time as fraction of year (0-1, 0 = 1jan.h00)
        real, intent(in) :: csol       !! The solar constant [W/m^2]
        real, intent(out) :: topsr(il) !! Daily-average insolation at the top of the atmosphere as a
                                       !! function of latitude

        integer :: j
        real :: ca1, ca2, ca3, cdecl, ch0, csolp, decl, fdis, h0, alpha, pigr, sa1
        real :: sa2, sa3, sdecl, sh0, tdecl

        ! 1. Compute declination angle and Earth-Sun distance factor
        pigr  = 2.0*asin(1.0)
        alpha = 2.0*pigr*tyear

        ca1 = cos(alpha)
        sa1 = sin(alpha)
        ca2 = ca1*ca1-sa1*sa1
        sa2 = 2.*sa1*ca1
        ca3 = ca1*ca2-sa1*sa2
        sa3 = sa1*ca2+sa2*ca1

        decl = 0.006918 - 0.399912*ca1 + 0.070257*sa1 - 0.006758*ca2 + 0.000907*sa2&
            & - 0.002697*ca3 + 0.001480*sa3

        fdis = 1.000110 + 0.034221*ca1 + 0.001280*sa1 + 0.000719*ca2 + 0.000077*sa2

        cdecl = cos(decl)
        sdecl = sin(decl)
        tdecl = sdecl/cdecl

        ! 2. Compute daily-average insolation at the atm. top
        csolp=csol/pigr

        do j = 1, il
            ch0 = min(1.0, max(-1.0, -tdecl*sia(j)/coa(j)))
            h0  = acos(ch0)
            sh0 = sin(h0)

            topsr(j) = csolp*fdis*(h0*sia(j)*sdecl + sh0*coa(j)*cdecl)
        end do
    end

    !>  Compute cloud-top level and cloud cover
    subroutine clouds(qa,rh,precnv,precls,iptop,gse,fmask,icltop,cloudc,clstr)
        integer :: iptop(ix,il)
        real, intent(in) :: qa(ix,il,kx)      !! Specific humidity [g/kg]
        real, intent(in) :: rh(ix,il,kx)      !! Relative humidity
        real, intent(in) :: precnv(ix,il)     !! Convection precipitation
        real, intent(in) :: precls(ix,il)     !! Large-scale condensational precipitation
        real, intent(in) :: gse(ix,il)        !! Vertical gradient of dry static energy
        real, intent(in) :: fmask(ix,il)      !! Fraction land-sea mask
        integer, intent(out) :: icltop(ix,il) !! Cloud top level
        real, intent(out) :: cloudc(ix,il)    !! Total cloud cover
        real, intent(out) :: clstr(ix,il)     !! Stratiform cloud cover

        integer :: i, j, k, nl1, nlp
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
    end
end module
