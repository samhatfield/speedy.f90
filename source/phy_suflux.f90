! Compute surface fluxes of momentum, energy and moisture, and define surface
! skin temperature from energy balance
! Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
!          UA     = u-wind                          (3-dim)
!          VA     = v-wind                          (3-dim)
!          TA     = temperature                     (3-dim)
!          QA     = specific humidity [g/kg]        (3-dim)
!          RH     = relative humidity [0-1]         (3-dim)
!          PHI    = geopotential                    (3-dim)
!          PHI0   = surface geopotential            (2-dim)
!          FMASK  = fractional land-sea mask        (2-dim)
!          TSEA   =  sea-surface temperature        (2-dim)
!          SSRD   = sfc sw radiation (downw. flux)  (2-dim)
!          SLRD   = sfc lw radiation (downw. flux)  (2-dim)
!          LFLUXLAND   = Logical related ti flux-correction
! Output:  USTR   = u stress                        (2-dim)
!          VSTR   = v stress                        (2-dim)
!          SHF    = sensible heat flux              (2-dim)
!          EVAP   = evaporation [g/(m^2 s)]         (2-dim)
!          SLRU   = sfc lw radiation (upward flux)  (2-dim)
!          HFLUXN = net heat flux into land/sea     (2-dim)
!          TSFC   = surface temperature (clim.)     (2-dim)
!          TSKIN  = skin surface temperature        (2-dim)
!          U0     = near-surface u-wind             (2-dim)
!          V0     = near-surface v-wind             (2-dim)
!          T0     = near-surface air temperature    (2-dim)
!          Q0     = near-surface sp. humidity [g/kg](2-dim)
subroutine suflux (psa, ua, va, ta, qa, rh, phi, phi0, fmask, tsea, ssrd, slrd, &
        & ustr, vstr, shf, evap, slru, hfluxn, tsfc, tskin, u0, v0, t0, q0, lfluxland)
    use mod_atparam
    use mod_sflcon
    use mod_physcon, only: p0, rd, cp, alhc, sbc, sigl, wvi, clat
    use mod_radcon, only: emisfc, alb_l, alb_s, snowc
	use mod_cpl_land_model, only: stl_am, soilw_am

    implicit none

    real, dimension(ix,il,kx), intent(in) :: ua, va, ta, qa, rh, phi
    real, dimension(ix,il), intent(in) :: phi0, fmask, tsea, ssrd,&
        & slrd

    real, dimension(ix,il,3), intent(inout) :: ustr, vstr, shf, evap, slru
    real, intent(inout) :: hfluxn(ix,il,2)
    real, dimension(ix,il), intent(inout) :: tsfc, tskin, u0, v0, t0, q0

    integer :: i, j, ks, nl1
    real, dimension(ix,il,2), save :: t1, q1
    real, dimension(ix,il,2) :: t2, qsat0
    real, save :: denvvs(ix,il,0:2)
    real :: dslr(ix,il), dtskin(ix,il), clamb(ix,il), astab, cdldv, cdsdv(ix,il), chlcp
    real :: dlambda, dt1, dthl, dths, esbc, ghum0, gtemp0
    real :: prd, qdummy, rcp, rdphi0, rdth, rdummy, tsk3(ix,il), vg2

    logical lscasym, lskineb
    logical lfluxland

    real :: psa(ix,il)

    lscasym = .true.   ! true : use an asymmetric stability coefficient
    lskineb = .true.   ! true : redefine skin temp. from energy balance

    esbc  = emisfc*sbc

    ghum0 = 1.0 - fhum0

    dlambda = clambsn - clambda

    ! =========================================================================
    ! Land surface
    ! =========================================================================

    if (lfluxland)  then
        ! 1. Extrapolation of wind, temp, hum. and density to the surface

        ! 1.1 Wind components
        u0 = fwind0 * ua(:,:,kx)
        v0 = fwind0 * va(:,:,kx)

        ! 1.2 Temperature
        gtemp0 = 1.0 - ftemp0
        rcp = 1.0/cp
        rdphi0 = -1.0/(rd*288.0*sigl(kx))
        nl1 = kx-1

        do i = 1, ix
            do j = 1, il
                ! Temperature difference between lowest level and sfc
                dt1 = wvi(kx,2)*(ta(i,j,kx) - ta(i,j,nl1))

                ! Extrapolated temperature using actual lapse rate (1:land, 2:sea)
                t1(i,j,1) = ta(i,j,kx) + dt1
                t1(i,j,2) = t1(i,j,1) + phi0(i,j)*dt1*rdphi0

                ! Extrapolated temperature using dry-adiab. lapse rate (1:land, 2:sea)
                t2(i,j,2) = ta(i,j,kx) + rcp*phi(i,j,kx)
                t2(i,j,1) = t2(i,j,2) - rcp*phi0(i,j)
            end do
        end do

        do i = 1, ix
            do j = 1, il
                if (ta(i,j,kx) > ta(i,j,nl1)) then
                    ! Use extrapolated temp. if dT/dz < 0
                    t1(i,j,1) = ftemp0*t1(i,j,1) + gtemp0*t2(i,j,1)
                    t1(i,j,2) = ftemp0*t1(i,j,2) + gtemp0*t2(i,j,2)
                else
                    ! Use temp. at lowest level if dT/dz > 0
                    t1(i,j,1) = ta(i,j,kx)
                    t1(i,j,2) = ta(i,j,kx)
                endif
                t0(i,j) = t1(i,j,2) + fmask(i,j)*(t1(i,j,1) - t1(i,j,2))
            end do
        end do

        ! 1.3 Density * wind speed (including gustiness factor)
        prd = p0/rd
        vg2 = vgust**2.0

        denvvs(:,:,0) = (prd*psa/t0)*sqrt(u0**2.0 + v0**2.0 + vg2)

        ! 2. Compute land-sfc. fluxes using prescribed skin temperature

        ! 2.1 Define effective skin temperature to compensate for
        !     non-linearity of heat/moisture fluxes during the daily cycle
        do j = 1, il
            tskin(:,j) = stl_am(:,j) + ctday*sqrt(clat(j))*ssrd(:,j)*(1.0 - alb_l(:,j))*psa(:,j)
        end do

        ! 2.2 Stability correction = f[pot.temp.(sfc)-pot.temp.(air)]
        rdth  = fstab/dtheta
        astab = 1.0
        if (lscasym) astab = 0.5   ! to get smaller ds/dt in stable conditions

        do i = 1, ix
            do j = 1,il
                ! Potential temp. difference (land+sea average)
                if (tskin(i,j) > t2(i,j,1)) then
                    dthl = min(dtheta, tskin(i,j) - t2(i,j,1))
                else
                    dthl = max(-dtheta, astab*(tskin(i,j) - t2(i,j,1)))
                end if
                denvvs(i,j,1) = denvvs(i,j,0)*(1.0 + dthl*rdth)
            end do
        end do

        ! 2.3 Wind stress
        do i = 1, ix
            do j = 1, il
                cdldv       =  cdl*denvvs(i,j,0)*forog(i,j)
                ustr(i,j,1) = -cdldv*ua(i,j,kx)
                vstr(i,j,1) = -cdldv*va(i,j,kx)
            end do
        end do

        ! 2.4 Sensible heat flux
        chlcp = chl*cp
        shf(:,:,1) = chlcp*denvvs(:,:,1)*(tskin - t1(:,:,1))

        ! 2.5 Evaporation
        if (fhum0 > 0.0) then
            call shtorh(-1, t1, psa, 1.0, q1, rh(:,:,kx), qsat0(:,:,1))

            q1(:,:,1) = fhum0*q1(:,:,1) + ghum0*qa(:,:,kx)
        else
            q1(:,:,1) = qa(:,:,kx)
        end if

        call shtorh(0, tskin, psa, 1.0, qdummy, rdummy, qsat0)
        evap(:,:,1) = chl*denvvs(:,:,1)*max(0.0, soilw_am*qsat0(:,:,1) - q1(:,:,1))

        ! 3. Compute land-surface energy balance;
        !    adjust skin temperature and heat fluxes

        ! 3.1. Emission of lw radiation from the surface
        !      and net heat fluxes into land surface
        tsk3 = tskin**3.0
        dslr = 4.0*esbc*tsk3
        slru(:,:,1) = esbc*tsk3*tskin
        hfluxn(:,:,1) = ssrd*(1.0 - alb_l) + slrd - (slru(:,:,1) + shf(:,:,1) + alhc*evap(:,:,1))

        ! 3.2 Re-definition of skin temperature from energy balance
        if (lskineb) then
            ! Compute net heat flux including flux into ground
            clamb = clambda + snowc*dlambda
            hfluxn(:,:,1) = hfluxn(:,:,1) - clamb*(tskin - stl_am)
            dtskin = tskin + 1.0

            ! Compute d(Evap) for a 1-degree increment of Tskin
            call shtorh(0, dtskin, psa, 1.0, qdummy, rdummy, qsat0(:,:,2))

            do i = 1, ix
                do j = 1, il
                    if (evap(i,j,1) > 0.0) then
                        qsat0(i,j,2) = soilw_am(i,j)*(qsat0(i,j,2) - qsat0(i,j,1))
                    else
                        qsat0(i,j,2) = 0.0
                    endif
                end do
            end do

            ! Redefine skin temperature to balance the heat budget
            dtskin = hfluxn(:,:,1)/(clamb + dslr + chl*denvvs(:,:,1)*(cp + alhc*qsat0(:,:,2)))
            tskin = tskin + dtskin

            ! Add linear corrections to heat fluxes
            shf(:,:,1) = shf(:,:,1) + chlcp*denvvs(:,:,1)*dtskin
            evap(:,:,1) = evap(:,:,1) + chl*denvvs(:,:,1)*qsat0(:,:,2)*dtskin
            slru(:,:,1) = slru(:,:,1) + dslr*dtskin
            hfluxn(:,:,1) = clamb*(tskin - stl_am)
        end if

        rdth  = fstab/dtheta
        astab = 1.0
        if (lscasym) astab = 0.5   ! to get smaller dS/dT in stable conditions

        do i = 1, ix
            do j = 1, il
                if (tsea(i,j) > t2(i,j,2)) then
                   dths = min(dtheta, tsea(i,j) - t2(i,j,2))
                else
                   dths = max(-dtheta, astab*(tsea(i,j) - t2(i,j,2)))
                end if
                denvvs(i,j,2) = denvvs(i,j,0)*(1.0 + dths*rdth)
            end do
        end do

        if (fhum0 > 0.0) then
            call shtorh(-1, t1(:,:,2), psa, 1.0, q1(:,:,2), rh(:,:,kx), qsat0(:,:,2))

            q1(:,:,2) = fhum0*q1(:,:,2) + ghum0*qa(:,:,kx)
        else
            q1(:,:,2) = qa(:,:,kx)
        end if

        ! 4.2 Wind stress
        ks = 2

        cdsdv = cds*denvvs(:,:,ks)
        ustr(:,:,2) = -cdsdv*ua(:,:,kx)
        vstr(:,:,2) = -cdsdv*va(:,:,kx)
    end if

    ! =========================================================================
    ! Sea surface
    ! =========================================================================

    ! 4.3 Sensible heat flux
    shf(:,:,2) = chs*cp*denvvs(:,:,ks)*(tsea - t1(:,:,2))

    ! 4.4 Evaporation
    call shtorh(0, tsea, psa, 1.0, qdummy, rdummy, qsat0(:,:,2))
    evap(:,:,2) = chs*denvvs(:,:,ks)*(qsat0(:,:,2) - q1(:,:,2))

    ! 4.5 Emission of lw radiation from the surface
    !     and net heat fluxes into sea surface
    slru(:,:,2) = esbc*tsea**4.0
    hfluxn(:,:,2) = ssrd*(1.0 - alb_s) + slrd - slru(:,:,2) + shf(:,:,2) + alhc*evap(:,:,2)

    ! =========================================================================
    ! Weighted average of surface fluxes and temperatures according to land-sea
    ! mask
    ! =========================================================================

    if (lfluxland)  then
        ustr(:,:,3) = ustr(:,:,2) + fmask*(ustr(:,:,1) - ustr(:,:,2))
        vstr(:,:,3) = vstr(:,:,2) + fmask*(vstr(:,:,1) - vstr(:,:,2))
        shf(:,:,3)  = shf(:,:,2)  + fmask*(shf(:,:,1)  - shf(:,:,2))
        evap(:,:,3) = evap(:,:,2) + fmask*(evap(:,:,1) - evap(:,:,2))
        slru(:,:,3) = slru(:,:,2) + fmask*(slru(:,:,1) - slru(:,:,2))

        tsfc  = tsea      + fmask*(stl_am - tsea)
        tskin = tsea      + fmask*(tskin  - tsea)
        t0    = t1(:,:,2) + fmask*(t1(:,:,1) - t1(:,:,2))
    end if
end

! Compute orographic factor for land surface drag
! Input:   phi0 = surface geopotential
subroutine sflset(phi0)
    use mod_atparam
    use mod_sflcon
    use mod_dyncon1, only: grav

    implicit none

    real, intent(in) :: phi0(ix,il)
    real :: rhdrag

    rhdrag = 1.0/(grav*hdrag)

    forog = 1.0 + rhdrag*(1.0 - exp(-max(phi0, 0.0)*rhdrag))
end
