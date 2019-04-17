subroutine phypar(vor1,div1,t1,q1,phi1,psl1,utend,vtend,ttend,qtend)
    !  subroutine phypar(ug1,vg1,tg1,q1,phi1,psl1,
    ! &                   utend,vtend,ttend,qtend)
    !
    !  Purpose: compute physical parametrization tendencies for u, v, t, q
    !  and add them to dynamical grid-point tendencies
    !  Input-only  arguments:   vor1   : vorticity (sp)
    !                           div1   : divergence (sp)
    !                           t1     : temperature (sp)
    !                           q1     : specific humidity (sp)
    !                           phi1   : geopotential (sp)
    !                           psl1   : log of sfc pressure (sp)
    !  Input-output arguments:  utend  : u-wind tendency (gp)
    !                           vtend  : v-wind tendency (gp)
    !                           ttend  : temp. tendency (gp)
    !                           qtend  : spec. hum. tendency (gp)

    use mod_cpl_flags, only: icsea
    use mod_lflags, only: lradsw
    use mod_atparam
    use mod_physcon, only: sig, sigh, grdsig, grdscp, cp
    use mod_surfcon, only: phis0
	use mod_cli_land, only: fmask_l
    use mod_var_land, only: stl_am, soilw_am
    use mod_var_sea, only: sst_am, ssti_om
    use mod_physvar
    use mod_sppt, only: mu, gen_sppt
    use mod_tsteps, only: sppt_on

    implicit none

    integer, parameter :: ngp=ix*il

    complex, dimension(mx,nx,kx), intent(in) :: vor1, div1, t1, q1, phi1
    complex, dimension(mx,nx), intent(in) :: psl1
    real, dimension(ngp,kx), intent(inout) :: utend, vtend, ttend, qtend

    complex, dimension(mx,nx) :: ucos, vcos
    real, dimension(ngp) :: pslg1, rps, gse
    real, dimension(ngp,kx) :: ug1, vg1, tg1, qg1, phig1, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn
    real, dimension(ngp,kx) :: se, rh, qsat
    real, dimension(ngp) :: psg, ts, tskin, u0, v0, t0, q0, cloudc, clstr, cltop, prtop
    real, dimension(ngp,kx) :: tt_cnv, qt_cnv, tt_lsc, qt_lsc, tt_rsw, tt_rlw, ut_pbl, vt_pbl,&
        & tt_pbl, qt_pbl
    integer :: iptop(ngp), icltop(ngp,2), icnv(ngp), j, k
    real :: sppt(ngp,kx)

    ! Keep a copy of the original (dynamics only) tendencies
    utend_dyn = utend
    vtend_dyn = vtend
    ttend_dyn = ttend
    qtend_dyn = qtend

	! =========================================================================
    ! Compute grid-point fields
    ! =========================================================================

    ! Convert model spectral variables to grid-point variables
    do k = 1, kx
		call uvspec(vor1(:,:,k), div1(:,:,k), ucos, vcos)
		call grid(ucos, ug1(:,k), 2)
		call grid(vcos, vg1(:,k), 2)
		call grid(t1(:,:,k), tg1(:,k), 1)
		call grid(q1(:,:,k), qg1(:,k), 1)
      	call grid(phi1(:,:,k),phig1(:,k),1)
    end do

    call grid(psl1,pslg1,1)

	! =========================================================================
    ! Compute thermodynamic variables
    ! =========================================================================

    do j = 1, ngp
     psg(j) = exp(pslg1(j))
     rps(j) = 1.0/psg(j)
    end do

    do k = 1, kx
        do j = 1, ngp
	        qg1(j,k) = max(qg1(j,k),0.0)
            se(j,k) = cp*tg1(j,k) + phig1(j,k)
        end do
    end do

    do k = 1, kx
        call shtorh(1, ngp, tg1(:,k), psg, sig(k), qg1(:,k), rh(:,k), qsat(:,k))
    end do

	! =========================================================================
    ! Precipitation
    ! =========================================================================

    ! Deep convection
    call convmf(psg, se, qg1, qsat, iptop, cbmf, precnv, tt_cnv, qt_cnv)

    do k = 2, kx
    	do j = 1, ngp
    		tt_cnv(j,k) = tt_cnv(j,k)*rps(j)*grdscp(k)
        	qt_cnv(j,k) = qt_cnv(j,k)*rps(j)*grdsig(k)
       	end do
    end do

    do j = 1, ngp
        icnv(j) = kx - iptop(j)
    end do

    ! Large-scale condensation
    call lscond(psg, qg1, qsat, iptop, precls, tt_lsc, qt_lsc)

    ttend = ttend + tt_cnv + tt_lsc
    qtend = qtend + qt_cnv + qt_lsc

	! =========================================================================
    ! Radiation (shortwave and longwave) and surface fluxes
    ! =========================================================================

    ! Compute shortwave tendencies and initialize lw transmissivity
    ! The shortwave radiation may be called at selected time steps
    if (lradsw) then
        do j = 1, ngp
            gse(j) = (se(j,kx-1) - se(j,kx))/(phig1(j,kx-1) - phig1(j,kx))
        end do

        call cloud(qg1, rh, precnv, precls, iptop, gse, fmask_l, icltop, cloudc, clstr)

        do j = 1, ngp
            cltop(j) = sigh(icltop(j,1) - 1)*psg(j)
            prtop(j) = float(iptop(j))
        end do

        call radsw(psg, qg1, icltop, cloudc, clstr, ssrd, ssr, tsr, tt_rsw)

        do k = 1, kx
            do j = 1, ngp
                tt_rsw(j,k) = tt_rsw(j,k)*rps(j)*grdscp(k)
            end do
        end do
    end if

    ! Compute downward longwave fluxes
    call radlw(-1, tg1, ts, slrd, slru(:,3), slr, olr, tt_rlw)

    ! Compute surface fluxes and land skin temperature
    call suflux(psg, ug1, vg1, tg1, qg1, rh, phig1, phis0, fmask_l, stl_am, sst_am, soilw_am, &
		& ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, q0, .true.)

    ! Recompute sea fluxes in case of anomaly coupling
    if (icsea > 0) then
       call suflux(psg, ug1, vg1, tg1, qg1, rh, phig1, phis0, fmask_l, stl_am, ssti_om, soilw_am, &
	   	& ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, q0, .false.)
    end if

    ! Compute upward longwave fluxes, convert them to tendencies and add
	! shortwave tendencies
    call radlw(1, tg1, ts, slrd, slru(:,3), slr, olr, tt_rlw)

    do k = 1, kx
        do j = 1, ngp
            tt_rlw(j,k) = tt_rlw(j,k)*rps(j)*grdscp(k)
            ttend(j,k) = ttend(j,k) + tt_rsw(j,k) + tt_rlw(j,k)
        end do
    end do

	! =========================================================================
    ! Planetary boundary later interactions with lower troposphere
    ! =========================================================================

    ! Vertical diffusion and shallow convection
    call vdifsc(ug1, vg1, se, rh, qg1, qsat, phig1, icnv, ut_pbl, vt_pbl, tt_pbl, qt_pbl)

    ! Add tendencies due to surface fluxes
    do j=1,ngp
        ut_pbl(j,kx) = ut_pbl(j,kx) + ustr(j,3)*rps(j)*grdsig(kx)
        vt_pbl(j,kx) = vt_pbl(j,kx) + vstr(j,3)*rps(j)*grdsig(kx)
        tt_pbl(j,kx) = tt_pbl(j,kx) + shf(j,3)*rps(j)*grdscp(kx)
        qt_pbl(j,kx) = qt_pbl(j,kx) + evap(j,3)*rps(j)*grdsig(kx)
    end do

    utend = utend + ut_pbl
    vtend = vtend + vt_pbl
    ttend = ttend + tt_pbl
    qtend = qtend + qt_pbl

    ! 5. Store all fluxes for coupling and daily-mean output
    call dmflux(1)

    ! Add SPPT noise
    if (sppt_on) then
        sppt = gen_sppt()

        ! The physical contribution to the tendency is *tend - *tend_dyn, where * is u, v, t, q
        do k = 1,kx
            utend(:,k) = (1 + sppt(:,k)*mu(k)) * (utend(:,k) - utend_dyn(:,k)) + utend_dyn(:,k)
            vtend(:,k) = (1 + sppt(:,k)*mu(k)) * (vtend(:,k) - vtend_dyn(:,k)) + vtend_dyn(:,k)
            ttend(:,k) = (1 + sppt(:,k)*mu(k)) * (ttend(:,k) - ttend_dyn(:,k)) + ttend_dyn(:,k)
            qtend(:,k) = (1 + sppt(:,k)*mu(k)) * (qtend(:,k) - qtend_dyn(:,k)) + qtend_dyn(:,k)
        end do
    end if
end
