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
    use mod_cpl_land_model, only: fmask_l
    use mod_cpl_sea_model, only: sst_am, ssti_om
    use mod_physvar
    use mod_sppt, only: mu, gen_sppt
    use mod_tsteps, only: sppt_on

    implicit none

    complex, dimension(mx,nx,kx), intent(in) :: vor1, div1, t1, q1, phi1
    complex, dimension(mx,nx), intent(in) :: psl1
    real, dimension(ix,il,kx), intent(inout) :: utend, vtend, ttend, qtend

    complex, dimension(mx,nx) :: ucos, vcos
    real, dimension(ix,il) :: pslg1, rps, gse
    real, dimension(ix,il,kx) :: ug1, vg1, tg1, qg1, phig1, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn
    real, dimension(ix,il,kx) :: se, rh, qsat
    real, dimension(ix,il) :: psg, ts, tskin, u0, v0, t0, q0, cloudc, clstr, cltop, prtop
    real, dimension(ix,il,kx) :: tt_cnv, qt_cnv, tt_lsc, qt_lsc, tt_rsw, tt_rlw, ut_pbl, vt_pbl,&
        & tt_pbl, qt_pbl
    integer :: iptop(ix,il), icltop(ix,il,2), icnv(ix,il), i, j, k
    real :: sppt(ix,il,kx)

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
		call grid(ucos, ug1(:,:,k), 2)
		call grid(vcos, vg1(:,:,k), 2)
		call grid(t1(:,:,k), tg1(:,:,k), 1)
		call grid(q1(:,:,k), qg1(:,:,k), 1)
      	call grid(phi1(:,:,k),phig1(:,:,k),1)
    end do

    call grid(psl1,pslg1,1)

    ! =========================================================================
    ! Compute thermodynamic variables
    ! =========================================================================

	do i = 1, ix
    	do j = 1, il
    		psg(i,j) = exp(pslg1(i,j))
     		rps(i,j) = 1.0/psg(i,j)
		end do
    end do

    do k = 1, kx
		do i = 1, ix
	        do j = 1, il
		        qg1(i,j,k) = max(qg1(i,j,k), 0.0)
	            se(i,j,k) = cp*tg1(i,j,k) + phig1(i,j,k)
	        end do
		end do
    end do

    do k = 1, kx
        call shtorh(1, tg1(:,:,k), psg, sig(k), qg1(:,:,k), rh(:,:,k), qsat(:,:,k))
    end do

    ! =========================================================================
    ! Precipitation
    ! =========================================================================

    ! Deep convection
    call convmf(psg, se, qg1, qsat, iptop, cbmf, precnv, tt_cnv, qt_cnv)

    do k = 2, kx
		tt_cnv(:,:,k) = tt_cnv(:,:,k)*rps*grdscp(k)
		qt_cnv(:,:,k) = qt_cnv(:,:,k)*rps*grdsig(k)
    end do

    icnv = kx - iptop

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
		do i = 1, ix
			do j = 1, il
	            gse(i,j) = (se(i,j,kx-1) - se(i,j,kx))/(phig1(i,j,kx-1) - phig1(i,j,kx))
			end do
        end do

        call cloud(qg1, rh, precnv, precls, iptop, gse, fmask_l, icltop, cloudc, clstr)

		do i = 1, ix
	        do j = 1, il
	            cltop(i,j) = sigh(icltop(i,j,1) - 1)*psg(i,j)
	            prtop(i,j) = float(iptop(i,j))
	        end do
		end do

        call radsw(psg, qg1, icltop, cloudc, clstr, ssrd, ssr, tsr, tt_rsw)

        do k = 1, kx
			tt_rsw(:,:,k) = tt_rsw(:,:,k)*rps*grdscp(k)
        end do
    end if

    ! Compute downward longwave fluxes
    call radlw(-1, tg1, ts, slrd, slru(:,:,3), slr, olr, tt_rlw)

    ! Compute surface fluxes and land skin temperature
    call suflux(psg, ug1, vg1, tg1, qg1, rh, phig1, phis0, fmask_l, sst_am, &
		& ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, q0, .true.)

    ! Recompute sea fluxes in case of anomaly coupling
    if (icsea > 0) then
       call suflux(psg, ug1, vg1, tg1, qg1, rh, phig1, phis0, fmask_l, ssti_om, &
	   	& ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, q0, .false.)
    end if

    ! Compute upward longwave fluxes, convert them to tendencies and add
	! shortwave tendencies
    call radlw(1, tg1, ts, slrd, slru(:,:,3), slr, olr, tt_rlw)

    do k = 1, kx
		tt_rlw(:,:,k) = tt_rlw(:,:,k)*rps*grdscp(k)
    end do

	ttend = ttend + tt_rsw + tt_rlw

    ! =========================================================================
    ! Planetary boundary later interactions with lower troposphere
    ! =========================================================================

    ! Vertical diffusion and shallow convection
    call vdifsc(ug1, vg1, se, rh, qg1, qsat, phig1, icnv, ut_pbl, vt_pbl, tt_pbl, qt_pbl)

    ! Add tendencies due to surface fluxes
	ut_pbl(:,:,kx) = ut_pbl(:,:,kx) + ustr(:,:,3)*rps*grdsig(kx)
	vt_pbl(:,:,kx) = vt_pbl(:,:,kx) + vstr(:,:,3)*rps*grdsig(kx)
	tt_pbl(:,:,kx) = tt_pbl(:,:,kx)  + shf(:,:,3)*rps*grdscp(kx)
	qt_pbl(:,:,kx) = qt_pbl(:,:,kx) + evap(:,:,3)*rps*grdsig(kx)

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
            utend(:,:,k) = (1 + sppt(:,:,k)*mu(k))*(utend(:,:,k) - utend_dyn(:,:,k)) &
				& + utend_dyn(:,:,k)
            vtend(:,:,k) = (1 + sppt(:,:,k)*mu(k))*(vtend(:,:,k) - vtend_dyn(:,:,k)) &
				& + vtend_dyn(:,:,k)
            ttend(:,:,k) = (1 + sppt(:,:,k)*mu(k))*(ttend(:,:,k) - ttend_dyn(:,:,k)) &
				& + ttend_dyn(:,:,k)
            qtend(:,:,k) = (1 + sppt(:,:,k)*mu(k))*(qtend(:,:,k) - qtend_dyn(:,:,k)) &
				& + qtend_dyn(:,:,k)
        end do
    end if
end
