module physics
    use mod_atparam

    implicit none

    private
    public precnv, precls, snowcv, snowls, cbmf, tsr, ssrd, ssr, slrd, slr,&
        & olr, slru, ustr, vstr, shf, evap, hfluxn
    public initialize_physics, get_physical_tendencies

    ! Physical variables shared among all physics schemes
    real, dimension(ix,il)   :: precnv ! Convective precipitation  [g/(m^2 s)], total
    real, dimension(ix,il)   :: precls ! Large-scale precipitation [g/(m^2 s)], total
    real, dimension(ix,il)   :: snowcv ! Convective precipitation  [g/(m^2 s)], snow only
    real, dimension(ix,il)   :: snowls ! Large-scale precipitation [g/(m^2 s)], snow only
    real, dimension(ix,il)   :: cbmf   ! Cloud-base mass flux
    real, dimension(ix,il)   :: tsr    ! Top-of-atmosphere shortwave radiation (downward)
    real, dimension(ix,il)   :: ssrd   ! Surface shortwave radiation (downward-only)
    real, dimension(ix,il)   :: ssr    ! Surface shortwave radiation (net downward)
    real, dimension(ix,il)   :: slrd   ! Surface longwave radiation (downward-only)
    real, dimension(ix,il)   :: slr    ! Surface longwave radiation (net upward)
    real, dimension(ix,il)   :: olr    ! Outgoing longwave radiation (upward)
    real, dimension(ix,il,3) :: slru   ! Surface longwave emission (upward)

    ! Third dimension -> 1:land, 2:sea, 3: weighted average
    real, dimension(ix,il,3) :: ustr   ! U-stress
    real, dimension(ix,il,3) :: vstr   ! V-stress
    real, dimension(ix,il,3) :: shf    ! Sensible heat flux
    real, dimension(ix,il,3) :: evap   ! Evaporation [g/(m^2 s)]
    real, dimension(ix,il,3) :: hfluxn ! Net heat flux into surface

contains
    ! Initialize physical parametrization routines
    subroutine initialize_physics
        use physical_constants
        use mod_dyncon1, only: grav, hsg, radang

        integer :: j, k

        ! 1.2 Functions of sigma and latitude
        sigh(0) = hsg(1)

        do k = 1, kx
            sig(k)  = 0.5*(hsg(k+1)+hsg(k))
            sigl(k) = log(sig(k))
            sigh(k) = hsg(k+1)
            dsig(k) = hsg(k+1)-hsg(k)
            grdsig(k) = grav/(dsig(k)*p0)
            grdscp(k) = grdsig(k)/cp
        end do

        ! Weights for vertical interpolation at half-levels(1,kx) and surface
        ! Note that for phys.par. half-lev(k) is between full-lev k and k+1
        ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
        ! Fsurf = Ffull(kx)+WVI(kx,2)*(Ffull(kx)-Ffull(kx-1))
        do k = 1, kx-1
            wvi(k,1) = 1./(sigl(k+1)-sigl(k))
            wvi(k,2) = (log(sigh(k))-sigl(k))*wvi(k,1)
        end do

        wvi(kx,1) = 0.
        wvi(kx,2) = (log(0.99)-sigl(kx))*wvi(kx-1,1)

        do j = 1, il
            slat(j) = sin(radang(j))
            clat(j) = cos(radang(j))
        end do
    end

    ! Compute physical parametrization tendencies for u, v, t, q and add them to
    ! dynamical grid-point tendencies
    ! Input-only  arguments:   vor   : vorticity (sp)
    !                          div   : divergence (sp)
    !                          t     : temperature (sp)
    !                          q     : specific humidity (sp)
    !                          phi   : geopotential (sp)
    !                          psl   : log of sfc pressure (sp)
    ! Input-output arguments:  utend  : u-wind tendency (gp)
    !                          vtend  : v-wind tendency (gp)
    !                          ttend  : temp. tendency (gp)
    !                          qtend  : spec. hum. tendency (gp)
    subroutine get_physical_tendencies(vor, div, t, q, phi, psl, utend, vtend, ttend, qtend)
        use physical_constants, only: sig, sigh, grdsig, grdscp, cp
        use boundaries, only: phis0
        use land_model, only: fmask_l
        use sea_model, only: sst_am, ssti_om, sea_coupling_flag
        use sppt, only: mu, gen_sppt
        use mod_tsteps, only: sppt_on
        use precipitation, only: convective_precipitation, large_scale_precipitation
        use shortwave_radiation, only: get_shortwave_rad_fluxes, clouds, compute_shortwave
        use longwave_radiation, only: get_longwave_rad_fluxes
        use surface_fluxes, only: get_surface_fluxes
        use vertical_diffusion, only: get_vertical_diffusion_tend
        use humidity, only: spec_hum_to_rel_hum

        complex, dimension(mx,nx,kx), intent(in) :: vor, div, t, q, phi
        complex, dimension(mx,nx), intent(in) :: psl
        real, dimension(ix,il,kx), intent(inout) :: utend, vtend, ttend, qtend

        complex, dimension(mx,nx) :: ucos, vcos
        real, dimension(ix,il) :: pslg, rps, gse
        real, dimension(ix,il,kx) :: ug, vg, tg, qg, phig, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn
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
    		call uvspec(vor(:,:,k), div(:,:,k), ucos, vcos)
    		call grid(ucos, ug(:,:,k), 2)
    		call grid(vcos, vg(:,:,k), 2)
    		call grid(t(:,:,k), tg(:,:,k), 1)
    		call grid(q(:,:,k), qg(:,:,k), 1)
          	call grid(phi(:,:,k),phig(:,:,k),1)
        end do

        call grid(psl,pslg,1)

        ! =========================================================================
        ! Compute thermodynamic variables
        ! =========================================================================

    	psg = exp(pslg)
    	rps = 1.0/psg

    	qg = max(qg, 0.0)
    	se = cp*tg + phig

        do k = 1, kx
            call spec_hum_to_rel_hum(tg(:,:,k), psg, sig(k), qg(:,:,k), rh(:,:,k), qsat(:,:,k))
        end do

        ! =========================================================================
        ! Precipitation
        ! =========================================================================

        ! Deep convection
        call convective_precipitation(psg, se, qg, qsat, iptop, cbmf, precnv, tt_cnv, qt_cnv)

        do k = 2, kx
    		tt_cnv(:,:,k) = tt_cnv(:,:,k)*rps*grdscp(k)
    		qt_cnv(:,:,k) = qt_cnv(:,:,k)*rps*grdsig(k)
        end do

        icnv = kx - iptop

        ! Large-scale condensation
        call large_scale_precipitation(psg, qg, qsat, iptop, precls, tt_lsc, qt_lsc)

        ttend = ttend + tt_cnv + tt_lsc
        qtend = qtend + qt_cnv + qt_lsc

        ! =========================================================================
        ! Radiation (shortwave and longwave) and surface fluxes
        ! =========================================================================

        ! Compute shortwave tendencies and initialize lw transmissivity
        ! The shortwave radiation may be called at selected time steps
        if (compute_shortwave) then
    		gse = (se(:,:,kx-1) - se(:,:,kx))/(phig(:,:,kx-1) - phig(:,:,kx))

            call clouds(qg, rh, precnv, precls, iptop, gse, fmask_l, icltop, cloudc, clstr)

    		do i = 1, ix
    	        do j = 1, il
    	            cltop(i,j) = sigh(icltop(i,j,1) - 1)*psg(i,j)
    	            prtop(i,j) = float(iptop(i,j))
    	        end do
    		end do

            call get_shortwave_rad_fluxes(psg, qg, icltop, cloudc, clstr, ssrd, ssr, tsr, tt_rsw)

            do k = 1, kx
    			tt_rsw(:,:,k) = tt_rsw(:,:,k)*rps*grdscp(k)
            end do
        end if

        ! Compute downward longwave fluxes
        call get_longwave_rad_fluxes(-1, tg, ts, slrd, slru(:,:,3), slr, olr, tt_rlw)

        ! Compute surface fluxes and land skin temperature
        call get_surface_fluxes(psg, ug, vg, tg, qg, rh, phig, phis0, fmask_l, sst_am, &
    		& ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, q0, .true.)

        ! Recompute sea fluxes in case of anomaly coupling
        if (sea_coupling_flag > 0) then
           call get_surface_fluxes(psg, ug, vg, tg, qg, rh, phig, phis0, fmask_l, ssti_om, &
    	   	& ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, q0, .false.)
        end if

        ! Compute upward longwave fluxes, convert them to tendencies and add
    	! shortwave tendencies
        call get_longwave_rad_fluxes(1, tg, ts, slrd, slru(:,:,3), slr, olr, tt_rlw)

        do k = 1, kx
    		tt_rlw(:,:,k) = tt_rlw(:,:,k)*rps*grdscp(k)
        end do

    	ttend = ttend + tt_rsw + tt_rlw

        ! =========================================================================
        ! Planetary boundary later interactions with lower troposphere
        ! =========================================================================

        ! Vertical diffusion and shallow convection
        call get_vertical_diffusion_tend(ug, vg, se, rh, qg, qsat, phig, icnv, ut_pbl, vt_pbl, &
            & tt_pbl, qt_pbl)

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
end module
