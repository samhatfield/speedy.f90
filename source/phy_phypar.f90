subroutine phypar(vor1,div1,t1,q1,phi1,psl1,utend,vtend,ttend,qtend)
    !  subroutine phypar(vor1,div1,t1,q1,phi1,psl1,
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
    !  Modified common blocks:  phygr1, phygr2, phygr3, phyten, fluxes

    use mod_cpl_flags, only: icsea
    use mod_lflags, only: lradsw, lrandf
    use mod_atparam
    use mod_physcon, only: sig, sigh, grdsig, grdscp, cp
    use mod_surfcon, only: fmask1, phis0
    use mod_var_land, only: stl_am, soilw_am
    use mod_var_sea, only: sst_am, ssti_om
    use mod_physvar
    use mod_sppt, only: mu, gen_sppt
    use mod_tsteps, only: sppt_on

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    complex, dimension(mx,nx,nlev) :: vor1, div1, t1, q1, phi1
    complex, dimension(mx,nx) :: psl1, ucos, vcos

    real, dimension(ngp,nlev) :: utend, vtend, ttend, qtend
    real, dimension(ngp,nlev) :: utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn

    integer :: iptop(ngp), icltop(ngp,2), icnv(ngp), iitest=0, j, k
    real, dimension(ngp) :: rps, gse
    real :: sppt(ngp,kx)

    ! Keep a copy of the original (dynamics only) tendencies
    utend_dyn = utend
    vtend_dyn = vtend
    ttend_dyn = ttend
    qtend_dyn = qtend

    ! 1. Compute grid-point fields
    ! 1.1 Convert model spectral variables to grid-point variables
    if (iitest.eq.1) print *, ' 1.1 in phypar'

    do k=1,nlev
      call uvspec(vor1(1,1,k),div1(1,1,k),ucos,vcos)
      call grid(ucos,ug1(1,k),2)
      call grid(vcos,vg1(1,k),2)
    end do

    do k=1,nlev
      call grid(t1(1,1,k),  tg1(1,k),  1)
      call grid(q1(1,1,k),  qg1(1,k),  1)
      call grid(phi1(1,1,k),phig1(1,k),1)
    end do

    call grid(psl1,pslg1,1)

    ! Remove negative humidity values
    !call qneg (qg1)

    ! 1.2 Compute thermodynamic variables
    if (iitest.eq.1) print *, ' 1.2 in phypar'

    do j=1,ngp
     psg(j)=exp(pslg1(j))
     rps(j)=1./psg(j)
    end do

    do k=1,nlev
        do j=1,ngp
            ! Remove when qneg is implemented
	        qg1(j,k)=max(qg1(j,k),0.)
            se(j,k)=cp*tg1(j,k)+phig1(j,k)
        end do
    end do

    do k=1,nlev
        call shtorh(1,ngp,tg1(1,k),psg,sig(k),qg1(1,k),rh(1,k),qsat(1,k))
    end do

    ! 2. Precipitation 
    ! 2.1 Deep convection
    call convmf(psg,se,qg1,qsat,iptop,cbmf,precnv,tt_cnv,qt_cnv)

    do k=2,nlev
       do j=1,ngp
        tt_cnv(j,k) = tt_cnv(j,k)*rps(j)*grdscp(k)
        qt_cnv(j,k) = qt_cnv(j,k)*rps(j)*grdsig(k)
       end do
    end do

    do j=1,ngp
        icnv(j)=nlev-iptop(j)
    end do

    ! 2.2 Large-scale condensation
!fk#if !defined(KNMI)
    call lscond(psg,qg1,qsat,iptop,precls,tt_lsc,qt_lsc)
!fk#else
!fk      call lscond (psg,qg1,qsat,ts,iptop,precls,snowls,tt_lsc,qt_lsc)
!fk#end if

    ttend = ttend + tt_cnv + tt_lsc
    qtend = qtend + qt_cnv + qt_lsc

    ! 3. Radiation (shortwave and longwave) and surface fluxes
    ! 3.1 Compute shortwave tendencies and initialize lw transmissivity
    if (iitest.eq.1) print *, ' 3.1 in PHYPAR'

    ! The sw radiation may be called at selected time steps
    if (lradsw) then
        do j=1,ngp
            gse(j) = (se(j,nlev-1)-se(j,nlev))/(phig1(j,nlev-1)-phig1(j,nlev))
        end do
  
        call cloud(qg1,rh,precnv,precls,iptop,gse,fmask1,icltop,cloudc,clstr)
  
        do j=1,ngp
            cltop(j)=sigh(icltop(j,1)-1)*psg(j)
            prtop(j)=float(iptop(j))
        end do
  
        call radsw(psg,qg1,icltop,cloudc,clstr,ssrd,ssr,tsr,tt_rsw)
  
        do k=1,nlev
            do j=1,ngp
                tt_rsw(j,k)=tt_rsw(j,k)*rps(j)*grdscp(k)
            end do
        end do
    end if

    ! 3.2 Compute downward longwave fluxes 
    call radlw(-1,tg1,ts,slrd,slru(1,3),slr,olr,tt_rlw)

    ! 3.3. Compute surface fluxes and land skin temperature
    if (iitest.eq.1) then 
        print *, ' 3.3 in PHYPAR'
        print *, 'mean(STL_AM) =', sum(STL_AM(:))/ngp
        print *, 'mean(SST_AM) =', sum(SST_AM(:))/ngp
    end if

    call suflux(psg,ug1,vg1,tg1,qg1,rh,phig1,phis0,fmask1,stl_am,sst_am,&
        & soilw_am,ssrd,slrd,ustr,vstr,shf,evap,slru,hfluxn,ts,tskin,u0,v0,t0,&
        & q0,.true.)
     
    ! 3.3.1. Recompute sea fluxes in case of anomaly coupling
    if (icsea .gt. 0) then 
       call suflux(psg,ug1,vg1,tg1,qg1,rh,phig1,phis0,fmask1,stl_am,ssti_om,&
           & soilw_am,ssrd,slrd,ustr,vstr,shf,evap,slru,hfluxn,ts,tskin,u0,v0,&
           & t0,q0,.false.)
    end if   

    ! 3.4 Compute upward longwave fluxes, convert them to tendencies 
    !     and add shortwave tendencies
    if (iitest.eq.1) print *, ' 3.4 in PHYPAR'

    call radlw (1,tg1,ts,slrd,slru(1,3),slr,olr,tt_rlw)

    do k=1,nlev
        do j=1,ngp
            tt_rlw(j,k) = tt_rlw(j,k)*rps(j)*grdscp(k)
            ttend(j,k) = ttend(j,k)+tt_rsw(j,k)+tt_rlw(j,k)
        end do
    end do

    ! 4. PBL interactions with lower troposphere
    ! 4.1 Vertical diffusion and shallow convection
    call vdifsc(ug1,vg1,se,rh,qg1,qsat,phig1,icnv,ut_pbl,vt_pbl,tt_pbl,qt_pbl)

    ! 4.2 Add tendencies due to surface fluxes 
    do j=1,ngp
        ut_pbl(j,nlev)=ut_pbl(j,nlev)+ustr(j,3)*rps(j)*grdsig(nlev)
        vt_pbl(j,nlev)=vt_pbl(j,nlev)+vstr(j,3)*rps(j)*grdsig(nlev)
        tt_pbl(j,nlev)=tt_pbl(j,nlev)+ shf(j,3)*rps(j)*grdscp(nlev)
        qt_pbl(j,nlev)=qt_pbl(j,nlev)+evap(j,3)*rps(j)*grdsig(nlev)
    end do

    utend = utend + ut_pbl
    vtend = vtend + vt_pbl
    ttend = ttend + tt_pbl
    qtend = qtend + qt_pbl

    ! 5. Store all fluxes for coupling and daily-mean output
    call dmflux(1)

    ! 6. Random diabatic forcing 
    if (lrandf) then
        ! 6.1 Compute zonal-mean cross sections of diabatic forcing
        if (lradsw) then
          call xs_rdf(tt_lsc,tt_cnv,1)
          call xs_rdf(tt_rsw,tt_rlw,2)
        end if

        ! 6.2 Compute and store 3-D pattern of random diabatic forcing
        tt_cnv = tt_cnv + tt_lsc

        call setrdf(tt_lsc)

        ttend = ttend + tt_lsc
    end if

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

subroutine xs_rdf(tt1,tt2,ivm)
    !  subroutine xs_rdf (tt1,tt2,ivm)
    !
    !  Purpose: compute zonal-mean cross-sec. of random diabatic forcing
    !  Input: tt1, tt2 = diabatic heating fields
    !         ivm      = index of vertical mode (1 or 2)

    use mod_atparam
    use mod_physcon, only: sig
    use mod_randfor, only: randfv

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real, dimension(nlon,nlat,nlev), intent(in) :: tt1, tt2
    integer, intent(in) :: ivm

    real :: rand1(0:nlat+1), pigr2, rnlon, rnsig
    integer :: i, j, k, nsmooth

    rnlon = 1./float(nlon)
    pigr2 = 4.*asin(1.)

    ! 1. Compute cross sections
    do k=1,nlev
         if (ivm.eq.1) then
            rnsig = rnlon
         else
            rnsig = rnlon*sin(pigr2*sig(k))
         endif

         do j=1,nlat
           randfv(j,k,ivm) = 0.
           do i=1,nlon
              randfv(j,k,ivm) = randfv(j,k,ivm)+tt1(i,j,k)+tt2(i,j,k)
           end do
           randfv(j,k,ivm) = randfv(j,k,ivm)*rnsig 
         end do
    end do

    ! 2. Perform smoothing in latitude
    do nsmooth=1,2
        do k=1,nlev

          do j=1,nlat
             rand1(j) = randfv(j,k,ivm)
          end do
          rand1(0) = rand1(2)
          rand1(nlat+1) = rand1(nlat-1)
             
          do j=1,nlat
             randfv(j,k,ivm) = 0.5*rand1(j)+0.25*(rand1(j-1)+rand1(j+1))
          end do
        end do
    end do
end

subroutine setrdf(tt_rdf)
    !  subroutine setrdf (tt_rdf)
    !
    !  Purpose: compute 3-D pattern of random diabatic forcing
    !  Output: tt_rdf = random diabatic forcing

    use mod_atparam
    use mod_randfor

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real :: tt_rdf(nlon,nlat,nlev)
    integer :: i, j, k

    do k=1,nlev
        do j=1,nlat
           do i=1,nlon
               tt_rdf(i,j,k) = randfh(i,j,1)*randfv(j,k,1)&
                   & +randfh(i,j,2)*randfv(j,k,2)
           end do
         end do
     end do
end
