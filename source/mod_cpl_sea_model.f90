module mod_cpl_sea_model
    use mod_atparam

    implicit none

    private
    public sea_model_init, couple_sea_atm
    public fmask_s, bmask_s, deglat_s, sst12, sice12, sstan3, hfseacl, sstom12
    public sstcl_ob, sst_am, sice_am, tice_am, ssti_om
    public hflux_s, hflux_i

    ! Constant parameters and fields in sea/ice model
    real :: rhcaps(ix,il) ! 1./heat_capacity (sea)
    real :: rhcapi(ix,il) ! 1./heat_capacity (ice)
    real :: cdsea(ix,il) ! 1./dissip_time (sea)
    real :: cdice(ix,il) ! 1./dissip_time (ice)
    real :: beta = 1.0 ! Heat flux coef. at sea/ice int.

    ! Sea masks
    real :: fmask_s(ix,il) ! Fraction of sea
    real :: bmask_s(ix,il) ! Binary sea mask
    real :: deglat_s(il) ! Grid latitudes

    ! Monthly-mean climatological fields over sea
    real :: sst12(ix,il,12) ! Sea/ice surface temperature
    real :: sice12(ix,il,12) ! Sea ice fraction

    ! SST anomaly fields
    real :: sstan3(ix,il,3) ! SST anomaly in 3 consecutive months

    ! Climatological fields from model output
    real :: hfseacl(ix,il) ! Annual-mean heat flux into sea sfc.
    real :: sstom12(ix,il,12) ! Ocean model SST climatology

    ! Daily observed climatological fields over sea
    real :: sstcl_ob(ix*il) ! Observed clim. SST
    real :: sicecl_ob(ix*il) ! Clim. sea ice fraction
    real :: ticecl_ob(ix*il) ! Clim. sea ice temperature
    real :: sstan_ob(ix*il) ! Daily observed SST anomaly

    ! Daily climatological fields from ocean model
    real :: sstcl_om(ix*il) ! Ocean model clim. SST

    ! Sea sfc. fields used by atmospheric model
    real :: sst_am(ix*il) ! SST (full-field)
    real :: sstan_am(ix*il) ! SST anomaly
    real :: sice_am(ix*il) ! Sea ice fraction
    real :: tice_am(ix*il) ! Sea ice temperature

    ! Sea sfc. fields from ocean/sea-ice model
    real :: sst_om(ix*il) ! Ocean model SST
    real :: sice_om(ix*il) ! Model sea ice fraction
    real :: tice_om(ix*il) ! Model sea ice temperature
    real :: ssti_om(ix*il) ! Model SST + sea ice temp.

    ! Weight for obs. SST anomaly in coupled runs
    real :: wsst_ob(ix*il)

    ! Fluxes at sea surface (all downward, except evaporation)
    real :: hflux_s(ix*il) ! Net heat flux into sea surface
    real :: hflux_i(ix*il) ! Net heat flux into sea-ice surface

contains
    ! Initialization of sea model
    subroutine sea_model_init
        integer, parameter :: ngp=ix*il

        ! Domain mask
        real :: dmask(ix,il)

        ! Domain flags
        logical :: l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian

        ! Heat capacities of mixed-layer and sea-ice
        real :: hcaps(il)
        real :: hcapi(il)

        integer :: i, j
        real :: coslat, crad

        ! 1. Set geographical domain, heat capacities and dissipation times
        !    for sea (mixed layer) and sea-ice

        ! Model parameters (default values)

        ! ocean mixed layer depth: d + (d0-d)*(cos_lat)^3
        real :: depth_ml = 60.               ! High-latitude depth
        real :: dept0_ml = 40.               ! Minimum depth (tropics)

        ! sea-ice depth : d + (d0-d)*(cos_lat)^2
        real :: depth_ice = 2.5              ! High-latitude depth
        real :: dept0_ice = 1.5              ! Minimum depth

        ! Dissipation time (days) for sea-surface temp. anomalies
        real :: tdsst  = 90.

        ! Dissipation time (days) for sea-ice temp. anomalies
        real :: dice = 30.

        ! Minimum fraction of sea for the definition of anomalies
        real :: fseamin = 1./3.

        ! Dissipation time (days) for sea-ice temp. anomalies
        real :: tdice

        ! Geographical domain
        ! note : more than one regional domain may be set .true.
        l_globe  =  .true.         ! global domain
        l_northe = .false.         ! Northern hem. oceans (lat > 20N)
        l_natlan = .false.         ! N. Atlantic (lat 20-80N, lon 100W-45E)
        l_npacif = .false.         ! N. Pacific  (lat 20-80N, lon 100E-100W)
        l_tropic = .false.         ! Tropics (lat 30S-30N)
        l_indian = .false.         ! Indian Ocean (lat 30S-30N, lon 30-120E)

        ! Reset model parameters
        include "cls_insea.h"

        ! Heat capacities per m^2 (depth*heat_cap/m^3)
        crad=asin(1.)/90.
        do j=1,il
            coslat   = cos(crad*deglat_s(j))
            hcaps(j) = 4.18e+6*(depth_ml +(dept0_ml -depth_ml) *coslat**3)
            hcapi(j) = 1.93e+6*(depth_ice+(dept0_ice-depth_ice)*coslat**2)
        end do

        ! 3. Compute constant parameters and fields

        ! Set domain mask
        if (l_globe) then
            dmask(:,:) = 1.
        else
            dmask(:,:) = 0.
            if (l_northe) call sea_domain('northe',dmask)
            if (l_natlan) call sea_domain('natlan',dmask)
            if (l_npacif) call sea_domain('npacif',dmask)
            if (l_tropic) call sea_domain('tropic',dmask)
            if (l_indian) call sea_domain('indian',dmask)
        end if

        ! Smooth latitudinal boundaries and blank out land points
        do j=2,il-1
            rhcaps(:,j) = 0.25*(dmask(:,j-1)+2*dmask(:,j)+dmask(:,j+1))
        end do
        dmask(:,2:il-1) = rhcaps(:,2:il-1)

        do j=1,il
            do i=1,ix
                if (fmask_s(i,j).lt.fseamin) dmask(i,j) = 0
            end do
        end do

        ! Set heat capacity and dissipation time over selected domain
        do j=1,il
            rhcaps(:,j) = 86400./hcaps(j)
            rhcapi(:,j) = 86400./hcapi(j)
        end do

        cdsea = dmask*tdsst/(1.+dmask*tdsst)
        cdice = dmask*tdice/(1.+dmask*tdice)
    end

    subroutine couple_sea_atm(day)
        use mod_cpl_flags, only: icsea, icice, isstan
        use mod_date, only: model_datetime, imont1
        use mod_cpl_bcinterp, only: forin5, forint

        integer, intent(in) :: day
        integer, parameter :: ngp=ix*il

        integer :: j
        real :: sstcl0, sstfr

        ! 1. Interpolate climatological fields and obs. SST anomaly
        !    to actual date

        ! Climatological SST
        call forin5(imont1,sst12,sstcl_ob)

        ! Climatological sea ice fraction
        call forint(imont1,sice12,sicecl_ob)

        ! SST anomaly
        if (isstan.gt.0) then
            if (model_datetime%day.eq.1.and.day.gt.0) call obs_ssta
            call forint (2,sstan3,sstan_ob)
        end if

        ! Ocean model climatological SST
        if (icsea.ge.3) then
            call forin5 (imont1,sstom12,sstcl_om)
        end if

        ! Adjust climatological fields over sea ice

        ! SST at freezing point
        sstfr = 273.2-1.8

        do j=1,ngp
            sstcl0 = sstcl_ob(j)

            if (sstcl_ob(j).gt.sstfr) then
                sicecl_ob(j) = min(0.5,sicecl_ob(j))
                ticecl_ob(j) = sstfr
                if (sicecl_ob(j).gt.0.) then
                    sstcl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/(1.-sicecl_ob(j))
                end if
            else
                sicecl_ob(j) = max(0.5,sicecl_ob(j))
                ticecl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/sicecl_ob(j)
                !ticecl_ob(j) = sstcl_ob(j)
                sstcl_ob(j)  = sstfr
            end if

            if (icsea.ge.3) sstcl_om(j) = sstcl_om(j)+(sstcl_ob(j)-sstcl0)
        end do

        if (day == 0) then
            ! 2. Initialize prognostic variables of ocean/ice model
            !    in case of no restart or no coupling
            sst_om(:)  = sstcl_ob(:)      ! SST
            tice_om(:) = ticecl_ob(:)     ! sea ice temperature
            sice_om(:) = sicecl_ob(:)     ! sea ice fraction

            if (icsea.le.0) sst_om(:) = 0.

            ! 3. Compute additional sea/ice variables
            wsst_ob(:) = 0.
            if (icsea.ge.4) call sea_domain('elnino',wsst_ob)
        else
            if (icsea > 0 .or. icice > 0) then
                ! 1. Run ocean mixed layer or
                !    call message-passing routines to receive data from ocean model
                call sea_model
            end if
        end if

        ! 3. Compute sea-sfc. anomalies and full fields for atm. model
        ! 3.1 SST
        sstan_am(:) = 0.

        if (icsea.le.1) then
            if (isstan.gt.0) sstan_am(:) = sstan_ob(:)

            ! Use observed SST (climatological or full field)
            sst_am(:) = sstcl_ob(:) + sstan_am(:)
        else if (icsea.eq.2) then
            ! Use full ocean model SST
            sst_am(:) = sst_om(:)
        else if (icsea.ge.3) then
            ! Define SST anomaly from ocean model ouput and climatology
            sstan_am(:) = sst_om(:) - sstcl_om(:)

            ! Merge with observed SST anomaly in selected area
            if (icsea.ge.4) then
                sstan_am(:) = sstan_am(:) + wsst_ob(:)*(sstan_ob(:)-sstan_am(:))
            end if

            ! Add observed SST climatology to model SST anomaly
            sst_am(:) = sstcl_ob(:) + sstan_am(:)
        end if

        ! 3.2 Sea ice fraction and temperature
        if (icice.gt.0) then
            sice_am(:) = sice_om(:)
            tice_am(:) = tice_om(:)
        else
            sice_am(:) = sicecl_ob(:)
            tice_am(:) = ticecl_ob(:)
        end if

        sst_am(:)  = sst_am(:)+sice_am(:)*(tice_am(:)-sst_am(:))
        ssti_om(:) = sst_om(:)+sice_am(:)*(tice_am(:)-sst_om(:))
    end subroutine

    ! Update observed SST anomaly array
    subroutine obs_ssta
        use mod_date, only: model_datetime, start_datetime
        use mod_tsteps, only: issty0
        use mod_input, only: load_boundary_file

        integer :: i, j, next_month

        sstan3(:,:,1) = sstan3(:,:,2)
        sstan3(:,:,2) = sstan3(:,:,3)

        ! Compute next month given initial SST year
        next_month = (start_datetime%year - issty0) * 12 + model_datetime%month

        ! Read next month SST anomalies
        sstan3(:,:,3) = load_boundary_file(30,next_month-1)

        call forchk(bmask_s, 1, -50.0, 50.0, 0.0, sstan3(:,:,3))
    end

    ! Purpose : Integrate slab ocean and sea-ice models for one day
    subroutine sea_model
        integer, parameter :: ngp=ix*il

        ! Input variables:
        real ::  sst0(ix,il)     ! SST at initial time
        real :: tice0(ix,il)     ! sea ice temp. at initial time
        real :: sice0(ix,il)     ! sea ice fraction at initial time
        real :: hfsea(ix,il)     ! sea+ice  sfc. heat flux between t0 and t1
        real :: hfice(ix,il)     ! ice-only sfc. heat flux between t0 and t1

        real ::  sstcl1(ix,il)   ! clim. SST at final time
        real :: ticecl1(ix,il)   ! clim. sea ice temp. at final time

        ! Output variables
        real ::  sst1(ix,il)     ! SST at final time
        real :: tice1(ix,il)     ! sea ice temp. at final time
        real :: sice1(ix,il)     ! sea ice fraction at final time

        ! Auxiliary variables
        real :: hflux(ix,il)   ! net sfc. heat flux
        real :: tanom(ix,il)   ! sfc. temperature anomaly
        real :: cdis(ix,il)    ! dissipation ceofficient

        real :: anom0, sstfr

        sst0 = reshape(sst_om, (/ix, il/))
        tice0 = reshape(tice_om, (/ix, il/))
        sice0 = reshape(sicecl_ob, (/ix, il/))
        hfsea = reshape(hflux_s, (/ix, il/))
        hfice = reshape(hflux_i, (/ix, il/))
        sstcl1 = reshape(sstcl_ob, (/ix, il/))
        ticecl1 = reshape(ticecl_ob, (/ix, il/))

        sstfr = 273.2-1.8       ! SST at freezing point

        !beta = 1.               ! heat flux coef. at sea-ice bottom

        ! 1. Ocean mixed layer
        ! Net heat flux
        hflux = hfsea-hfseacl-sice0*(hfice+beta*(sstfr-tice0))

        ! Anomaly at t0 minus climatological temp. tendency
        tanom = sst0 - sstcl1

        ! Time evoloution of temp. anomaly
        tanom = cdsea*(tanom+rhcaps*hflux)

        ! Full SST at final time
        sst1 = tanom + sstcl1

        ! 2. Sea-ice slab model

        ! Net heat flux
        hflux = hfice + beta*(sstfr-tice0)

        ! Anomaly w.r.t final-time climatological temp.
        tanom = tice0 - ticecl1

        ! Definition of non-linear damping coefficient
        anom0     = 20.
        cdis = cdice*(anom0/(anom0+abs(tanom)))
        !cdis(:,:) = cdice(:,:)

        ! Time evoloution of temp. anomaly
        tanom = cdis*(tanom+rhcapi*hflux)

        ! Full ice temperature at final time
        tice1 = tanom + ticecl1

        ! Persistence of sea ice fraction
        sice1 = sice0

        sst_om = reshape(sst1, (/ngp/))
        tice_om = reshape(tice1, (/ngp/))
        sice_om = reshape(sice1, (/ngp/))
    end

    ! Definition of ocean domains
    subroutine sea_domain(cdomain,dmask)
        character(len=6), intent(in) :: cdomain           ! domain name

        ! Output variables (initialized by calling routine)
        real, intent(inout) :: dmask(ix,il)         ! domain mask

        integer :: i, j
        real :: arlat, dlon, rlon, rlonw, wlat

        print *, 'sea domain : ', cdomain

        dlon = 360./float(ix)

        if (cdomain.eq.'northe') then
            do j=1,il
                if (deglat_s(j).gt.20.0) dmask(:,j) = 1.
            end do
        end if

        if (cdomain.eq.'natlan') then
             do j=1,il
               if (deglat_s(j).gt.20.0.and.deglat_s(j).lt.80.0) then
                 do i=1,ix
                   rlon = (i-1)*dlon
                   if (rlon.lt.45.0.or.rlon.gt.260.0) dmask(i,j) = 1.
                 end do
               end if
             end do
        end if

        if (cdomain.eq.'npacif') then
            do j=1,il
                if (deglat_s(j).gt.20.0.and.deglat_s(j).lt.65.0) then
                    do i=1,ix
                        rlon = (i-1)*dlon
                        if (rlon.gt.120.0.and.rlon.lt.260.0) dmask(i,j) = 1.
                    end do
                end if
            end do
        end if

        if (cdomain.eq.'tropic') then
            do j=1,il
                if (deglat_s(j).gt.-30.0.and.deglat_s(j).lt.30.0) dmask(:,j) = 1.
            end do
        end if

        if (cdomain.eq.'indian') then
            do j=1,il
                if (deglat_s(j).gt.-30.0.and.deglat_s(j).lt.30.0) then
                    do i=1,ix
                        rlon = (i-1)*dlon
                        if (rlon.gt.30.0.and.rlon.lt.120.0) dmask(i,j) = 1.
                    end do
                end if
            end do
        end if

        if (cdomain.eq.'elnino') then
            do j=1,il
                arlat = abs(deglat_s(j))
                if (arlat.lt.25.0) then
                    wlat = 1.
                    if (arlat.gt.15.0) wlat = (0.1*(25.-arlat))**2
                    rlonw = 300.-2*max(deglat_s(j),0.)
                    do i=1,ix
                        rlon = (i-1)*dlon
                        if (rlon.gt.165.0.and.rlon.lt.rlonw) then
                            dmask(i,j) = wlat
                        else if (rlon.gt.155.0.and.rlon.lt.165.0) then
                            dmask(i,j) = wlat*0.1*(rlon-155.)
                        end if
                    end do
                end if
            end do
        end if
    end
end module
