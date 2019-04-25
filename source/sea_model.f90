module sea_model
    use params

    implicit none

    private
    public sea_model_init, couple_sea_atm
    public fmask_s
    public sstcl_ob, sst_am, sice_am, tice_am, ssti_om
    public hflux_s, hflux_i
    public sea_coupling_flag, sst_anomaly_coupling_flag

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
    real :: sstcl_ob(ix,il) ! Observed clim. SST
    real :: sicecl_ob(ix,il) ! Clim. sea ice fraction
    real :: ticecl_ob(ix,il) ! Clim. sea ice temperature
    real :: sstan_ob(ix,il) ! Daily observed SST anomaly

    ! Daily climatological fields from ocean model
    real :: sstcl_om(ix,il) ! Ocean model clim. SST

    ! Sea sfc. fields used by atmospheric model
    real :: sst_am(ix,il) ! SST (full-field)
    real :: sstan_am(ix,il) ! SST anomaly
    real :: sice_am(ix,il) ! Sea ice fraction
    real :: tice_am(ix,il) ! Sea ice temperature

    ! Sea sfc. fields from ocean/sea-ice model
    real :: sst_om(ix,il) ! Ocean model SST
    real :: sice_om(ix,il) ! Model sea ice fraction
    real :: tice_om(ix,il) ! Model sea ice temperature
    real :: ssti_om(ix,il) ! Model SST + sea ice temp.

    ! Weight for obs. SST anomaly in coupled runs
    real :: wsst_ob(ix,il)

    ! Fluxes at sea surface (all downward, except evaporation)
    real :: hflux_s(ix,il) ! Net heat flux into sea surface
    real :: hflux_i(ix,il) ! Net heat flux into sea-ice surface

    ! Flag for sea-surface temperature coupling
    ! 0 = precribed SST, no coupling
    ! 1 = precribed SST, ocean model forced by atmosphere
    ! 2 = full (uncorrected) SST from coupled ocean model
    ! 3 = SST anomaly from coupled ocean model + observed SST climatology
    ! 4 = as 3 with prescribed SST anomaly in ElNino region
    integer :: sea_coupling_flag  = 0

    ! Flag for sea-ice coupling
    integer :: ice_coupling_flag  = 1

    ! Flag for observed SST anomaly
    ! 0 = climatological SST
    ! 1 = observed anomaly
    ! (active if sea_coupling_flag = 0, 1; set to 1 if sea_coupling_flag = 4)
    integer :: sst_anomaly_coupling_flag = 1

contains
    ! Initialization of sea model
    subroutine sea_model_init
        use boundaries, only: fmask, fillsf, forchk
        use date, only: isst0
        use geometry, only: radang
        use input_output, only: load_boundary_file

        ! Domain mask
        real :: dmask(ix,il)

        ! Domain flags
        logical :: l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian

        ! Heat capacities of mixed-layer and sea-ice
        real :: hcaps(il)
        real :: hcapi(il)

        integer :: i, j, month
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

        ! Minimum fraction of sea for the definition of anomalies
        real :: fseamin = 1./3.

        ! Dissipation time (days) for sea-ice temp. anomalies
        real :: tdice = 30.0

        ! Threshold for land-sea mask definition (i.e. minimum fraction of
        ! either land or sea)
        real :: thrsh = 0.1

        ! Geographical domain
        ! note : more than one regional domain may be set .true.
        l_globe  =  .true.         ! global domain
        l_northe = .false.         ! Northern hem. oceans (lat > 20N)
        l_natlan = .false.         ! N. Atlantic (lat 20-80N, lon 100W-45E)
        l_npacif = .false.         ! N. Pacific  (lat 20-80N, lon 100E-100W)
        l_tropic = .false.         ! Tropics (lat 30S-30N)
        l_indian = .false.         ! Indian Ocean (lat 30S-30N, lon 30-120E)

        ! =========================================================================
        ! Initialize sea-surface boundary conditions
        ! =========================================================================

        ! Fractional and binary sea masks
        do j = 1, il
            do i = 1, ix
                fmask_s(i,j) = 1.0 - fmask(i,j)

                if (fmask_s(i,j) >= thrsh) then
                    bmask_s(i,j) = 1.0
                    if (fmask_s(i,j) > (1.0 - thrsh)) fmask_s(i,j) = 1.0
                else
                    bmask_s(i,j) = 0.0
                    fmask_s(i,j) = 0.0
                end if
            end do
        end do

        ! Grid latitudes for sea-surface variables
        deglat_s =  radang*90.0/asin(1.0)

        ! SST
        do month = 1, 12
            sst12(:,:,month) = load_boundary_file("sea_surface_temperature.nc", "sst", month)

            call fillsf(sst12(:,:,month), 0.0)
        end do

        call forchk(bmask_s, 12, 100.0, 400.0, 273.0, sst12)

        ! Sea ice concentration
        do month = 1, 12
            sice12(:,:,month) = max(load_boundary_file("sea_ice.nc", "icec", month), 0.0)
        end do

        call forchk(bmask_s, 12, 0.0, 1.0, 0.0, sice12)

        ! SST anomalies for initial and preceding/following months
        if (sst_anomaly_coupling_flag > 0) then
            write (*,'(A,I0.2)') 'SST anomalies are read starting from month ', isst0
            do month = 1, 3
                if ((isst0 <= 1 .and. month /= 2) .or. isst0 > 1) then
                    sstan3(:,:,month) = load_boundary_file("sea_surface_temperature_anomaly.nc", &
                        & "ssta", isst0-2+month, 420)
                end if
            end do

            call forchk(bmask_s, 3, -50.0, 50.0, 0.0, sstan3)
        end if

        ! Climatological fields for the ocean model (TO BE RECODED)
        ! Annual-mean heat flux into sea-surface
        hfseacl = 0.0

        if (sea_coupling_flag >= 1) then
            stop "Model behaviour when sea_coupling_flag >= 1 not implemented yet"
        end if

        ! Ocean model SST climatology:
        ! defined by adding SST model bias to observed climatology
        ! (bias may be defined in a different period from climatology)

        if (sea_coupling_flag >= 3) then
            stop "Model behaviour when sea_coupling_flag >= 3 not implemented yet"
        end if

        ! =========================================================================
        ! Compute heat capacities
        ! =========================================================================

        ! Heat flux coefficient at sea/ice interface [(W/m^2)/deg]
        beta = 1.

        ! Heat capacities per m^2 (depth*heat_cap/m^3)
        crad=asin(1.)/90.
        do j=1,il
            coslat   = cos(crad*deglat_s(j))
            hcaps(j) = 4.18e+6*(depth_ml +(dept0_ml -depth_ml) *coslat**3)
            hcapi(j) = 1.93e+6*(depth_ice+(dept0_ice-depth_ice)*coslat**2)
        end do

        ! =========================================================================
        ! Compute constant parameters and fields
        ! =========================================================================

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
        use date, only: model_datetime, imont1
        use mod_cpl_bcinterp, only: forin5, forint

        integer, intent(in) :: day

        integer :: i, j
        real :: sstcl0, sstfr

        ! 1. Interpolate climatological fields and obs. SST anomaly
        !    to actual date

        ! Climatological SST
        call forin5(imont1,sst12,sstcl_ob)

        ! Climatological sea ice fraction
        call forint(imont1,sice12,sicecl_ob)

        ! SST anomaly
        if (sst_anomaly_coupling_flag.gt.0) then
            if (model_datetime%day.eq.1.and.day.gt.0) call obs_ssta
            call forint (2,sstan3,sstan_ob)
        end if

        ! Ocean model climatological SST
        if (sea_coupling_flag.ge.3) then
            call forin5 (imont1,sstom12,sstcl_om)
        end if

        ! Adjust climatological fields over sea ice

        ! SST at freezing point
        sstfr = 273.2-1.8

        do i = 1, ix
            do j = 1, il
                sstcl0 = sstcl_ob(i,j)

                if (sstcl_ob(i,j) > sstfr) then
                    sicecl_ob(i,j) = min(0.5, sicecl_ob(i,j))
                    ticecl_ob(i,j) = sstfr
                    if (sicecl_ob(i,j) .gt. 0.0) then
                        sstcl_ob(i,j) = sstfr + (sstcl_ob(i,j) - sstfr)/(1.0 - sicecl_ob(i,j))
                    end if
                else
                    sicecl_ob(i,j) = max(0.5, sicecl_ob(i,j))
                    ticecl_ob(i,j) = sstfr + (sstcl_ob(i,j) - sstfr)/sicecl_ob(i,j)
                    sstcl_ob(i,j)  = sstfr
                end if

                if (sea_coupling_flag >= 3) sstcl_om(i,j) = sstcl_om(i,j) + (sstcl_ob(i,j) - sstcl0)
            end do
        end do

        if (day == 0) then
            ! 2. Initialize prognostic variables of ocean/ice model
            !    in case of no restart or no coupling
            sst_om  = sstcl_ob      ! SST
            tice_om = ticecl_ob     ! sea ice temperature
            sice_om = sicecl_ob     ! sea ice fraction

            if (sea_coupling_flag <= 0) sst_om = 0.0

            ! 3. Compute additional sea/ice variables
            wsst_ob = 0.
            if (sea_coupling_flag >= 4) call sea_domain('elnino',wsst_ob)
        else
            if (sea_coupling_flag > 0 .or. ice_coupling_flag > 0) then
                ! 1. Run ocean mixed layer or
                !    call message-passing routines to receive data from ocean model
                call run_sea_model
            end if
        end if

        ! 3. Compute sea-sfc. anomalies and full fields for atm. model
        ! 3.1 SST
        sstan_am = 0.0

        if (sea_coupling_flag <= 1) then
            if (sst_anomaly_coupling_flag > 0) sstan_am = sstan_ob

            ! Use observed SST (climatological or full field)
            sst_am = sstcl_ob + sstan_am
        else if (sea_coupling_flag.eq.2) then
            ! Use full ocean model SST
            sst_am = sst_om
        else if (sea_coupling_flag >= 3) then
            ! Define SST anomaly from ocean model ouput and climatology
            sstan_am = sst_om - sstcl_om

            ! Merge with observed SST anomaly in selected area
            if (sea_coupling_flag >= 4) then
                sstan_am = sstan_am + wsst_ob*(sstan_ob - sstan_am)
            end if

            ! Add observed SST climatology to model SST anomaly
            sst_am = sstcl_ob + sstan_am
        end if

        ! 3.2 Sea ice fraction and temperature
        if (ice_coupling_flag > 0) then
            sice_am = sice_om
            tice_am = tice_om
        else
            sice_am = sicecl_ob
            tice_am = ticecl_ob
        end if

        sst_am  = sst_am + sice_am*(tice_am - sst_am)
        ssti_om = sst_om + sice_am*(tice_am - sst_om)
    end subroutine

    ! Update observed SST anomaly array
    subroutine obs_ssta
        use date, only: model_datetime, start_datetime
        use input_output, only: load_boundary_file
        use boundaries, only: forchk

        integer :: i, j, next_month

        sstan3(:,:,1) = sstan3(:,:,2)
        sstan3(:,:,2) = sstan3(:,:,3)

        ! Compute next month given initial SST year
        next_month = (start_datetime%year - issty0) * 12 + model_datetime%month

        ! Read next month SST anomalies
        sstan3(:,:,3) = load_boundary_file("sea_surface_temperature_anomaly.nc", "ssta", &
            & next_month, 420)

        call forchk(bmask_s, 1, -50.0, 50.0, 0.0, sstan3(:,:,3))
    end

    ! Purpose : Integrate slab ocean and sea-ice models for one day
    subroutine run_sea_model
        ! Auxiliary variables
        real :: hflux(ix,il)   ! net sfc. heat flux
        real :: tanom(ix,il)   ! sfc. temperature anomaly
        real :: cdis(ix,il)    ! dissipation ceofficient

        real :: anom0, sstfr

        sstfr = 273.2-1.8       ! SST at freezing point

        ! 1. Ocean mixed layer
        ! Net heat flux
        hflux = hflux_s-hfseacl-sicecl_ob*(hflux_i+beta*(sstfr-tice_om))

        ! Anomaly at t0 minus climatological temp. tendency
        tanom = sst_om - sstcl_ob

        ! Time evoloution of temp. anomaly
        tanom = cdsea*(tanom+rhcaps*hflux)

        ! Full SST at final time
        sst_om = tanom + sstcl_ob

        ! 2. Sea-ice slab model

        ! Net heat flux
        hflux = hflux_i + beta*(sstfr-tice_om)

        ! Anomaly w.r.t final-time climatological temp.
        tanom = tice_om - ticecl_ob

        ! Definition of non-linear damping coefficient
        anom0     = 20.
        cdis = cdice*(anom0/(anom0+abs(tanom)))
        !cdis(:,:) = cdice(:,:)

        ! Time evoloution of temp. anomaly
        tanom = cdis*(tanom+rhcapi*hflux)

        ! Full ice temperature at final time
        tice_om = tanom + ticecl_ob

        ! Persistence of sea ice fraction
        sice_om = sicecl_ob
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
