subroutine sea_model_init(fmask_s,rlat)
    !  subroutine sea_model_init (fmask_s,rlat)
    !
    !  Purpose : Initialization of sea model
    !  Initialized common blocks: sea_mc

    use mod_atparam
    use mod_cplcon_sea

    implicit none
							
    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    ! Input variables
    real fmask_s(nlon,nlat)            ! sea mask (fraction of sea)
    real rlat(nlat)                    ! latitudes in degrees

    ! Auxiliary variables

    ! Domain mask
    real :: dmask(nlon,nlat)              

    ! Domain flags
    logical :: l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian            

    ! Heat capacity of mixed-l
    real :: hcaps(nlat)

    ! Heat capacity of sea-ice
    real :: hcapi(nlat)                      

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
    do j=1,nlat
        coslat   = cos(crad*rlat(j))
        hcaps(j) = 4.18e+6*(depth_ml +(dept0_ml -depth_ml) *coslat**3)
        hcapi(j) = 1.93e+6*(depth_ice+(dept0_ice-depth_ice)*coslat**2)
    end do

    ! 3. Compute constant parameters and fields

    ! Set domain mask
    if (l_globe) then
        dmask(:,:) = 1.
    else
        dmask(:,:) = 0.
        if (l_northe) call SEA_DOMAIN ('northe',rlat,dmask)
        !fkif (l_arctic) call SEA_DOMAIN ('arctic',rlat,dmask)
        if (l_natlan) call SEA_DOMAIN ('natlan',rlat,dmask)
        if (l_npacif) call SEA_DOMAIN ('npacif',rlat,dmask)
        if (l_tropic) call SEA_DOMAIN ('tropic',rlat,dmask)
        if (l_indian) call SEA_DOMAIN ('indian',rlat,dmask)
    end if

    ! Smooth latitudinal boundaries and blank out land points
    do j=2,nlat-1
        rhcaps(:,j) = 0.25*(dmask(:,j-1)+2*dmask(:,j)+dmask(:,j+1))
    end do
    dmask(:,2:nlat-1) = rhcaps(:,2:nlat-1)

    do j=1,nlat
        do i=1,nlon
            if (fmask_s(i,j).lt.fseamin) dmask(i,j) = 0
        end do
    end do

    ! Set heat capacity and dissipation time over selected domain
    do j=1,nlat
        rhcaps(:,j) = 86400./hcaps(j)
        rhcapi(:,j) = 86400./hcapi(j)
    end do

    cdsea = dmask*tdsst/(1.+dmask*tdsst)
    cdice = dmask*tdice/(1.+dmask*tdice)
end

subroutine sea_model 
    ! subroutine sea_model

    ! Purpose : Integrate slab ocean and sea-ice models for one day

    use mod_atparam
    use mod_cplcon_sea
    use mod_cplvar_sea

    implicit none

    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    !real vsea_input(nlon,nlat,8), vsea_output(nlon,nlat,3)

    ! Input variables:
    real ::  sst0(nlon,nlat)     ! SST at initial time
    real :: tice0(nlon,nlat)     ! sea ice temp. at initial time
    real :: sice0(nlon,nlat)     ! sea ice fraction at initial time
    real :: hfsea(nlon,nlat)     ! sea+ice  sfc. heat flux between t0 and t1
    real :: hfice(nlon,nlat)     ! ice-only sfc. heat flux between t0 and t1

    real ::  sstcl1(nlon,nlat)   ! clim. SST at final time 
    real :: ticecl1(nlon,nlat)   ! clim. sea ice temp. at final time
    real :: hfseacl(nlon,nlat)   ! clim. heat flux due to advection/upwelling

    ! Output variables
    real ::  sst1(nlon,nlat)     ! SST at final time
    real :: tice1(nlon,nlat)     ! sea ice temp. at final time 
    real :: sice1(nlon,nlat)     ! sea ice fraction at final time 

    ! Auxiliary variables
    real :: hflux(nlon,nlat)   ! net sfc. heat flux
    real :: tanom(nlon,nlat)   ! sfc. temperature anomaly
    real :: cdis(nlon,nlat)    ! dissipation ceofficient

    real :: anom0, sstfr

    sst0 = reshape(vsea_input(:,1), (/nlon, nlat/))
    tice0 = reshape(vsea_input(:,2), (/nlon, nlat/))
    sice0 = reshape(vsea_input(:,3), (/nlon, nlat/))
    hfsea = reshape(vsea_input(:,4), (/nlon, nlat/))
    hfice = reshape(vsea_input(:,5), (/nlon, nlat/))
    sstcl1 = reshape(vsea_input(:,6), (/nlon, nlat/))
    ticecl1 = reshape(vsea_input(:,7), (/nlon, nlat/))
    hfseacl = reshape(vsea_input(:,8), (/nlon, nlat/))

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

    vsea_output(:,1) = reshape(sst1, (/ngp/))
    vsea_output(:,2) = reshape(tice1, (/ngp/))
    vsea_output(:,3) = reshape(sice1, (/ngp/))
end

subroutine sea_domain(cdomain,rlat,dmask) 
    ! subroutine sea_domain (cdomain,rlat,dmask)

    ! Purpose : Definition of ocean domains

    use mod_atparam

    implicit none

    integer, parameter :: nlon=ix, nlat=il

    ! Input variables

    character(len=6), intent(in) :: cdomain           ! domain name
    real, intent(in) :: rlat(nlat)               ! latitudes in degrees

    ! Output variables (initialized by calling routine) 
    real, intent(inout) :: dmask(nlon,nlat)         ! domain mask 

    integer :: i, j
    real :: arlat, dlon, rlon, rlonw, wlat

    print *, 'sea domain : ', cdomain

    dlon = 360./float(nlon)
                                   
    if (cdomain.eq.'northe') then
        do j=1,nlat
            if (rlat(j).gt.20.0) dmask(:,j) = 1.
        end do
    end if

    if (cdomain.eq.'natlan') then
         do j=1,nlat
           if (rlat(j).gt.20.0.and.rlat(j).lt.80.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.lt.45.0.or.rlon.gt.260.0) dmask(i,j) = 1.
             end do
           end if
         end do
    end if

    if (cdomain.eq.'npacif') then
        do j=1,nlat
            if (rlat(j).gt.20.0.and.rlat(j).lt.65.0) then
                do i=1,nlon
                    rlon = (i-1)*dlon
                    if (rlon.gt.120.0.and.rlon.lt.260.0) dmask(i,j) = 1.
                end do
            end if
        end do
    end if

    if (cdomain.eq.'tropic') then
        do j=1,nlat
            if (rlat(j).gt.-30.0.and.rlat(j).lt.30.0) dmask(:,j) = 1.
        end do
    end if

    if (cdomain.eq.'indian') then
        do j=1,nlat
            if (rlat(j).gt.-30.0.and.rlat(j).lt.30.0) then
                do i=1,nlon
                    rlon = (i-1)*dlon
                    if (rlon.gt.30.0.and.rlon.lt.120.0) dmask(i,j) = 1.
                end do
            end if
        end do
    end if

    if (cdomain.eq.'elnino') then
        do j=1,nlat
            arlat = abs(rlat(j))
            if (arlat.lt.25.0) then
                wlat = 1.
                if (arlat.gt.15.0) wlat = (0.1*(25.-arlat))**2
                rlonw = 300.-2*max(rlat(j),0.)
                do i=1,nlon
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

    !do j=1,nlat
    !    print *, 'lat = ',rlat(j),' sea model domain  = ',dmask(nlon/2,j)
    !end do
end
