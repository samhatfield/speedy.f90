

      SUBROUTINE SEA_MODEL_INIT (fmask_s,rlat) 
C--
C--   SUBROUTINE SEA_MODEL_INIT (fmask_s,rlat)
C--
C--   Purpose : Initialization of sea model
C--   Initialized common blocks: SEA_MC
C--	
							
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

C     Input variables

      real fmask_s(nlon,nlat)            ! sea mask (fraction of sea)
      real rlat(nlat)                    ! latitudes in degrees

      include "com_cplcon_sea.h"

C     Auxiliary variables

      real    dmask(nlon,nlat)              ! domain mask
      logical l_globe,
     &        l_northe, l_natlan, l_npacif, 
     &        l_tropic, l_indian            ! domain flags

      real hcaps(nlat)                      ! heat capacity of mixed-l.
      real hcapi(nlat)                      ! heat capacity of sea-ice

C--  
C--   1. Set geographical domain, heat capacities and dissipation times
C--      for sea (mixed layer) and sea-ice 

C     Model parameters (default values)

C     ocean mixed layer depth: d + (d0-d)*(cos_lat)^3
      depth_ml = 60.               ! High-latitude depth
      dept0_ml = 40.               ! Minimum depth (tropics)

C     sea-ice depth : d + (d0-d)*(cos_lat)^2
      depth_ice = 2.5              ! High-latitude depth
      dept0_ice = 1.5              ! Minimum depth 

C     Dissipation time (days) for sea-surface temp. anomalies
      tdsst  = 90.

C     Dissipation time (days) for sea-ice temp. anomalies
      tdice = 30.

C     Heat flux coefficient at sea/ice interface [(W/m^2)/deg]
      beta = 1.

C     Minimum fraction of sea for the definition of anomalies
      fseamin = 1./3.

C     Geographical domain
C     note : more than one regional domain may be set .true.

      l_globe  =  .true.         ! global domain
      l_northe = .false.         ! Northern hem. oceans (lat > 20N)
      l_natlan = .false.         ! N. Atlantic (lat 20-80N, lon 100W-45E)
      l_npacif = .false.         ! N. Pacific  (lat 20-80N, lon 100E-100W)
      l_tropic = .false.         ! Tropics (lat 30S-30N)
      l_indian = .false.         ! Indian Ocean (lat 30S-30N, lon 30-120E)

C     Reset model parameters

      include "cls_insea.h"

C     Heat capacities per m^2 (depth*heat_cap/m^3)

      crad=asin(1.)/90.
      do j=1,nlat
        coslat   = cos(crad*rlat(j))
        hcaps(j) = 4.18e+6*(depth_ml +(dept0_ml -depth_ml) *coslat**3)
        hcapi(j) = 1.93e+6*(depth_ice+(dept0_ice-depth_ice)*coslat**2)
      enddo
C--
C--   3. Compute constant parameters and fields

C     Set domain mask

      if (l_globe) then
        dmask(:,:) = 1.
      else
        dmask(:,:) = 0.
        if (l_northe) call SEA_DOMAIN ('northe',rlat,dmask)
cfk        if (l_arctic) call SEA_DOMAIN ('arctic',rlat,dmask)
        if (l_natlan) call SEA_DOMAIN ('natlan',rlat,dmask)
        if (l_npacif) call SEA_DOMAIN ('npacif',rlat,dmask)
        if (l_tropic) call SEA_DOMAIN ('tropic',rlat,dmask)
        if (l_indian) call SEA_DOMAIN ('indian',rlat,dmask)
      endif

C     Smooth latitudinal boundaries and blank out land points

      do j=2,nlat-1
         rhcaps(:,j) = 0.25*(dmask(:,j-1)+2*dmask(:,j)+dmask(:,j+1))
      enddo
      dmask(:,2:nlat-1) = rhcaps(:,2:nlat-1)

      do j=1,nlat
        do i=1,nlon
           if (fmask_s(i,j).lt.fseamin) dmask(i,j) = 0
        enddo
      enddo

C     Set heat capacity and dissipation time over selected domain

      do j=1,nlat
         rhcaps(:,j) = 86400./hcaps(j)
         rhcapi(:,j) = 86400./hcapi(j)
      enddo

      cdsea(:,:) = dmask(:,:)*tdsst/(1.+dmask(:,:)*tdsst)
      cdice(:,:) = dmask(:,:)*tdice/(1.+dmask(:,:)*tdice)

      return
      end

      SUBROUTINE SEA_MODEL 
C--
C--   SUBROUTINE SEA_MODEL
C--
C--   Purpose : Integrate slab ocean and sea-ice models for one day
							
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

c      real vsea_input(nlon,nlat,8), vsea_output(nlon,nlat,3)

      include "com_cplcon_sea.h"

      include "com_cplvar_sea.h"

C     Input variables:

      real  sst0(nlon,nlat)     ! SST at initial time
      real tice0(nlon,nlat)     ! sea ice temp. at initial time
      real sice0(nlon,nlat)     ! sea ice fraction at initial time
      real hfsea(nlon,nlat)     ! sea+ice  sfc. heat flux between t0 and t1
      real hfice(nlon,nlat)     ! ice-only sfc. heat flux between t0 and t1

      real  sstcl1(nlon,nlat)   ! clim. SST at final time 
      real ticecl1(nlon,nlat)   ! clim. sea ice temp. at final time
      real hfseacl(nlon,nlat)   ! clim. heat flux due to advection/upwelling

      equivalence    (sst0,vsea_input(1,1))
      equivalence   (tice0,vsea_input(1,2))
      equivalence   (sice0,vsea_input(1,3))
      equivalence   (hfsea,vsea_input(1,4))
      equivalence   (hfice,vsea_input(1,5))
      equivalence  (sstcl1,vsea_input(1,6))
      equivalence (ticecl1,vsea_input(1,7))
      equivalence (hfseacl,vsea_input(1,8))

C     Output variables
 
      real  sst1(nlon,nlat)     ! SST at final time
      real tice1(nlon,nlat)     ! sea ice temp. at final time 
      real sice1(nlon,nlat)     ! sea ice fraction at final time 

      equivalence  (sst1,vsea_output(1,1))
      equivalence (tice1,vsea_output(1,2))
      equivalence (sice1,vsea_output(1,3))

C     Auxiliary variables

      real hflux(nlon,nlat)   ! net sfc. heat flux
      real tanom(nlon,nlat)   ! sfc. temperature anomaly
      real  cdis(nlon,nlat)   ! dissipation ceofficient

      sstfr = 273.2-1.8       ! SST at freezing point

c      beta = 1.               ! heat flux coef. at sea-ice bottom
C--
C--   1. Ocean mixed layer

C     Net heat flux
      hflux(:,:) = hfsea(:,:)-hfseacl(:,:)
     &             -sice0(:,:)*(hfice(:,:)+beta*(sstfr-tice0(:,:)))

C     Anomaly at t0 minus climatological temp. tendency
      tanom(:,:) = sst0(:,:) - sstcl1(:,:)

C     Time evoloution of temp. anomaly 
      tanom(:,:) = cdsea(:,:)*(tanom(:,:)+rhcaps(:,:)*hflux(:,:))

C     Full SST at final time
      sst1(:,:) = tanom(:,:) + sstcl1(:,:)

C--
C--   2. Sea-ice slab model

C     Net heat flux
      hflux(:,:) = hfice(:,:)
     &           + beta*(sstfr-tice0(:,:))

C     Anomaly w.r.t final-time climatological temp.
      tanom(:,:) = tice0(:,:) - ticecl1(:,:)

C     Definition of non-linear damping coefficient
      anom0     = 20. 
      cdis(:,:) = cdice(:,:)*(anom0/(anom0+abs(tanom(:,:))))
c     cdis(:,:) = cdice(:,:)

C     Time evoloution of temp. anomaly 
      tanom(:,:) = cdis(:,:)*(tanom(:,:)+rhcapi(:,:)*hflux(:,:))

C     Full ice temperature at final time
      tice1(:,:) = tanom(:,:) + ticecl1(:,:)

C     Persistence of sea ice fraction
      sice1(:,:) = sice0(:,:)

      return
      end

      SUBROUTINE SEA_DOMAIN (cdomain,rlat,dmask) 
C--
C--   SUBROUTINE SEA_DOMAIN (cdomain,rlat,dmask)
C--
C--   Purpose : Definition of ocean domains
C--	
							
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

C     Input variables

      character*6 cdomain           ! domain name
      real rlat(nlat)               ! latitudes in degrees

C     Output variables (initialized by calling routine) 

      real dmask(nlon,nlat)         ! domain mask 

      print *, 'sea domain : ', cdomain

      dlon = 360./float(nlon)
                                   
      if (cdomain.eq.'northe') then
        do j=1,nlat
          if (rlat(j).gt.20.0) dmask(:,j) = 1.
        enddo
      endif

      if (cdomain.eq.'natlan') then
         do j=1,nlat
           if (rlat(j).gt.20.0.and.rlat(j).lt.80.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.lt.45.0.or.rlon.gt.260.0) dmask(i,j) = 1.
             enddo
           endif
         enddo
      endif

      if (cdomain.eq.'npacif') then
         do j=1,nlat
           if (rlat(j).gt.20.0.and.rlat(j).lt.65.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.gt.120.0.and.rlon.lt.260.0) dmask(i,j) = 1.
             enddo
           endif
         enddo
      endif

      if (cdomain.eq.'tropic') then
        do j=1,nlat
          if (rlat(j).gt.-30.0.and.rlat(j).lt.30.0) dmask(:,j) = 1.
        enddo
      endif

      if (cdomain.eq.'indian') then
         do j=1,nlat
           if (rlat(j).gt.-30.0.and.rlat(j).lt.30.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.gt.30.0.and.rlon.lt.120.0) dmask(i,j) = 1.
             enddo
           endif
         enddo
      endif

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
               endif
             enddo
           endif
         enddo
      endif

c      do j=1,nlat
c       print *, 'lat = ',rlat(j),' sea model domain  = ',dmask(nlon/2,j)
c      enddo

      return
      end

