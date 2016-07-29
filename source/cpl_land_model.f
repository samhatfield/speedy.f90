

      SUBROUTINE LAND_MODEL_INIT (fmask_l,alb0) 
C--
C--   SUBROUTINE LAND_MODEL_INIT (fmask_l,alb0)
C--
C--   Purpose : Initialization of land model
C--   Initialized common blocks: LAND_MC
C--	
							
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

C     Input variables

      real fmask_l(nlon,nlat)            ! land mask (fraction of land)
      real    alb0(nlon,nlat)            ! annual-mean albedo

      include "com_cplcon_land.h"

C     Auxiliary variables

      real    dmask(nlon,nlat)           ! domain mask

C--  
C--   1. Set heat capacities and dissipation times for 
C--      soil and ice-sheet layers 

C     Model parameters (default values)

C     soil layer depth (m)
      depth_soil = 1.0

C     land-ice depth (m)
      depth_lice = 5.0

C     Dissipation time (days) for land-surface temp. anomalies
      tdland  = 40.

C     Minimum fraction of land for the definition of anomalies
      flandmin = 1./3.

C     Reset model parameters

      include "cls_inland.h"

C     Heat capacities per m^2 (depth*heat_cap/m^3)
      hcapl  = depth_soil*2.50e+6
      hcapli = depth_lice*1.93e+6

C--
C--   2. Compute constant fields

C     Set domain mask (blank out sea points)

      dmask(:,:) = 1.

      do j=1,nlat
        do i=1,nlon
           if (fmask_l(i,j).lt.flandmin) dmask(i,j) = 0
        enddo
      enddo

C     Set time_step/heat_capacity and dissipation fields

      do j=1,nlat
        do i=1,nlon
           if (alb0(i,j).lt.0.4) then
              rhcapl(i,j) = 86400./hcapl
           else
              rhcapl(i,j) = 86400./hcapli
           endif
        enddo
      enddo

      cdland(:,:) = dmask(:,:)*tdland/(1.+dmask(:,:)*tdland)

      return
      end

      SUBROUTINE LAND_MODEL 
C--
C--   SUBROUTINE LAND_MODEL
C--
C--   Purpose : Integrate slab land-surface model for one day
							
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

c      real vland_input(nlon,nlat,3), vland_output(nlon,nlat,2)

      include "com_cplcon_land.h"

      include "com_cplvar_land.h"

C     Input variables:

      real   stl0(nlon,nlat)    ! land temp. at initial time
      real hfland(nlon,nlat)    ! land sfc. heat flux between t0 and t1
      real stlcl1(nlon,nlat)    ! clim. land temp. at final time 

      equivalence    (stl0,vland_input(1,1))
      equivalence  (hfland,vland_input(1,2))
      equivalence  (stlcl1,vland_input(1,3))

C     Output variables
 
      real  stl1(nlon,nlat)     ! land temp. at final time

      equivalence  (stl1,vland_output(1,1))

C     Auxiliary variables

      real hflux(nlon,nlat)   ! net sfc. heat flux
      real tanom(nlon,nlat)   ! sfc. temperature anomaly

C--
C--   1. Land-surface (soil/ice-sheet) layer

C     Net heat flux
c     (snow correction to be added?)
      hflux(:,:) = hfland(:,:)

C     Anomaly w.r.t final-time climatological temp.
      tanom(:,:) = stl0(:,:) - stlcl1(:,:)

C     Time evoloution of temp. anomaly 
      tanom(:,:) = cdland(:,:)*(tanom(:,:)+rhcapl(:,:)*hflux(:,:))

C     Full SST at final time
      stl1(:,:) = tanom(:,:) + stlcl1(:,:)

      return
      end

