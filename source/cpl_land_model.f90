subroutine land_model_init(fmask_l,alb0) 
    ! subroutine land_model_init (fmask_l,alb0)
    !
    ! purpose : initialization of land model
    ! initialized common blocks: land_mc
    
    use mod_atparam

    implicit none
							
    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    ! Input variables
    real, intent(in) :: fmask_l(nlon,nlat)            ! land mask (fraction of land)
    real, intent(in) :: alb0(nlon,nlat)            ! annual-mean albedo

    include "com_cplcon_land.h"

    ! Auxiliary variables
    integer :: i, j
    real dmask(nlon,nlat)           ! domain mask
    real :: depth_soil, depth_lice, tdland, hcapl, hcapli, flandmin

    ! 1. Set heat capacities and dissipation times for 
    !    soil and ice-sheet layers 

    ! Model parameters (default values)

    ! Soil layer depth (m)
    depth_soil = 1.0

    ! Land-ice depth (m)
    depth_lice = 5.0

    ! Dissipation time (days) for land-surface temp. anomalies
    tdland  = 40.

    ! Minimum fraction of land for the definition of anomalies
    flandmin = 1./3.

    ! Reset model parameters
    include "cls_inland.h"

    ! Heat capacities per m^2 (depth*heat_cap/m^3)
    hcapl  = depth_soil*2.50e+6
    hcapli = depth_lice*1.93e+6

    ! 2. Compute constant fields
    ! Set domain mask (blank out sea points)
    dmask(:,:) = 1.

    do j=1,nlat
        do i=1,nlon
            if (fmask_l(i,j).lt.flandmin) dmask(i,j) = 0
        end do
    end do

    ! Set time_step/heat_capacity and dissipation fields
    do j=1,nlat
        do i=1,nlon
            if (alb0(i,j).lt.0.4) then
                rhcapl(i,j) = 86400./hcapl
            else
                rhcapl(i,j) = 86400./hcapli
            endif
        end do
    end do

    cdland(:,:) = dmask(:,:)*tdland/(1.+dmask(:,:)*tdland)
end

subroutine land_model 
    ! subroutine land_model
    !
    ! purpose : integrate slab land-surface model for one day
							
    use mod_atparam

    implicit none

    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    !real vland_input(nlon,nlat,3), vland_output(nlon,nlat,2)
    include "com_cplcon_land.h"
    include "com_cplvar_land.h"

    ! Input variables:
    real :: stl0(nlon,nlat)    ! land temp. at initial time
    real :: hfland(nlon,nlat)    ! land sfc. heat flux between t0 and t1
    real :: stlcl1(nlon,nlat)    ! clim. land temp. at final time 

    equivalence   (stl0,vland_input(1,1))
    equivalence (hfland,vland_input(1,2))
    equivalence (stlcl1,vland_input(1,3))

    ! Output variables
    real :: stl1(nlon,nlat)     ! land temp. at final time

    equivalence (stl1,vland_output(1,1))

    ! Auxiliary variables
    real :: hflux(nlon,nlat)   ! net sfc. heat flux
    real :: tanom(nlon,nlat)   ! sfc. temperature anomaly

    ! 1. Land-surface (soil/ice-sheet) layer

    ! Net heat flux
    ! (snow correction to be added?)
    hflux(:,:) = hfland(:,:)

    ! Anomaly w.r.t final-time climatological temp.
    tanom(:,:) = stl0(:,:) - stlcl1(:,:)

    ! Time evoloution of temp. anomaly 
    tanom(:,:) = cdland(:,:)*(tanom(:,:)+rhcapl(:,:)*hflux(:,:))

    ! Full SST at final time
    stl1(:,:) = tanom(:,:) + stlcl1(:,:)
end
