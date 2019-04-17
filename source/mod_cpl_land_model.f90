module mod_cpl_land_model
    use mod_atparam

    implicit none

    private
    public rhcapl, cdland, vland_input, vland_output
    public stlcl_ob, snowdcl_ob, soilwcl_ob, stl_am, snowd_am, soilw_am, stl_lm
    public land_model_init, ini_land, land_model
    public atm2land, land2atm

    ! 1./heat_capacity (land)
    real :: rhcapl(ix,il)           

    ! 1./dissip_time (land)
    real :: cdland(ix,il)           

    ! Input and output land variables exchanged by coupler
    ! Land model input variables
    real :: vland_input(ix*il,4)            

    ! Land model output variables
    real :: vland_output(ix*il,2)           

    ! Daily observed climatological fields over land
    real :: stlcl_ob(ix*il)              ! clim. land sfc. temperature
    real :: snowdcl_ob(ix*il)              ! clim. snow depth (water equiv)
    real :: soilwcl_ob(ix*il)              ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    real :: stl_am(ix*il)                   ! land sfc. temperature
    real :: snowd_am(ix*il)                 ! snow depth (water equiv)
    real :: soilw_am(ix*il)                 ! soil water availability

    ! Land sfc. fields from land model
    real :: stl_lm(ix*il)                 ! land-model sfc. temperature

    contains
        subroutine land_model_init(fmask_l,alb0) 
            ! subroutine land_model_init (fmask_l,alb0)
            !
            ! purpose : initialization of land model
            ! initialized common blocks: land_mc
            
            ! Input variables
            ! Land mask (fraction of land)
            real, intent(in) :: fmask_l(ix,il)            
            ! Annual-mean albedo
            real, intent(in) :: alb0(ix,il)            
        
            ! Auxiliary variables
            integer :: i, j
            real :: dmask(ix,il)           ! domain mask
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
        
            do j=1,il
                do i=1,ix
                    if (fmask_l(i,j).lt.flandmin) dmask(i,j) = 0
                end do
            end do
        
            ! Set time_step/heat_capacity and dissipation fields
            do j=1,il
                do i=1,ix
                    if (alb0(i,j).lt.0.4) then
                        rhcapl(i,j) = 86400./hcapl
                    else
                        rhcapl(i,j) = 86400./hcapli
                    endif
                end do
            end do
        
            cdland(:,:) = dmask(:,:)*tdland/(1.+dmask(:,:)*tdland)
        end

        subroutine ini_land()
            ! 1. Compute climatological fields for initial date
            call atm2land(0)

            ! 2. Initialize prognostic variables of land model
            stl_lm(:)  = stlcl_ob(:)      ! land sfc. temperature

            ! 3. Compute additional land variables
            call land2atm(0)
        end

        subroutine atm2land(jday)
            use mod_cpl_flags, only: icland
            use mod_flx_land, only: hflux_l
            use mod_cli_land, only: stl12, snowd12, soilw12
            use mod_date, only: imont1, tmonth

            integer, intent(in) :: jday
            integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

            ! 1. Interpolate climatological fields to actual date

            ! Climatological land sfc. temperature
            call forin5(ngp,imont1,tmonth,stl12,stlcl_ob)

            ! Climatological snow depth
            call forint(ngp,imont1,tmonth,snowd12,snowdcl_ob)

            ! Climatological soil water availability
            call forint(ngp,imont1,tmonth,soilw12,soilwcl_ob)

            if (jday.le.0) return

            ! 2. Set input variables for mixed-layer/ocean model
            if (icland.gt.0) then
                vland_input(:,1) = stl_lm(:)
                vland_input(:,2) = hflux_l(:)
                vland_input(:,3) = stlcl_ob(:)
            end if

            ! 3. Call message-passing routines to send data (if needed)
        end

        subroutine land2atm(jday)
            use mod_cpl_flags, only: icland

            integer, intent(in) :: jday

            if (jday.gt.0.and.icland.gt.0) then
                ! 1. Run ocean mixed layer or
                !    call message-passing routines to receive data from ocean model
                call land_model

                ! 2. Get updated variables for mixed-layer/ocean model
                stl_lm(:) = vland_output(:,1)      ! land sfc. temperature
            end if

            ! 3. Compute land-sfc. fields for atm. model
            ! 3.1 Land sfc. temperature
            if (icland.le.0) then
                ! Use observed climatological field
                stl_am(:) = stlcl_ob(:)
            else
                ! Use land model sfc. temperature
                stl_am(:) = stl_lm(:)
            end if

            ! 3.2 Snow depth and soil water availability
            snowd_am(:) = snowdcl_ob(:)
            soilw_am(:) = soilwcl_ob(:)
        end

        subroutine land_model
            ! subroutine land_model
            !
            ! purpose : integrate slab land-surface model for one day
        							
            !real vland_input(ix,il,3), vland_output(ix,il,2)
        
            ! Input variables:
            real :: stl0(ix*il)    ! land temp. at initial time
            real :: hfland(ix*il)    ! land sfc. heat flux between t0 and t1
            real :: stlcl1(ix*il)    ! clim. land temp. at final time 
        
            ! Output variables
            real :: stl1(ix*il)     ! land temp. at final time
        
            ! Auxiliary variables
            real :: hflux(ix*il)   ! net sfc. heat flux
            real :: tanom(ix*il)   ! sfc. temperature anomaly

            ! Initialise variables
            stl0 = vland_input(:,1)
            hfland = vland_input(:,2)
            stlcl1 = vland_input(:,3)
 
            ! 1. Land-surface (soil/ice-sheet) layer
 
            ! Net heat flux
            ! (snow correction to be added?)
            hflux = hfland
 
            ! Anomaly w.r.t final-time climatological temp.
            tanom = stl0 - stlcl1
        
            ! Time evoloution of temp. anomaly 
            tanom = reshape(cdland, (/ ix*il /))*&
                & (tanom+reshape(rhcapl, (/ ix*il /))*hflux)
 
            ! Full SST at final time
            stl1 = tanom + stlcl1

            vland_output(:,1) = stl1
        end
end
