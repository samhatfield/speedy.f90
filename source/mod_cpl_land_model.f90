module mod_cpl_land_model
    use mod_atparam

    implicit none

    private
    public stlcl_ob, stl_am, snowd_am, soilw_am
    public land_model_init, ini_land, land_model
    public atm2land, land2atm
    public fmask_l, bmask_l, stl12, snowd12, soilw12

    ! 1./heat_capacity (land)
    real :: rhcapl(ix,il)

    ! 1./dissip_time (land)
    real :: cdland(ix,il)

    ! Daily observed climatological fields over land
    real :: stlcl_ob(ix*il)                ! clim. land sfc. temperature
    real :: snowdcl_ob(ix*il)              ! clim. snow depth (water equiv)
    real :: soilwcl_ob(ix*il)              ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    real :: stl_am(ix*il)                   ! land sfc. temperature
    real :: snowd_am(ix*il)                 ! snow depth (water equiv)
    real :: soilw_am(ix*il)                 ! soil water availability

    ! Land sfc. fields from land model
    real :: stl_lm(ix*il)                 ! land-model sfc. temperature

    ! Land masks
    ! Fraction of land
    real :: fmask_l(ix,il)

    ! Binary land mask
    real :: bmask_l(ix,il)

    ! Monthly-mean climatological fields over land
    ! Land surface temperature
    real :: stl12(ix,il,12)

    ! Snow depth (water equiv.)
    real :: snowd12(ix,il,12)

    ! Soil water availabilityend module
    real :: soilw12(ix,il,12)

    contains
        subroutine land_model_init(alb0)
            ! purpose : initialization of land model
            ! initialized common blocks: land_mc

            ! Input variable
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
            call atm2land

            ! 2. Initialize prognostic variables of land model
            stl_lm(:)  = stlcl_ob(:)      ! land sfc. temperature

            ! 3. Compute additional land variables
            call land2atm(0)
        end

        subroutine atm2land
            use mod_date, only: imont1
            use mod_cpl_bcinterp, only: forint, forin5

            ! Interpolate climatological fields to actual date

            ! Climatological land sfc. temperature
            call forin5(imont1, stl12, stlcl_ob)

            ! Climatological snow depth
            call forint(imont1, snowd12, snowdcl_ob)

            ! Climatological soil water availability
            call forint(imont1, soilw12, soilwcl_ob)
        end

        subroutine land2atm(jday)
            use mod_cpl_flags, only: icland

            integer, intent(in) :: jday

            if (jday > 0.and. icland > 0) then
                ! Run land model
                call land_model
            end if

            ! 3. Compute land-sfc. fields for atm. model
            ! 3.1 Land sfc. temperature
            if (icland <= 0) then
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

        ! Integrate slab land-surface model for one day
        subroutine land_model
            use mod_flx_land, only: hflux_l

            ! Surface temperature anomaly
            real :: tanom(ix*il)

            ! Land-surface (soil/ice-sheet) layer
            ! Anomaly w.r.t. final-time climatological temperature
            tanom = stl_lm - stlcl_ob

            ! Time evolution of temperature anomaly
            tanom = reshape(cdland, (/ ix*il /))*(tanom + reshape(rhcapl, (/ ix*il /))*hflux_l)

            ! Full surface temperature at final time
            stl_lm = tanom + stlcl_ob
        end
end
