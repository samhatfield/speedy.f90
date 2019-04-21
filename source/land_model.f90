module land_model
    use mod_atparam

    implicit none

    private
    public stlcl_ob, stl_am, snowd_am, soilw_am
    public land_model_init, couple_land_atm
    public fmask_l, bmask_l, stl12, snowd12, soilw12
    public hflux_l
    public land_coupling_flag

    ! 1./heat_capacity (land)
    real :: rhcapl(ix,il)

    ! 1./dissip_time (land)
    real :: cdland(ix,il)

    ! Daily observed climatological fields over land
    real :: stlcl_ob(ix,il)                ! clim. land sfc. temperature
    real :: snowdcl_ob(ix,il)              ! clim. snow depth (water equiv)
    real :: soilwcl_ob(ix,il)              ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    real :: stl_am(ix,il)                   ! land sfc. temperature
    real :: snowd_am(ix,il)                 ! snow depth (water equiv)
    real :: soilw_am(ix,il)                 ! soil water availability

    ! Land sfc. fields from land model
    real :: stl_lm(ix,il)                 ! land-model sfc. temperature

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

    ! Net heat flux into land surface
    real :: hflux_l(ix,il)

    ! Flag for land-coupling
    ! 0 = no coupling
    ! 1 = land-model)
    integer :: land_coupling_flag = 1

    contains
        subroutine land_model_init(alb0)
            ! purpose : initialization of land model

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

        subroutine couple_land_atm(day)
            use mod_date, only: imont1
            use mod_cpl_bcinterp, only: forin5, forint

            integer, intent(in) :: day

            ! Interpolate climatological fields to actual date

            ! Climatological land surface temperature
            call forin5(imont1, stl12, stlcl_ob)

            ! Climatological snow depth
            call forint(imont1, snowd12, snowdcl_ob)

            ! Climatological soil water availability
            call forint(imont1, soilw12, soilwcl_ob)

            ! If it's the first day then initialise the land surface
            ! temperature from climatology
            if (day == 0) then
                stl_lm = stlcl_ob
                stl_am = stlcl_ob
            else
                ! Run the land model if the land model flags is switched on
                if (land_coupling_flag == 1) then
                    call run_land_model

                    stl_am = stl_lm
                ! Otherwise get the land surface from climatology
                else
                    stl_am = stlcl_ob
                end if
            end if

            ! Always get snow depth and soil water availability from climatology
            snowd_am = snowdcl_ob
            soilw_am = soilwcl_ob
        end subroutine

        ! Integrate slab land-surface model for one day
        subroutine run_land_model
            ! Surface temperature anomaly
            real :: tanom(ix,il)

            ! Land-surface (soil/ice-sheet) layer
            ! Anomaly w.r.t. final-time climatological temperature
            tanom = stl_lm - stlcl_ob

            ! Time evolution of temperature anomaly
            tanom = cdland*(tanom + rhcapl*hflux_l)

            ! Full surface temperature at final time
            stl_lm = tanom + stlcl_ob
        end
end
