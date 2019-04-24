module land_model
    use mod_atparam

    implicit none

    private
    public stl_am, snowd_am, soilw_am
    public land_model_init, couple_land_atm
    public fmask_l, hflux_l
    public land_coupling_flag
    public sd2sc

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

    ! Snow depth (mm water) corresponding to snow cover = 1
    real :: sd2sc = 60.0

    contains
        subroutine land_model_init
            ! purpose : initialization of land model
            use input_output, only: load_boundary_file
            use boundaries, only: forchk, fmask, alb0, fillsf

            ! Auxiliary variables
            integer :: i, j, month
            real :: dmask(ix,il)           ! domain mask
            real :: depth_soil, depth_lice, tdland, hcapl, hcapli, flandmin

            ! Soil moisture parameters
            ! Soil wetness at field capacity (volume fraction)
            real :: swcap = 0.30

            ! Soil wetness at wilting point  (volume fraction)
            real :: swwil = 0.17

            ! Threshold for land-sea mask definition (i.e. minimum fraction of
            ! either land or sea)
            real :: thrsh = 0.1

            real :: rsw, sdep1, sdep2, swroot, swwil2
            real, dimension(ix,il) :: veg_low, veg_high, veg, swl1, swl2
            integer :: idep2

            ! =========================================================================
            ! Initialize land-surface boundary conditions
            ! =========================================================================

            ! Fractional and binary land masks
            fmask_l = fmask
            do j = 1, il
                do i = 1, ix
                    if (fmask_l(i,j) >= thrsh) then
                        bmask_l(i,j) = 1.0
                        if (fmask(i,j) > (1.0 - thrsh)) fmask_l(i,j) = 1.0
                    else
                        bmask_l(i,j) = 0.0
                        fmask_l(i,j) = 0.0
                    end if
                end do
            end do

            ! Land-surface temperature
            do month = 1, 12
                stl12(:,:,month) = load_boundary_file("land.nc", "stl", month)

                call fillsf(stl12(:,:,month), 0.0)
            end do

            call forchk(bmask_l, 12, 0.0, 400.0, 273.0, stl12)

            ! Snow depth
            do month = 1, 12
                snowd12(:,:,month) = load_boundary_file("snow.nc", "snowd", month)
            end do

            call forchk(bmask_l, 12, 0.0, 20000.0, 0.0, snowd12)

            ! Read soil moisture and compute soil water availability using vegetation fraction
            ! Read vegetation fraction
            veg_high = load_boundary_file("surface.nc", "vegh")
            veg_low  = load_boundary_file("surface.nc", "vegl")

            ! Combine high and low vegetation fractions
            veg = max(0.0, veg_high + 0.8*veg_low)

            ! Read soil moisture
            sdep1 = 70.0
            idep2 = 3
            sdep2 = idep2*sdep1

            swwil2 = idep2*swwil
            rsw    = 1.0/(swcap + idep2*(swcap - swwil))

            do month = 1, 12
                ! Combine soil water content from two top layers
                swl1 = load_boundary_file("soil.nc", "swl1", month)
                swl2 = load_boundary_file("soil.nc", "swl2", month)

                do j = 1, il
                    do i = 1, ix
                        swroot = idep2*swl2(i,j)
                        soilw12(i,j,month) = min(1.0, rsw*(swl1(i,j) + veg(i,j) &
                            & *max(0.0, swroot - swwil2)))
                    end do
                end do
            end do

            call forchk(bmask_l, 12, 0.0, 10.0, 0.0, soilw12)

            ! =========================================================================
            ! Set heat capacities and dissipation times for soil and ice-sheet layers
            ! =========================================================================

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
            use date, only: imont1
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
