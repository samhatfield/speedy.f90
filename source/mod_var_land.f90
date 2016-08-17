module mod_var_land
    use mod_atparam

    implicit none

    private
    public stlcl_ob, snowdcl_ob, soilwcl_ob, stl_am, snowd_am, soilw_am, stl_lm

    ! Daily observed climatological fields over land
    real :: stlcl_ob(ix*il)              ! clim. land sfc. temperature 
    real :: snowdcl_ob(ix*il)              ! clim. snow depth (water equiv)
    real :: soilwcl_ob(ix*il)              ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    real :: stl_am(ix*il)                 ! land sfc. temperature
    real :: snowd_am(ix*il)                 ! snow depth (water equiv)
    real :: soilw_am(ix*il)                 ! soil water availability

    ! Land sfc. fields from land model
    real :: stl_lm(ix*il)                 ! land-model sfc. temperature 
end module
