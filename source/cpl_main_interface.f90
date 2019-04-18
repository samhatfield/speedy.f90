subroutine ini_coupler()
    use mod_atparam
    use mod_cpl_land_model, only: land_model_init, ini_land
    use mod_cpl_sea_model, only: sea_model_init
    use mod_surfcon, only: fmask, alb0
    use mod_cli_sea, only: fmask_s, deglat_s

    implicit none

    ! 1.1 initialize land model constants
    call land_model_init(alb0)

    ! 1.2 initialize land model variables
    call ini_land()

    ! 2.1 initialize sea and ice model constants
    call sea_model_init(fmask_s,deglat_s)

    ! 2.2 initialize sea and ice model variables
    call ini_sea()
end

subroutine agcm_to_coupler(jday)
    use mod_cpl_land_model, only: atm2land

    implicit none

    integer, intent(in) :: jday

    ! 1. send fields to land model
    call atm2land

    ! 2. send fields to sea and ice model
    call atm2sea(jday)
end

subroutine coupler_to_agcm(jday)
    use mod_cpl_land_model, only: land2atm

    implicit none

    integer, intent(in) :: jday

    ! 1. get updated fields from land model
    call land2atm(jday)

    ! 2. get updated fields from sea and ice model
    call sea2atm(jday)
end
