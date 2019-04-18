subroutine ini_coupler()
    use mod_atparam
    use mod_cpl_land_model, only: land_model_init, couple_land_atm
    use mod_cpl_sea_model, only: sea_model_init, ini_sea
    use mod_surfcon, only: alb0

    implicit none

    ! 1.1 initialize land model constants
    call land_model_init(alb0)

    ! 1.2 initialize land model variables
    call couple_land_atm(0)

    ! 2.1 initialize sea and ice model constants
    call sea_model_init

    ! 2.2 initialize sea and ice model variables
    call ini_sea
end

subroutine coupler(day)
    use mod_cpl_land_model, only: couple_land_atm
    use mod_cpl_sea_model, only: atm2sea, sea2atm

    implicit none

    integer, intent(in) :: day

    call couple_land_atm(day)

    ! 2. send fields to sea and ice model
    call atm2sea(day)

    ! 2. get updated fields from sea and ice model
    call sea2atm(day)
end