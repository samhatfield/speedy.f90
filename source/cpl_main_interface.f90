subroutine ini_coupler()
    use mod_atparam
    use land_model, only: land_model_init, couple_land_atm
    use sea_model, only: sea_model_init, couple_sea_atm
    use mod_surfcon, only: alb0

    implicit none

    ! 1.1 initialize land model constants
    call land_model_init(alb0)

    ! 1.2 initialize land model variables
    call couple_land_atm(0)

    ! 2.1 initialize sea and ice model constants
    call sea_model_init

    ! 2.2 initialize sea and ice model variables
    call couple_sea_atm(0)
end

subroutine coupler(day)
    use land_model, only: couple_land_atm
    use sea_model, only: couple_sea_atm

    implicit none

    integer, intent(in) :: day

    call couple_land_atm(day)
    call couple_sea_atm(day)
end
