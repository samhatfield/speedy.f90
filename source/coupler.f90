module coupler
    implicit none

    private
    public initialize_coupler, couple_sea_land

contains
    subroutine initialize_coupler
        use land_model, only: land_model_init, couple_land_atm
        use sea_model, only: sea_model_init, couple_sea_atm

        ! Initialize land model constants
        call land_model_init

        ! Initialize land model variables
        call couple_land_atm(0)

        ! Initialize sea and ice model constants
        call sea_model_init

        ! Initialize sea and ice model variables
        call couple_sea_atm(0)
    end

    subroutine couple_sea_land(day)
        use land_model, only: couple_land_atm
        use sea_model, only: couple_sea_atm

        integer, intent(in) :: day

        call couple_land_atm(day)
        call couple_sea_atm(day)
    end
end module
