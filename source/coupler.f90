!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 29/04/2019
!> Highest level interface to land and sea models.
module coupler
    implicit none

    private
    public initialize_coupler, couple_sea_land

contains
    !> Initialize both land and sea models.
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
    end subroutine

    !> Exchange fluxes between atmosphere and land/sea.
    subroutine couple_sea_land(day)
        use land_model, only: couple_land_atm
        use sea_model, only: couple_sea_atm

        integer, intent(in) :: day !! The current day of the model integration (starting from 0)

        call couple_land_atm(day)
        call couple_sea_atm(day)
    end subroutine
end module
