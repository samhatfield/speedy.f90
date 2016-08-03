!> @brief
!> Logical flags to control certain behaviour.
module mod_lflags
    implicit none

    private
    public lppres, lco2, lradsw, lrandf

    ! Logical flags to activate processes throughout the integration 

    ! Flag to post-process upper-air fields on pressure levels (.false. for
    ! model level p.p.)
    logical, parameter :: lppres = .true.

    ! Flag for CO2 optical thickness increase
    logical, parameter :: lco2 = .false.

    ! Logical flags to activate processes in selected time steps (updated in
    ! STLOOP)

    ! Flag for shortwave radiation routine
    logical :: lradsw  = .true.

    ! Flag for random diabatic forcing
    logical :: lrandf = .false.
end module
