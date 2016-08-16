module mod_cplvar_sea
    use mod_atparam

    implicit none

    private
    public vsea_input, vsea_output

    ! Input and output sea variables exchanged by coupler
    ! Ocean model input variables
    real :: vsea_input(ix*il,8)

    ! Ocean model output variablesend module
    real :: vsea_output(ix*il,3)
end module
