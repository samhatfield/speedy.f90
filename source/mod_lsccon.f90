!> @brief
!> Constants for large-scale condensation.
module mod_lsccon
    implicit none

    private
    public trlsc, rhlsc, drhlsc, rhblsc

    ! Relaxation time (in hours) for specific humidity 
    real, parameter :: trlsc  = 4.0

    ! Maximum relative humidity threshold (at sigma=1)
    real, parameter :: rhlsc  = 0.9

    ! Vertical range of relative humidity threshold
    real, parameter :: drhlsc = 0.1

    ! Relative humidity threshold for boundary layer
    real, parameter :: rhblsc = 0.95
end module
