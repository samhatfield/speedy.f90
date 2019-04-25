!> @brief
!> Length of the integration and time stepping constants.
module mod_tsteps
    implicit none

    private
    public isst0

    ! Record in SST anomaly file corr. to the initial month
    ! Initialized in agcm_init
    integer :: isst0


end module
