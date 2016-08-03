!> @brief
!> Flags to set coupling options (see doc_instep.txt).
module mod_cpl_flags
    implicit none

    private
    public icland, icsea, icice, isstan

    ! Flag for land-coupling
    integer :: icland = 1

    ! Flag for sea (SST) coupling
    integer :: icsea  = 0

    ! Flag for sea-ice coupling
    integer :: icice  = 1

    ! Flag for observed SST anomaly
    integer :: isstan = 1
end module
