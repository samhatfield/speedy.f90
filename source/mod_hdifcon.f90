module mod_hdifcon
    use mod_atparam

    implicit none

    private
    public dmp, dmpd, dmps, dmp1, dmp1d, dmp1s, tcorv, qcorv, tcorh, qcorh

    ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
    real, dimension(mx,nx) :: dmp, dmpd, dmps

    ! Damping coef. for horizontal diffusion (implicit) (initial. in indyns)
    real, dimension(mx,nx) :: dmp1, dmp1d, dmp1s

    ! Vertical comp. of orographic correction (initial. in INDYNS)
    real, dimension(kx) :: tcorv, qcorv

    ! Horizontal component of orographic correction (updated in FORDATE)
    complex, dimension(mx,nx) :: tcorh, qcorh
end module
