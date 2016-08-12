module mod_tmean
    use mod_atparam

    implicit none

    private
    public ns2d_1, ns2d_2, ns2d_d1, ns2d_d2, ns3d1, ns3d2, ns3d3, ns3d, save3d,&
        & save2d_1, save2d_2, rnsave, save2d_d1, save2d_d2

    ! Post-processing and output parameters (time-mean and variance fields) 
    ! No. of 2-d time-mean fields, incremented at post-proc. steps 
    integer, parameter :: ns2d_1 = 18

    ! No. of 2-d time-mean fields, incremented every step (fluxes)
    integer, parameter :: ns2d_2 = 12

    ! No. of 2-d daily-mean fields, incremented at post-proc. steps 
    integer, parameter :: ns2d_d1 = 8

    ! No. of 2-d daily-mean fields, incremented every step (fluxes)
    integer, parameter :: ns2d_d2 = 7

    ! No. of 3-d time-mean model variables
    integer, parameter :: ns3d1 = 9

    ! No. of 3-d time-mean variances and covariances
    integer, parameter :: ns3d2 = 6

    ! No. of 3-d time-mean diabatic heating fields
    integer, parameter :: ns3d3 = 5

    integer, parameter :: ns3d = ns3d1 + ns3d2 + ns3d3

    ! Arrays for the computation of time-means
    ! (initial./used by tmout, updated in tminc and dmflux)
    real :: save3d  (ix*il,kx,ns3d)     ! 3-D fields saved every post-proc. step
    real :: save2d_1(ix*il,ns2d_1)        ! 2-D fields saved every post-proc. step
    real :: save2d_2(ix*il,ns2d_2)        ! 2-D fields saved every step (fluxes)
    real :: rnsave                      ! post-processing counter

    ! Arrays for the computation of daily-means
    ! (initial./used by DMOUT, updated in TMINC and DMFLUX)
    real :: save2d_d1(ix*il,ns2d_d1)     ! Daily output saved every post-proc step
    real :: save2d_d2(ix*il,ns2d_d2)     ! Daily output saved every step (fluxes)
end module
