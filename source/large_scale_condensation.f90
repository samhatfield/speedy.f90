!> Parametrization of large-scale condensation
!
!  Large-scale condensation is modelled as a relaxation of humidity to a
!  sigma-dependent threshold value. The temperature tendency is computed as the
!  resultant latent heating. Precipitation is diagnosed as all the moisture lost
!  to condensation falling out of the atmospheric column in the timestep.
module large_scale_condensation
    use params

    implicit none

    private
    public get_large_scale_condensation_tendencies

    ! Constants for large-scale condensation
    real, parameter :: trlsc  = 4.0  !! Relaxation time (in hours) for specific humidity
    real, parameter :: rhlsc  = 0.9  !! Maximum relative humidity threshold (at sigma=1)
    real, parameter :: drhlsc = 0.1  !! Vertical range of relative humidity threshold
    real, parameter :: rhblsc = 0.95 !! Relative humidity threshold for boundary layer

contains
    !> Compute large-scale condensation and associated tendencies of temperature
    !  and moisture
    subroutine get_large_scale_condensation_tendencies(psa, qa, qsat, itop, precls, dtlsc, dqlsc)
        use physical_constants, only: p0, cp, alhc, alhs, grav
        use geometry, only: fsg, dhs

        real, intent(in) :: psa(ix,il)        !! Normalised surface pressure [p/p0]
        real, intent(in) :: qa(ix,il,kx)      !! Specific humidity [g/kg]
        real, intent(in) :: qsat(ix,il,kx)    !! Saturation specific humidity [g/kg]
        integer, intent(inout) :: itop(ix,il) !! Cloud top diagnosed from precipitation due to
                                              !! convection and large-scale condensation
        real, intent(out) :: precls(ix,il)    !! Precipitation due to large-scale condensation
        real, intent(out) :: dtlsc(ix,il,kx)  !! Temperature tendency due to large-scale
                                              !! condensation
        real, intent(out) :: dqlsc(ix,il,kx)  !! Specific humidity tendency due to large-scale
                                              !! condensation

        integer :: i, j, k
        real :: psa2(ix,il), dqa, dqmax, pfact, prg, qsmax, rhref, rtlsc, sig2, tfact

        ! 1. Initialization
        qsmax = 10.0

        rtlsc = 1.0/(trlsc*3600.0)
        tfact = alhc/cp
        prg = p0/grav

        dtlsc(:,:,1) = 0.0
        dqlsc(:,:,1) = 0.0
        precls  = 0.0

        psa2 = psa**2.0

        ! Tendencies of temperature and moisture
        ! NB. A maximum heating rate is imposed to avoid grid-point-storm
        ! instability
        do k = 2, kx
            sig2 = fsg(k)**2.0
            rhref = rhlsc + drhlsc*(sig2 - 1.0)
            if (k == kx) rhref = max(rhref, rhblsc)
            dqmax = qsmax*sig2*rtlsc

            do i = 1, ix
                do j = 1, il
                    dqa = rhref*qsat(i,j,k) - qa(i,j,k)
                    if (dqa < 0.0) then
                        itop(i,j)    = min(k,itop(i,j))
                        dqlsc(i,j,k) = dqa*rtlsc
                        dtlsc(i,j,k) = tfact*min(-dqlsc(i,j,k), dqmax*psa2(i,j))
                    else
                        dqlsc(i,j,k) = 0.0
                        dtlsc(i,j,k) = 0.0
                    end if
                end do
            end do
        end do

        ! Large-scale precipitation
        do k = 2, kx
            pfact = dhs(k)*prg
            precls = precls - pfact*dqlsc(:,:,k)
        end do

        precls = precls*psa
    end
end module
