! Compute large-scale precipitation and associated tendencies of temperature
! and moisture
! Input:   psa    = norm. surface pressure [p/p0]           (2-dim)
!          qa     = specific humidity [g/kg]                (3-dim)
!          qsat   = saturation spec. hum. [g/kg]            (3-dim)
!          itop   = top of convection (layer index)         (2-dim)
! Output:  itop   = top of conv+l.s.condensat.(layer index) (2-dim)
!          precls = large-scale precipitation [g/(m^2 s)]   (2-dim)
!          dtlsc  = temperature tendency from l.s. cond     (3-dim)
!          dqlsc  = hum. tendency [g/(kg s)] from l.s. cond (3-dim)
subroutine lscond(psa, qa, qsat, itop, precls, dtlsc, dqlsc)
    use mod_lsccon
    use mod_atparam
    use mod_physcon, only: p0, cp, alhc, alhs, sig, dsig
    use mod_dyncon1, only: grav

    implicit none

    real, intent(in) :: psa(ix,il), qa(ix,il,kx), qsat(ix,il,kx)

    integer, intent(inout) :: itop(ix,il)
    real, intent(inout) :: precls(ix,il), dtlsc(ix,il,kx), dqlsc(ix,il,kx)

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
        sig2 = sig(k)**2.0
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
        pfact = dsig(k)*prg
        precls = precls - pfact*dqlsc(:,:,k)
    end do

    precls = precls*psa
end
