!fk#if !defined(KNMI)
subroutine lscond(psa,qa,qsat,itop,precls,dtlsc,dqlsc)
!fk#else
!fksubroutine lscond (psa,qa,qsat,ts,itop,precls,snowls,dtlsc,dqlsc)
!fk#endif
    !  subroutine lscond (psa,qa,qsat,
    ! *                   itop,precls,dtlsc,dqlsc) 
    !
    !  Purpose: Compute large-scale precipitation and
    !           associated tendencies of temperature and moisture
    !  Input:   psa    = norm. surface pressure [p/p0]           (2-dim)
    !           qa     = specific humidity [g/kg]                (3-dim)
    !           qsat   = saturation spec. hum. [g/kg]            (3-dim)
    !           itop   = top of convection (layer index)         (2-dim)
    !  Output:  itop   = top of conv+l.s.condensat.(layer index) (2-dim)
    !           precls = large-scale precipitation [g/(m^2 s)]   (2-dim)
    !           dtlsc  = temperature tendency from l.s. cond     (3-dim)
    !           dqlsc  = hum. tendency [g/(kg s)] from l.s. cond (3-dim)

    use mod_lsccon
    use mod_atparam
    use mod_physcon, only: p0, gg, cp, alhc, alhs, sig, dsig

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real, intent(in) :: psa(ngp), qa(ngp,nlev), qsat(ngp,nlev)

    integer, intent(inout) :: itop(ngp)
    real, intent(inout) :: precls(ngp), dtlsc(ngp,nlev), dqlsc(ngp,nlev)
    !fk#if defined(KNMI)
    !fkreal, intent(in) :: ts(ngp)
    !fkreal, intent(inout) :: snowls(ngp)
    !fk#endif

    integer :: j, k
    real :: psa2(ngp), dqa, dqmax, pfact, prg, qsmax, rhref, rtlsc, sig2, tfact

    ! 1. Initialization
    qsmax = 10.

    rtlsc = 1./(trlsc*3600.)
    tfact = alhc/cp
    !fk#if defined(KNMI)
    !fktfacts= alhs/cp
    !fk#endif
    prg = p0/gg

    dtlsc(:,1) = 0.
    dqlsc(:,1) = 0.
    precls  = 0.
   !fk#if defined(KNMI)
   !fksnowls = 0.
   !fk#endif
    do j=1,ngp
        psa2(j) = psa(j)*psa(j)
    end do

    ! 2. Tendencies of temperature and moisture
    !    NB. A maximum heating rate is imposed to avoid 
    !        grid-point-storm instability 
    do k=2,nlev
        sig2=sig(k)*sig(k)
        rhref = rhlsc+drhlsc*(sig2-1.)
        if (k.eq.nlev) rhref = max(rhref,rhblsc)
        dqmax = qsmax*sig2*rtlsc

        do j=1,ngp
            dqa = rhref*qsat(j,k)-qa(j,k)
            if (dqa.lt.0.0) then
                itop(j)    = min(k,itop(j))
                dqlsc(j,k) = dqa*rtlsc
                !fk#if !defined(KNMI)
                dtlsc(j,k) = tfact*min(-dqlsc(j,k),dqmax*psa2(j))
                !fk#else
                !fkif (ts(j).gt.273.15) then
                !fk    dtlsc(j,k) = tfact*min(-dqlsc(j,k),dqmax*psa2(j))
                !fkelse
                !fk    dtlsc(j,k) = tfacts*min(-dqlsc(j,k),dqmax*psa2(j))
                !fkend if
                !fk#endif
            else
                dqlsc(j,k) = 0.
                dtlsc(j,k) = 0.
            endif
        end do
    end do

    ! 3. Large-scale precipitation
    do k=2,nlev
        pfact = dsig(k)*prg
        do j=1,ngp
            precls(j) = precls(j)-pfact*dqlsc(j,k)
            !fk#if defined(KNMI)
            !fkif (ts(j).lt.273.15) then
            !fk    snowls(j) = snowls(j)-pfact*dqlsc(j,k)
            !fkendif
            !fk#endif
        end do
    end do

    do j=1,ngp
        precls(j) = precls(j)*psa(j)
        !fk#if defined(KNMI)
        !fksnowls(j) = snowls(j)*psa(j)
        !fk#endif
    end do
end
