subroutine indyns
    ! subroutine indyns
    !
    ! Purpose : set time-stepping constants and initialize coefficients
    !           and spectral operators for model dynamics

    use mod_tsteps, only: nsteps
    use mod_dyncon0
    use mod_dyncon1
    use mod_atparam
    use mod_hdifcon, only: dmp, dmpd, dmps, tcorv, qcorv
    use geometry, only: fsg, hsg, sia, cosg

    implicit none

    integer :: j, k, jj, npowhd
    real :: elap, elapn, hdifd, hdiff, hdifs, qexp, rgam, rlap, twn

    ! 1. Definition of constants
    if (mod(nsteps,2) /= 0) stop ' Invalid no. of time steps'

    ! Power of Laplacian in horizontal diffusion
    npowhd = 4

    ! Coefficients to compute geopotential
    do k = 1, kx
      xgeop1(k) = rgas*log(hsg(k+1)/fsg(k))
      if (k /= kx) xgeop2(k+1) = rgas*log(fsg(k+1)/hsg(k+1))
    end do

    ! Coefficients for horizontal diffusion
    ! Spectral damping coefficients
    hdiff = 1./(thd *3600.)
    hdifd = 1./(thdd*3600.)
    hdifs = 1./(thds*3600.)
    rlap  = 1./float(mtrun*(mtrun+1))

    do j = 1, nx
        do k = 1, mx
            twn = float(isc*(k-1)+j-1)
            elap = (twn*(twn+1.)*rlap)
            elapn = elap**npowhd
            dmp(k,j)  = hdiff*elapn
            dmpd(k,j) = hdifd*elapn
            dmps(k,j) = hdifs*elap
        end do
        ! dmps(1,j)=0.
    end do

    ! 5.2 Orographic correction terms for temperature and humidity
    !     (vertical component)
    rgam = rgas*gamma/(1000.*grav)
    qexp = hscale/hshum

    tcorv(1)=0.
    qcorv(1)=0.
    qcorv(2)=0.

    do k = 2, kx
        tcorv(k) = fsg(k)**rgam
        if (k.gt.2) qcorv(k) = fsg(k)**qexp
    end do
end
