subroutine indyns
    ! subroutine indyns
    !
    ! Purpose : set time-stepping constants and initialize coefficients
    !           and spectral operators for model dynamics
    ! Initialized common blocks: dync0, dync1,  dync2,  dync3,  dync4, 
    !                            hdifc1, hdifc3,
    !                            common blocks for spectral transforms 
    !                            (through routine parmtr) 
    !

    use mod_tsteps, only: nsteps, alph
    use mod_dyncon0
    use mod_dyncon1
    use mod_atparam
    use mod_hdifcon, only: dmp, dmpd, dmps, tcorv, qcorv
    use mod_spectral, only: sia, cosg

    implicit none

    integer :: j, k, jj, npowhd
    real :: elap, elapn, hdifd, hdiff, hdifs, qexp, rad1, rgam, rlap, twn

    ! 1. Definition of constants
    if (mod(nsteps,2) /= 0) stop ' Invalid no. of time steps'

    ! alph = 0 ---- forward step for gravity wave terms
    ! alph = 1 ---- backward implicit -----------------
    ! alph = 0.5 -- centered implicit -----------------
    alph = 0.5 

    ! Power of Laplacian in horizontal diffusion
    npowhd = 4

    ! 2. Definition of model levels  

    ! 2.1 Half (vertical velocity) levels
    if (kx == 5) then
        hsg(:6) = (/ 0.000, 0.150, 0.350, 0.650, 0.900, 1.000 /)
    else if (kx == 7) then
        hsg(:8) = (/ 0.020, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000 /)
    else if (kx == 8) then
        hsg(:9) = (/ 0.000, 0.050, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000 /)
    end if

    do k = 1, kxp
        print *, ' Model half-level (*1000)', k, nint(HSG(k)*1000)
    end do

    ! 2.2 Layer thicknesses and full (u,v,T) levels
    do k = 1, kx
        dhs(k) = hsg(k+1)-hsg(k)
        fsg(k) = 0.5*(hsg(k+1)+hsg(k))
    end do

    do k = 1, kx
        print *, ' Model full-level (*1000)', k, nint(FSG(k)*1000)
    end do

    ! 2.3 Additional functions of sigma
    do k = 1, kx
        dhsr(k) = 0.5/dhs(k)
        fsgr(k) = akap/(2.*fsg(k))
    end do

    ! 3. Horizontal functions and spectral operators

    ! 3.1 Initialization of spectral operators 
    call parmtr(rearth)

    ! 3.2 Latitudes and functions of latitude
    !     NB: J=1 is Southernmost point!
    do j = 1, iy
        jj = il + 1 - j
        rad1 = asin(sia(j))
        radang(j)  = -rad1
        radang(jj) =  rad1
        gsin(j)    = -sia(j)
        gsin(jj)   =  sia(j)
    end do

    do j = 1, il
        gcos(j) = cosg(j)
        coriol(j) = 2.*omega*gsin(j)
    end do


    ! 4. Coefficients to compute geopotential
    do k = 1, kx
      xgeop1(k) = rgas*log(hsg(k+1)/fsg(k))
      if (k /= kx) xgeop2(k+1) = rgas*log(fsg(k+1)/hsg(k+1))
    end do

    ! 5. Coefficients for horizontal diffusion

    ! 5.1 Spectral damping coefficients
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
        print *, ' temp/hum correction at level ', k, tcorv(k), qcorv(k)
    end do
end
