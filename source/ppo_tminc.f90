subroutine tminc
    ! subroutine tminc
    !
    ! Purpose : perform post-processing on pressure levels
    !           and increment time-mean arrays
    ! Modified common blocks : TMSAVE
    !

    use mod_lflags, only: lppres
    use mod_atparam
    use mod_tmean, only: ns3d1, ns3d2, ns3d3, save3d, save2d_1, rnsave,&
        & save2d_d1
    use mod_physcon, only: gg, rd, sigl, pout
    use mod_surfcon, only: phis0
    use mod_cli_land, only: bmask_l
    use mod_var_land, only: stl_am, soilw_am
    use mod_cli_sea, only: bmask_s
    use mod_var_sea, only: sst_am, sstan_am, sst_om, ssti_om
    use mod_physvar
    use mod_radcon, only: albsfc

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real :: adsave(ngp,6), phisg(ngp), pmsl(ngp), qsatpl(ngp), st0(ngp)

    ! Fields for vertical interpolation
    integer :: k0(ngp), iitest, nv, nvt, nuv, n0, n, kmid, klow, kj1, j, k, kj
    integer :: khigh
    real :: w0(ngp), zout(ngp), zinp(nlev), rdzinp(nlev)
    real :: rg, rdr2, gam0, rgam, zmin, tsg, wref, tref, textr, rrgam, plog
    real :: phi1, phi2, fwind, aref

    phisg = reshape(phis0, (/ngp/))

    ! Level indices for daily-mean storage of upper-air fields
    ! Upper tropos. u, v :
    khigh = 3
    ! Mid-tropos. geopotential height :
    kmid  = 5
    ! Lower tropos. u, v, q :
    klow  = 7

    iitest=0
    if (iitest.eq.1) print *, ' inside tminc'

    rg    = 1./gg
    rdr2  = 0.5*rd
    gam0  = 0.006*rg
    rgam  = rd*gam0
    rrgam = 1./rgam

    ! 0. Increment post-processing counter
    rnsave = rnsave +1

    ! 1. Store 2-d time-mean fields

    ! 1.1 Compute additional surface fields
    if (iitest.eq.1) print*,' store 2d fields'

    ! Mean-sea-level pressure
    do j=1,ngp
        tsg=0.5*(t0(j)+max(255.,min(295.,t0(j))))
        pmsl(j)=psg(j)*(1.+gam0*phisg(j)/tsg)**rrgam
    end do

    ! 1.2 Increment time-mean arrays
    n0=0
    call add1f(save2d_1,psg,       ngp,n0,1000)
    call add1f(save2d_1,pmsl,      ngp,n0,1000)
    call add1f(save2d_1,ts,        ngp,n0,1)
    call add1f(save2d_1,tskin,     ngp,n0,1)
    call add1f(save2d_1,soilw_am,  ngp,n0,100)
    call add1f(save2d_1,albsfc,    ngp,n0,100)
    call add1f(save2d_1,u0,        ngp,n0,1)
    call add1f(save2d_1,v0,        ngp,n0,1)
    call add1f(save2d_1,t0,        ngp,n0,1)
    call add1f(save2d_1,rh(1,nlev),ngp,n0,100)
    call add1f(save2d_1,cloudc,    ngp,n0,100)
    call add1f(save2d_1,clstr,     ngp,n0,100)
    call add1f(save2d_1,cltop,     ngp,n0,1000)
    call add1f(save2d_1,prtop,     ngp,n0,1)

    ! Land and sea surface temperatures
    call maskout(stl_am,st0,bmask_l,ngp)
    call add1f(save2d_1,st0,ngp,n0,1)
    call maskout(sst_am,st0,bmask_s,ngp)
    call add1f(save2d_1,st0,ngp,n0,1)
    ! Ocean model SST
    !call maskout(sst_om,st0,bmask_s,ngp)
    ! Ocean model SST + T_ice
    call maskout(ssti_om,st0,bmask_s,ngp)
    call add1f(save2d_1,st0,ngp,n0,1)
    ! SST anomaly (wrt obs. clim.)
    call maskout(sstan_am,st0,bmask_s,ngp)
    call add1f(save2d_1,st0,ngp,n0,1)

    ! NB Fluxes of water, energy and momentum are stored 
    ! every time step by subroutine fluxinc

    ! 1.3 Increment daily-mean arrays
    do j=1,ngp
         save2d_d1(j,1)=save2d_d1(j,1)+pmsl(j)*1000
         save2d_d1(j,2)=save2d_d1(j,2)+t0(j)
    end do

    ! 2. Perform vertical interpolation from sigma to pressure levels
    !    and increment time-mean fields
    if (iitest.eq.1) print*, ' store 3d fields'

    zinp(1)  =-sigl(1)
    do k=2,nlev
        zinp(k)  =-sigl(k)
        rdzinp(k)= 1./(zinp(k-1)-zinp(k))
    end do

    zmin = zinp(nlev)

    do k=1,kx
        ! 2.1 Set coefficients for vertical interpolation
        ! using coordinate Z = log (p_s/p) 
        if (lppres) then
            plog=log(pout(k))
            do j=1,ngp
                zout(j)=pslg1(j)-plog
            end do
        else
            ! Set zout=zinp(k) to do post-proc. on sigma levels
            do j=1,ngp
                zout(j)=zinp(k)
            end do
        end if

        call setvin(zinp,rdzinp,zout,ngp,kx,k0,w0)

        ! 2.2 Interpolate 3-d fields
        ! Temperature (extrapolated below the lowest level when W0(j)<0)
        call verint(adsave(1,2),tg1,ngp,kx,k0,w0) 

        ! Remove extrapolation of temperature inversions 
        ! and correct extrap. values using a reference lapse rate
        wref = 0.7

        do j=1,ngp
            if (zout(j).lt.zmin) then
                textr = max(adsave(j,2),tg1(j,nlev))
                aref = rgam*(zmin-zout(j))
                tref = tg1(j,nlev)*(1.+aref+0.5*aref*aref)
                adsave(j,2) = textr+wref*(tref-textr)
            end if
        end do

        ! Geopotential (computed from the closest levels 
        ! using the hydrostatic equation)
        do j=1,ngp
            w0(j)=max(w0(j),0.)
        end do

        if (lppres) then
            do j=1,ngp
                kj=k0(j)
                kj1=kj-1
                phi1=phig1(j,kj)+rdr2*(adsave(j,2)+tg1(j,kj ))&
                    & *(zout(j)-zinp(kj))
                phi2=phig1(j,kj1)+rdr2*(adsave(j,2)+tg1(j,kj1))&
                    & *(zout(j)-zinp(kj1))
                adsave(j,1)=phi1+w0(j)*(phi2-phi1)
            end do
        else
            call verint(adsave(1,1),phig1,ngp,kx,k0,w0) 
        end if

        ! Wind and relative humidity 
        ! a) Interpolate above the lowest level
        call verint(adsave(1,3),ug1,ngp,kx,k0,w0) 
        call verint(adsave(1,4),vg1,ngp,kx,k0,w0) 
        call verint(adsave(1,6),rh, ngp,kx,k0,w0) 

        ! b) Decrease wind speed below the lowest level
        do j=1,ngp
            if (zout(j).lt.zmin) then
                fwind=adsave(j,1)/phig1(j,nlev)
                adsave(j,3)=adsave(j,3)*fwind
                adsave(j,4)=adsave(j,4)*fwind  
            end if
        end do

        ! Estimate specific humidity using interpolated rel.hum. and
        ! sat. spec.hum. at interpolated temperature
        if (lppres) then
            call shtorh(-1,ngp,adsave(1,2),pout(k),-1.,adsave(1,5),adsave(1,6),&
                & qsatpl)

            ! Below the surface, set spec.hum. = near-surface value 
            do j=1,ngp
                if (zout(j).lt.0.0) then
                    adsave(j,5)=q0(j)
                    adsave(j,6)=q0(j)/qsatpl(j)
                end if
            end do
        else
            call verint(adsave(1,5),qg1,ngp,kx,k0,w0) 
        end if

        ! Rescale geopotential and rel. humidity
        do j=1,ngp
            adsave(j,1)=adsave(j,1)*rg
            adsave(j,6)=adsave(j,6)*100.
        end do

        ! 2.3 Save upper-air fields

        ! 2.3.1 Add 3d upper-air fields to time-mean arrays
        do n=1,6
            do j=1,ngp
                save3d(j,k,n)=save3d(j,k,n)+adsave(j,n)
            end do
        end do

        ! 2.3.2 Add upper-air fields fields at selected levels to daily-means
        ! arrays
        if (k.eq.kmid) then
            do j=1,ngp
                save2d_d1(j,3)=save2d_d1(j,3)+adsave(j,1)
            end do
        end if

        if (k.eq.klow) then
            do j=1,ngp
                save2d_d1(j,4)=save2d_d1(j,4)+adsave(j,3)
                save2d_d1(j,5)=save2d_d1(j,5)+adsave(j,4)
                save2d_d1(j,6)=save2d_d1(j,6)+adsave(j,5)
            end do
        else if (k.eq.khigh) then
            do j=1,ngp
                save2d_d1(j,7)=save2d_d1(j,7)+adsave(j,3)
                save2d_d1(j,8)=save2d_d1(j,8)+adsave(j,4)
            end do
        end if

        ! 2.4 Store variances on pressure levels
        if (ns3d2.gt.0) then
            do n=1,4
                nv=ns3d1+n
                do j=1,ngp
                    save3d(j,k,nv)=save3d(j,k,nv)+adsave(j,n)*adsave(j,n)
                end do
            end do

            nuv=ns3d1+5
            nvt=ns3d1+6
            do j=1,ngp
                save3d(j,k,nuv)=save3d(j,k,nuv)+adsave(j,3)*adsave(j,4)
                save3d(j,k,nvt)=save3d(j,k,nvt)+adsave(j,2)*adsave(j,4)
            end do
        end if
       ! end-of-loop over pressure levels
    end do

    ! 3. Store diabatic forcing terms on model levels
    if (ns3d3.gt.0) then
        n0=ns3d1+ns3d2
        call add1f(save3d,tt_lsc,ngp*nlev,n0,1)
        call add1f(save3d,tt_cnv,ngp*nlev,n0,1)
        call add1f(save3d,tt_rsw,ngp*nlev,n0,1) 
        call add1f(save3d,tt_rlw,ngp*nlev,n0,1) 
        call add1f(save3d,tt_pbl,ngp*nlev,n0,1) 
    end if

    if (iitest.eq.1) print *, 'end of tminc'
end

subroutine setvin(zinp,rdzinp,zout,ngp,nlev,k0,w0)
    implicit none

    real, intent(in) :: zinp(nlev), rdzinp(nlev), zout(ngp)
    integer, intent(in) :: ngp, nlev
    real, intent(inout) :: w0(ngp)
    integer, intent(inout) :: k0(ngp)
    integer :: j, k
 
    ! *** 1. Select closest vertical levels
    do j=1,ngp
        k0(j)=2
    end do
 
    do k=2,nlev-1
        do j=1,ngp
            if (zout(j).lt.zinp(k)) k0(j)=k+1
        end do
    end do
 
    ! *** 2. Compute interpolation weight
    do j=1,ngp
        w0(j)=(zout(j)-zinp(k0(j)))*rdzinp(k0(j))
    end do
end

subroutine verint(f2d,f3d,ngp,nlev,k0,w0)
    implicit none

    ! *** 1. Perform vertical interpolation 
    integer, intent(in) :: ngp, nlev, k0(ngp)
    real, intent(in) :: f3d(ngp,nlev), w0(ngp)
    real, intent(inout) :: f2d(ngp)
    integer :: j

    do j=1,ngp
        f2d(j)=f3d(j,k0(j))+w0(j)*(f3d(j,k0(j)-1)-f3d(j,k0(j)))
    end do
end

subroutine add1f(fsave,fadd,ngp,nf,ifact)
    implicit none

    ! *** Add one field to storage array 
    real, intent(in) :: fadd(ngp)
    integer, intent(in) :: ngp, ifact
    real, intent(inout) :: fsave(ngp,*)
    integer, intent(inout) :: nf
    integer :: j

    nf=nf+1

    if (ifact.eq.1) then 
        do j=1,ngp
            fsave(j,nf)=fsave(j,nf)+fadd(j)
        end do
    else
        do j=1,ngp
            fsave(j,nf)=fsave(j,nf)+fadd(j)*ifact
        end do
    end if
end

subroutine maskout(finp,fout,fmask,ngp)
    implicit none

    ! *** Set undefined values according to binary land-sea mask
    real, intent(in) :: finp(ngp), fmask(ngp) 
    integer, intent(in) :: ngp
    real, intent(inout) :: fout(ngp)
    integer :: j
    real :: xundef

    xundef = 9.999e+19

    do j=1,ngp
        if (fmask(j).le.0.) then
            fout(j) = xundef
        else
            fout(j) = finp(j) 
        end if
    end do
end
