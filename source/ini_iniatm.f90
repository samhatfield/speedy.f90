subroutine ini_atm(cexp)
    !   subroutine ini_atm (cexp)
    !
    !   purpose : call initialization routines for all model common blocks 

    use mod_tsteps, only: nmonts, nsteps, nstout, idout, iyear0, imont0, indrdf, ipout, ihout
    use mod_atparam
    use mod_dyncon1, only: grav, hsg, fsg, radang
    use mod_tmean
    use mod_date, only: ndaytot

    implicit none

    character(len=3) :: cexp        ! experiment identifier
    real :: ppl(kx)            ! post-processing levels (hpa/1000)
    integer :: iitest = 1, is3d = 1, k, nddm, ndm, ndtm, ntm

    ! 1. initialize ffts
    if (iitest == 1) print *, 'calling inifft'
    call inifft

    ! 2. initialize dynamical constants and operators
    if (iitest == 1) print *, 'calling indyns'
    call indyns

    ! 3. set post-processing levels
    do k = 1, kx
        ppl(k) = prlev(fsg(k))
    end do

    ! 4. initialize constants for physical parametrization
    if (iitest == 1) print *, 'calling inphys'
    call inphys(hsg, ppl, radang)

    ! 5. initialize forcing fields (boundary cond. + random forcing)
    if (iitest == 1) print *, 'calling inbcon'
    call inbcon(grav,radang)

    if (iitest == 1) print *, 'calling inirdf'
    call inirdf(indrdf)

    ! 6. initialize model variables
    if (iitest == 1) print *, 'calling invars'
    call invars

    ! 7. initialize time-mean arrays for surface fluxes and output fields
    if (iitest == 1) print *, 'calling dmflux'
    call dmflux(0)

    if (iitest == 1) print *, 'calling dmout'
    call dmout(0)

    if (ihout .eqv. .false.) then ! Do not call time-mean procedures if IHOUT = true 
        if (iitest == 1) print *, 'calling tmout'
        call tmout(0)
    
        ! 8. set up the time-mean and daily-mean output (grads format)
        ! 8.1 control files for time-means
        if (nstout <= 0) then
            ntm  = nmonts
            ndtm = - 1
        else
            ntm  = ndaytot*nsteps/nstout
            ndtm = 1440*nstout/nsteps
        end if
    
        if (iitest == 1) print *, 'calling setctl'
    
        call setctl(12, ix, il, kx, ntm, ndtm, is3d, ns3d1, ns2d_1, ns2d_2,&
            & radang, ppl, 'attm', cexp, iyear0, imont0)
        is3d = is3d + ns3d1
        call setctl(14, ix, il, kx, ntm, ndtm, is3d, ns3d2, 0, 0, radang, ppl,&
            & 'atva', cexp, iyear0, imont0)
    
        is3d = is3d + ns3d2
        call setctl(16, ix, il, kx, ntm, ndtm, is3d, ns3d3, 0, 0, radang, ppl,&
            & 'atdf', cexp, iyear0, imont0)
    
        ! 8.2 control files for daily means
        ndm = ndaytot
        nddm = 1
    
        if (iitest == 1) print *, 'calling setctl_d'
    
        if (idout == 1) then
            call setctl_d(18, ix, il, kx, ndm, nddm, 3, 1, radang, ppl, 'daytm',&
                & cexp, iyear0, imont0)
        else if (idout == 2) then
            call setctl_d(18, ix, il, kx, ndm, nddm, ns2d_d1, 1, radang, ppl,&
                & 'daytm', cexp, iyear0, imont0)
        else if (idout >= 3) then
            call setctl_d(18, ix, il, kx, ndm, nddm, ns2d_d1, ns2d_d2, radang, ppl,&
                & 'daytm', cexp, iyear0, imont0)
        end if
    else
        call iogrid(5) ! create control file for 6-hourly output
        if (ipout) then
            call iogrid(3)
            call geop(1)
        end if
    end if

    ! 8.3 output files for grid-point fields
    if (ihout .eqv. .false.) call setgrd(0,cexp)

    ! Write initial data
    if (ihout .and. ipout) call iogrid(2)
    if (ihout) call iogrid(4)

    contains
        function prlev(siglev)
            ! function prlev (siglev)
            ! purpose : select the closest standard pressure level for post-proc.
            ! input :   siglev = sigma level
            implicit none
        
            real, intent(in) :: siglev
            real :: plev(14) = (/ 0.925, 0.850, 0.775, 0.700, 0.600, 0.500, 0.400,&
                & 0.300, 0.250, 0.200, 0.150, 0.100, 0.050, 0.030 /)
            real :: prlev, dif, adif
            integer :: k
        
            dif = 1.0 - siglev
        
            prlev = 1.0
        
            do k = 1, 14
                adif = abs(plev(k) - siglev)
                if (adif <= dif) then
                    dif = adif
                    prlev = plev(k)
                end if
            end do
        end
end
