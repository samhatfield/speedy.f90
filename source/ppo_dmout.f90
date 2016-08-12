subroutine dmout(imode)
    ! subroutine dmout(imode)

    ! Purpose : write daily-means into output files
    ! Input :   imode = 0 initialize daily-mean arrays to 0
    !           imode > 0 write daily-means and reset arrays to 0
    ! Modified common blocks : tmsave 

    use mod_tsteps, only: nsteps, nstppr, idout
    use mod_atparam
    use mod_tmean, only: ns2d_d1, ns2d_d2, save2d_d1, save2d_d2

    implicit none

    integer, intent(in) :: imode
    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real*4 :: r4out(ngp), fmean
    integer :: iitest=0, n, nout_d1, nout_d2

    if (iitest.eq.1) print *, 'inside DMOUT'

    if (imode /= 0) then
        ! 1. Divide the accumulated fields to get the means
        ! Fields saved at post-processing steps
        fmean = real(nstppr)/real(nsteps)
    
        save2d_d1(:,:) = save2d_d1(:,:)*fmean
    
        ! Fields saved at every step (fluxes)
        fmean = 1./real(nsteps)
    
        save2d_d2(:,:) = save2d_d2(:,:)*fmean
    
        ! 2. Write daily-mean output file 
        if (idout.eq.1) then
             nout_d1 = 3
             nout_d2 = 1
        else if (idout.eq.2) then
             nout_d1 = ns2d_d1
             nout_d2 = 1
        else
             nout_d1 = ns2d_d1
             nout_d2 = ns2d_d2
        end if
    
        do n=1,nout_d1
            r4out(:) = save2d_d1(:,n)
            write (17) r4out
        end do
      
        do n=1,nout_d2
            r4out(:) = save2d_d2(:,n)
            write (17) r4out
        end do
    end if

    ! 3. Reset arrays to zero for the next daily-mean

    if (iitest.eq.1) print*,' reset to zero'

    save2d_d1(:,:) = 0.

    save2d_d2(:,:) = 0.

    if (iitest.eq.1) print *, 'end of DMOUT'
end
