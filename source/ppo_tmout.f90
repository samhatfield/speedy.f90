subroutine tmout(imode)
    !  subroutine tmout (imode)
    !
    !  Purpose : write time-means and variances into output files
    !  Input :   imode = 0 initialize time-mean arrays to 0
    !            imode > 0 write time-means and reset arrays to 0
    !  Modified common blocks : TMSAVE 
    !

    use mod_tsteps, only: nstppr
    use mod_atparam
    use mod_tmean
    use mod_physcon, only: p0, pout

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer, intent(in) :: imode

    ! Fields used to compute omega, psi and chi
    complex :: vorsp(mx,nx), divsp(mx,nx), psisp(mx,nx)

    real :: div3d(ngp,nlev), dpr2, fmean
    real*4 :: r4out(ngp)

    integer :: iitest=1, k, n, nuv, nv, nvt
    if (iitest.eq.1) print *, 'inside TMOUT'

    if (imode.eq.0) go to 700

    ! 1. Divide the accumulated fields to get the means
    ! Fields saved at post-processing steps
    fmean = 1./rnsave

    save3d(:,:,:) = save3d(:,:,:)*fmean
    save2d_1(:,:) = save2d_1(:,:)*fmean

    ! Fields saved at every step (fluxes)
    fmean = fmean/real(nstppr)

    save2d_2(:,:) = save2d_2(:,:)*fmean

    ! 2. Compute omega, psi and chi on p surfaces from wind 
    do k=1,kx
        call vdspec (save3d(1,k,3),save3d(1,k,4),vorsp,divsp,2)
        if (ix.eq.iy*4) then
            call TRUNCT (VORSP)
            call TRUNCT (DIVSP)
        end if

        call invlap (vorsp,psisp)
        call grid (psisp,save3d(1,k,8),1)
        call invlap (divsp,psisp)
        call grid (psisp,save3d(1,k,9),1)
        call grid (divsp,div3d(1,k),1)
    end do

    dpr2 = 0.5*pout(1)*p0
    save3d(:,1,7) = -div3d(:,1)*dpr2

    do k=2,kx
        dpr2 = 0.5*(POUT(k)-POUT(k-1))*p0
        save3d(:,k,7) = save3d(:,k-1,7)-(div3d(:,k)+div3d(:,k-1))*dpr2
    end do

    save3d(:,:,8) = save3d(:,:,8)*1.e-6
    save3d(:,:,9) = save3d(:,:,9)*1.e-6

    ! 3. Write time-mean output file including 3-d and 2-d fields
    if (iitest.eq.1) print*,' write model output'

    do n=1,ns3d1
        do k=kx,1,-1
            r4out(:) = save3d(:,k,n)
            write (11) r4out
        end do
    end do

    do n=1,ns2d_1
        r4out(:) = save2d_1(:,n)
        write (11) r4out
    end do
  
    do n=1,ns2d_2
        r4out(:) = save2d_2(:,n)
        write (11) r4out
    end do

    !----------------------------------------------------------------

    if (ns3d2.gt.0) then
        ! 4. Compute variances and covariances
        do n=1,4
            nv=n+ns3d1
            save3d(:,:,nv) = save3d(:,:,nv)-save3d(:,:,n)**2
        end do

        nuv=ns3d1+5
        save3d(:,:,nuv) = save3d(:,:,nuv)-save3d(:,:,3)*save3d(:,:,4)

        nvt=ns3d1+6
        save3d(:,:,nvt) = save3d(:,:,nvt)-save3d(:,:,2)*save3d(:,:,4)

        ! 5. Write 2-nd order moments 
        do n=ns3d1+1,ns3d1+ns3d2
            do k=kx,1,-1
                r4out(:) = save3d(:,k,n)
                write (13) r4out
            end do
        end do
    end if

    !----------------------------------------------------------------

    if (ns3d3.gt.0) then
        ! 6. Write diabatic forcing fields (in degK/day)
        do n=ns3d1+ns3d2+1,ns3d
            do k=kx,1,-1
                r4out(:) = save3d(:,k,n)*86400.
                write (15) r4out
            end do
        end do
    end if

    !----------------------------------------------------------------

    ! 7. Reset arrays to zero for the next time-mean
    700 continue

    if (iitest.eq.1) print*,' reset to zero'

    rnsave = 0.

    SAVE3D(:,:,:) = 0.
    SAVE2D_1(:,:) = 0.
    SAVE2D_2(:,:) = 0.

    if (iitest.eq.1) print *, 'end of TMOUT'
end
