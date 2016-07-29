      SUBROUTINE TMOUT (IMODE)
C--
C--   SUBROUTINE TMOUT (IMODE)
C--
C--   Purpose : write time-means and variances into output files
C--   Input :   IMODE = 0 initialize time-mean arrays to 0
C--             IMODE > 0 write time-means and reset arrays to 0
C--   Modified common blocks : TMSAVE 
C--
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Parameters for post-processing arrays
      include "par_tmean.h"

C     Post-processing arrays (time means)
      include "com_tmean.h"
c      include "com_tmean_daily.h"

C     Time stepping constants
      include "com_tsteps.h"

C     Physical constants
      include "com_physcon.h"

C     Fields used to compute omega, psi and chi
      complex VORSP(mx,nx), DIVSP(mx,nx), PSISP(mx,nx), CHISP(mx,nx)
      equivalence (PSISP,CHISP)

      real    DIV3D(ngp,nlev)
      real*4  R4OUT(ngp)

      iitest=1
      if (iitest.eq.1) print *, 'inside TMOUT'

      if (imode.eq.0) go to 700


C--   1. Divide the accumulated fields to get the means

C     Fields saved at post-processing steps
      fmean = 1./rnsave

      SAVE3D(:,:,:) = SAVE3D(:,:,:)*fmean
      SAVE2D_1(:,:) = SAVE2D_1(:,:)*fmean

C     Fields saved at every step (fluxes)
      fmean = fmean/real(nstppr)

      SAVE2D_2(:,:) = SAVE2D_2(:,:)*fmean

C--   2. Compute omega, psi and chi on p surfaces from wind 

      do k=1,kx

        call VDSPEC (SAVE3D(1,k,3),SAVE3D(1,k,4),
     &               VORSP,DIVSP,2)
        if (ix.eq.iy*4) then
          call TRUNCT (VORSP)
          call TRUNCT (DIVSP)
        endif
        call INVLAP (VORSP,PSISP)
        call GRID (PSISP,SAVE3D(1,k,8),1)
        call INVLAP (DIVSP,CHISP)
        call GRID (CHISP,SAVE3D(1,k,9),1)
        call GRID (DIVSP,DIV3D(1,k),1)

      enddo

      dpr2 = 0.5*pout(1)*p0
      SAVE3D(:,1,7) = -DIV3D(:,1)*dpr2

      do k=2,kx
        dpr2 = 0.5*(POUT(k)-POUT(k-1))*p0
        SAVE3D(:,k,7) = SAVE3D(:,k-1,7)-
     &                  (DIV3D(:,k)+DIV3D(:,k-1))*dpr2
      enddo

      SAVE3D(:,:,8) = SAVE3D(:,:,8)*1.e-6
      SAVE3D(:,:,9) = SAVE3D(:,:,9)*1.e-6

C--   3. Write time-mean output file including 3-d and 2-d fields

      if (iitest.eq.1) print*,' write model output'

      do n=1,ns3d1
        do k=kx,1,-1
          R4OUT(:) = SAVE3D(:,k,n)
          write (11) R4OUT
        enddo
      enddo

      do n=1,ns2d_1
        R4OUT(:) = SAVE2D_1(:,n)
        write (11) R4OUT
      enddo
  
      do n=1,ns2d_2
        R4OUT(:) = SAVE2D_2(:,n)
        write (11) R4OUT
      enddo

C     ----------------------------------------------------------------

      if (ns3d2.gt.0) then

C--   4. Compute variances and covariances

        do n=1,4
          nv=n+ns3d1
          SAVE3D(:,:,nv) = SAVE3D(:,:,nv)-SAVE3D(:,:,n)**2
        enddo

        nuv=ns3d1+5
        SAVE3D(:,:,nuv) = SAVE3D(:,:,nuv)-SAVE3D(:,:,3)*SAVE3D(:,:,4)

        nvt=ns3d1+6
        SAVE3D(:,:,nvt) = SAVE3D(:,:,nvt)-SAVE3D(:,:,2)*SAVE3D(:,:,4)

C--   5. Write 2-nd order moments 

        do n=ns3d1+1,ns3d1+ns3d2
          do k=kx,1,-1
            R4OUT(:) = SAVE3D(:,k,n)
            write (13) R4OUT
          enddo
        enddo

      endif

C     ----------------------------------------------------------------

      if (ns3d3.gt.0) then

C--   6. Write diabatic forcing fields (in degK/day)

        do n=ns3d1+ns3d2+1,ns3d
          do k=kx,1,-1
            R4OUT(:) = SAVE3D(:,k,n)*86400.
            write (15) R4OUT
          enddo
        enddo

      endif

C     ----------------------------------------------------------------

C--   7. Reset arrays to zero for the next time-mean

  700 continue

      if (iitest.eq.1) print*,' reset to zero'

      rnsave = 0.

      SAVE3D(:,:,:) = 0.

      SAVE2D_1(:,:) = 0.

      SAVE2D_2(:,:) = 0.

C--
      if (iitest.eq.1) print *, 'end of TMOUT'

      RETURN
      END
