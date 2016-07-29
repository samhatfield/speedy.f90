      SUBROUTINE DMOUT (IMODE)
C--
C--   SUBROUTINE DMOUT (IMODE)
C--
C--   Purpose : write daily-means into output files
C--   Input :   IMODE = 0 initialize daily-mean arrays to 0
C--             IMODE > 0 write daily-means and reset arrays to 0
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

C     Time stepping constants
      include "com_tsteps.h"

      real*4  R4OUT(ngp)

      iitest=0
      if (iitest.eq.1) print *, 'inside DMOUT'

      if (imode.eq.0) go to 300


C--   1. Divide the accumulated fields to get the means

C     Fields saved at post-processing steps
      fmean = real(nstppr)/real(nsteps)

      SAVE2D_D1(:,:) = SAVE2D_D1(:,:)*fmean

C     Fields saved at every step (fluxes)
      fmean = 1./real(nsteps)

      SAVE2D_D2(:,:) = SAVE2D_D2(:,:)*fmean

C--   2. Write daily-mean output file 

      if (idout.eq.1) then
         nout_d1 = 3
         nout_d2 = 1
      else if (idout.eq.2) then
         nout_d1 = ns2d_d1
         nout_d2 = 1
      else
         nout_d1 = ns2d_d1
         nout_d2 = ns2d_d2
      endif

      do n=1,nout_d1
        R4OUT(:) = SAVE2D_D1(:,n)
        write (17) R4OUT
      enddo
  
      do n=1,nout_d2
        R4OUT(:) = SAVE2D_D2(:,n)
        write (17) R4OUT
      enddo

C     ----------------------------------------------------------------

C--   3. Reset arrays to zero for the next daily-mean

  300 continue

      if (iitest.eq.1) print*,' reset to zero'

      SAVE2D_D1(:,:) = 0.

      SAVE2D_D2(:,:) = 0.

C--
      if (iitest.eq.1) print *, 'end of DMOUT'

      RETURN
      END
