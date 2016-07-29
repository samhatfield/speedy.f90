      SUBROUTINE TMINC
C--
C--   SUBROUTINE TMINC
C--
C--   Purpose : perform post-processing on pressure levels
C--             and increment time-mean arrays
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

C     Constants and conversion factors
      include "com_physcon.h"

C     Model variables, tendencies and fluxes on gaussian grid
      include "com_physvar.h"

C     Surface variables on gaussian grid
      include "com_cli_sea.h"
      include "com_cli_land.h"
      include "com_var_sea.h"
      include "com_var_land.h"
      include "com_radcon.h"


C     Surface properties (time-inv.)
      include "com_surfcon.h"

C     Logical flags
      include "com_lflags.h"

      real ADSAVE(ngp,6), PHISG(ngp), PMSL(ngp), QSATPL(ngp), ST0(ngp)
      equivalence (PHISG,PHIS0)

C     fields for vertical interpolation
      integer K0(ngp)
      real    W0(ngp), ZOUT(ngp), ZINP(nlev), RDZINP(nlev)

C     Level indices for daily-mean storage of upper-air fields
C     upper tropos. u, v :
      khigh = 3
C     mid-tropos. geopotential height :
      kmid  = 5
C     lower tropos. u, v, q :
      klow  = 7

      iitest=0
      if (iitest.eq.1) print *, ' inside TMINC'

      rg    = 1./gg
      rdr2  = 0.5*rd
      gam0  = 0.006*rg
      rgam  = rd*gam0
      rrgam = 1./rgam

C     0. Increment post-processing counter

      rnsave = rnsave +1

C--   1. Store 2-d time-mean fields

C     1.1 Compute additional surface fields

      if (iitest.eq.1) print*,' store 2d fields'

c     mean-sea-level pressure
      do j=1,ngp
        tsg=0.5*(T0(j)+max(255.,min(295.,T0(j))))
        PMSL(j)=PSG(j)*(1.+gam0*PHISG(j)/tsg)**rrgam
      enddo

C     1.2 Increment time-mean arrays

      n0=0
      call ADD1F (SAVE2D_1,PSG,       ngp,n0,1000)
      call ADD1F (SAVE2D_1,PMSL,      ngp,n0,1000)
      call ADD1F (SAVE2D_1,TS,        ngp,n0,1)
      call ADD1F (SAVE2D_1,TSKIN,     ngp,n0,1)
      call ADD1F (SAVE2D_1,SOILW_AM,  ngp,n0,100)
      call ADD1F (SAVE2D_1,ALBSFC,    ngp,n0,100)
      call ADD1F (SAVE2D_1,U0,        ngp,n0,1)
      call ADD1F (SAVE2D_1,V0,        ngp,n0,1)
      call ADD1F (SAVE2D_1,T0,        ngp,n0,1)
      call ADD1F (SAVE2D_1,RH(1,NLEV),ngp,n0,100)
      call ADD1F (SAVE2D_1,CLOUDC,    ngp,n0,100)
      call ADD1F (SAVE2D_1,CLSTR,     ngp,n0,100)
      call ADD1F (SAVE2D_1,CLTOP,     ngp,n0,1000)
      call ADD1F (SAVE2D_1,PRTOP,     ngp,n0,1)

c     land and sea surface temperatures
      call MASKOUT (STL_AM,ST0,BMASK_L,ngp)
      call ADD1F   (SAVE2D_1,ST0,ngp,n0,1)
      call MASKOUT (SST_AM,ST0,BMASK_S,ngp)
      call ADD1F   (SAVE2D_1,ST0,ngp,n0,1)
c     ocean model SST
c      call MASKOUT (SST_OM,ST0,BMASK_S,ngp)
c     ocean model SST + T_ice
      call MASKOUT (SSTI_OM,ST0,BMASK_S,ngp)
      call ADD1F   (SAVE2D_1,ST0,ngp,n0,1)
c     SST anomaly (wrt obs. clim.)
      call MASKOUT (SSTAN_AM,ST0,BMASK_S,ngp)
      call ADD1F   (SAVE2D_1,ST0,ngp,n0,1)

c     NB Fluxes of water, energy and momentum are stored 
c        every time step by subroutine FLUXINC

C     1.3 Increment daily-mean arrays

      do j=1,ngp
         SAVE2D_D1(j,1)=SAVE2D_D1(j,1)+PMSL(j)*1000
         SAVE2D_D1(j,2)=SAVE2D_D1(j,2)+T0(j)
      enddo

C--   2. Perform vertical interpolation from sigma to pressure levels
C--      and increment time-mean fields

      if (iitest.eq.1) print*, ' store 3d fields'

      ZINP(1)  =-SIGL(1)
      do k=2,nlev
        ZINP(k)  =-SIGL(k)
        RDZINP(k)= 1./(ZINP(k-1)-ZINP(k))
      enddo

      zmin  = ZINP(nlev)

      do k=1,kx

C       2.1 Set coefficients for vertical interpolation
C           using coordinate Z = log (p_s/p) 

        if (lppres) then
          plog=log(POUT(k))
          do j=1,ngp
            ZOUT(j)=PSLG1(j)-plog
          enddo
        else
c         Set zout=zinp(k) to do post-proc. on sigma levels
          do j=1,ngp
            ZOUT(j)=ZINP(k)
          enddo
        endif

        call SETVIN (ZINP,RDZINP,ZOUT,ngp,kx,K0,W0)

C       2.2 Interpolate 3-d fields

c       Temperature (extrapolated below the lowest level when W0(j)<0)

        call VERINT (ADSAVE(1,2),TG1,ngp,kx,K0,W0) 

c       Remove extrapolation of temperature inversions 
c       and correct extrap. values using a reference lapse rate

        wref = 0.7

        do j=1,ngp
          if (ZOUT(j).lt.zmin) then
            textr = max(ADSAVE(j,2),TG1(j,nlev))
            aref = rgam*(zmin-ZOUT(j))
            tref = TG1(j,nlev)*(1.+aref+0.5*aref*aref)
            ADSAVE(j,2) = textr+wref*(tref-textr)
          endif
        enddo

c       Geopotential (computed from the closest levels 
c                     using the hydrostatic equation)

        do j=1,ngp
          W0(j)=max(W0(j),0.)
        enddo

        if (lppres) then

          do j=1,ngp
            kj=K0(j)
            kj1=kj-1
            phi1=PHIG1(j,kj)
     &           +rdr2*(ADSAVE(j,2)+TG1(j,kj ))*(ZOUT(j)-ZINP(kj ))
            phi2=PHIG1(j,kj1)
     &           +rdr2*(ADSAVE(j,2)+TG1(j,kj1))*(ZOUT(j)-ZINP(kj1))
            ADSAVE(j,1)=phi1+W0(j)*(phi2-phi1)
          enddo

        else

          call VERINT (ADSAVE(1,1),PHIG1,ngp,kx,K0,W0) 

        endif

c       Wind and relative humidity 

c       a) Interpolate above the lowest level

        call VERINT (ADSAVE(1,3),UG1,ngp,kx,K0,W0) 
        call VERINT (ADSAVE(1,4),VG1,ngp,kx,K0,W0) 
        call VERINT (ADSAVE(1,6),RH, ngp,kx,K0,W0) 

c       b) Decrease wind speed below the lowest level

        do j=1,ngp
          if (ZOUT(j).lt.zmin) then
            fwind=ADSAVE(j,1)/PHIG1(j,nlev)
            ADSAVE(j,3)=ADSAVE(j,3)*fwind
            ADSAVE(j,4)=ADSAVE(j,4)*fwind  
          endif
        enddo

c       Estimate specific humidity using interpolated rel.hum. and
c       sat. spec.hum. at interpolated temperature

        if (lppres) then

          call SHTORH (-1,ngp,ADSAVE(1,2),POUT(k),-1.,
     *                 ADSAVE(1,5),ADSAVE(1,6),QSATPL)

c         Below the surface, set spec.hum. = near-surface value 

          do j=1,ngp
            if (ZOUT(j).lt.0.0) then
              ADSAVE(j,5)=Q0(j)
              ADSAVE(j,6)=Q0(j)/QSATPL(j)
            endif
          enddo

        else

          call VERINT (ADSAVE(1,5),QG1,ngp,kx,K0,W0) 

        endif

c        rescale geopotential and rel. humidity

        do j=1,ngp
          ADSAVE(j,1)=ADSAVE(j,1)*rg
          ADSAVE(j,6)=ADSAVE(j,6)*100.
        enddo

C       2.3 Save upper-air fields

C       2.3.1 Add 3d upper-air fields to time-mean arrays

        do n=1,6
         do j=1,ngp
           SAVE3D(j,k,n)=SAVE3D(j,k,n)+ADSAVE(j,n)
         enddo
        enddo

C       2.3.2 Add upper-air fields fields at selected levels 
C             to daily-means arrays

        if (k.eq.kmid) then

          do j=1,ngp
           SAVE2D_D1(j,3)=SAVE2D_D1(j,3)+ADSAVE(j,1)
          enddo

        endif

        if (k.eq.klow) then

         do j=1,ngp
           SAVE2D_D1(j,4)=SAVE2D_D1(j,4)+ADSAVE(j,3)
           SAVE2D_D1(j,5)=SAVE2D_D1(j,5)+ADSAVE(j,4)
           SAVE2D_D1(j,6)=SAVE2D_D1(j,6)+ADSAVE(j,5)
         enddo

        else if (k.eq.khigh) then

          do j=1,ngp
           SAVE2D_D1(j,7)=SAVE2D_D1(j,7)+ADSAVE(j,3)
           SAVE2D_D1(j,8)=SAVE2D_D1(j,8)+ADSAVE(j,4)
          enddo
        endif


C       2.4 Store variances on pressure levels

        if (ns3d2.gt.0) then

         do n=1,4
          nv=ns3d1+n
          do j=1,ngp
            SAVE3D(j,k,nv)=SAVE3D(j,k,nv)+ADSAVE(j,n)*ADSAVE(j,n)
          enddo
         enddo

         nuv=ns3d1+5
         nvt=ns3d1+6
         do j=1,ngp
           SAVE3D(j,k,nuv)=SAVE3D(j,k,nuv)+ADSAVE(j,3)*ADSAVE(j,4)
           SAVE3D(j,k,nvt)=SAVE3D(j,k,nvt)+ADSAVE(j,2)*ADSAVE(j,4)
         enddo

        endif

C       end-of-loop over pressure levels

      enddo

C--   3. Store diabatic forcing terms on model levels

      if (ns3d3.gt.0) then

        n0=ns3d1+ns3d2
        call ADD1F (SAVE3D,TT_LSC,ngp*nlev,n0,1)
        call ADD1F (SAVE3D,TT_CNV,ngp*nlev,n0,1)
        call ADD1F (SAVE3D,TT_RSW,ngp*nlev,n0,1) 
        call ADD1F (SAVE3D,TT_RLW,ngp*nlev,n0,1) 
        call ADD1F (SAVE3D,TT_PBL,ngp*nlev,n0,1) 

      endif

      if (iitest.eq.1) print *, 'end of TMINC'

      RETURN
      END

      SUBROUTINE SETVIN (ZINP,RDZINP,ZOUT,NGP,NLEV,K0,W0)
C
      INTEGER K0(NGP)
      REAL    ZINP(NLEV), RDZINP(NLEV), ZOUT(NGP), W0(NGP)
C
C *** 1. Select closest vertical levels
C
      DO J=1,NGP
        K0(J)=2
      ENDDO
C
      DO K=2,NLEV-1
       DO J=1,NGP
         IF (ZOUT(J).LT.ZINP(K)) K0(J)=K+1
       ENDDO
      ENDDO
C
C *** 2. Compute interpolation weight
C
      DO J=1,NGP
        W0(J)=(ZOUT(J)-ZINP(K0(J)))*RDZINP(K0(J))
      ENDDO
C
      RETURN
      END

      SUBROUTINE VERINT (F2D,F3D,NGP,NLEV,K0,W0)

C *** 1. Perform vertical interpolation 

      INTEGER K0(NGP)
      REAL    F2D(NGP), F3D(NGP,NLEV), W0(NGP)

      DO J=1,NGP
        F2D(J)=F3D(J,K0(J))+W0(J)*(F3D(J,K0(J)-1)-F3D(J,K0(J)))
      ENDDO
C
      RETURN
      END

      SUBROUTINE ADD1F (FSAVE,FADD,NGP,NF,IFACT)

C *** Add one field to storage array 

      REAL FSAVE(NGP,*), FADD(NGP)

      NF=NF+1

      if (ifact.eq.1) then 
        do j=1,NGP
           FSAVE(j,nf)=FSAVE(j,nf)+FADD(j)
         enddo
      else
        do j=1,NGP
           FSAVE(j,nf)=FSAVE(j,nf)+FADD(j)*ifact
         enddo
      endif
C
      RETURN
      END

      SUBROUTINE MASKOUT (FINP,FOUT,FMASK,NGP)

C *** Set undefined values according to binary land-sea mask

      REAL FINP(NGP), FOUT(NGP), FMASK(NGP) 

      xundef = 9.999E+19

      DO j=1,NGP
         if (FMASK(j).le.0.) then
           FOUT(j) = xundef
         else
           FOUT(j) = FINP(j) 
         endif
      ENDDO

      RETURN
      END

