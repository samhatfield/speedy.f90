cfk#if !defined(KNMI)
      SUBROUTINE CONVMF (PSA,SE,QA,QSAT,
     *                   ITOP,CBMF,PRECNV,DFSE,DFQA)
cfk#else
cfk      SUBROUTINE CONVMF (PSA,SE,QA,QSAT,TS,
cfk     *                   ITOP,CBMF,PRECNV,SNOWCV,DFSE,DFQA)
cfk#endif
C--
C--   SUBROUTINE CONVMF (PSA,SE,QA,QSAT,
C--  *                   ITOP,CBMF,PRECNV,DFSE,DFQA)
C--
C--   Purpose: Compute convective fluxes of dry static energy and moisture
C--            using a simplified mass-flux scheme
C--   Input:   PSA    = norm. surface pressure [p/p0]            (2-dim)
C--            SE     = dry static energy                        (3-dim)
C--            QA     = specific humidity [g/kg]                 (3-dim)
C--            QSAT   = saturation spec. hum. [g/kg]             (3-dim)
C--   Output:  ITOP   = top of convection (layer index)          (2-dim)
C--            CBMF   = cloud-base mass flux                     (2-dim)
C--            PRECNV = convective precipitation [g/(m^2 s)]     (2-dim)
C--            DFSE   = net flux of d.s.en. into each atm. layer (3-dim)
C--            DFQA   = net flux of sp.hum. into each atm. layer (3-dim)
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

      include "com_physcon.h"

C     Convection constants

      include "com_cnvcon.h"

      REAL PSA(NGP), SE(NGP,NLEV), QA(NGP,NLEV), QSAT(NGP,NLEV)

      INTEGER ITOP(NGP)
      REAL CBMF(NGP), PRECNV(NGP), DFSE(NGP,NLEV), DFQA(NGP,NLEV)
cfk#if defined(KNMI)
cfk      REAL ts(ngp),snowcv(ngp)
cfk#endif

      REAL MSS(NGP,2:NLEV), MSE0, MSE1, MSS0, MSS2, MSTHR, 
     &     QDIF(NGP), ENTR(2:NLEV-1)
      LOGICAL LQTHR
C
C--   1. Initialization of output and workspace arrays
C
      NL1=NLEV-1
      NLP=NLEV+1
      FQMAX=5.
 
      FM0=P0*DSIG(NLEV)/(GG*TRCNV*3600)
      RDPS=2./(1.-PSMIN)

C      Used in exp 566 to 604:
C      psmin=0.8
C      rdps=1./(1.-psmin)
   
      DO K=1,NLEV
        DO J=1,NGP
          DFSE(J,K)=0.0
          DFQA(J,K)=0.0
        ENDDO
      ENDDO

      DO J=1,NGP
        CBMF(J)  =0.0
        PRECNV(J)=0.0
cfk#if defined(KNMI)
cfk        SNOWCV(J)=0.0
cfk#endif
      ENDDO

C     Saturation moist static energy
      DO K=2,NLEV
        DO J=1,NGP
cfk#if !defined(KNMI)
          MSS(J,K)=SE(J,K)+ALHC*QSAT(J,K)
cfk#else
cfk          IF (TS(J).GT.273.15) THEN
cfk            MSS(J,K)=SE(J,K)+ALHC*QSAT(J,K)
cfk          ELSE
cfk            MSS(J,K)=SE(J,K)+ALHS*QSAT(J,K)
cfk          ENDIF
cfk#endif
        ENDDO
      ENDDO

C     Entrainment profile (up to sigma = 0.5)
      SENTR=0.
      DO K=2,NL1
        ENTR(K)=(MAX(0.,SIG(K)-0.5))**2
        SENTR=SENTR+ENTR(K)
      ENDDO

      SENTR=ENTMAX/SENTR
      DO K=2,NL1
        ENTR(K)=ENTR(K)*SENTR
      ENDDO
C
C--   2. Check of conditions for convection
C
      RLHC=1./ALHC

      DO J=1,NGP

        ITOP(J)=NLP

        IF (PSA(J).GT.PSMIN) THEN

C         Minimum of moist static energy in the lowest two levels
          MSE0=SE(J,NLEV)+ALHC*QA(J,NLEV)
          MSE1=SE(J,NL1) +ALHC*QA(J,NL1)
          MSE1=MIN(MSE0,MSE1)

C         Saturation (or super-saturated) moist static energy in PBL 
          MSS0=MAX(MSE0,MSS(J,NLEV))

          KTOP1=NLEV
          KTOP2=NLEV

          DO K=NLEV-3,3,-1

            MSS2=MSS(J,K)+WVI(K,2)*(MSS(J,K+1)-MSS(J,K))

C           Check 1: conditional instability 
C                    (MSS in PBL > MSS at top level)
            IF (MSS0.GT.MSS2) THEN
               KTOP1=K
            ENDIF

C           Check 2: gradient of actual moist static energy 
C                    between lower and upper troposphere                     
            IF (MSE1.GT.MSS2) THEN
               KTOP2=K
               MSTHR=MSS2
            ENDIF

          ENDDO

          IF (KTOP1.LT.NLEV) THEN

C           Check 3: RH > RH_c at both k=NLEV and k=NL1
            QTHR0=RHBL*QSAT(J,NLEV)
            QTHR1=RHBL*QSAT(J,NL1)
            LQTHR=(QA(J,NLEV).GT.QTHR0.AND.QA(J,NL1).GT.QTHR1)

            IF (KTOP2.LT.NLEV) THEN
               ITOP(J)=KTOP1
               QDIF(J)=MAX(QA(J,NLEV)-QTHR0,(MSE0-MSTHR)*RLHC)
            ELSE IF (LQTHR) THEN
               ITOP(J)=KTOP1
               QDIF(J)=QA(J,NLEV)-QTHR0
            ENDIF

          ENDIF

        ENDIF

      ENDDO
C
C--   3. Convection over selected grid-points
C
      DO 300 J=1,NGP
      IF (ITOP(J).EQ.NLP) GO TO 300

C       3.1 Boundary layer (cloud base)

        K =NLEV
        K1=K-1

C       Maximum specific humidity in the PBL
        QMAX=MAX(1.01*QA(J,K),QSAT(J,K))

C       Dry static energy and moisture at upper boundary
        SB=SE(J,K1)+WVI(K1,2)*(SE(J,K)-SE(J,K1))
        QB=QA(J,K1)+WVI(K1,2)*(QA(J,K)-QA(J,K1))
        QB=MIN(QB,QA(J,K))

C       Cloud-base mass flux, computed to satisfy:
C       fmass*(qmax-qb)*(g/dp)=qdif/trcnv
        FPSA=PSA(J)*MIN(1.,(PSA(J)-PSMIN)*RDPS)
        FMASS=FM0*FPSA*MIN(FQMAX,QDIF(J)/(QMAX-QB))
        CBMF(J)=FMASS

C       Upward fluxes at upper boundary
        FUS=FMASS*SE(J,K)
        FUQ=FMASS*QMAX

C       Downward fluxes at upper boundary
        FDS=FMASS*SB
        FDQ=FMASS*QB

C       Net flux of dry static energy and moisture
        DFSE(J,K)=FDS-FUS
        DFQA(J,K)=FDQ-FUQ

C       3.2 Intermediate layers (entrainment)

        DO K=NLEV-1,ITOP(J)+1,-1
        K1=K-1

C         Fluxes at lower boundary
          DFSE(J,K)=FUS-FDS
          DFQA(J,K)=FUQ-FDQ

C         Mass entrainment
          ENMASS=ENTR(K)*PSA(J)*CBMF(J)
          FMASS=FMASS+ENMASS

C         Upward fluxes at upper boundary
          FUS=FUS+ENMASS*SE(J,K)
          FUQ=FUQ+ENMASS*QA(J,K)

C         Downward fluxes at upper boundary
          SB=SE(J,K1)+WVI(K1,2)*(SE(J,K)-SE(J,K1))
          QB=QA(J,K1)+WVI(K1,2)*(QA(J,K)-QA(J,K1))
          FDS=FMASS*SB
          FDQ=FMASS*QB

C         Net flux of dry static energy and moisture
          DFSE(J,K)=DFSE(J,K)+FDS-FUS
          DFQA(J,K)=DFQA(J,K)+FDQ-FUQ

C         Secondary moisture flux
          DELQ=RHIL*QSAT(J,K)-QA(J,K)
          IF (DELQ.GT.0.0) THEN
            FSQ=SMF*CBMF(J)*DELQ
            DFQA(J,K)   =DFQA(J,K)   +FSQ 
            DFQA(J,NLEV)=DFQA(J,NLEV)-FSQ
          ENDIF

        ENDDO

C       3.3 Top layer (condensation and detrainment)

        K=ITOP(J)
C
C       Flux of convective precipitation
        QSATB=QSAT(J,K)+WVI(K,2)*(QSAT(J,K+1)-QSAT(J,K))
cfk#if !defined(KNMI)
        PRECNV(J)=MAX(FUQ-FMASS*QSATB,0.0)

C       Net flux of dry static energy and moisture
        DFSE(J,K)=FUS-FDS+ALHC*PRECNV(J)
        DFQA(J,K)=FUQ-FDQ-PRECNV(J)
cfk#else
cfk        IF (TS(J).GT.273.15) THEN
cfk           PRECNV(J)=MAX(FUQ-FMASS*QSATB,0.0)
cfkC          Net flux of dry static energy and moisture
cfk           DFSE(J,K)=FUS-FDS+ALHC*PRECNV(J)
cfk           DFQA(J,K)=FUQ-FDQ-PRECNV(J)
cfk        ELSE
cfk           SNOWCV(J)=MAX(FUQ-FMASS*QSATB,0.0)
cfk           PRECNV(J)=SNOWCV(J)
cfkC          Net flux of dry static energy and moisture
cfk           DFSE(J,K)=FUS-FDS+ALHS*SNOWCV(J)
cfk           DFQA(J,K)=FUQ-FDQ-SNOWCV(J)
cfk        ENDIF
cfk#endif

 300  CONTINUE

      RETURN
      END
