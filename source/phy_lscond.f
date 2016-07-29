cfk#if !defined(KNMI)
      SUBROUTINE LSCOND (PSA,QA,QSAT,
     *                   ITOP,PRECLS,DTLSC,DQLSC)
cfk#else
cfk      SUBROUTINE LSCOND (PSA,QA,QSAT,TS,
cfk     *                   ITOP,PRECLS,SNOWLS,DTLSC,DQLSC)
cfk#endif
C--
C--   SUBROUTINE LSCOND (PSA,QA,QSAT,
C--  *                   ITOP,PRECLS,DTLSC,DQLSC) 
C--
C--   Purpose: Compute large-scale precipitation and
C--            associated tendencies of temperature and moisture
C--   Input:   PSA    = norm. surface pressure [p/p0]           (2-dim)
C--            QA     = specific humidity [g/kg]                (3-dim)
C--            QSAT   = saturation spec. hum. [g/kg]            (3-dim)
C--            ITOP   = top of convection (layer index)         (2-dim)
C--   Output:  ITOP   = top of conv+l.s.condensat.(layer index) (2-dim)
C--            PRECLS = large-scale precipitation [g/(m^2 s)]   (2-dim)
C--            DTLSC  = temperature tendency from l.s. cond     (3-dim)
C--            DQLSC  = hum. tendency [g/(kg s)] from l.s. cond (3-dim)
C--
C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

      include "com_physcon.h"

C     Large-scale condensation constants

      include "com_lsccon.h"

      REAL PSA(NGP), QA(NGP,NLEV), QSAT(NGP,NLEV)

      INTEGER ITOP(NGP)
      REAL PRECLS(NGP), DTLSC(NGP,NLEV), DQLSC(NGP,NLEV)
cfk#if defined(KNMI)
cfk      REAL ts(ngp),snowls(ngp)
cfk#endif

      REAL PSA2(NGP)

C--   1. Initialization

      QSMAX = 10.

      RTLSC = 1./(TRLSC*3600.)
      TFACT = ALHC/CP
cfk#if defined(KNMI)
cfk      TFACTS= ALHS/CP
cfk#endif
      PRG = P0/GG

      DO J=1,NGP
        DTLSC(J,1) = 0.
        DQLSC(J,1) = 0.
        PRECLS(J)  = 0.
cfk#if defined(KNMI)
cfk        SNOWLS(J)  = 0.
cfk#endif
        PSA2(J)    = PSA(J)*PSA(J)
      ENDDO

C--   2. Tendencies of temperature and moisture
C        NB. A maximum heating rate is imposed to avoid 
C            grid-point-storm instability 

      DO K=2,NLEV

        SIG2=SIG(K)*SIG(K)
        RHREF = RHLSC+DRHLSC*(SIG2-1.)
        IF (K.EQ.NLEV) RHREF = MAX(RHREF,RHBLSC)
        DQMAX = QSMAX*SIG2*RTLSC

        DO J=1,NGP
          DQA = RHREF*QSAT(J,K)-QA(J,K)
          IF (DQA.LT.0.0) THEN
            ITOP(J)    = MIN(K,ITOP(J))
            DQLSC(J,K) = DQA*RTLSC
cfk#if !defined(KNMI)
            DTLSC(J,K) = TFACT*MIN(-DQLSC(J,K),DQMAX*PSA2(J))
cfk#else
cfk            IF (TS(J).GT.273.15) THEN
cfk              DTLSC(J,K) = TFACT*MIN(-DQLSC(J,K),DQMAX*PSA2(J))
cfk            ELSE
cfk              DTLSC(J,K) = TFACTS*MIN(-DQLSC(J,K),DQMAX*PSA2(J))
cfk            ENDIF
cfk#endif
          ELSE
            DQLSC(J,K) = 0.
            DTLSC(J,K) = 0.
          ENDIF
        ENDDO

      ENDDO

C--   3. Large-scale precipitation

      DO K=2,NLEV
        PFACT = DSIG(K)*PRG
        DO J=1,NGP
          PRECLS(J) = PRECLS(J)-PFACT*DQLSC(J,K)
cfk#if defined(KNMI)
cfk          IF (TS(J).LT.273.15) THEN
cfk             SNOWLS(J) = SNOWLS(J)-PFACT*DQLSC(J,K)
cfk          ENDIF
cfk#endif
        ENDDO
      ENDDO

      DO J=1,NGP
        PRECLS(J) = PRECLS(J)*PSA(J)
cfk#if defined(KNMI)
cfk        SNOWLS(J) = SNOWLS(J)*PSA(J)
cfk#endif
      ENDDO

C--
      RETURN
      END
