
      SUBROUTINE VDIFSC (UA,VA,SE,RH,QA,QSAT,PHI,ICNV,
     &                   UTENVD,VTENVD,TTENVD,QTENVD)
C--
C--   SUBROUTINE VDIFSC (UA,VA,SE,RH,QA,QSAT,PHI,ICNV,
C--  &                   UTENVD,VTENVD,TTENVD,QTENVD)
C--
C--   Purpose: Compute tendencies of momentum, energy and moisture
C--            due to vertical diffusion and shallow convection
C--   Input:   UA     = u-wind                           (3-dim)
C--            VA     = v-wind                           (3-dim)
C--            SE     = dry static energy                (3-dim)
C--            RH     = relative humidity [0-1]          (3-dim)
C--            QA     = specific humidity [g/kg]         (3-dim)
C--            QSAT   = saturation sp. humidity [g/kg]   (3-dim)
C--            PHI    = geopotential                     (3-dim)
C--            ICNV   = index of deep convection         (2-dim)
C--   Output:  UTENVD = u-wind tendency                  (3-dim)
C--            VTENVD = v-wind tendency                  (3-dim)
C--            TTENVD = temperature tendency             (3-dim)
C--            QTENVD = sp. humidity tendency [g/(kg s)] (3-dim)
C--
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

      include "com_physcon.h"

C     Vertical diffusion constants

      include "com_vdicon.h"

      REAL UA(NGP,NLEV), VA(NGP,NLEV), SE(NGP,NLEV),
     &     RH(NGP,NLEV), QA(NGP,NLEV), QSAT(NGP,NLEV),
     &     PHI(NGP,NLEV)

      INTEGER ICNV(NGP)

      REAL UTENVD(NGP,NLEV), VTENVD(NGP,NLEV),
     &     TTENVD(NGP,NLEV), QTENVD(NGP,NLEV)

      REAL RSIG(NLEV), RSIG1(NLEV)


C--   1. Initalization

C     N.B. In this routine, fluxes of dry static energy and humidity
C          are scaled in such a way that:
C          d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma

      NL1  = NLEV-1
      CSHC = DSIG(NLEV)/3600.
      CVDI = (SIGH(NL1)-SIGH(1))/((NL1-1)*3600.)

      FSHCQ  = CSHC/TRSHC
      FSHCSE = CSHC/(TRSHC*CP)

      FVDIQ  = CVDI/TRVDI
      FVDISE = CVDI/(TRVDS*CP)

      DO K=1,NL1
        RSIG(K)=1./DSIG(K)
	RSIG1(K)=1./(1.-SIGH(K))
      ENDDO
      RSIG(NLEV)=1./DSIG(NLEV)
   
      DO K=1,NLEV
        DO J=1,NGP
          UTENVD(J,K) = 0.
          VTENVD(J,K) = 0.
          TTENVD(J,K) = 0.
          QTENVD(J,K) = 0.
        ENDDO
      ENDDO


C--   2. Shallow convection

      DRH0   = RHGRAD*(SIG(NLEV)-SIG(NL1))
      FVDIQ2 = FVDIQ*SIGH(NL1)

      DO J=1,NGP

        DMSE = (SE(J,NLEV)-SE(J,NL1))+ALHC*(QA(J,NLEV)-QSAT(J,NL1))
        DRH  = RH(J,NLEV)-RH(J,NL1)
        FCNV = 1.

        IF (DMSE.GE.0.0) THEN

          IF (ICNV(J).GT.0) FCNV = REDSHC

          FLUXSE         = FCNV*FSHCSE*DMSE
          TTENVD(J,NL1)  = FLUXSE*RSIG(NL1)
          TTENVD(J,NLEV) =-FLUXSE*RSIG(NLEV)

          IF (DRH.GE.0.0) THEN
            FLUXQ          = FCNV*FSHCQ*QSAT(J,NLEV)*DRH
            QTENVD(J,NL1)  = FLUXQ*RSIG(NL1) 
            QTENVD(J,NLEV) =-FLUXQ*RSIG(NLEV)
          ENDIF

        ELSE IF (DRH.GE.DRH0) THEN

          FLUXQ          = FVDIQ2*QSAT(J,NL1)*DRH
          QTENVD(J,NL1)  = FLUXQ*RSIG(NL1) 
          QTENVD(J,NLEV) =-FLUXQ*RSIG(NLEV)

        ENDIF

      ENDDO

C--   3. Vertical diffusion of moisture above the PBL

      DO K=3,NLEV-2

        IF (SIGH(K).GT.0.5) THEN

          DRH0   = RHGRAD*(SIG(K+1)-SIG(K))
          FVDIQ2 = FVDIQ*SIGH(K)

          DO J=1,NGP

            DRH=RH(J,K+1)-RH(J,K)

            IF (DRH.GE.DRH0) THEN
              FLUXQ        = FVDIQ2*QSAT(J,K)*DRH
              QTENVD(J,K)  = QTENVD(J,K)  +FLUXQ*RSIG(K)
              QTENVD(J,K+1)= QTENVD(J,K+1)-FLUXQ*RSIG(K+1)
            ENDIF

          ENDDO

        ENDIF

      ENDDO

C--   4. Damping of super-adiabatic lapse rate

      DO K=1,NL1
       DO J=1,NGP

         SE0 = SE(J,K+1)+SEGRAD*(PHI(J,K)-PHI(J,K+1))

         IF (SE(J,K).LT.SE0) THEN
           FLUXSE      = FVDISE*(SE0-SE(J,K))
           TTENVD(J,K) = TTENVD(J,K)+FLUXSE*RSIG(K)
           DO K1=K+1,NLEV
             TTENVD(J,K1) = TTENVD(J,K1)-FLUXSE*RSIG1(K)
           ENDDO
         ENDIF

       ENDDO
      ENDDO

C--

      RETURN
      END
