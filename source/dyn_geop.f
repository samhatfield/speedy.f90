      SUBROUTINE GEOP (JJ)
C--
C--   SUBROUTINE GEOP (JJ)
C--
C--   Purpose : compute spectral geopotential from spectral temperature T
C--             and spectral topography PHIS, as in GFDL Climate Group GCM
C--   Input :   JJ = time level index (1 or 2)
C--   Modified common blocks : DYNSP2
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_dyncon1.h"
      include "com_dynvar.h"

C--   1. Bottom layer (integration over half a layer)

      DO M=1,MX
       DO N=1,NX
         PHI(M,N,KX)=PHIS(M,N)+XGEOP1(KX)*T(M,N,KX,JJ)
       ENDDO
      ENDDO

C--   2. Other layers (integration two half-layers)

      DO K=KX-1,1,-1
        DO M=1,MX
         DO N=1,NX
           PHI(M,N,K)=PHI(M,N,K+1)+XGEOP2(K+1)*T(M,N,K+1,JJ)
     &                           +XGEOP1(K)  *T(M,N,K,  JJ)
         ENDDO
        ENDDO
      ENDDO

C     3. lapse-rate correction in the free troposphere

      DO K=2,KX-1
        CORF=XGEOP1(K)*0.5*LOG(HSG(K+1)/FSG(K))/LOG(FSG(K+1)/FSG(K-1))
        DO N=1,NX
          PHI(1,N,K)=PHI(1,N,K)+CORF*(T(1,N,K+1,JJ)-T(1,N,K-1,JJ))
        ENDDO
      ENDDO

C--
      RETURN
      END
