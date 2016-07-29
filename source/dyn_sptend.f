      SUBROUTINE SPTEND (DIVDT,TDT,PSDT,J4)
C--
C--   SUBROUTINE SPTEND (DIVDT,TDT,PSDT,J4)
C--
C--   Purpose : compute spectral tendencies of divergence, temperature
C--             and log_surf.pressure)
C--   Input/output : DIVDT = divergence tendency (spec.)
C--                  TDT   = temperature tendency (spec.)
C--                  PSDT  = tendency of log_surf.pressure (spec.)
C--                  J4    = time level index (1 or 2)
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_dyncon1.h"
      include "com_dyncon2.h"
      include "com_dynvar.h"

      COMPLEX DIVDT(MX,NX,KX), TDT(MX,NX,KX), PSDT(MX,NX)

      COMPLEX DUMK(MX,NX,KXP),DMEANC(MX,NX),SIGDTC(MX,NX,KXP)
      COMPLEX TEMPC(MX,NX,3)
      COMPLEX DUMC(MX,NX,2),ZERO

      ZERO=(0.,0.)

c   vertical mean div and pressure tendency

      DO M=1,MX 
       DO N=1,NX
        DMEANC(M,N)=ZERO
       ENDDO
      ENDDO

      DO 53 K=1,KX
        DO M=1,MX
         DO N=1,NX
           DMEANC(M,N)=DMEANC(M,N)+DIV(M,N,K,J4)*DHS(K)
         ENDDO
        ENDDO
   53 CONTINUE

      DO M=1,MX
       DO N=1,NX
        PSDT(M,N)=PSDT(M,N)-DMEANC(M,N)
       ENDDO
      ENDDO
      PSDT(1,1)=ZERO

c  sigma-dot "velocity" and temperature tendency

      DO M=1,MX 
       DO N=1,NX
        SIGDTC(M,N,1)=ZERO
        SIGDTC(M,N,KXP)=ZERO
       ENDDO
      ENDDO

      DO 217 K=1,KXM
        DO M=1,MX
         DO N=1,NX
          SIGDTC(M,N,K+1)=SIGDTC(M,N,K)
     *     -DHS(K)*(DIV(M,N,K,J4)-DMEANC(M,N))
         ENDDO
        ENDDO
  217 CONTINUE

      DO M=1,MX
       DO N=1,NX
        DUMK(M,N,1)=ZERO
        DUMK(M,N,KXP)=ZERO
       ENDDO
      ENDDO

      DO 5702 K=2,KX
      DO M=1,MX
        DO N=1,NX
         DUMK(M,N,K)=SIGDTC(M,N,K)*(TREF(K)-TREF(K-1))
        ENDDO
      ENDDO
 5702 CONTINUE

      DO 5502 K=1,KX
      DO M=1,MX
       DO N=1,NX
        TDT(M,N,K)= TDT(M,N,K)-(DUMK(M,N,K+1)+DUMK(M,N,K))*DHSR(K)
     *     +TREF3(K)*(SIGDTC(M,N,K+1)+SIGDTC(M,N,K))
     *     -TREF2(K)*DMEANC(M,N)
       ENDDO
      ENDDO
 5502 CONTINUE

c   geopotential and divergence tendency

      CALL GEOP(J4)

      DO 18 K=1,KX
        DO M=1,MX
         DO N=1,NX
          DUMC(M,N,1)=PHI(M,N,K)+RGAS*TREF(K)*PS(M,N,J4)
         ENDDO
        ENDDO
        CALL LAP(DUMC(1,1,1),DUMC(1,1,2))
        DO  M=1,MX
          DO N=1,NX
           DIVDT(M,N,K)=DIVDT(M,N,K)-DUMC(M,N,2)
          ENDDO
        ENDDO
   18 CONTINUE

      RETURN
      END   
