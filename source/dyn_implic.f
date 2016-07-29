      SUBROUTINE IMPLIC (DIVDT,TDT,PSDT)
C--
C--   SUBROUTINE IMPLIC (DIVDT,TDT,PSDT)
C--
C--   Purpose : Correct tendencies for implicit gravity wave model
C--   Input/output : DIVDT = divergence tendency
C--                  TDT   = temperature tendency
C--                  PSDT  = tendency of log(surf.pressure)
C--
      include "atparam.h"
      include "atparam1.h"
      PARAMETER (MXNXKX=MX*NX*KX)

      include "com_dyncon1.h"
      include "com_dyncon2.h"

      COMPLEX DIVDT(MX,NX,KX),TDT(MX,NX,KX),PSDT(MX,NX)
      COMPLEX YE(MX,NX,KX),YF(MX,NX,KX),ZERO

      ZERO=(0.,0.)

      DO 1 K=1,KX
      DO 1 N=1,NX
      DO 1 M=1,MX
        YE(M,N,K)=ZERO
    1 CONTINUE

      DO 2 K1=1,KX
      DO 2 K=1,KX
        DO  M=1,MX
          DO  N=1,NX
            YE(M,N,K)=YE(M,N,K)+XD(K,K1)*TDT(M,N,K1)
          ENDDO
        ENDDO
    2 CONTINUE

      DO 21 K=1,KX
      DO M=1,MX
       DO N=1,NX
        YE(M,N,K)=YE(M,N,K)+TREF1(K)*PSDT(M,N)
       ENDDO
      ENDDO
   21 CONTINUE

      DO 4 K=1,KX
        DO M=1,MX
          DO N=1,NX
           YF(M,N,K)=DIVDT(M,N,K)+ELZ(M,N)*YE(M,N,K)
          ENDDO
        ENDDO
    4 CONTINUE

      DO M=1,MX
       DO N=1,NX 
        DO K=1,KX
          DIVDT(M,N,K)=ZERO
        ENDDO
       ENDDO
      ENDDO

      DO 6 N=1,NX
      DO 6 M=1,MX
        MM=ISC*(M-1)+1
        LL=MM+N-2
        IF(LL.NE.0) THEN
          DO 66 K1=1,KX
          DO 66 K=1,KX
            DIVDT(M,N,K)=DIVDT(M,N,K)+XJ(K,K1,LL)*YF(M,N,K1)
   66     CONTINUE
        ENDIF
    6 CONTINUE

      DO 46 K=1,KX
        DO  M=1,MX
          DO N=1,NX
           PSDT(M,N)=PSDT(M,N)-DIVDT(M,N,K)*DHSX(K)
          ENDDO
        ENDDO
   46 CONTINUE

      DO 7 K=1,KX
      DO 7 K1=1,KX
        DO M=1,MX
         DO N=1,NX
           TDT(M,N,K)=TDT(M,N,K)+XC(K,K1)*DIVDT(M,N,K1)
         ENDDO
        ENDDO
    7 CONTINUE

      RETURN
      END
