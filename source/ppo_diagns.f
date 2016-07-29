      SUBROUTINE DIAGNS (JJ,ISTEP)
C--
C--   SUBROUTINE DIAGNS (JJ,ISTEP)
C--
C--   Purpose: print global means of eddy kinetic energy and temperature 
C--   Input : JJ    = time level index (1 or 2)
C--           ISTEP = time step index
C--
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_tsteps.h"

      include "com_dynvar.h"

      COMPLEX TEMP(MX,NX)
      REAL DIAG (KX,3)

C--   1. Get global-mean temperature and compute eddy kinetic energy 

      SQHALF=SQRT(0.5)

      DO K=1,KX

        DIAG(K,1)=0.
        DIAG(K,2)=0.
        DIAG(K,3)=SQHALF*REAL(T(1,1,K,JJ))

        CALL INVLAP (VOR(1,1,K,JJ),TEMP)

        DO M=2,MX
         DO N=1,NX
          DIAG(K,1)=DIAG(K,1)-
     &              REAL(TEMP(M,N)*CONJG(VOR(M,N,K,JJ)))
         ENDDO
        ENDDO

        CALL INVLAP (DIV(1,1,K,JJ),TEMP)

        DO M=2,MX
         DO N=1,NX
          DIAG(K,2)=DIAG(K,2)-
     &              REAL(TEMP(M,N)*CONJG(DIV(M,N,K,JJ)))
         ENDDO
        ENDDO

      ENDDO

C--   2. Print results to screen

      IF (MOD(ISTEP,NSTDIA).EQ.0) THEN
        PRINT 2001, ISTEP, (DIAG(K,1),K=1,KX)
        PRINT 2002,        (DIAG(K,2),K=1,KX)
        PRINT 2003,        (DIAG(K,3),K=1,KX)
      ENDIF

C--   3. Stop integration if model variables are out of range

      DO K=1,KX
        IF (DIAG(K,1).GT.500.OR.
     &      DIAG(K,2).GT.500.OR.
     &      DIAG(K,3).LT.180.OR.
     &      DIAG(K,3).GT.320.) THEN

            PRINT 2001, ISTEP, (DIAG(KK,1),KK=1,KX)
            PRINT 2002,        (DIAG(KK,2),KK=1,KX)
            PRINT 2003,        (DIAG(KK,3),KK=1,KX)

C           Write model fields at t-1 on output file 
            CALL TMOUT (0)
            CALL TMINC
            NSTOUT=NSTPPR
            CALL TMOUT (1)

            STOP '*** Model variables out of accepted range ***'

        ENDIF
      ENDDO

 2001 FORMAT(' step =',i6,' reke =',(10F8.2))
 2002 FORMAT         (13X,' deke =',(10F8.2))
 2003 FORMAT         (13X,' temp =',(10F8.2))

C--
      RETURN
      END
