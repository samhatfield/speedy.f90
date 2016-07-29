
      SUBROUTINE INVARS 
C--
C--   SUBROUTINE INVARS (ISTART)
C--
C--   Purpose : initialize all spectral variables starting from
C--             either a reference atmosphere or a restart file
C--   Input :   ISTART = 0 : reference atmosphere (at rest)
C--                    = 1 : restart file
C--   Initialized common blocks : DATE1, DYNSP1, DYNSP2 (PHIS only),
C--                               SFCANOM, SFCFLUX
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_date.h"
      include "com_tsteps.h"

      include "com_surfcon.h"

      include "com_dyncon0.h"
      include "com_dyncon1.h"

      include "com_dynvar.h"


      COMPLEX ZERO, CCON, SURFS(MX,NX)
      REAL  SURFG(IX,IL)

      GAM1 = GAMMA/(1000.*GRAV)
      ZERO = (0.,0.)
      CCON = (1.,0.)*SQRT(2.)


C--   1. Compute spectral surface geopotential

      CALL SPEC (PHI0,PHIS)
      IF (IX.EQ.IY*4) CALL TRUNCT (PHIS)

      CALL GRID (PHIS,PHIS0,1)


      IF (ISTART.EQ.0) THEN

C--   2. Start from reference atmosphere (at rest) 

        print*, ' starting from rest'

        IYEAR  = IYEAR0
        IMONTH = IMONT0

C       2.1 Set vorticity, divergence and tracers to zero

        DO K=1,KX
         DO N=1,NX
          DO M=1,MX
            VOR(M,N,K,1)=ZERO
            DIV(M,N,K,1)=ZERO
          ENDDO
         ENDDO
        ENDDO

        DO ITR=1,NTR
         DO K=1,KX
          DO N=1,NX
           DO M=1,MX
             TR(M,N,K,1,ITR)=ZERO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C       2.2 Set reference temperature :
C           tropos:  T = 288 degK at z = 0, constant lapse rate
C           stratos: T = 216 degK, lapse rate = 0

        TREF  = 288.
        TTOP  = 216.
        GAM2  = GAM1/TREF
        RGAM  = RGAS*GAM1
        RGAMR = 1./RGAM

C       Surface and stratospheric air temperature

        DO N=1,NX
         DO M=1,MX
           T(M,N,1,1)=ZERO
	   T(M,N,2,1)=ZERO
           SURFS(M,N)=-GAM1*PHIS(M,N)
         ENDDO
        ENDDO

        T(1,1,1,1)=CCON*TTOP
        T(1,1,2,1)=CCON*TTOP
        SURFS(1,1)=CCON*TREF-GAM1*PHIS(1,1)

C       Temperature at tropospheric levels
        DO K=3,KX
          FACTK=FSG(K)**RGAM
          DO N=1,NX
           DO M=1,MX
             T(M,N,K,1)=SURFS(M,N)*FACTK
           ENDDO
          ENDDO
        ENDDO

C       2.3 Set log(ps) consistent with temperature profile
C           p_ref = 1013 hPa at z = 0   

        RLOG0=LOG(1.013)
        DO J=1,IL
         DO I=1,IX
           SURFG(I,J)=RLOG0+RGAMR*LOG(1.-GAM2*PHIS0(I,J))
         ENDDO
        ENDDO

        CALL SPEC (SURFG,PS)
        IF (IX.EQ.IY*4) CALL TRUNCT (PS)

C       2.4 Set tropospheric spec. humidity in g/kg
C           Qref = RHref * Qsat(288K, 1013hPa)

        ESREF=17.
        QREF=REFRH1*0.622*ESREF
        QEXP=HSCALE/HSHUM
        
C       Spec. humidity at the surface 

        DO J=1,IL
         DO I=1,IX
           SURFG(I,J)=QREF*EXP(QEXP*SURFG(I,J))
         ENDDO
         print *, ' Q0 jlat = ', j, surfg(1,j)
        ENDDO

        CALL SPEC (SURFG,SURFS)
        IF (IX.EQ.IY*4) CALL TRUNCT (SURFS)

C       Spec. humidity at tropospheric levels      

        DO K=3,KX
          FACTK=FSG(K)**QEXP
          print *, 'vertical scale factor at level ', k, factk
          DO N=1,NX
           DO M=1,MX
             TR(M,N,K,1,1)=SURFS(M,N)*FACTK
           ENDDO
          ENDDO
        ENDDO

C       Print diagnostics from initial conditions
  
        CALL DIAGNS (1,0)

      ELSE

C--   3. Start from restart file 

        print*,' reading a restart file'

        CALL RESTART (0)

C       Print diagnostics from initial conditions
  
        CALL DIAGNS (2,0)

      ENDIF

C--
      RETURN
      END

