      SUBROUTINE STLOOP (ISTEP)
C--
C--   SUBROUTINE STLOOP (ISTEP)
C--
C--   Purpose: Perform a series of time steps calling 
C--            post-processing/output routines at selected steps
C--   Input/output : ISTEP = time step index
C--   Updated common block : LFLAG2
C-- 

      include "com_tsteps.h"
      include "com_lflags.h"
 
      iitest=0

      DO J=1,NSTEPS

        if (iitest.eq.1) print*, 'STLOOP: calling step ', istep

C       Set logical flags

        LRADSW = (MOD(ISTEP,NSTRAD).EQ.1)
        LRANDF = ((ISTEP.LE.NSTRDF).OR.(NSTRDF.LT.0))

C       Perform one leapfrog time step

        CALL STEP (2,2,DELT2,ALPH,ROB,WIL)   

C       Do diagnostic, post-processing and I/O tasks 
 
        CALL DIAGNS (2,ISTEP)

        IF (MOD(ISTEP,NSTPPR).EQ.0) CALL TMINC

        IF (NSTOUT.GT.0.AND.MOD(ISTEP,NSTOUT).EQ.0) CALL TMOUT (1)

        ISTEP=ISTEP+1

      ENDDO

      RETURN
      END
  
