      SUBROUTINE STEPONE
C--
C--   SUBROUTINE STEPONE
C--
C--   Purpose : call initialization of semi-implicit scheme 
C--             and perform initial time step
C--   Initialized common block : LFLAG2
C--

      include "com_tsteps.h"
      include "com_lflags.h"
  
      iitest=1
      if (iitest.eq.1) print*, ' instep: initial time step'

      IF (ISTART.EQ.0) THEN

        DELTH = 0.5*DELT
        LRADSW = .TRUE.
        LRANDF = .FALSE.

        if (iitest.eq.1) print*, ' semi-impl. initialization'
        CALL IMPINT (DELTH,ALPH)

        if (iitest.eq.1) print*, ' forward half-step'
        CALL STEP (1,1,DELTH,ALPH,ROB)

        if (iitest.eq.1) print*, ' semi-impl. initialization'
        CALL IMPINT (DELT,ALPH)

        if (iitest.eq.1) print*, ' leapfrog half-step'
        CALL STEP (1,2,DELT,ALPH,ROB)

      ENDIF

      if (iitest.eq.1) print*, ' semi-impl. initialization'
      CALL IMPINT (DELT2,ALPH)

C-- 
      RETURN
      END
  

