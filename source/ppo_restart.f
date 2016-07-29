
      SUBROUTINE RESTART (JDAY)
C--
C--   SUBROUTINE RESTART (JDAY)
C--
C--   Purpose : read or write a restart file
C--   Input :   JDAY  = 0 : read model variables from a restart file
C--                   > 0 : write model variables  to a restart file
C--                         at selected dates and at the end of run 
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_date.h"
      include "com_tsteps.h"

      include "com_dynvar.h"

      IF (JDAY.EQ.0) THEN

C--   1. Read the restart dataset corresponding to the specified initial date

  100    CONTINUE
         READ (3,END=200) IYEAR, IMONTH

         IF (IYEAR.EQ.IYEAR0.AND.IMONTH.EQ.IMONT0) THEN

           print*, 'Read restart dataset for year/month: ', IYEAR,IMONTH
           
           READ (3) VOR
           READ (3) DIV
           READ (3) T
           READ (3) PS
           READ (3) TR

           CALL REST_LAND (0)

           CALL REST_SEA (0)

         ELSE
									
           print*, 'Skip restart dataset for year/month: ', IYEAR,IMONTH
           
           DO JREC=1,5
             READ (3) ADUMMY
           ENDDO

           CALL REST_LAND (0)

           CALL REST_SEA (0)

           GO TO 100

         ENDIF

C     Check for write-up dates

      ELSE IF ( (IDAY.eq.1) .and.
     &          (mod(IMONTH-1,NMONRS).eq.0.or.JDAY.eq.NDAYTOT) ) THEN

C--   2. Write date and model variables to the restart file

         print*, 'Write restart dataset for year/month: ', IYEAR,IMONTH

         WRITE (10) IYEAR, IMONTH

         WRITE (10) VOR
         WRITE (10) DIV
         WRITE (10) T
         WRITE (10) PS
         WRITE (10) TR

         CALL REST_LAND (1)

         CALL REST_SEA (1)

      ENDIF
C--
      RETURN

C--   4. Stop integration if restart file is not found

  200 CONTINUE

      print*, ' No restart dataset for the specified initial date'

      STOP 'invalid restart'

C--
      END

