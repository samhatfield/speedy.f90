
      SUBROUTINE AGCM_INIT (cexp,inidate,ntimes,irstart,
     &                      ndays)
C--
C--   SUBROUTINE AGCM_INIT (cexp,inidate,ntimes,irstart,
C--  &                      ndays)
C--
C--   Purpose: Initialization of atmos. model and coupling interface 
C--
c      include "atparam.h"
c      include "atparam1.h"

c      parameter ( NGP = IX*IL )

      include "com_tsteps.h"
      include "com_date.h"

      include "com_lflags.h"
      include "com_cpl_flags.h"


C     Input (reset by input/include files if inidate = 0):
      CHARACTER*3 cexp        ! experiment identifier
      INTEGER     inidate     ! initial date YYYYMM
      INTEGER     ntimes      ! integr. length in months (< 0) or days (> 0)
      INTEGER     irstart     ! restart flag: 0 = no, > 0 = yes

C     Output:
      INTEGER     ndays       ! total no. of integration days

      PRINT *, ' Hallo from Speedy_AGCM'

C--   1. Set run initial time, duration, time-stepping and coupling options

      if (inidate.le.0) then
         read (2,*) istart
         read (2,'(a3)') cexp
      else
         istart = irstart
      endif

      if (istart.ne.0) istart = 1

      include "cls_instep.h"

      if (inidate.gt.0) then

         iyear0 = inidate/100
         imont0 = mod(inidate,100)

         isst0 = (iyear0-issty0)*12+imont0

         if (ntimes.lt.0) then
            nmonts = -ntimes
            ndaysl = 0
         else
            nmonts = 0
            ndaysl = ntimes
         endif

      endif

      CALL NEWDATE (0)

      ndays = ndaytot

C     check consistency of coupling and prescribed SST anomaly flags
      if (icsea.ge.4) isstan = 1

C--   2. Initialization of atmospheric model constants and variables 

      CALL INI_ATM (cexp)

C--   3. Initialization of coupled modules (land, sea, ice)
 
      CALL INI_COUPLER (istart)

C--   4. Set up the forcing fields for the first time step

      CALL FORDATE (0)

C--   5. Do the initial (2nd-order) time step, initialize the semi-impl. scheme

      CALL STEPONE

C--
      RETURN
      END


      subroutine NEWDATE (imode)
C--
C--   subroutine NEWDATE (imode)
C--   Purpose:   initilialize and update date variables 
C--   Input :    imode = 0 for initialization, > 0 for update  

      parameter ( NCAL=365 )
      
      include "com_date.h"
      include "com_tsteps.h"

C     365-day calendar
      integer ncal365(12)
      data    ncal365/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      if (imode.le.0) then

C        Calendar

         if (ncal.eq.365) then
            ndaycal(:,1) = ncal365(:)
         else
            ndaycal(:,1) = 30
         endif
 
         ndaycal(1,2) = 0
         do jm=2,12
            ndaycal(jm,2) = ndaycal(jm-1,1)+ndaycal(jm-1,2)
         enddo

C        Total no. of integration days

         ndaytot = ndaysl
         im = imont0
      
         do jm=1,nmonts
            ndaytot = ndaytot+ndaycal(im,1)
            im = im+1
            if (im.eq.13) im=1
         enddo

C        Initial date

         iyear  = iyear0
         imonth = imont0
         iday   = 1

      else

C        Set new date

         iday = iday+1

         if (iday.gt.ndaycal(imonth,1)) then
            iday   = 1
            imonth = imonth+1
         endif

         if (imonth.gt.12) then
            imonth = 1
            iyear  = iyear+1
         endif

      endif

C     Additional variables to define forcing terms and boundary cond.

      if (iseasc.ge.1) then

         imont1 = imonth
         tmonth = (iday-0.5)/float(ndaycal(imonth,1))
         tyear  = (ndaycal(imonth,2)+iday-0.5)/float(ncal)

      else

         imont1 = imont0
         tmonth = 0.5
         tyear  = (ndaycal(imont1,2)+0.5*ndaycal(imont1,2))/float(ncal)

      endif

      return
      end



