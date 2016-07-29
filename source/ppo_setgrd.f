
      SUBROUTINE SETGRD (ind,norun)
C--
C--   SUBROUTINE SETGRD (ind)
C--   Purpose : open and close output files (.grd)
C--
C--   Input : ind = 0 for initialization, 1 otherwise
C--         : norun = run identifier
C--

      include "com_date.h"
      include "com_tsteps.h"

      character*3  norun
      character*16 ofile11, ofile13, ofile15
      character*17 ofile17

      save ofile11, ofile13, ofile15, ofile17

      if (ind.eq.0) then

         ofile11='attmNNN_YYYY.grd'
         ofile13='atvaNNN_YYYY.grd'
         ofile15='atdfNNN_YYYY.grd'

         ofile11(5:7)=norun
         ofile13(5:7)=norun
         ofile15(5:7)=norun

         if (IDOUT .gt. 0) then
            ofile17='daytmNNN_YYYY.grd'
            ofile17(6:8)=norun
         endif

      endif

      write (ofile11(9:12),'(i4)') iyear
      write (ofile13(9:12),'(i4)') iyear
      write (ofile15(9:12),'(i4)') iyear

      if (IDOUT.gt.0) write (ofile17(10:13),'(i4)') iyear

      if (ind.ne.0) then

         close( unit=11 )
         close( unit=13 )
         close( unit=15 )

         if (IDOUT.gt.0) close( unit=17 )

      endif
      
      open ( unit=11, file=ofile11, 
     &       status='new', form='unformatted', access='sequential' )
      open ( unit=13, file=ofile13, 
     &       status='new', form='unformatted', access='sequential' )
      open ( unit=15, file=ofile15, 
     &       status='new', form='unformatted', access='sequential' )

      if (IDOUT.gt.0) then
         open ( unit=17, file=ofile17,
     &          status='new', form='unformatted', access='sequential' )
      endif

      return
      end


