!
!   /DATE1/: date and time variables (updated in NEWDATE)
      INTEGER IYEAR, IMONTH, IDAY, IMONT1
      REAL TMONTH, TYEAR
      COMMON /DATE1/ IYEAR, IMONTH, IDAY, IMONT1, TMONTH, TYEAR

!   /DATE2/: calendar set-up (initialized in NEWDATE)
      INTEGER NDAYCAL, NDAYTOT
      COMMON /DATE2/ NDAYCAL(12,2), NDAYTOT
