!
!  /ISTEPS/: Integer time-stepping constants (initial. in INDYNS)
!   NMONTS = Integration length in months
!   NDAYSL = No. of days in the last month of int. (max=30)
!   NSTEPS = No. of time steps in one day
!   NSTDIA = Period (no. of steps) for diagnostic print-out
!   NSTPPR = Period (no. of steps) for post-processing 
!   NSTOUT = Period (no. of steps) for time-mean output
!   IDOUT  = daily output flag 
!            (0=no, 1=basic (Z500,PREC,MSLP,TEMP0), 2=full)
!   NMONRS = Period (no. of months) for restart file update
!   ISEASC = Seasonal cycle flag (0=no, 1=yes)
!   ISTART = Start flag (0: from rest, 1: from restart file)
!   IYEAR0 = Year of initial date (4-digit, eg 1900)
!   IMONT0 = Month of initial date (1 to 12)

      INTEGER NMONTS, NDAYSL, NSTEPS, NSTDIA, NSTPPR, NSTOUT, IDOUT
      INTEGER NMONRS, ISEASC, ISTART, IYEAR0, IMONT0
      COMMON /ISTEPS/ NMONTS, NDAYSL, NSTEPS, NSTDIA, NSTPPR, NSTOUT,   &
     & IDOUT, NMONRS, ISEASC, ISTART, IYEAR0, IMONT0

!
!  /ISTFOR/: Integer forcing indices (initial. in INDYNS)
!   NSTRAD = Period (no. of steps) for shortwave radiation 
!   NSTRDF = Duration of random diabatic forcing ( 0 : no forcing, 
!            > 0 : no. of initial steps, < 0 : whole integration)									  
!   INDRDF = Initialization index for random diabatic forcing
!   ISST0  = record in SST anomaly file corr. to the initial month     

      INTEGER NSTRAD, NSTRDF, INDRDF, ISST0
      COMMON /ISTFOR/ NSTRAD, NSTRDF, INDRDF, ISST0

!
!   /RSTEPS/: Real time-stepping constants (initial. in INDYNS)
!    DELT   = Time step in seconds
!    DELT2  = 2 * time step in seconds
!    ROB    = Damping factor in Robert time filter
!    WIL    = Parameter of Williams filter
!    ALPH   = Coefficient for semi-implicit computations

      REAL DELT, DELT2, ROB, WIL, ALPH
      COMMON /RSTEPS/ DELT, DELT2, ROB, WIL, ALPH

