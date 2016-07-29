
      SUBROUTINE SETCTL (IUNIT,NLON,NLAT,NLEV,NTM,NDTM,
     &                   I3D,N3D,N2D_1,N2D_2,
     &                   RLAT,RLEV,NAME,NORUN,IYEAR0,IMONT0)
C--
C--   Aux. routine SETCTL : write descriptor (.ctl) output file 
C--
      CHARACTER*80 LINE(10), LN3D(30), LN2D_1(20), LN2D_2(15)
      CHARACTER*4  LMON(12), NAME
      CHARACTER*3  NORUN
      CHARACTER*11 CTLNAME
      INTEGER ILEV(30)
      REAL RLAT(NLAT), RLEV(NLEV)
C
C *** 1. Initialization
C
      DATA LMON /'1jan','1feb','1mar','1apr','1may','1jun',
     &           '1jul','1aug','1sep','1oct','1nov','1dec'/

      DATA LN3D/
     &     'GH         n  99  geopotential height               [m]',
     &     'TEMP       n  99  abs. temperature               [degK]',
     &     'U          n  99  zonal (u) wind                  [m/s]',
     &     'V          n  99  meridional (v) wind             [m/s]',
     &     'Q          n  99  specific humidity              [g/Kg]',
     &     'RH         n  99  relative humidity                 [%]',
     &     'OMEGA      n  99  pressure vertical velocity     [Pa/s]',
     &     'PSI        n  99  streamfunction           [10^6 m^2/s]',
     &     'CHI        n  99  velocity potential       [10^6 m^2/s]',
     
     &     'VARGH      n  99  variance of geop. height        [m^2]',
     &     'VART       n  99  variance of temperature      [degK^2]',
     &     'VARU       n  99  variance of u-wind             [J/Kg]',
     &     'VARV       n  99  variance of v-wind             [J/Kg]',
     &     "COVUV      n  99  u'v' covariance (trans.)       [J/Kg]",
     &     "COVVT      n  99  v'T' covariance (trans.)   [degK m/s]",
 
     &     'DTLSC      n  99  dT/dt by large-scale cond. [degK/day]',
     &     'DTCNV      n  99  dT/dt by convection        [degK/day]',
     &     'DTRSW      n  99  dT/dt by shortwave rad.    [degK/day]',
     &     'DTRLW      n  99  dT/dt by longwave  rad.    [degK/day]',
     &     'DTPBL      n  99  dT/dt by PBL processes     [degK/day]',
     &      10*' '/

      DATA LN2D_1/
     &     'SP         0  99  surface pressure                [hPa]',
     &     'MSLP       0  99  mean-sea-level pressure         [hPa]',
     &     'ST         0  99  surface temperature            [degK]',
     &     'SKINT      0  99  skin temperature               [degK]',
     &     'SWAV       0  99  soil wetness availability         [%]',
     &     'ALB        0  99  surface albedo                    [%]',
     &     'U0         0  99  near-surface u-wind             [m/s]',
     &     'V0         0  99  near-surface v-wind             [m/s]',
     &     'TEMP0      0  99  near-surface air temperature   [degK]',
     &     'RH0        0  99  near-surface relative humidity    [%]',
     &     'CLC        0  99  cloud cover (deep clouds)         [%]',
     &     'CLSTR      0  99  cloud cover (strat. clouds)       [%]',
     &     'CLTOP      0  99  pressure at cloud top           [hPa]',
     &     'IPTOP      0  99  highest precipitation level index  []',
     &     'LST        0  99  land-surface temp.             [degK]',
     &     'SST        0  99   sea-surface temp.             [degK]',
     &     'SSTOM      0  99  ocean model sea-surface temp.  [degK]',
     &     'SSTA       0  99  SST anomaly w.r.t. obs. clim.  [degK]',
     &      2*' '/

      DATA LN2D_2/
     &     'PRECLS     0  99  large-scale precipitation    [mm/day]',
     &     'PRECNV     0  99  convective precipitation     [mm/day]',
     &     'EVAP       0  99  evaporation                  [mm/day]',
     &     'USTR       0  99  u-stress                (dw.) [N/m^2]',
     &     'VSTR       0  99  v-stress                (dw.) [N/m^2]',
     &     'TSR        0  99  top shortwave rad.      (dw.) [W/m^2]',
     &     'OLR        0  99  outgoing longwave rad.  (uw.) [W/m^2]',
     &     'SSR        0  99  surface shortwave rad.  (dw.) [W/m^2]',
     &     'SLR        0  99  surface longwave rad.   (uw.) [W/m^2]',
     &     'SHF        0  99  sensible heat flux      (uw.) [W/m^2]',
     &     'LSHF       0  99  heat flux into land sfc (dw.) [W/m^2]',
     &     'SSHF       0  99  heat flux into  sea sfc (dw.) [W/m^2]',
     &      3*' '/

      LINE( 1)='DSET   ^attmxxx_%y4.grd'
      LINE( 2)='TITLE   Means/variances from run no. xxx'                 
      LINE( 3)='UNDEF   9.999E+19'
      LINE( 4)='OPTIONS sequential template big_endian 365_day_calendar'
      LINE( 5)='XDEF     nnn  LINEAR     0.000     x.xxx'
      LINE( 6)='YDEF     nnn  LEVELS'
      LINE( 7)='ZDEF      nn  LEVELS       950'
      LINE( 8)='TDEF    nnnn  LINEAR  1jan1900    nnnndy'
      LINE( 9)='VARS      nn'
      LINE(10)='ENDVARS'

      CTLNAME=NAME//NORUN//'.ctl'
      OPEN ( UNIT=IUNIT, FILE=CTLNAME, FORM='FORMATTED' )

      C1=90./ASIN(1.)

      do K=1,NLEV
         ILEV(K)=NINT(1000.*RLEV(K))
      enddo
C
C *** 2. Insert parameters in strings
C
      LINE(1)( 9:12)= NAME(1:4)
      LINE(1)(13:15)=NORUN(1:3)
      LINE(2)(43:45)=NORUN(1:3)

      write (LINE(5)(10:12),'(I3)') NLON
      write (LINE(5)(31:40),'(F10.3)') (360./NLON)
      write (LINE(6)(10:12),'(I3)') NLAT
      write (LINE(7)(10:12),'(I3)') NLEV

      write (LINE(8) (7:12),'(I6)') NTM
             LINE(8)(23:26) =       LMON(IMONT0)(1:4)
      write (LINE(8)(27:30),'(I4)') IYEAR0

      if (NDTM.lt.0) then
        write (LINE(8)(35:38),'(I4)') -NDTM
               LINE(8)(39:40) =       'mo'
      else if (mod(NDTM,1440).eq.0) then
        write (LINE(8)(35:38),'(I4)') NDTM/1440
               LINE(8)(39:40) =       'dy'
      else if (mod(NDTM,60).eq.0) then
        write (LINE(8)(35:38),'(I4)') NDTM/60 
               LINE(8)(39:40) =       'hr'
      else 
        write (LINE(8)(35:38),'(I4)') NDTM 
               LINE(8)(39:40) =       'mn'
      endif

      write (LINE(9)(11:12),'(I2)') N3D+N2D_1+N2D_2
C
C *** 3. Write ASCII control file 
C
      do JLINE=1,5
        write (IUNIT,1000) LINE(JLINE)
      enddo
C
      write (IUNIT,1010) LINE(6)(1:20), (C1*RLAT(J),J=1,NLAT)
      write (IUNIT,1020) LINE(7)(1:20), (ILEV(K),K=NLEV,1,-1)
C
      do JLINE=8,9
        write (IUNIT,1000) LINE(JLINE)
      enddo
C
      do JLINE=I3D,I3D+N3D-1
        write (LN3D(JLINE)(10:12),'(I3)') NLEV
        write (IUNIT,1000) LN3D(JLINE)
      enddo
C
      do JLINE=1,N2D_1
        write (IUNIT,1000) LN2D_1(JLINE)
      enddo

      do JLINE=1,N2D_2
        write (IUNIT,1000) LN2D_2(JLINE)
      enddo

      write (IUNIT,1000) LINE(10)

      CLOSE ( UNIT=IUNIT )

 1000 FORMAT (A80)
 1010 FORMAT (A20,6F10.3/(8F10.3))
 1020 FORMAT (A20,10I6)

      RETURN
      END 


      SUBROUTINE SETCTL_D (IUNIT,NLON,NLAT,NLEV,NTM,NDTM,N2D_1,N2D_2,
     *                     RLAT,RLEV,NAME,NORUN,IYEAR0,IMONT0)
C--
C--   Aux. routine SETCTL_D : write descriptor (.ctl) output file 
C--

      CHARACTER*80 LINE(10), LN2D_1(10), LN2D_2(10)
      CHARACTER*4  LMON(12)
      CHARACTER*5  NAME
      CHARACTER*3  NORUN
      CHARACTER*12 CTLNAME
      INTEGER ILEV(30), NCOUNT
      REAL RLAT(NLAT), RLEV(NLEV)
C
C *** 1. Initialization
C
      DATA LMON /'1jan','1feb','1mar','1apr','1may','1jun',
     &           '1jul','1aug','1sep','1oct','1nov','1dec'/

      DATA LN2D_1/
     &     'MSLP       0  99  mean-sea-level pressure         [hPa]',
     &     'TEMP0      0  99  near-surface air temperature   [degK]',
     &     'GH_500     0  99  geopotential height at 500 hPa    [m]',
     &     'U_850      0  99  zonal (u) wind at 850 hPa       [m/s]',
     &     'V_850      0  99  meridional (v) wind at 850 hPa  [m/s]',
     &     'Q_850      0  99  specific humidity at 850 hPa   [g/Kg]',
     &     'U_200      0  99  zonal (u) wind at 200 hPa       [m/s]',
     &     'V_200      0  99  meridional (v) wind at 200 hPa  [m/s]',
     &      2*' '/

      DATA LN2D_2/
     &     'PREC       0  99  precipitation                [mm/day]',
     &     'EVAP       0  99  evaporation                  [mm/day]',
     &     'USTR       0  99  u-stress                (dw.) [N/m^2]',
     &     'VSTR       0  99  v-stress                (dw.) [N/m^2]',
     &     'OLR        0  99  outgoing longwave rad.  (uw.) [W/m^2]',
     &     'LSHF       0  99  heat flux into land sfc (dw.) [W/m^2]',
     &     'SSHF       0  99  heat flux into  sea sfc (dw.) [W/m^2]',
     &      3*' '/

      LINE( 1)='DSET   ^attmdxxx_%y4.grd'
      LINE( 2)='TITLE   Daily means from run no. xxx'                 
      LINE( 3)='UNDEF   9.999E+19'
      LINE( 4)='OPTIONS sequential template big_endian 365_day_calendar'
      LINE( 5)='XDEF     nnn  LINEAR     0.000     x.xxx'
      LINE( 6)='YDEF     nnn  LEVELS'
      LINE( 7)='ZDEF      nn  LEVELS       950'
      LINE( 8)='TDEF  nnnnnn  LINEAR  1jan1900      nndy'
      LINE( 9)='VARS      nn'
      LINE(10)='ENDVARS'

      CTLNAME=NAME//NORUN//'.ctl'
      OPEN ( UNIT=IUNIT, FILE=CTLNAME, FORM='FORMATTED' )
C
      C1=90./ASIN(1.)
C
      do K=1,NLEV
        ILEV(K)=NINT(1000.*RLEV(K))
      enddo
C
C--   2. Insert parameters in strings
C
      LINE(1)( 9:13)= NAME(1:5)
      LINE(1)(14:16)=NORUN(1:3)
      LINE(2)(39:41)=NORUN(1:3)

      write (LINE(5)(10:12),'(I3)') NLON
      write (LINE(5)(31:40),'(F10.3)') (360./NLON)
      write (LINE(6)(10:12),'(I3)') NLAT
      write (LINE(7)(10:12),'(I3)') NLEV

      write (LINE(8) (7:12),'(I6)') NTM
             LINE(8)(23:26) =       LMON(IMONT0)(1:4)
      write (LINE(8)(27:30),'(I4)') IYEAR0
      write (LINE(8)(37:38),'(I2)') NDTM
             LINE(8)(39:40) =       'dy'

      write (LINE(9)(11:12),'(I2)') N2D_1+N2D_2
C
C--   3. Write ASCII control file 
C
      do JLINE=1,5
        write (IUNIT,1000) LINE(JLINE)
      enddo

      write (IUNIT,1010) LINE(6)(1:20), (C1*RLAT(J),J=1,NLAT)
      write (IUNIT,1020) LINE(7)(1:20), (ILEV(K),K=NLEV,1,-1)

      do JLINE=8,9
        write (IUNIT,1000) LINE(JLINE)
      enddo
      
      do JLINE=1,N2D_1
         write (IUNIT,1000) LN2D_1(JLINE)
      enddo
      
      do JLINE=1,N2D_2
         write (IUNIT,1000) LN2D_2(JLINE)
      enddo

      write (IUNIT,1000) LINE(10)

      CLOSE ( UNIT=IUNIT )

 1000 FORMAT (A80)
 1010 FORMAT (A20,6F10.3/(8F10.3))
 1020 FORMAT (A20,10I6)

      RETURN
      END 
