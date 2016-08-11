subroutine setctl(iunit,nlon,nlat,nlev,ntm,ndtm,i3d,n3d,n2d_1,n2d_2,rlat,rlev,&
        & name,norun,iyear0,imont0)
    !  Aux. routine setctl : write descriptor (.ctl) output file 
    !

    implicit none

    integer, intent(in) :: iunit, nlon, nlat, nlev, ntm, ndtm, i3d, n3d, n2d_1,&
        & n2d_2, iyear0, imont0
    real, intent(in) :: rlat(nlat), rlev(nlev)
    character(len=80) :: line(10), ln3d(30), ln2d_1(20), ln2d_2(15)
    character(len=4)  :: lmon(12), name
    character(len=3)  :: norun
    character(len=11) :: ctlname
    integer :: ilev(30), j, jline, k
    real :: c1
 
    ! 1. Initialization
    lmon = (/'1jan','1feb','1mar','1apr','1may','1jun',&
        & '1jul','1aug','1sep','1oct','1nov','1dec'/)

    ln3d(:20) = (/&
     &     'GH         n  99  geopotential height               [m]',&
     &     'TEMP       n  99  abs. temperature               [degK]',&
     &     'U          n  99  zonal (u) wind                  [m/s]',&
     &     'V          n  99  meridional (v) wind             [m/s]',&
     &     'Q          n  99  specific humidity              [g/Kg]',&
     &     'RH         n  99  relative humidity                 [%]',&
     &     'OMEGA      n  99  pressure vertical velocity     [Pa/s]',&
     &     'PSI        n  99  streamfunction           [10^6 m^2/s]',&
     &     'CHI        n  99  velocity potential       [10^6 m^2/s]',&
                                                                     
     &     'VARGH      n  99  variance of geop. height        [m^2]',&
     &     'VART       n  99  variance of temperature      [degK^2]',&
     &     'VARU       n  99  variance of u-wind             [J/Kg]',&
     &     'VARV       n  99  variance of v-wind             [J/Kg]',&
     &     "COVUV      n  99  u'v' covariance (trans.)       [J/Kg]",&
     &     "COVVT      n  99  v't' covariance (trans.)   [degK m/s]",&
                                                                     
     &     'DTLSC      n  99  dt/dt by large-scale cond. [degK/day]',&
     &     'DTCNV      n  99  dt/dt by convection        [degK/day]',&
     &     'DTRSW      n  99  dt/dt by shortwave rad.    [degK/day]',&
     &     'DTRLW      n  99  dt/dt by longwave  rad.    [degK/day]',&
     &     'DTPBL      n  99  dt/dt by pbl processes     [degK/day]'/)
    ln3d(21:) = ' '

    ln2d_1(:18) = (/&
     &     'SP         0  99  surface pressure                [hPa]',&
     &     'MSLP       0  99  mean-sea-level pressure         [hPa]',&
     &     'ST         0  99  surface temperature            [degK]',&
     &     'SKINT      0  99  skin temperature               [degK]',&
     &     'SWAV       0  99  soil wetness availability         [%]',&
     &     'ALB        0  99  surface albedo                    [%]',&
     &     'U0         0  99  near-surface u-wind             [m/s]',&
     &     'V0         0  99  near-surface v-wind             [m/s]',&
     &     'TEMP0      0  99  near-surface air temperature   [degK]',&
     &     'RH0        0  99  near-surface relative humidity    [%]',&
     &     'CLC        0  99  cloud cover (deep clouds)         [%]',&
     &     'CLSTR      0  99  cloud cover (strat. clouds)       [%]',&
     &     'CLTOP      0  99  pressure at cloud top           [hPa]',&
     &     'IPTOP      0  99  highest precipitation level index  []',&
     &     'LST        0  99  land-surface temp.             [degK]',&
     &     'SST        0  99   sea-surface temp.             [degK]',&
     &     'SSTOM      0  99  ocean model sea-surface temp.  [degK]',&
     &     'SSTA       0  99  sst anomaly w.r.t. obs. clim.  [degK]'/)
    ln2d_1(19:) = ' '

    ln2d_2(:12) = (/&
     &     'PRECLS     0  99  large-scale precipitation    [mm/day]',&
     &     'PRECNV     0  99  convective precipitation     [mm/day]',&
     &     'EVAP       0  99  evaporation                  [mm/day]',&
     &     'USTR       0  99  u-stress                (dw.) [N/m^2]',&
     &     'VSTR       0  99  v-stress                (dw.) [M/m^2]',&
     &     'TSR        0  99  top shortwave rad.      (dw.) [W/m^2]',&
     &     'OLR        0  99  outgoing longwave rad.  (uw.) [W/m^2]',&
     &     'SSR        0  99  surface shortwave rad.  (dw.) [W/m^2]',&
     &     'SLR        0  99  surface longwave rad.   (uw.) [W/m^2]',&
     &     'SHF        0  99  sensible heat flux      (uw.) [W/m^2]',&
     &     'LSHF       0  99  heat flux into land sfc (dw.) [W/m^2]',&
     &     'SSHF       0  99  heat flux into  sea sfc (dw.) [W/m^2]'/)
     ln2d_2(13:) = ' '

    line( 1)='dset   ^attmxxx_%y4.grd'
    line( 2)='title   means/variances from run no. xxx'                 
    line( 3)='undef   9.999e+19'
    line( 4)='options sequential template big_endian 365_day_calendar'
    line( 5)='xdef     nnn  linear     0.000     x.xxx'
    line( 6)='ydef     nnn  levels'
    line( 7)='zdef      nn  levels       950'
    line( 8)='tdef    nnnn  linear  1jan1900    nnnndy'
    line( 9)='vars      nn'
    line(10)='endvars'

    ctlname=name//norun//'.ctl'
    open ( unit=iunit, file=ctlname, form='formatted' )
        c1=90./asin(1.)

        do k=1,nlev
           ilev(k)=nint(1000.*rlev(k))
        end do
 
        ! 2. Insert parameters in strings
        line(1)( 9:12)= name(1:4)
        line(1)(13:15)=norun(1:3)
        line(2)(43:45)=norun(1:3)

        write (line(5)(10:12),'(i3)') nlon
        write (line(5)(31:40),'(f10.3)') (360./nlon)
        write (line(6)(10:12),'(i3)') nlat
        write (line(7)(10:12),'(i3)') nlev

        write (line(8) (7:12),'(i6)') ntm
               line(8)(23:26) =       lmon(imont0)(1:4)
        write (line(8)(27:30),'(i4)') iyear0

        if (ndtm.lt.0) then
            write (line(8)(35:38),'(i4)') -ndtm
                 line(8)(39:40) =       'mo'
        else if (mod(ndtm,1440).eq.0) then
            write (line(8)(35:38),'(i4)') ndtm/1440
                 line(8)(39:40) =       'dy'
        else if (mod(ndtm,60).eq.0) then
            write (line(8)(35:38),'(i4)') ndtm/60 
                 line(8)(39:40) =       'hr'
        else 
            write (line(8)(35:38),'(i4)') ndtm 
                 line(8)(39:40) =       'mn'
        endif

        write (line(9)(11:12),'(i2)') n3d+n2d_1+n2d_2
 
        ! 3. Write ASCII control file 
        do jline=1,5
            write (iunit,1000) line(jline)
        end do
 
        write (iunit,1010) line(6)(1:20), (c1*rlat(j),j=1,nlat)
        write (iunit,1020) line(7)(1:20), (ilev(k),k=nlev,1,-1)
 
        do jline=8,9
            write (iunit,1000) line(jline)
        end do
 
        do jline=i3d,i3d+n3d-1
            write (ln3d(jline)(10:12),'(i3)') nlev
            write (iunit,1000) ln3d(jline)
        end do
 
        do jline=1,n2d_1
            write (iunit,1000) ln2d_1(jline)
        end do

        do jline=1,n2d_2
            write (iunit,1000) ln2d_2(jline)
        end do

        write (iunit,1000) line(10)

    close ( unit=iunit )

    1000 format (a80)
    1010 format (a20,6f10.3/(8f10.3))
    1020 format (a20,10i6)
end 

subroutine setctl_d(iunit,nlon,nlat,nlev,ntm,ndtm,n2d_1,n2d_2,rlat,rlev,name,&
        & norun,iyear0,imont0)
    !  Aux. routine setctl_d : write descriptor (.ctl) output file 

    implicit none

    integer, intent(in) :: iunit, nlon, nlat, nlev, ntm, ndtm, n2d_1, n2d_2,&
        & iyear0, imont0
    real, intent(in) :: rlat(nlat), rlev(nlev)
    character(len=80) :: line(10), ln2d_1(10), ln2d_2(10)
    character(len=4)  :: lmon(12)
    character(len=5)  :: name
    character(len=3)  :: norun
    character(len=12) :: ctlname
    integer :: ilev(30), ncount, j, jline, k
    real :: c1
 
    ! 1. Initialization
    lmon  = (/'1jan','1feb','1mar','1apr','1may','1jun',&
        & '1jul','1aug','1sep','1oct','1nov','1dec'/)

    ln2d_1(:8) = (/&
     &     'MSLP       0  99  mean-sea-level pressure         [hPa]',&
     &     'TEMP0      0  99  near-surface air temperature   [degK]',&
     &     'GH_500     0  99  geopotential height at 500 hPa    [m]',&
     &     'U_850      0  99  zonal (u) wind at 850 hPa       [m/s]',&
     &     'V_850      0  99  meridional (v) wind at 850 hPa  [m/s]',&
     &     'Q_850      0  99  specific humidity at 850 hPa   [g/Kg]',&
     &     'U_200      0  99  zonal (u) wind at 200 hPa       [m/s]',&
     &     'V_200      0  99  meridional (v) wind at 200 hPa  [m/s]'/)
    ln2d_1(9:) = ' '

    ln2d_2(:7) = (/&
     &     'PREC       0  99  precipitation                [mm/day]',&
     &     'EVAP       0  99  evaporation                  [mm/day]',&
     &     'USTR       0  99  u-stress                (dw.) [N/m^2]',&
     &     'VSTR       0  99  v-stress                (dw.) [N/m^2]',&
     &     'OLR        0  99  outgoing longwave rad.  (uw.) [W/m^2]',&
     &     'LSHF       0  99  heat flux into land sfc (dw.) [W/m^2]',&
     &     'SSHF       0  99  heat flux into  sea sfc (dw.) [W/m^2]'/)
    ln2d_2(8:) = ' '

    line( 1)='dset   ^attmdxxx_%y4.grd'
    line( 2)='title   daily means from run no. xxx'                 
    line( 3)='undef   9.999e+19'
    line( 4)='options sequential template big_endian 365_day_calendar'
    line( 5)='xdef     nnn  linear     0.000     x.xxx'
    line( 6)='ydef     nnn  levels'
    line( 7)='zdef      nn  levels       950'
    line( 8)='tdef  nnnnnn  linear  1jan1900      nndy'
    line( 9)='vars      nn'
    line(10)='endvars'

    ctlname=name//norun//'.ctl'
    open ( unit=iunit, file=ctlname, form='formatted' )
        c1=90./asin(1.)
 
        do k=1,nlev
            ilev(k)=nint(1000.*rlev(k))
        end do
 
        ! 2. Insert parameters in strings
        line(1)( 9:13)= name(1:5)
        line(1)(14:16)=norun(1:3)
        line(2)(39:41)=norun(1:3)

        write (line(5)(10:12),'(i3)') nlon
        write (line(5)(31:40),'(f10.3)') (360./nlon)
        write (line(6)(10:12),'(i3)') nlat
        write (line(7)(10:12),'(i3)') nlev

        write (line(8) (7:12),'(i6)') ntm
               line(8)(23:26) =       lmon(imont0)(1:4)
        write (line(8)(27:30),'(i4)') iyear0
        write (line(8)(37:38),'(i2)') ndtm
               line(8)(39:40) =       'dy'

        write (line(9)(11:12),'(i2)') n2d_1+n2d_2
 
        ! 3. Write ASCII control file 
        do jline=1,5
            write (iunit,1000) line(jline)
        end do

        write (iunit,1010) line(6)(1:20), (c1*rlat(j),j=1,nlat)
        write (iunit,1020) line(7)(1:20), (ilev(k),k=nlev,1,-1)

        do jline=8,9
            write (iunit,1000) line(jline)
        end do
          
        do jline=1,n2d_1
             write (iunit,1000) ln2d_1(jline)
        end do
          
        do jline=1,n2d_2
             write (iunit,1000) ln2d_2(jline)
        end do

        write (iunit,1000) line(10)
    close ( unit=iunit )

    1000 format (a80)
    1010 format (a20,6f10.3/(8f10.3))
    1020 format (a20,10i6)
end 
