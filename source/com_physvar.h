!
!   /PHYGR1/ : Model variables on gaussian grid (updated in PHYPAR)
!    UG1    = u-wind
!    VG1    = v-wind
!    TG1    = abs. temperature
!    QG1    = specific humidity (g/kg)
!    PHIG1  = geopotential
!    PSLG1  = log. of surface pressure

      REAL UG1, VG1, TG1, QG1, PHIG1, PSLG1
      COMMON /PHYGR1/ UG1(NGP,NLEV), VG1(NGP,NLEV), TG1(NGP,NLEV),      &
     &                QG1(NGP,NLEV), PHIG1(NGP,NLEV), PSLG1(NGP)

!   
!   /PHYGR2/ : Diagnosed upper-air variables (updated in PHYPAR)
!    SE     = dry static energy
!    RH     = relative humidity
!    QSAT   = saturation specific humidity (g/kg)

      REAL SE, RH, QSAT
      COMMON /PHYGR2/ SE(NGP,NLEV), RH(NGP,NLEV), QSAT(NGP,NLEV)

!
!   /PHYGR3/ : Diagnosed surface variables (updated in PHYPAR)
!    PSG    = surface pressure
!    TS     = surface temperature
!    TSKIN  = skin temperature
!    U0     = near-surface u-wind
!    V0     = near-surface v-wind
!    T0     = near-surface air temperature
!    Q0     = near-surface specific humidity (g/kg)
!    CLOUDC = total cloud cover (fraction)
!    CLSTR  = stratiform cloud cover (fraction)
!    CLTOP  = norm. pressure at cloud top
!    PRTOP  = top of precipitation (level index)
	

      REAL PSG, TS, TSKIN, U0, V0, T0, Q0, CLOUDC, CLSTR, CLTOP, PRTOP
      COMMON /PHYGR3/ PSG(NGP), TS(NGP), TSKIN(NGP),                    &
     &                U0(NGP), V0(NGP), T0(NGP), Q0(NGP),               &
     &                CLOUDC(NGP), CLSTR(NGP),                          &
     &                CLTOP(NGP),  PRTOP(NGP)

!
!   /PHYTEN/ : Physical param. tendencies (updated in PHYPAR)
!    TT_CNV  =  temperature tendency due to convection
!    QT_CNV  = sp. humidity tendency due to convection
!    TT_LSC  =  temperature tendency due to large-scale condensation
!    QT_LSC  = sp. humidity tendency due to large-scale condensation
!    TT_RSW  =  temperature tendency due to short-wave radiation
!    TT_RLW  =  temperature tendency due to long-wave radiation
!    UT_PBL  =       u-wind tendency due to PBL and diffusive processes
!    VT_PBL  =       v-wind tendency due to PBL and diffusive processes
!    TT_PBL  =  temperature tendency due to PBL and diffusive processes
!    QT_PBL  = sp. humidity tendency due to PBL and diffusive processes

      REAL TT_CNV, QT_CNV, TT_LSC, QT_LSC, TT_RSW, TT_RLW, UT_PBL,      &
     &                VT_PBL, TT_PBL, QT_PBL
      COMMON /PHYTEN/ TT_CNV(NGP,NLEV), QT_CNV(NGP,NLEV),               &
     &                TT_LSC(NGP,NLEV), QT_LSC(NGP,NLEV),               &
     &                TT_RSW(NGP,NLEV), TT_RLW(NGP,NLEV),               &
     &                UT_PBL(NGP,NLEV), VT_PBL(NGP,NLEV),               &
     &                TT_PBL(NGP,NLEV), QT_PBL(NGP,NLEV)

!
!   /FLUXES/ : Surface and upper boundary fluxes (updated in PHYPAR)
!    PRECNV = convective precipitation  [g/(m^2 s)], total
!    PRECLS = large-scale precipitation [g/(m^2 s)], total
!    SNOWCV = convective precipitation  [g/(m^2 s)], snow only
!    SNOWLS = large-scale precipitation [g/(m^2 s)], snow only
!    CBMF   = cloud-base mass flux 
!    TSR    = top-of-atm. shortwave radiation (downward)
!    SSRD   = surface shortwave radiation (downward-only)
!    SSR    = surface shortwave radiation (net downward)
!    SLRD   = surface longwave radiation  (downward-only)
!    SLR    = surface longwave radiation  (net upward) 
!    OLR    = outgoing longwave radiation (upward)
!    SLRU   = surface longwave emission   (upward)
!                                      (1:land, 2:sea, 3: wgt. average)
!    USTR   = u-stress                 (1:land, 2:sea, 3: wgt. average)
!    VSTR   = v-stress                 (1:land, 2:sea, 3: wgt. average)
!    SHF    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
!    EVAP   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
!    HFLUXN = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)

      REAL PRECNV, PRECLS, SNOWCV, SNOWLS, CBMF, TSR, SSRD, SSR, SLRD,  &
     &                SLR, OLR, SLRU, USTR, VSTR, SHF, EVAP, HFLUXN
      COMMON /FLUXES/ PRECNV(NGP), PRECLS(NGP),                         &
     &                SNOWCV(NGP), SNOWLS(NGP), CBMF(NGP),              &
     &                TSR(NGP), SSRD(NGP), SSR(NGP),                    &
     &                SLRD(NGP), SLR(NGP), OLR(NGP), SLRU(NGP,3),       &
     &                USTR(NGP,3), VSTR(NGP,3),                         &
     &                SHF(NGP,3), EVAP(NGP,3), HFLUXN(NGP,3)
