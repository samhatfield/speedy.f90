!
!:  /RADCON/: Radiation and cloud constants (initial. in INPHYS)

!:   SOLC   = Solar constant (area averaged) in W/m^2
!:   ALBSEA = Albedo over sea 
!:   ALBICE = Albedo over sea ice (for ice fraction = 1)
!:   ALBSN  = Albedo over snow (for snow cover = 1)

!:   RHCL1  = relative hum. threshold corr. to cloud cover = 0
!:   RHCL2  = relative hum. corr. to cloud cover = 1
!:   QACL   = specific hum. threshold for cloud cover
!:   WPCL   = cloud c. weight for the sq. root of precip. (for p = 1 mm/day)
!:   PMAXCL = max. value of precip. (mm/day) contributing to cloud cover 

!:   CLSMAX = maximum stratiform cloud cover
!:   CLSMINL= minimum stratiform cloud cover over land (for RH = 1)
!:   GSE_S0 = gradient of dry static energy corresp. to strat.c.c. = 0
!:   GSE_S1 = gradient of dry static energy corresp. to strat.c.c. = 1

!:   ALBCL  = cloud albedo (for cloud cover = 1)
!:   ALBCLS = stratiform cloud albedo (for st. cloud cover = 1)
!:   EPSSW  = fraction of incoming solar radiation absorbed by ozone
!:   EPSLW  = fraction of blackbody spectrum absorbed/emitted by PBL only
!:   EMISFC = longwave surface emissivity

!:            shortwave absorptivities (for dp = 10^5 Pa) :
!:   ABSDRY = abs. of dry air      (visible band)
!:   ABSAER = abs. of aerosols     (visible band)
!:   ABSWV1 = abs. of water vapour (visible band, for dq = 1 g/kg)
!:   ABSWV2 = abs. of water vapour (near IR band, for dq = 1 g/kg)
!:   ABSCL2 = abs. of clouds       (visible band, for dq_base = 1 g/kg)
!:   ABSCL1 = abs. of clouds       (visible band, maximum value)

!:            longwave absorptivities (per dp = 10^5 Pa) :
!:   ABLWIN = abs. of air in "window" band
!:   ABLCO2 = abs. of air in CO2 band
!:   ABLWV1 = abs. of water vapour in H2O band 1 (weak),   for dq = 1 g/kg
!:   ABLWV2 = abs. of water vapour in H2O band 2 (strong), for dq = 1 g/kg
!:   ABLCL1 = abs. of "thick" clouds in window band (below cloud top) 
!:   ABLCL2 = abs. of "thin" upper clouds in window and H2O bands
!
      REAL SOLC, ALBSEA, ALBICE, ALBSN, RHCL1, RHCL2, QACL, WPCL,       &
     &                PMAXCL, CLSMAX, CLSMINL,GSE_S0, GSE_S1, ALBCL,    &
     &                ALBCLS, EPSSW,  EPSLW,  EMISFC, ABSDRY, ABSAER,   &
     &                ABSWV1, ABSWV2, ABSCL1, ABSCL2, ABLWIN, ABLCO2,   &
     &                ABLCO2_ref, ABLWV1, ABLWV2, ABLCL1, ABLCL2

      COMMON /RADCON/ SOLC,   ALBSEA, ALBICE, ALBSN,                    &
     &                RHCL1,  RHCL2,  QACL,   WPCL,   PMAXCL,           &
     &                CLSMAX, CLSMINL,GSE_S0, GSE_S1,                   &
     &                ALBCL,  ALBCLS, EPSSW,  EPSLW,  EMISFC,           &
     &                ABSDRY, ABSAER, ABSWV1, ABSWV2, ABSCL1, ABSCL2,   &
     &                ABLWIN, ABLCO2, ABLCO2_ref,     ABLWV1, ABLWV2,   &
     &                ABLCL1, ABLCL2

								
!   /RADFIX/: Time-invariant fields (initial. in RADSET)
!    FBAND  = energy fraction emitted in each LW band = f(T)
!
      REAL FBAND(100:400,4)
      COMMON /RADFIX/ FBAND
								
!   /RADZON/: Zonally-averaged fields for SW/LW scheme (updated in SOL_OZ)
!    FSOL   = flux of incoming solar radiation
!    OZONE  = flux absorbed by ozone (lower stratos.)
!    OZUPP  = flux absorbed by ozone (upper stratos.)
!    ZENIT  = optical depth ratio (function of solar zenith angle)
!    STRATZ = stratospheric correction for polar night
!
      REAL FSOL(NGP), OZONE(NGP), OZUPP(NGP), ZENIT(NGP), STRATZ(NGP)
      COMMON /RADZON/ FSOL, OZONE, OZUPP, ZENIT, STRATZ

!   /RADSFC/: Radiative properties of the surface (updated in FORDATE)
!    ALB_L  = daily-mean albedo over land (bare-land + snow)
!    ALB_S  = daily-mean albedo over sea  (open sea + sea ice)
!    ALBSFC = combined surface albedo (land + sea)
!    SNOWC  = effective snow cover (fraction)
!
      REAL ALB_L(NGP), ALB_S(NGP), ALBSFC(NGP), SNOWC(NGP)
      COMMON /RADSFC/ ALB_L, ALB_S, ALBSFC, SNOWC

!   /RADFLD/: Transmissivity and blackbody rad. (updated in RADSW/RADLW)
!    TAU2   = transmissivity of atmospheric layers
!    ST4A   = blackbody emission from full and half atmospheric levels
!    STRATC = stratospheric correction term 
!    FLUX   = radiative flux in different spectral bands
!
      REAL TAU2(NGP,NLEV,4), ST4A(NGP,NLEV,2), STRATC(NGP,2),           &
     &                FLUX(NGP,4)
      COMMON /RADFLD/ TAU2, ST4A, STRATC, FLUX
								
!   /RADCLD/: Radiative properties of clouds (updated in CLOUD)
!    QCLOUD = Equivalent specific humidity of clouds 
!
      REAL QCLOUD(NGP), IRHTOP(NGP)
      COMMON /RADCLD/ QCLOUD, IRHTOP
