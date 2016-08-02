!
!   /PHYCON/: Physical constants (initial. in INPHYS)
!    P0   = reference pressure
!    GG   = gravity accel.
!    RD   = gas constant for dry air
!    CP   = specific heat at constant pressure
!    ALHC = latent heat of condensation
!    ALHS = latent heat of sublimation
!    SBC  = Stefan-Boltzmann constant

      REAL P0, GG, RD, CP, ALHC, ALHS, SBC
      COMMON /PHYCON/ P0, GG, RD, CP, ALHC, ALHS, SBC
!
!   /FSIGLT/: Functions of sigma and latitude (initial. in INPHYS)
!    SIG    = full-level sigma 
!    SIGL   = logarithm of full-level sigma
!    SIGH   = half-level sigma
!    DSIG   = layer depth in sigma
!    POUT   = norm. pressure level [p/p0] for post-processing
!    GRDSIG = g/(d_sigma p0) : to convert fluxes of u,v,q into d(u,v,q)/dt
!    GRDSCP = g/(d_sigma p0 c_p): to convert energy fluxes into dT/dt
!    WVI    = weights for vertical interpolation
!    SLAT   = sin(lat)
!    CLAT   = cos(lat)

      REAL SIG(NLEV), SIGL(NLEV), SIGH(0:NLEV), DSIG(NLEV), POUT(NLEV), &
     &                GRDSIG(NLEV), GRDSCP(NLEV), WVI(NLEV,2),          &
     &                SLAT(NLAT), CLAT(NLAT)
      COMMON /FSIGLT/ SIG, SIGL, SIGH, DSIG, POUT, GRDSIG, GRDSCP, WVI, &
     &                SLAT, CLAT
