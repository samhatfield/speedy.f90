C--
C--   /PHYCON/: Physical constants (initial. in INPHYS)
C--    P0   = reference pressure
C--    GG   = gravity accel.
C--    RD   = gas constant for dry air
C--    CP   = specific heat at constant pressure
C--    ALHC = latent heat of condensation
C--    ALHS = latent heat of sublimation
C--    SBC  = Stefan-Boltzmann constant

      COMMON /PHYCON/ P0, GG, RD, CP, ALHC, ALHS, SBC
C--
C--   /FSIGLT/: Functions of sigma and latitude (initial. in INPHYS)
C--    SIG    = full-level sigma 
C--    SIGL   = logarithm of full-level sigma
C--    SIGH   = half-level sigma
C--    DSIG   = layer depth in sigma
C--    POUT   = norm. pressure level [p/p0] for post-processing
C--    GRDSIG = g/(d_sigma p0) : to convert fluxes of u,v,q into d(u,v,q)/dt
C--    GRDSCP = g/(d_sigma p0 c_p): to convert energy fluxes into dT/dt
C--    WVI    = weights for vertical interpolation
C--    SLAT   = sin(lat)
C--    CLAT   = cos(lat)

      COMMON /FSIGLT/ SIG(NLEV), SIGL(NLEV), SIGH(0:NLEV), DSIG(NLEV),
     *                POUT(NLEV), GRDSIG(NLEV), GRDSCP(NLEV), 
     *                WVI(NLEV,2), SLAT(NLAT), CLAT(NLAT)
