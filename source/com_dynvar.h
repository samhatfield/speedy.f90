C--
C--   /DYNSP1/ : prognostic spectral variables for model dynamics 
C--              (initial. in INVARS, updated in STEP)
C--    VOR     = vorticity
C--    DIV     = divergence
C--    T       = abs. temperature
C--    PS      = log of (norm.) sfc pressure (p_s/p0)
C--    TR      = tracers (tr.1: spec. humidity in g/kg)

      COMMON /DYNSP1/ VOR(MX,NX,KX,2), DIV(MX,NX,KX,2), T(MX,NX,KX,2),
     *                PS(MX,NX,2), TR(MX,NX,KX,2,NTR)
      COMPLEX         VOR, DIV, T, PS, TR
C--
C--   /DYNSP2/ : geopotential 
C--    PHI     = atmos. geopotential  (updated in GEOP)
C--    PHIS    = surface geopotential (initial. in INVARS)

      COMMON /DYNSP2/ PHI(MX,NX,KX), PHIS(MX,NX)
      COMPLEX         PHI, PHIS
