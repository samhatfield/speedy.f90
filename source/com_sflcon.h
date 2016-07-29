C--
C--:  /SFLCON/: Constants for surface fluxes (initial. in INPHYS)
C--:   FWIND0 = ratio of near-sfc wind to lowest-level wind
C--:   FTEMP0 = weight for near-sfc temperature extrapolation (0-1) :
C--:            1 : linear extrapolation from two lowest levels
C--:            0 : constant potential temperature ( = lowest level)
C--:   FHUM0  = weight for near-sfc specific humidity extrapolation (0-1) :
C--:            1 : extrap. with constant relative hum. ( = lowest level)
C--:            0 : constant specific hum. ( = lowest level)
C--:   CDL    = drag coefficient for momentum over land
C--:   CDS    = drag coefficient for momentum over sea
C--:   CHL    = heat exchange coefficient over land
C--:   CHS    = heat exchange coefficient over sea
C--:   VGUST  = wind speed for sub-grid-scale gusts
C--:   CTDAY  = daily-cycle correction (dTskin/dSSRad)
C--:   DTHETA = Potential temp. gradient for stability correction
C--:   FSTAB  = Amplitude of stability correction (fraction)
C--:   HDRAG  = Height scale for orographic correction
C--:   FHDRAG = Amplitude of orographic correction (fraction)
C--:   CLAMBDA = Heat conductivity in skin-to-root soil layer
C--:   CLAMBSN = Heat conductivity in soil for snow cover = 1

      COMMON /SFLCON/ FWIND0, FTEMP0, FHUM0,
     &                CDL, CDS, CHL, CHS, VGUST, CTDAY, 
     &                DTHETA, FSTAB, HDRAG, FHDRAG,
     &                CLAMBDA, CLAMBSN

C--
C--   /SFLFIX/: Time-invariant fields (initial. in SFLSET)

      COMMON /SFLFIX/ FOROG(NGP)
