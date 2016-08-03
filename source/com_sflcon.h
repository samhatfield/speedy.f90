!
!  /SFLCON/: Constants for surface fluxes (initial. in INPHYS)
!   FWIND0 = ratio of near-sfc wind to lowest-level wind
!   FTEMP0 = weight for near-sfc temperature extrapolation (0-1) :
!            1 : linear extrapolation from two lowest levels
!            0 : constant potential temperature ( = lowest level)
!   FHUM0  = weight for near-sfc specific humidity extrapolation (0-1) :
!            1 : extrap. with constant relative hum. ( = lowest level)
!            0 : constant specific hum. ( = lowest level)
!   CDL    = drag coefficient for momentum over land
!   CDS    = drag coefficient for momentum over sea
!   CHL    = heat exchange coefficient over land
!   CHS    = heat exchange coefficient over sea
!   VGUST  = wind speed for sub-grid-scale gusts
!   CTDAY  = daily-cycle correction (dTskin/dSSRad)
!   DTHETA = Potential temp. gradient for stability correction
!   FSTAB  = Amplitude of stability correction (fraction)
!   HDRAG  = Height scale for orographic correction
!   FHDRAG = Amplitude of orographic correction (fraction)
!   CLAMBDA = Heat conductivity in skin-to-root soil layer
!   CLAMBSN = Heat conductivity in soil for snow cover = 1

      REAL FWIND0, FTEMP0, FHUM0, CDL, CDS, CHL, CHS, VGUST, CTDAY,     &
     &                DTHETA, FSTAB, HDRAG, FHDRAG, CLAMBDA, CLAMBSN
      COMMON /SFLCON/ FWIND0, FTEMP0, FHUM0,                            &
     &                CDL, CDS, CHL, CHS, VGUST, CTDAY,                 &
     &                DTHETA, FSTAB, HDRAG, FHDRAG,                     &
     &                CLAMBDA, CLAMBSN

!
!   /SFLFIX/: Time-invariant fields (initial. in SFLSET)

      REAL FOROG
      COMMON /SFLFIX/ FOROG(NGP)
