C--
C--   /DYNC1/ : Physical constants for dynamics (initial. in INDYNS)
      COMMON /DYNC1/ REARTH, OMEGA, GRAV, AKAP, RGAS, 
     &               PI,A,G
C--
C--   /DYNC2/ : Vertical level parameters (initial. in INDYNS)
      COMMON /DYNC2/ HSG(KXP), DHS(KX), FSG(KX), DHSR(KX), FSGR(KX)
C--
C--   /DYNC3/ : Functions of lat. and lon. (initial. in INDYNS)
      COMMON /DYNC3/ RADANG(IL), GSIN(IL), GCOS(IL), CORIOL(IL)
C--
C--   /DYNC4/ : Constants for hydrostatic eq. (initial. in INDYNS)
      COMMON /DYNC4/ XGEOP1(KX), XGEOP2(KX)
