!   /DYNC1/ : Physical constants for dynamics (initial. in INDYNS)
      REAL REARTH, OMEGA, GRAV, AKAP, RGAS, PI, A, G
      COMMON /DYNC1/ REARTH, OMEGA, GRAV, AKAP, RGAS, PI, A, G
!
!   /DYNC2/ : Vertical level parameters (initial. in INDYNS)
      REAL HSG(KXP), DHS(KX), FSG(KX), DHSR(KX), FSGR(KX)
      COMMON /DYNC2/ HSG, DHS, FSG, DHSR, FSGR
!
!   /DYNC3/ : Functions of lat. and lon. (initial. in INDYNS)
      REAL RADANG(IL), GSIN(IL), GCOS(IL), CORIOL(IL)
      COMMON /DYNC3/ RADANG, GSIN, GCOS, CORIOL
!
!   /DYNC4/ : Constants for hydrostatic eq. (initial. in INDYNS)
      REAL XGEOP1(KX), XGEOP2(KX)
      COMMON /DYNC4/ XGEOP1, XGEOP2
