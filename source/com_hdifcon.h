C--
C--   /HDIFC1/ : Damping coef. for horizontal diffusion (explicit)
C--              (initial. in INDYNS)  
      COMMON /HDIFC1/ DMP(MX,NX), DMPD(MX,NX), DMPS(MX,NX)
C--
C--   /HDIFC2/ : Damping coef. for horizontal diffusion (implicit)
C--              (initial. in INDYNS)
      COMMON /HDIFC2/ DMP1(MX,NX), DMP1D(MX,NX), DMP1S(MX,NX)
C--
C--   /HDIFC3/ : Vertical comp. of orographic correction (initial. in INDYNS)
      COMMON /HDIFC3/ TCORV(KX), QCORV(KX) 
C--
C--   /HDIFC4/ : Horizontal component of orographic correction 
C--              (updated in FORDATE)
      COMMON /HDIFC4/ TCORH(MX,NX), QCORH(MX,NX) 
      COMPLEX         TCORH, QCORH
