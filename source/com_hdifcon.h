!
!   /HDIFC1/ : Damping coef. for horizontal diffusion (explicit)
!              (initial. in INDYNS)  
      REAL DMP(MX,NX), DMPD(MX,NX), DMPS(MX,NX)
      COMMON /HDIFC1/ DMP, DMPD, DMPS
!
!   /HDIFC2/ : Damping coef. for horizontal diffusion (implicit)
!              (initial. in INDYNS)
      REAL DMP1(MX,NX), DMP1D(MX,NX), DMP1S(MX,NX)
      COMMON /HDIFC2/ DMP1, DMP1D, DMP1S
!
!   /HDIFC3/ : Vertical comp. of orographic correction (initial. in INDYNS)
      REAL TCORV, QCORV
      COMMON /HDIFC3/ TCORV(KX), QCORV(KX) 
!
!   /HDIFC4/ : Horizontal component of orographic correction 
!              (updated in FORDATE)
      COMPLEX         TCORH, QCORH
      COMMON /HDIFC4/ TCORH(MX,NX), QCORH(MX,NX) 
