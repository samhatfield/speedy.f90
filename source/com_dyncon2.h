C--
C--   /DYNC5/ : Temp. profile for semi-imp. scheme (initial. in IMPINT)
      COMMON /DYNC5/ TREF(KX), TREF1(KX), TREF2(KX), TREF3(KX)
C--
C--   /DYNC6/ : Arrays for semi-implicit scheme (initial. in IMPINT)
      COMMON /DYNC6/ XA(KX,KX), XB(KX,KX), XC(KX,KX),
     &               XD(KX,KX), XE(KX,KX),
     &               XF(KX,KX,LMAX), XG(KX,KX,LMAX),
     &               XH(KX,KX,LMAX), XJ(KX,KX,LMAX),
     &               DHSX(KX), ELZ(MX,NX)
