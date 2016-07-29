 
      SUBROUTINE SHTORH (IMODE,NGP,TA,PS,SIG,QA,RH,QSAT)
C--
C--   SUBROUTINE SHTORH (IMODE,NGP,TA,PS,SIG,QA,RH,QSAT)
C--
C--   Purpose: compute saturation specific humidity and 
C--            relative hum. from specific hum. (or viceversa)
C--   Input:   IMODE  : mode of operation
C--            NGP    : no. of grid-points
C--            TA     : abs. temperature
C--            PS     : normalized pressure   (=  p/1000_hPa) [if SIG < 0]
C--                   : normalized sfc. pres. (= ps/1000_hPa) [if SIG > 0]
C--            SIG    : sigma level
C--            QA     : specific humidity in g/kg [if IMODE > 0]
C--            RH     : relative humidity         [if IMODE < 0]
C--            QSAT   : saturation spec. hum. in g/kg 
C--   Output:  RH     : relative humidity         [if IMODE > 0] 
C--            QA     : specific humidity in g/kg [if IMODE < 0]
C--        
      REAL TA(NGP), PS(*), QA(NGP), RH(NGP), QSAT(NGP)
C
C---  1. Compute Qsat (g/kg) from T (degK) and normalized pres. P (= p/1000_hPa)
C        If SIG > 0, P = Ps * sigma, otherwise P = Ps(1) = const. 
C
      E0=  6.108E-3
      C1= 17.269
      C2= 21.875
      T0=273.16
      T1= 35.86
      T2=  7.66
      
      DO 110 J=1,NGP
        IF (TA(J).GE.T0) THEN
          QSAT(J)=E0*EXP(C1*(TA(J)-T0)/(TA(J)-T1))
        ELSE
          QSAT(J)=E0*EXP(C2*(TA(J)-T0)/(TA(J)-T2))
        ENDIF
  110 CONTINUE
C
      IF (SIG.LE.0.0) THEN
        DO 120 J=1,NGP
          QSAT(J)=622.*QSAT(J)/(PS(1)-0.378*QSAT(J))
  120   CONTINUE
      ELSE
        DO 130 J=1,NGP
          QSAT(J)=622.*QSAT(J)/(SIG*PS(J)-0.378*QSAT(J))
  130   CONTINUE
      ENDIF
C
C---  2. Compute rel.hum. RH=Q/Qsat (IMODE>0), or Q=RH*Qsat (IMODE<0)
C
      IF (IMODE.GT.0) THEN
        DO 210 J=1,NGP
          RH(J)=QA(J)/QSAT(J)
  210   CONTINUE
      ELSE IF (IMODE.LT.0) THEN
        DO 220 J=1,NGP
          QA(J)=RH(J)*QSAT(J)
  220   CONTINUE
      ENDIF
C				
      RETURN
      END

      SUBROUTINE ZMEDDY (NLON,NLAT,FF,ZM,EDDY)
C
C *** Decompose a field into zonal-mean and eddy component
C
      REAL FF(NLON,NLAT), ZM(NLAT), EDDY(NLON,NLAT)
C
      RNLON=1./NLON
C
      DO 130 J=1,NLAT
C
        ZM(J)=0.
        DO 110 I=1,NLON
          ZM(J)=ZM(J)+FF(I,J)
 110    CONTINUE
        ZM(J)=ZM(J)*RNLON
C
        DO 120 I=1,NLON
          EDDY(I,J)=FF(I,J)-ZM(J)
 120    CONTINUE
C
 130  CONTINUE
C
C--
      RETURN
      END




