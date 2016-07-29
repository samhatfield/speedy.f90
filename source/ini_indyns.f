
      SUBROUTINE INDYNS
C--
C--   SUBROUTINE INDYNS
C--
C--   Purpose : set time-stepping constants and initialize coefficients
C--             and spectral operators for model dynamics
C--   Initialized common blocks: DYNC0, DYNC1,  DYNC2,  DYNC3,  DYNC4, 
C--                              HDIFC1, HDIFC3,
C--                              common blocks for spectral transforms 
C--                              (through routine PARMTR) 
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_tsteps.h"
      include "com_lflags.h"

      include "com_dyncon0.h"
      include "com_dyncon1.h"
      include "com_hdifcon.h"
      include "com_spectral.h"

C--   1. Definition of constants
 
      IF (MOD(NSTEPS,2).NE.0) STOP ' Invalid no. of time steps'
      DELT=86400./NSTEPS
      DELT2=2.*DELT

c     ROB = Robert filter parameter
      ROB  = 0.05

c     WIL = Williams filter parameter
c      WIL = 1.0
      WIL = 0.53

c     ALPH = 0 ---- forward step for gravity wave terms
c     ALPH = 1 ---- backward implicit -----------------
c     ALPH = 0.5 -- centered implicit -----------------
      ALPH = 0.5 

C     1.2 Physical constants required by the dynamical core

      REARTH = 6.371E+6
      OMEGA  = 7.292E-05
      GRAV   = 9.81
      AKAP   = 2./7.
      RGAS   = AKAP*1004.

      PI = 4.*ATAN(1.)
      A  = REARTH
      G  = GRAV

C     1.3 Reference vertical profiles of temperature and humidity
C         and horizontal diffusion constants

      include "cls_indyns.h"

c     Power of Laplacian in horizontal diffusion
      NPOWHD = 4

C--   2. Definition of model levels  

C     2.1 Half (vertical velocity) levels

      IF (KX.EQ.5) THEN
        HSG(1)=0.000
        HSG(2)=0.150
        HSG(3)=0.350
        HSG(4)=0.650   
        HSG(5)=0.900
        HSG(6)=1.000
      ELSE IF (KX.EQ.7) THEN
        HSG(1)=0.020
        HSG(2)=0.140
        HSG(3)=0.260
        HSG(4)=0.420   
        HSG(5)=0.600
        HSG(6)=0.770
        HSG(7)=0.900
        HSG(8)=1.000
       ELSE IF (KX.EQ.8) THEN
        HSG(1)=0.000
        HSG(2)=0.050
        HSG(3)=0.140
        HSG(4)=0.260
        HSG(5)=0.420
        HSG(6)=0.600
        HSG(7)=0.770
        HSG(8)=0.900
        HSG(9)=1.000
      ENDIF

      DO K=1,KXP
        PRINT *, ' Model half-level (*1000)', k, nint(HSG(k)*1000)
      ENDDO

c     2.2 Layer thicknesses and full (u,v,T) levels

      DO K=1,KX
        DHS(K)=HSG(K+1)-HSG(K)
        FSG(K)=0.5*(HSG(K+1)+HSG(K))
      ENDDO

      DO K=1,KX
        PRINT *, ' Model full-level (*1000)', k, nint(FSG(k)*1000)
      ENDDO

c    2.3 Additional functions of sigma

      DO K=1,KX
        DHSR(K)=0.5/DHS(K)
        FSGR(K)=AKAP/(2.*FSG(K))
      ENDDO

C--   3. Horizontal functions and spectral operators

C     3.1 Initialization of spectral operators 

      CALL PARMTR (REARTH)

C     3.2 Latitudes and functions of latitude
C         NB: J=1 is Southernmost point!

      DO J=1,IY
        JJ=IL+1-J
        RAD1=ASIN(SIA(J))
        RADANG(J) =-RAD1
        RADANG(JJ)= RAD1
        GSIN(J)   =-SIA(J)
        GSIN(JJ)  = SIA(J)
      ENDDO

      DO J=1,IL
        GCOS(J)=COSG(J)
        CORIOL(J)=2.*OMEGA*GSIN(J)
      ENDDO


C--   4. Coefficients to compute geopotential

      DO K=1,KX
        XGEOP1(K)=RGAS*LOG(HSG(K+1)/FSG(K))
        IF(K.NE.KX) XGEOP2(K+1)=RGAS*LOG(FSG(K+1)/HSG(K+1))
      ENDDO

C--   5. Coefficients for horizontal diffusion

C     5.1 Spectral damping coefficients

      HDIFF = 1./(THD *3600.)
      HDIFD = 1./(THDD*3600.)
      HDIFS = 1./(THDS*3600.)
      RLAP  = 1./FLOAT(MTRUN*(MTRUN+1))

      DO N=1,NX
        DO M=1,MX
          TWN =FLOAT(ISC*(M-1)+N-1)
          ELAP=(TWN*(TWN+1.)*RLAP)
          ELAPN=ELAP**NPOWHD
          DMP (M,N)=HDIFF*ELAPN
          DMPD(M,N)=HDIFD*ELAPN
          DMPS(M,N)=HDIFS*ELAP
        ENDDO
C       DMPS(1,N)=0.
      ENDDO

C     5.2 Orographic correction terms for temperature and humidity
C         (vertical component) 

      RGAM = RGAS*GAMMA/(1000.*GRAV)
      QEXP = HSCALE/HSHUM
 
      TCORV(1)=0.
      QCORV(1)=0.
      QCORV(2)=0.

      DO K=2,KX
        TCORV(K)=FSG(K)**RGAM
        IF (K.GT.2) QCORV(K)=FSG(K)**QEXP
        print *, ' temp/hum correction at level ', k, TCORV(k), QCORV(k)
      ENDDO

C--   
      RETURN
      END


