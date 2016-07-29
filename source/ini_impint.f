      SUBROUTINE IMPINT (DT,ALPH)
C--
C--   SUBROUTINE IMPINT (DT,ALPH)
C--
C--   Purpose : initialize constants for implicit computation of
C--             horizontal diffusion and gravity waves
C--   Input :   DT   = time step
C--             ALPH = stepping coefficient for gravity wave scheme
C--                    (0.0 = forward, 0.5 = centred, 1.0 = backward)
C--   Initialized common blocks : DYNC5, DYNC6, HDIFC2
C--

c     IMPINT initializes constants for the implicit gravity wave computation.
c     It is assumed that that all implicit steps are of length DELT2 and use 
c     the forward/backward parameter ALPH.  IMPINT has to be re-called 
c     whenever either of these two parameters is changed. IMPINT should
c     be called even if the explicit option is chosen for the gravity wave
c     terms (the reference state temperature TREF is subtracted from some
c     terms anyway to reduce roundoff error; also the constants needed for
c     the biharmonic diffusion, which is assumed always to be backwards 
c     implicit, are defined in IMPINT)
										
      include "atparam.h"
      include "atparam1.h"

      include "com_dyncon0.h"
      include "com_dyncon1.h"
      include "com_dyncon2.h"
      include "com_hdifcon.h"

      DIMENSION DSUM(KX),INDX(KX),YA(KX,KX)

C--   1. Constants for backwards implicit biharmonic diffusion

      DO M=1,MX
        DO N=1,NX 
          DMP1 (M,N)=1./(1.+DMP (M,N)*DT)
          DMP1D(M,N)=1./(1.+DMPD(M,N)*DT)
          DMP1S(M,N)=1./(1.+DMPS(M,N)*DT)
        ENDDO
      ENDDO

C--   1. Constants for implicit gravity wave computation

c     reference atmosphere, function of sigma only

      RGAM = RGAS*GAMMA/(1000.*GRAV)

      DO 101 K=1,KX
        TREF(K)=288.*MAX(0.2,FSG(K))**RGAM
        print *, '  Tref = ', TREF(K)
        TREF1(K)=RGAS*TREF(K)
        TREF2(K)=AKAP*TREF(K)
        TREF3(K)=FSGR(K)*TREF(K)
  101 CONTINUE

c     Other constants 

      XI=DT*ALPH
      XXI = XI/(A*A)

      DO 201 K=1,KX
        DHSX(K)=XI*DHS(K)
  201 CONTINUE

      DO 60 N=1,NX
      DO 60 M=1,MX      
        MM=ISC*(M-1)+1
        LL=MM+N-2
        ELZ(M,N)=FLOAT(LL)*FLOAT(LL+1)*XXI
   60 CONTINUE

 
c   T(K) = TEX(K)+YA(K,K')*D(K') + XA(K,K')*SIG(K')

      DO 1 K=1,KX
      DO 1 K1=1,KXM
        XA(K,K1)=0.
    1 CONTINUE

      DO 2 K=1,KX
      DO 2 K1=1,KX
        YA(K,K1)=-AKAP*TREF(K)*DHS(K1)
    2 CONTINUE

      DO 3 K=2,KX
        XA(K,K-1)=0.5*(AKAP*TREF(K)/FSG(K)
     *      -(TREF(K)-TREF(K-1))/DHS(K))
    3 CONTINUE

      DO 4 K=1,KXM
        XA(K,K)=0.5*(AKAP*TREF(K)/FSG(K)
     *    -(TREF(K+1)-TREF(K))/DHS(K))
    4 CONTINUE

c     SIG(K)=XB(K,K')*D(K')

      DSUM(1)=DHS(1)
      DO 6 K=2,KX
        DSUM(K)=DSUM(K-1)+DHS(K)
    6 CONTINUE

      DO 5 K=1,KXM
      DO 5 K1=1,KX
        XB(K,K1)=DHS(K1)*DSUM(K)
        IF(K1.LE.K) XB(K,K1)=XB(K,K1)-DHS(K1)
    5 CONTINUE

c     T(K)=TEX(K)+XC(K,K')*D(K')

      DO 7 K=1,KX
      DO 7 K1=1,KX
        XC(K,K1)=YA(K,K1)
        DO 77 K2=1,KXM
          XC(K,K1)=XC(K,K1)+XA(K,K2)*XB(K2,K1)
   77   CONTINUE
    7 CONTINUE

c      P(K)=XD(K,K')*T(K') 

      DO 8 K=1,KX
      DO 8 K1=1,KX
        XD(K,K1)=0.
    8 CONTINUE

      DO 91 K=1,KX
      DO 91 K1=K+1,KX
        XD(K,K1)=RGAS*LOG(HSG(K1+1)/HSG(K1))
   91 CONTINUE
      DO 92 K=1,KX
        XD(K,K)=RGAS*LOG(HSG(K+1)/FSG(K))
   92 CONTINUE

c      P(K)=YE(K)+XE(K,K')*D(K')

      DO 10 K=1,KX
      DO 10 K1=1,KX
        XE(K,K1)=0.
      DO 10 K2=1,KX
        XE(K,K1)=XE(K,K1)+XD(K,K2)*XC(K2,K1)
   10 CONTINUE

      DO 11 L=1,LMAX
        XXX=(FLOAT(L)*FLOAT(L+1))/(A*A)
        DO 12 K=1,KX
        DO 12 K1=1,KX
          XF(K,K1,L)=XI*XI*XXX*(RGAS*TREF(K)*DHS(K1)-XE(K,K1))
   12   CONTINUE
        DO 13 K=1,KX
          XF(K,K,L)=XF(K,K,L)+1.
   13   CONTINUE
   11 CONTINUE

      DO 30 L=1,LMAX
        CALL INV(XF(1,1,L),XJ(1,1,L),INDX,KX)
   30 CONTINUE

      DO 50 K=1,KX
      DO 50 K1=1,KX
        XC(K,K1)=XC(K,K1)*XI
   50 CONTINUE

      RETURN
      END
