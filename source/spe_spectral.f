******************************************************************
      SUBROUTINE GAUSSL(X,W,M)
c   a slightly modified version of a program in Numerical Recipes 
c       (Cambridge Univ. Press, 1989)
c   input:
c      m    = number of gaussian latitudes between pole and equator
c   output:
c      x(m) = sin(gaussian latitude) 
c      w(m) = weights in gaussian quadrature (sum should equal 1.0)
      DOUBLE PRECISION Z,Z1,P1,P2,P3,PP,EPS
      DIMENSION X(M),W(M)
      PARAMETER (EPS=3.D-14)
      N=2*M
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
    1   CONTINUE
        P1=1.D0
        P2=0.D0
        DO 11 J=1,N
          P3=P2
          P2=P1
          P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   11   CONTINUE
        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS) GO TO 1
        X(I)=Z
        W(I)=2.D0/((1.D0-Z*Z)*PP*PP)
   12 CONTINUE
      RETURN
      END
****************************************************************
      SUBROUTINE PARMTR(A)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

c   initializes Legendre transforms and constants used for other
c   subroutines that manipulate spherical harmonics
c
c   input:  A = radius of the sphere
c   first compute Gaussian latitudes and weights at the IY points from 
c       pole to equator
c   SIA(IY) is sin of latitude, WT(IY) are Gaussian weights for quadratures,
c     saved in COMMON/GAUSS (com_spectral.h)
      CALL GAUSSL(SIA,WT,IY)
      AM1 = 1./A
      AM2=  1./(A*A)
c   COA(IY) = cos(lat); WGHT needed for transforms, 
c             saved in COMMON/GAUSS (com_spectral.h)
      DO 1 J=1,IY
        COSQR = 1.0-SIA(J)**2
        COA(J)=SQRT(COSQR)
        WGHT(J)=WT(J)/(A*COSQR)
    1 CONTINUE
c   expand cosine and its reciprocal to cover both hemispheres, 
c      saved in COMMON/GAUSS (com_spectral.h)
      DO 11 J=1,IY
        JJ=IL+1-J
        COSG(J)=COA(J)
        COSG(JJ)=COA(J)
        COSGR(J)=1./COA(J)
        COSGR(JJ)=1./COA(J)
        COSGR2(J)=1./(COA(J)*coa(j))
        COSGR2(JJ)=1./(COA(J)*coa(j))
   11 CONTINUE
c  MM = zonal wavenumber = m
c     ISC=3 implies that only wavenumber 0,3,6,9,etc are included in model
c  LL = total wavenumber of spherical harmonic = l
c  L2 = l*(l+1)
c  EL2 = l*(l+1)/(a**2)
c  EL4 = EL2*EL2 ; for biharmonic diffusion
c  ELM2 = 1./EL2
c  TRFILT used to filter out "non-triangular" part of rhomboidal truncation
c   saved in COMMON/CSPEC (com_spectral.h)

      DO 2 N=1,NX
      nsh2(n)=0
       DO 2 M=1,MX
         MM(M)=ISC*(M-1)
         LL(M,N)=MM(M)+N-1
         L2(M,N)=LL(M,N)*(LL(M,N)+1)
         EL2(M,N)=FLOAT(L2(M,N))*AM2
         EL4(M,N)=EL2(M,N)*EL2(M,N)
         if (ll(m,n).le.ntrun1.or.ix.ne.4*iy) nsh2(n)=nsh2(n)+2
         if (ll(m,n).le.ntrun) then
           trfilt(m,n)=1.
         else
           trfilt(m,n)=0.
         endif
    2 CONTINUE
      ELM2(1,1)=0.
      DO 3 M=2,MX
      DO 3 N=1,NX
        ELM2(M,N)=1./EL2(M,N)
    3 CONTINUE
      DO 4 N=2,NX
        ELM2(1,N)=1./EL2(1,N)
    4 CONTINUE

c    quantities needed to generate and differentiate Legendre polynomials
c    all m values up to MXP = ISC*MTRUN+1 are needed by recursion relation 
c    saved in COMMON/LGND (com_spectral.h)
      DO 10 M=1,MXP
      DO 10 N=1,NXP
        EMM(M)=FLOAT(M-1)
        ELL(M,N)=FLOAT(N+M-2)
        EMM2=EMM(M)**2
        ELL2=ELL(M,N)**2
        IF(N.EQ.NXP) THEN
          EPSI(M,N)=0.0
        ELSE IF(N.EQ.1.AND.M.EQ.1) THEN
          EPSI(M,N)=0.0
        ELSE
          EPSI(M,N)=SQRT((ELL2-EMM2)/(4.*ELL2-1.))
        ENDIF
        REPSI(M,N)=0.0
        IF(EPSI(M,N).GT.0.) REPSI(M,N)=1./EPSI(M,N)
   10 CONTINUE
      SQRHLF=SQRT(.5)
      DO 20 M=2,MXP
        CONSQ(M) = SQRT(.5*(2.*EMM(M)+1.)/EMM(M))
   20 CONTINUE
c  quantities required by subroutines GRAD, UVSPEC, and VDS
c  saved in COMMON/UV, COMMON/GRADC, and COMMON/VDS (com_spectral.h)
      DO 15 M=1,MX
      DO 15 N=1,NX
        M1=MM(M)
        M2=M1+1
        EL1=FLOAT(LL(M,N))
        IF(N.EQ.1) THEN
          GRADX(M)=FLOAT(M1)/A
          UVDX(M,1)=-A/FLOAT(M1+1)
          UVDYM(M,1)=0.0
          VDDYM(M,1)=0.0
        ELSE
          UVDX(M,N)=-A*FLOAT(M1)/(EL1*(EL1+1))
          GRADYM(M,N)=(EL1-1.)*EPSI(M2,N)/A
          UVDYM(M,N)=-A*EPSI(M2,N)/EL1
          VDDYM(M,N)=(EL1+1)*EPSI(M2,N)/A
        ENDIF
        GRADYP(M,N)=(EL1+2.)*EPSI(M2,N+1)/A
        UVDYP(M,N)=-A*EPSI(M2,N+1)/(EL1+1.)
        VDDYP(M,N)=EL1*EPSI(M2,N+1)/A
   15 CONTINUE
c  generate associated Legendre polynomial
c  LGNDRE computes the polynomials at a particular latitiude, POLY(MX,NX), and stores
c  them in COMMON/POL (com_spectral.h)
c  polynomials and 'clones' stored in COMMON/POL1 (com_spectral.h)
      DO 30 J=1,IY
        CALL LGNDRE(J)
        DO 51 N=1,NX
        DO 51 M=1,MX
          M1=2*M-1
          M2=2*M
          CPOL(M1,N,J)=POLY(M,N)
          CPOL(M2,N,J)=POLY(M,N)
   51   CONTINUE
   30 CONTINUE
      RETURN
      END
****************************************************************
      SUBROUTINE LGNDRE(J)
c  follows Leith Holloways code 

      include "atparam.h"
c      include "param1spec.h"
      PARAMETER (SMALL = 1.E-30)
      include "com_spectral.h"

      DIMENSION ALP(MXP,NX)
        Y = COA(J)
        X = SIA(J)
c  start recursion with N=1 (M=L) diagonal 
        ALP(1,1) = SQRHLF
        DO 31 M=2,MXP
          ALP(M,1) = CONSQ(M)*Y*ALP(M-1,1)
   31   CONTINUE
c  continue with other elements
        DO 32 M=1,MXP
          ALP(M,2)=(X*ALP(M,1))*REPSI(M,2)
   32   CONTINUE
        DO 33 N=3,NX
        DO 33 M=1,MXP
          ALP(M,N)=(X*ALP(M,N-1)-EPSI(M,N-1)*ALP(M,N-2))*REPSI(M,N)
   33   CONTINUE
c  zero polynomials with absolute values smaller than 10**(-30)
        DO 34 N=1,NX
        DO 34 M=1,MXP
          IF(ABS(ALP(M,N)) .LE. SMALL) ALP(M,N)=0.0
   34   CONTINUE
c   pick off the required polynomials
        DO 35 N=1,NX
        DO 35 M=1,MX
          MM2=ISC*(M-1)+1
          POLY(M,N)=ALP(MM2,N)
   35   CONTINUE
      RETURN
      END
***************************************************************
      SUBROUTINE LAP(STRM,VORM)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      COMPLEX VORM(MX,NX),STRM(MX,NX)
      DO 1 N=1,NX
      DO 1 M=1,MX
        VORM(M,N)=-STRM(M,N)*EL2(M,N)
    1 CONTINUE  

      RETURN
      END
*******************************************************************
      SUBROUTINE INVLAP(VORM,STRM)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      COMPLEX VORM(MX,NX),STRM(MX,NX)
      DO 1 N=1,NX
      DO 1 M=1,MX
        STRM(M,N)=-VORM(M,N)*ELM2(M,N)
    1 CONTINUE
c     DO 1 M=1,MXNX
c       STRM(M,1)=-VORM(M,1)*ELM2(M,1)
c   1 CONTINUE
      RETURN
      END
*********************************************************************
      SUBROUTINE GRAD(PSI,PSDX,PSDY)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      DIMENSION PSI(2,MX,NX),PSDX(2,MX,NX),PSDY(2,MX,NX)
      DO 1 N=1,NX
      DO 1 M=1,MX
        PSDX(2,M,N)=GRADX(M)*PSI(1,M,N)
        PSDX(1,M,N)=-GRADX(M)*PSI(2,M,N)
    1 CONTINUE
      DO 2 K=1,2
      DO 2 M=1,MX
        PSDY(K,M,1)=GRADYP(M,1)*PSI(K,M,2)
        PSDY(K,M,NX)=-GRADYM(M,NX)*PSI(K,M,NTRUN1)
    2 CONTINUE
      DO 3 K=1,2
      DO 3 N=2,NTRUN1
      DO 3 M=1,MX
        PSDY(K,M,N)=-GRADYM(M,N)*PSI(K,M,N-1)+GRADYP(M,N)*PSI(K,M,N+1)
    3 CONTINUE
      RETURN
      END
******************************************************************
      SUBROUTINE VDS(UCOSM,VCOSM,VORM,DIVM)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"
                                                        
      DIMENSION VORM(2,MX,NX),DIVM(2,MX,NX),UCOSM(2,MX,NX),
     * VCOSM(2,MX,NX),ZC(2,MX,NX),ZP(2,MX,NX)
      do 1 n=1,nx
      DO 1 M=1,MX
        ZP(2,M,n)=GRADX(M)*UCOSM(1,M,n)
        ZP(1,M,n)=-GRADX(M)*UCOSM(2,M,n)
        ZC(2,M,n)=GRADX(M)*VCOSM(1,M,n)
        ZC(1,M,n)=-GRADX(M)*VCOSM(2,M,n)
    1 CONTINUE
      DO 3 K=1,2
      DO 3 M=1,MX
        VORM(K,M,1)=ZC(K,M,1)-VDDYP(M,1)*UCOSM(K,M,2)
        VORM(K,M,NX)=VDDYM(M,NX)*UCOSM(K,M,NTRUN1)
        DIVM(K,M,1)=ZP(K,M,1)+VDDYP(M,1)*VCOSM(K,M,2)
        DIVM(K,M,NX)=-VDDYM(M,NX)*VCOSM(K,M,NTRUN1)
    3 CONTINUE
      DO 4 K=1,2
      DO 4 N=2,NTRUN1
      DO 4 M=1,MX
        VORM(K,M,N)=VDDYM(M,N)*UCOSM(K,M,N-1)-VDDYP(M,N)*
     *   UCOSM(K,M,N+1)+ZC(K,M,N)  
        DIVM(K,M,N)=-VDDYM(M,N)*VCOSM(K,M,N-1)+VDDYP(M,N)*
     *   VCOSM(K,M,N+1)+ZP(K,M,N)
    4 CONTINUE
      RETURN     
      END
******************************************************************
      SUBROUTINE UVSPEC(VORM,DIVM,UCOSM,VCOSM)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"
                                                        
      DIMENSION VORM(2,MX,NX),DIVM(2,MX,NX),UCOSM(2,MX,NX),
     * VCOSM(2,MX,NX),ZC(2,MX,NX),ZP(2,MX,NX)
c     DO 1 M=1,MXNX
      DO 1 M=1,MX
      DO 1 n=1,NX
        ZP(2,M,n)=UVDX(M,n)*VORM(1,M,n)
        ZP(1,M,n)=-UVDX(M,n)*VORM(2,M,n)
        ZC(2,M,n)=UVDX(M,n)*DIVM(1,M,n)
        ZC(1,M,n)=-UVDX(M,n)*DIVM(2,M,n)
    1 CONTINUE
      DO 3 K=1,2
      DO 3 M=1,MX
        UCOSM(K,M,1)=ZC(K,M,1)-UVDYP(M,1)*VORM(K,M,2)
        UCOSM(K,M,NX)=UVDYM(M,NX)*VORM(K,M,NTRUN1)
        VCOSM(K,M,1)=ZP(K,M,1)+UVDYP(M,1)*DIVM(K,M,2)
        VCOSM(K,M,NX)=-UVDYM(M,NX)*DIVM(K,M,NTRUN1)
    3 CONTINUE
      DO 4 K=1,2
      DO 4 N=2,NTRUN1
      DO 4 M=1,MX
        VCOSM(K,M,N)=-UVDYM(M,N)*DIVM(K,M,N-1)+UVDYP(M,N)*
     *   DIVM(K,M,N+1)+ZP(K,M,N)  
        UCOSM(K,M,N)= UVDYM(M,N)*VORM(K,M,N-1)-UVDYP(M,N)*
     *   VORM(K,M,N+1)+ZC(K,M,N)
    4 CONTINUE
      RETURN     
      END
*******************************************************************
      SUBROUTINE GRID(VORM,VORG,KCOS)

      include "atparam.h"
c      include "param1spec.h"

      DIMENSION VORG(IX,IL)
      DIMENSION VORM(MX2,NX),VARM(MX2,IL)
      CALL GRIDY(VORM,VARM)
      CALL GRIDX(VARM,VORG,KCOS)
      RETURN
      END
*********************************************************************
      SUBROUTINE SPEC(VORG,VORM)

      include "atparam.h"
c      include "param1spec.h"

      DIMENSION VORG(IX,IL)
      DIMENSION VORM(MX2,NX),VARM(MX2,IL)
      CALL SPECX(VORG,VARM)
      CALL SPECY(VARM,VORM)
      RETURN
      END
*********************************************************************
      SUBROUTINE VDSPEC(UG,VG,VORM,DIVM,KCOS)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      DIMENSION UG(IX,IL),VG(IX,IL)
      DIMENSION VORM(MX2,NX),DIVM(MX2,NX),UM(MX2,IL),VM(MX2,IL)
      DIMENSION UG1(IX,IL),VG1(IX,IL)
      DIMENSION DUMC1(MX2,NX),DUMC2(MX2,NX)
      IF (KCOS.EQ.2) THEN
        DO 7 J=1,IL
          DO 8 I=1,IX
            UG1(I,J)=UG(I,J)*COSGR(J)
            VG1(I,J)=VG(I,J)*COSGR(J)
    8     CONTINUE
    7   CONTINUE
      ELSE
        DO 77 J=1,IL
          DO 88 I=1,IX
            UG1(I,J)=UG(I,J)*COSGR2(J)
            VG1(I,J)=VG(I,J)*COSGR2(J)
   88     CONTINUE
   77   CONTINUE
      ENDIF
      CALL SPECX(UG1,UM)  
      CALL SPECX(VG1,VM)
      CALL SPECY(UM,DUMC1)
      CALL SPECY(VM,DUMC2)
      CALL VDS(DUMC1,DUMC2,VORM,DIVM)
      RETURN
      END
*********************************************************************
      SUBROUTINE GRIDY(V,VARM)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      DIMENSION VARM(MX2,IL),V(MX2,NX)
      DIMENSION VM1(MX2),VM2(MX2)

      DO 100 J=1,IY
        J1=IL+1-J

        DO M=1,MX2
          VM1(M)=0.
          VM2(M)=0.
        ENDDO

        DO N=1,NX,2
C        DO M=1,MX2
         DO M=1,nsh2(n)
           VM1(M)=VM1(M)+V(M,N)*CPOL(M,N,J)
         ENDDO
        ENDDO

        DO N=2,NX,2
C        DO M=1,MX2
         DO M=1,nsh2(n)
           VM2(M)=VM2(M)+V(M,N)*CPOL(M,N,J)
         ENDDO
        ENDDO

        DO M=1,MX2
          VARM(M,J1)=VM1(M)+VM2(M)
          VARM(M,J) =VM1(M)-VM2(M)
        ENDDO

  100 CONTINUE

      RETURN
      END
******************************************************************
      SUBROUTINE SPECY(VARM,VORM)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      DIMENSION VARM(MX2,IL), VORM(MX2,NX)
      dimension svarm(mx2,iy), dvarm(mx2,iy)

      DO N=1,NX
       DO M=1,MX2
         VORM(M,N)=0.
       ENDDO
      ENDDO

      do j=1,iy
        j1=il+1-j
        do m=1,mx2
          svarm(m,j)=(varm(m,j1)+varm(m,j))*wt(j)
          dvarm(m,j)=(varm(m,j1)-varm(m,j))*wt(j)
        ENDDO
      ENDDO

      DO 100 J=1,IY
        J1=IL+1-J

        DO N=1,NTRUN1,2
C        DO M=1,MX2
         DO M=1,nsh2(n)
           VORM(M,N) = VORM(M,N)+CPOL(M,N,J)*svarm(m,j)
         ENDDO
        ENDDO

        DO N=2,NTRUN1,2
C        DO M=1,MX2
         DO M=1,nsh2(n)
           VORM(M,N) = VORM(M,N)+CPOL(M,N,J)*dvarm(m,j)
         ENDDO
        ENDDO

  100 CONTINUE

      RETURN
      END
******************************************************************
      SUBROUTINE TRUNCT(VOR)

      include "atparam.h"
c      include "param1spec.h"
      include "com_spectral.h"

      COMPLEX VOR(MX,NX)

c       DO 20 M=1,MXNX
        DO 20 N=1,NX
        DO 20 M=1,MX
          VOR(M,N)=VOR(M,N)*TRFILT(M,N)
   20   CONTINUE

      RETURN
      END
