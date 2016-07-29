
      SUBROUTINE SOL_OZ (TYEAR)

C--
C--   SUBROUTINE SOL_OZ (TYEAR)
C--
C--   Purpose: Compute zonally-averaged fields to be used 
C--            in the computation of SW absorption:
C--            FSOL   = flux of incoming solar radiation
C--            OZONE  = flux absorbed by ozone (lower stratos.)
C--            OZUPP  = flux absorbed by ozone (upper stratos.)
C--            ZENIT  = function of solar zenith angle
C--   Input:   TYEAR  = time as fraction of year (0-1, 0 = 1jan.h00)
C--   Updated common blocks: RADZON
C--
C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Constants + functions of sigma and latitude
      include "com_physcon.h"

C     Radiation constants
      include "com_radcon.h"

      real topsr(nlat)

C     ALPHA = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
      ALPHA=4.*ASIN(1.)*(TYEAR+10./365.)
      DALPHA=0.
c      DALPHA=ASIN(0.5)

      COZ1= 1.0*MAX(0.,COS(ALPHA-DALPHA))
      COZ2= 1.8
C
      AZEN=1.0
      NZEN=2

      RZEN=-COS(ALPHA)*23.45*ASIN(1.)/90.
      CZEN=COS(RZEN)
      SZEN=SIN(RZEN)

      FS0=6.

C     solar radiation at the top
      call solar (tyear,4.*solc,nlat,clat,slat,
     &            topsr)

      DO J=1,NLAT

        J0=1+NLON*(J-1)
        FLAT2=1.5*SLAT(J)**2-0.5

C       solar radiation at the top
        FSOL(J0)=topsr(j)

C       ozone depth in upper and lower stratosphere 
        OZUPP(J0)=0.5*EPSSW
        OZONE(J0)=0.4*EPSSW*(1.0+COZ1*SLAT(J)+COZ2*FLAT2)

C       zenith angle correction to (downward) absorptivity 
        ZENIT(J0)=1.+AZEN*(1.-(CLAT(J)*CZEN+SLAT(J)*SZEN))**NZEN

C       ozone absorption in upper and lower stratosphere 
        OZUPP(J0)=FSOL(J0)*OZUPP(J0)*ZENIT(J0)
        OZONE(J0)=FSOL(J0)*OZONE(J0)*ZENIT(J0)

C       Polar night cooling in the stratosphere
        STRATZ(J0)=MAX(FS0-FSOL(J0),0.)

        DO I=1,NLON-1
          FSOL  (I+J0)=FSOL  (J0)
          OZONE (I+J0)=OZONE (J0)
          OZUPP (I+J0)=OZUPP (J0)
          ZENIT (I+J0)=ZENIT (J0)
          STRATZ(I+J0)=STRATZ(J0)
        ENDDO

      ENDDO
C
      RETURN
      END

      subroutine solar (tyear,csol,nlat,clat,slat,
     &                  topsr)

C--   Average daily flux of solar radiation, from Hartmann (1994)

      real clat(nlat), slat(nlat),
     &     topsr(nlat)

C--   1. Compute declination angle and Earth-Sun distance factor

      pigr  = 2.*asin(1.)
      alpha = 2.*pigr*tyear

      ca1 = cos(alpha)
      sa1 = sin(alpha)
      ca2 = ca1*ca1-sa1*sa1
      sa2 = 2.*sa1*ca1
      ca3 = ca1*ca2-sa1*sa2
      sa3 = sa1*ca2+sa2*ca1

      decl = 0.006918-0.399912*ca1+0.070257*sa1
     &               -0.006758*ca2+0.000907*sa2
     &               -0.002697*ca3+0.001480*sa3

      fdis = 1.000110+0.034221*ca1+0.001280*sa1
     &               +0.000719*ca2+0.000077*sa2

      cdecl = cos(decl)
      sdecl = sin(decl)
      tdecl = sdecl/cdecl

C--   2. Compute daily-average insolation at the atm. top

      csolp=csol/pigr

      do j=1,nlat

        ch0 = min(1.,max(-1.,-tdecl*slat(j)/clat(j)))
        h0  = acos(ch0)
        sh0 = sin(h0)

        topsr(j) = csolp*fdis*(h0*slat(j)*sdecl+sh0*clat(j)*cdecl)

      enddo

C--
      return
      end

      SUBROUTINE CLOUD (QA,RH,PRECNV,PRECLS,IPTOP,GSE,FMASK,
     &                  ICLTOP,CLOUDC,CLSTR)
C--
C--   SUBROUTINE CLOUD (QA,RH,PRECNV,PRECLS,IPTOP,GSE,FMASK,
C--  &                  ICLTOP,CLOUDC,CLSTR)
C--
C--   Purpose: Compute cloud-top level and cloud cover
C--   Input:   QA     = specific humidity [g/kg]                (3-dim)
C--            RH     = relative humidity                       (3-dim)
C--            PRECNV = convective precipitation                (2-dim)
C--            PRECLS = large-scale precipitation               (2-dim)
C--            IPTOP  = top level of precipitating cloud        (2-dim)
C--            GSE    = gradient of dry st. energy (dSE/dPHI)   (2-dim)
C--            FMASK  = fractional land-sea mask                (2-dim)
C--   Output:  ICLTOP = cloud top level (all clouds)            (2-dim)
C--            CLOUDC = total cloud cover                       (2-dim)
C--            CLSTR  = stratiform cloud cover                  (2-dim)

C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )
C
C     Constants + functions of sigma and latitude
C
      include "com_physcon.h"
C
C     Cloud and radiation parameters
C
      include "com_radcon.h"

      INTEGER IPTOP(NGP)
      REAL QA(NGP,NLEV), RH(NGP,NLEV),  
     &     PRECNV(NGP),  PRECLS(NGP), GSE(NGP), FMASK(NGP)

      INTEGER ICLTOP(NGP)
      REAL CLOUDC(NGP), CLSTR(NGP)
      
      NL1  = NLEV-1
      NLP  = NLEV+1
      RRCL = 1./(RHCL2-RHCL1)

C--   1.  Cloud cover, defined as the sum of:
C         - a term proportional to the square-root of precip. rate 
C         - a quadratic function of the max. relative humidity
C           in tropospheric layers above PBL where Q > QACL :
C           ( = 0 for RHmax < RHCL1, = 1 for RHmax > RHCL2 )
C         Cloud-top level: defined as the highest (i.e. least sigma)
C           between the top of convection/condensation and
C           the level of maximum relative humidity. 


      DO J=1,NGP
        IF (RH(J,NL1).GT.RHCL1) THEN
          CLOUDC(J) = RH(J,NL1)-RHCL1
          ICLTOP(J) = NL1
        ELSE
          CLOUDC(J) = 0.
          ICLTOP(J) = NLP
        ENDIF
      ENDDO

      DO K=3,NLEV-2
        DO J=1,NGP
          DRH = RH(J,K)-RHCL1
          IF (DRH.GT.CLOUDC(J).AND.QA(J,K).GT.QACL) THEN
            CLOUDC(J) = DRH
            ICLTOP(J) = K
          ENDIF
        ENDDO
      ENDDO

      DO J=1,NGP
        CL1 = MIN(1.,CLOUDC(J)*RRCL)
        PR1 = MIN(PMAXCL,86.4*(PRECNV(J)+PRECLS(J)))
        CLOUDC(J) = MIN(1.,WPCL*SQRT(PR1)+CL1*CL1)
        ICLTOP(J) = MIN(IPTOP(J),ICLTOP(J))
      ENDDO

C--   2.  Equivalent specific humidity of clouds 

      DO J=1,NGP
        QCLOUD(J) = QA(J,NL1)
      ENDDO

C--   3. Stratiform clouds at the top of PBL

      inew = 1

      if (inew.gt.0) then

c        CLSMAX  = 0.6
c        CLSMINL = 0.15
c        GSE_S0  = 0.25
c        GSE_S1  = 0.40

        CLFACT = 1.2
        RGSE   = 1./(GSE_S1-GSE_S0)

        DO J=1,NGP
c         Stratocumulus clouds over sea
          FSTAB    = MAX(0.,MIN(1.,RGSE*(GSE(J)-GSE_S0)))
          CLSTR(J) = FSTAB*MAX(CLSMAX-CLFACT*CLOUDC(J),0.)
c         Stratocumulus clouds over land
          CLSTRL   = MAX(CLSTR(J),CLSMINL)*RH(J,NLEV)
          CLSTR(J) = CLSTR(J)+FMASK(J)*(CLSTRL-CLSTR(J))
        ENDDO

      else

        CLSMAX  = 0.3
        CLSMINL = 0.1
        ALBCOR  = ALBCL/0.5
 
        DO J=1,NGP
c         Stratocumulus clouds over sea
          CLSTR(J) = MAX(CLSMAX-CLOUDC(J),0.)
c         Rescale for consistency with previous albedo values
          CLSTR(J) = CLSTR(J)*ALBCOR
c         Correction for aerosols over land
          CLSTR(J) = CLSTR(J)+FMASK(J)*(CLSMINL-CLSTR(J))
        ENDDO

      endif

      RETURN
      END


      SUBROUTINE RADSW (PSA,QA,ICLTOP,CLOUDC,CLSTR,
     &                  FSFCD,FSFC,FTOP,DFABS)
C--
C--   SUBROUTINE RADSW (PSA,QA,ICLTOP,CLOUDC,CLSTR,
C--  &                  FSFCD,FSFC,FTOP,DFABS)
C--
C--   Purpose: Compute the absorption of shortwave radiation and
C--            initialize arrays for longwave-radiation routines
C--   Input:   PSA    = norm. surface pressure [p/p0]           (2-dim)
C--            QA     = specific humidity [g/kg]                (3-dim)
C--            ICLTOP = cloud top level                         (2-dim)
C--            CLOUDC = total cloud cover                       (2-dim)
C--            CLSTR  = stratiform cloud cover                  (2-dim)
C--   Output:  FSFCD  = downward-only flux of sw rad. at the surface (2-dim)
C--            FSFC   = net (downw.) flux of sw rad. at the surface  (2-dim)
C--            FTOP   = net (downw.) flux of sw rad. at the atm. top (2-dim)
C--            DFABS  = flux of sw rad. absorbed by each atm. layer  (3-dim)
C--
C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Constants + functions of sigma and latitude

      include "com_physcon.h"

C     Radiation parameters

      include "com_radcon.h"
C
      INTEGER ICLTOP(NGP)
      REAL PSA(NGP), QA(NGP,NLEV), CLOUDC(NGP), CLSTR(NGP)

      REAL FTOP(NGP), FSFC(NGP), FSFCD(NGP), DFABS(NGP,NLEV)

      REAL ACLOUD(NGP), PSAZ(NGP), FREFL(NGP,NLEV)

      EQUIVALENCE (FREFL(1,1),TAU2(1,1,3))


      NL1 = NLEV-1

      FBAND2 = 0.05
      FBAND1 = 1.-FBAND2

c      ALBMINL=0.05
c      ALBCLS = 0.5
C
C--   1.  Initialization

      DO K=1,NLEV
       DO J=1,NGP
         FREFL(J,K)=0.
       ENDDO
      ENDDO

      DO J=1,NGP
cfk-- change to ensure only ICLTOP <= NLEV used
        IF(ICLTOP(J) .LE. NLEV) THEN
          FREFL(J,ICLTOP(J))= ALBCL*CLOUDC(J)
        ENDIF
cfk-- end change
        FREFL(J,NLEV)     = ALBCLS*CLSTR(J)
      ENDDO

C
C--   2. Shortwave transmissivity:
C        function of layer mass, ozone (in the statosphere),
C        abs. humidity and cloud cover (in the troposphere)

      DO J=1,NGP
        PSAZ(J)=PSA(J)*ZENIT(J)
        ACLOUD(J)=CLOUDC(J)*MIN(ABSCL1*QCLOUD(J),ABSCL2)
      ENDDO

      DO J=1,NGP
        DELTAP=PSAZ(J)*DSIG(1)
        TAU2(J,1,1)=EXP(-DELTAP*ABSDRY)
      ENDDO

      DO K=2,NL1
       ABS1=ABSDRY+ABSAER*SIG(K)**2
       DO J=1,NGP
         DELTAP=PSAZ(J)*DSIG(K)
         IF (K.GE.ICLTOP(J)) THEN
           TAU2(J,K,1)=EXP(-DELTAP*
     &                 (ABS1+ABSWV1*QA(J,K)+ACLOUD(J)))
         ELSE
           TAU2(J,K,1)=EXP(-DELTAP*(ABS1+ABSWV1*QA(J,K)))
         ENDIF
       ENDDO
      ENDDO

      ABS1=ABSDRY+ABSAER*SIG(NLEV)**2
      DO J=1,NGP
        DELTAP=PSAZ(J)*DSIG(NLEV)
        TAU2(J,NLEV,1)=EXP(-DELTAP*(ABS1+ABSWV1*QA(J,NLEV)))
      ENDDO

      DO K=2,NLEV
       DO J=1,NGP
         DELTAP=PSAZ(J)*DSIG(K)
         TAU2(J,K,2)=EXP(-DELTAP*ABSWV2*QA(J,K))
       ENDDO
      ENDDO

C
C---  3. Shortwave downward flux 

C     3.1 Initialization of fluxes 

      DO J=1,NGP
        FTOP(J)  =FSOL(J)
        FLUX(J,1)=FSOL(J)*FBAND1
        FLUX(J,2)=FSOL(J)*FBAND2
      ENDDO

C     3.2 Ozone and dry-air absorption in the stratosphere

      K=1
      DO J=1,NGP
        DFABS(J,K)=FLUX(J,1)
        FLUX (J,1)=TAU2(J,K,1)*(FLUX(J,1)-OZUPP(J)*PSA(J))
        DFABS(J,K)=DFABS(J,K)-FLUX(J,1)
      ENDDO

      K=2
      DO J=1,NGP
        DFABS(J,K)=FLUX(J,1)
        FLUX (J,1)=TAU2(J,K,1)*(FLUX(J,1)-OZONE(J)*PSA(J))
        DFABS(J,K)=DFABS(J,K)-FLUX(J,1)
      ENDDO

C     3.3  Absorption and reflection in the troposphere

      DO K=3,NLEV
       DO J=1,NGP
         FREFL(J,K)=FLUX(J,1)*FREFL(J,K)
         FLUX (J,1)=FLUX(J,1)-FREFL(J,K)
         DFABS(J,K)=FLUX(J,1)
         FLUX (J,1)=TAU2(J,K,1)*FLUX(J,1)
         DFABS(J,K)=DFABS(J,K)-FLUX(J,1)
       ENDDO
      ENDDO

      DO K=2,NLEV
       DO J=1,NGP
         DFABS(J,K)=DFABS(J,K)+FLUX(J,2)
         FLUX (J,2)=TAU2(J,K,2)*FLUX(J,2)
         DFABS(J,K)=DFABS(J,K)-FLUX(J,2)
       ENDDO
      ENDDO

C
C---  4. Shortwave upward flux 
	
C     4.1  Absorption and reflection at the surface

      DO J=1,NGP
        FSFCD(J)  = FLUX(J,1)+FLUX(J,2)
        FLUX(J,1) = FLUX(J,1)*ALBSFC(J)
        FSFC(J)   = FSFCD(J)-FLUX(J,1)
      ENDDO
	
C     4.2  Absorption of upward flux

      DO K=NLEV,1,-1
       DO J=1,NGP
         DFABS(J,K)=DFABS(J,K)+FLUX(J,1)
         FLUX (J,1)=TAU2(J,K,1)*FLUX(J,1)
         DFABS(J,K)=DFABS(J,K)-FLUX(J,1)
         FLUX (J,1)=FLUX(J,1)+FREFL(J,K)
       ENDDO
      ENDDO

C     4.3  Net solar radiation = incoming - outgoing

      DO J=1,NGP
        FTOP(J)=FTOP(J)-FLUX(J,1)
      ENDDO
C
C---  5.  Initialization of longwave radiation model

C     5.1  Longwave transmissivity:
C          function of layer mass, abs. humidity and cloud cover.

C     Cloud-free levels (stratosphere + PBL)

      K=1
        DO J=1,NGP
          DELTAP=PSA(J)*DSIG(K)
          TAU2(J,K,1)=EXP(-DELTAP*ABLWIN)
          TAU2(J,K,2)=EXP(-DELTAP*ABLCO2)
          TAU2(J,K,3)=1.
          TAU2(J,K,4)=1.
        ENDDO

      DO K=2,NLEV,NLEV-2
        DO J=1,NGP
          DELTAP=PSA(J)*DSIG(K)
          TAU2(J,K,1)=EXP(-DELTAP*ABLWIN)
          TAU2(J,K,2)=EXP(-DELTAP*ABLCO2)
          TAU2(J,K,3)=EXP(-DELTAP*ABLWV1*QA(J,K))
          TAU2(J,K,4)=EXP(-DELTAP*ABLWV2*QA(J,K))
        ENDDO
      ENDDO

C     Cloudy layers (free troposphere)

      DO J=1,NGP
        ACLOUD(J)=CLOUDC(J)*ABLCL2
      ENDDO

      DO K=3,NL1
       DO J=1,NGP
         DELTAP=PSA(J)*DSIG(K)
         IF (K.LT.ICLTOP(J)) THEN
           ACLOUD1=ACLOUD(J)
         ELSE
           ACLOUD1=ABLCL1*CLOUDC(J)
         ENDIF
         TAU2(J,K,1)=EXP(-DELTAP*(ABLWIN+ACLOUD1))
         TAU2(J,K,2)=EXP(-DELTAP*ABLCO2)
         TAU2(J,K,3)=EXP(-DELTAP*MAX(ABLWV1*QA(J,K),ACLOUD(J)))
         TAU2(J,K,4)=EXP(-DELTAP*MAX(ABLWV2*QA(J,K),ACLOUD(J)))
       ENDDO
      ENDDO

C     5.2  Stratospheric correction terms

      EPS1=EPSLW/(DSIG(1)+DSIG(2))
      DO J=1,NGP
        STRATC(J,1)=STRATZ(J)*PSA(J)
        STRATC(J,2)=EPS1*PSA(J)
      ENDDO
C
C---
      RETURN
      END


      SUBROUTINE RADLW (IMODE,TA,TS,
     &                  FSFCD,FSFCU,
     &                  FSFC,FTOP,DFABS)

C--
C--   SUBROUTINE RADLW (IMODE,TA,TS,
C--  &                  FSFCD,FSFCU,
C--  &                  FSFC,FTOP,DFABS)
C--
C--   Purpose: Compute the absorption of longwave radiation
C--   Input:   IMODE  = index for operation mode 
C--                     -1 : downward flux only
C--                      0 : downward + upward flux 
C--                     +1 : upward flux only
C--            TA     = absolute temperature (3-dim)
C--            TS     = surface temperature                    [if IMODE=0]
C--            FSFCD  = downward flux of lw rad. at the sfc.   [if IMODE=1]
C--            FSFCU  = surface blackbody emission (upward)    [if IMODE=1]
C--            DFABS  = DFABS output from RADLW(-1,... )       [if IMODE=1]
C--   Output:  FSFCD  = downward flux of lw rad. at the sfc.[if IMODE=-1,0]
C--            FSFCU  = surface blackbody emission (upward)  [if IMODE=  0]
C--            FSFC   = net upw. flux of lw rad. at the sfc. [if IMODE=0,1]
C--            FTOP   = outgoing flux of lw rad. at the top  [if IMODE=0,1]
C--            DFABS  = flux of lw rad. absorbed by each atm. layer (3-dim)
C--
C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Number of radiation bands with tau < 1
      PARAMETER ( NBAND=4 )

C     Constants + functions of sigma and latitude

      include "com_physcon.h"

C     Radiation parameters

      include "com_radcon.h"

      REAL TA(NGP,NLEV), TS(NGP)

      REAL FSFCD(NGP), FSFCU(NGP)

      REAL FTOP(NGP), FSFC(NGP), DFABS(NGP,NLEV)

C
      NL1=NLEV-1

      REFSFC=1.-EMISFC

      IF (IMODE.EQ.1) GO TO 410

C---  1. Blackbody emission from atmospheric levels.
C        The linearized gradient of the blakbody emission is computed
C        from temperatures at layer boundaries, which are interpolated 
C        assuming a linear dependence of T on log_sigma.
C        Above the first (top) level, the atmosphere is assumed isothermal.

C     Temperature at level boundaries
      DO K=1,NL1
       DO J=1,NGP
         ST4A(J,K,1)=TA(J,K)+WVI(K,2)*(TA(J,K+1)-TA(J,K))
       ENDDO
      ENDDO 

C     Mean temperature in stratospheric layers 
      DO J=1,NGP
        ST4A(J,1,2)=0.75*TA(J,1)+0.25* ST4A(J,1,1)
        ST4A(J,2,2)=0.50*TA(J,2)+0.25*(ST4A(J,1,1)+ST4A(J,2,1))
      ENDDO

C     Temperature gradient in tropospheric layers 

      ANIS =1.0
      ANISH=0.5*ANIS

      DO K=3,NL1
       DO J=1,NGP
         ST4A(J,K,2)=ANISH*MAX(ST4A(J,K,1)-ST4A(J,K-1,1),0.)
       ENDDO
      ENDDO

      DO J=1,NGP
        ST4A(J,NLEV,2)=ANIS*MAX(TA(J,NLEV)-ST4A(J,NL1,1),0.)
      ENDDO

C     Blackbody emission in the stratosphere
      DO K=1,2
       DO J=1,NGP
         ST4A(J,K,1)=SBC*ST4A(J,K,2)**4
         ST4A(J,K,2)=0.
       ENDDO
      ENDDO

C     Blackbody emission in the troposphere
      DO K=3,NLEV
       DO J=1,NGP
         ST3A=SBC*TA(J,K)**3
         ST4A(J,K,1)=ST3A*TA(J,K)
         ST4A(J,K,2)=4.*ST3A*ST4A(J,K,2)
       ENDDO
      ENDDO


C---  2. Initialization of fluxes

      DO J=1,NGP
        FSFCD(J)=0.
      ENDDO

      DO K=1,NLEV
       DO J=1,NGP
         DFABS(J,K)=0.
       ENDDO
      ENDDO

C---  3. Emission ad absorption of longwave downward flux.
C        For downward emission, a correction term depending on the      
C        local temperature gradient and on the layer transmissivity is  
C        added to the average (full-level) emission of each layer. 
	
C     3.1  Stratosphere

      K=1
      DO JB=1,2
       DO J=1,NGP
         EMIS=1.-TAU2(J,K,JB)
         BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)+EMIS*ST4A(J,K,2))
         FLUX(J,JB)=EMIS*BRAD
         DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
       ENDDO
      ENDDO

      DO JB=3,NBAND
       DO J=1,NGP
         FLUX(J,JB)=0.
       ENDDO
      ENDDO
	
C     3.2  Troposphere

      DO JB=1,NBAND
       DO K=2,NLEV
        DO J=1,NGP
          EMIS=1.-TAU2(J,K,JB)
          BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)+EMIS*ST4A(J,K,2))
          DFABS(J,K)=DFABS(J,K)+FLUX(J,JB)
          FLUX(J,JB)=TAU2(J,K,JB)*FLUX(J,JB)+EMIS*BRAD
          DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
        ENDDO
       ENDDO
      ENDDO

C     3.3 Surface downward flux

      DO JB=1,NBAND
       DO J=1,NGP
         FSFCD(J)=FSFCD(J)+EMISFC*FLUX(J,JB)
       ENDDO
      ENDDO

C     3.4 Correction for "black" band (incl. surface reflection)

      EPS1=EPSLW*EMISFC
      DO J=1,NGP
        CORLW=EPS1*ST4A(J,NLEV,1)
        DFABS(J,NLEV)=DFABS(J,NLEV)-CORLW
        FSFCD(J)     =FSFCD(J)     +CORLW
      ENDDO

      IF (IMODE.EQ.-1) RETURN

C---  4. Emission ad absorption of longwave upward flux. 
C        For upward emission, a correction term depending on the      
C        local temperature gradient and on the layer transmissivity is  
C        subtracted from the average (full-level) emission of each layer. 
	
C     4.1  Surface

C     Black-body (or grey-body) emission 
      ESBC=EMISFC*SBC
      DO J=1,NGP
        TSQ=TS(J)*TS(J)
        FSFCU(J)=ESBC*TSQ*TSQ
      ENDDO

C     Entry point for upward-only mode (IMODE=1)
 410  CONTINUE

      DO J=1,NGP
        FSFC(J)=FSFCU(J)-FSFCD(J)
      ENDDO

      DO JB=1,NBAND
       DO J=1,NGP
         FLUX(J,JB)=FBAND(NINT(TS(J)),JB)*FSFCU(J)
     &              +REFSFC*FLUX(J,JB)
       ENDDO
      ENDDO
	
C     4.2  Troposphere

C     Correction for "black" band
      DO J=1,NGP
        DFABS(J,NLEV)=DFABS(J,NLEV)+EPSLW*FSFCU(J)
      ENDDO

      DO JB=1,NBAND
       DO K=NLEV,2,-1
        DO J=1,NGP
          EMIS=1.-TAU2(J,K,JB)
          BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)-EMIS*ST4A(J,K,2))
          DFABS(J,K)=DFABS(J,K)+FLUX(J,JB)
          FLUX(J,JB)=TAU2(J,K,JB)*FLUX(J,JB)+EMIS*BRAD
          DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
        ENDDO
       ENDDO
      ENDDO
	
C     4.3  Stratosphere

      K=1
      DO JB=1,2
       DO J=1,NGP
         EMIS=1.-TAU2(J,K,JB)
         BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)-EMIS*ST4A(J,K,2))
         DFABS(J,K)=DFABS(J,K)+FLUX(J,JB)
         FLUX(J,JB)=TAU2(J,K,JB)*FLUX(J,JB)+EMIS*BRAD
         DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
       ENDDO
      ENDDO

C     Correction for "black" band and polar night cooling

      DO J=1,NGP
        CORLW1=DSIG(1)*STRATC(J,2)*ST4A(J,1,1)+STRATC(J,1)
        CORLW2=DSIG(2)*STRATC(J,2)*ST4A(J,2,1)
        DFABS(J,1)=DFABS(J,1)-CORLW1
        DFABS(J,2)=DFABS(J,2)-CORLW2
        FTOP(J)   =CORLW1+CORLW2
      ENDDO

C     4.4  Outgoing longwave radiation 

      DO JB=1,NBAND
       DO J=1,NGP
         FTOP(J)=FTOP(J)+FLUX(J,JB)
       ENDDO
      ENDDO

C---						
      RETURN
      END


      SUBROUTINE RADSET
C--
C--   SUBROUTINE RADSET
C--
C--   Purpose: compute energy fractions in LW bands
C--            as a function of temperature
C--   Initialized common blocks: RADFIX

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Radiation constants
      include "com_radcon.h"

      EPS1=1.-EPSLW

      DO JTEMP=200,320
        FBAND(JTEMP,2)=(0.148-3.0e-6*(JTEMP-247)**2)*EPS1
        FBAND(JTEMP,3)=(0.356-5.2e-6*(JTEMP-282)**2)*EPS1
        FBAND(JTEMP,4)=(0.314+1.0e-5*(JTEMP-315)**2)*EPS1
        FBAND(JTEMP,1)=EPS1-(FBAND(JTEMP,2)+
     &                       FBAND(JTEMP,3)+FBAND(JTEMP,4))
      ENDDO

      DO JB=1,4
        DO JTEMP=100,199
          FBAND(JTEMP,JB)=FBAND(200,JB)
        ENDDO
        DO JTEMP=321,400
          FBAND(JTEMP,JB)=FBAND(320,JB)
        ENDDO
      ENDDO

C--						
      RETURN
      END
