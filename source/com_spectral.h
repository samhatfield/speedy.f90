C
C--   /CSPEC/: initial. in PARMTR
      COMMON /CSPEC/ EL2(MX,NX),ELM2(MX,NX),EL4(MX,NX), 
     *               TRFILT(MX,NX),
     *               L2(MX,NX),LL(MX,NX),MM(MX),
     *               nsh2(nx)
C
C--   /GAUSS/: initial. in PARMTR
      COMMON /GAUSS/ SIA(IY),COA(IY),WT(IY),WGHT(IY),
     *               COSG(IL),COSGR(IL),cosgr2(il)
C
C--   /GRADC/: initial. in PARMTR
      COMMON /GRADC/ GRADX(MX),GRADYM(MX,NX),GRADYP(MX,NX)
C
C--   /LGND/: initial. in PARMTR
      COMMON /LGND/ SQRHLF,CONSQ(MXP),EPSI(MXP,NXP),REPSI(MXP,NXP),
     *              EMM(MXP),ELL(MXP,NXP)
C
C--   /POLS/: initial. in PARMTR
      COMMON /POLS/ POLY(MX,NX)
C
C--   /POLS1/: initial. in PARMTR
      COMMON /POLS1/ CPOL(MX2,NX,IY)
C
C--   /UV/: initial. in PARMTR
      COMMON /UV/ UVDX(MX,NX),UVDYM(MX,NX),UVDYP(MX,NX)
C
C--   /VD/: initial. in PARMTR
      COMMON /VD/ VDDYM(MX,NX),VDDYP(MX,NX)
