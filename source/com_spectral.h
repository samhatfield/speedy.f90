!
!--   /CSPEC/: initial. in PARMTR
      REAL EL2, ELM2, EL4, TRFILT
      INTEGER L2, LL, MM, NSH2
      COMMON /CSPEC/ EL2(MX,NX),ELM2(MX,NX),EL4(MX,NX),                 &
     &               TRFILT(MX,NX),                                     &
     &               L2(MX,NX),LL(MX,NX),MM(MX),                        &
     &               nsh2(nx)
!
!--   /GAUSS/: initial. in PARMTR
      REAL SIA, COA, WT, WGHT, COSG, COSGR, COSGR2
      COMMON /GAUSS/ SIA(IY),COA(IY),WT(IY),WGHT(IY),                   &
     &               COSG(IL),COSGR(IL),cosgr2(il)
!
!--   /GRADC/: initial. in PARMTR
      REAL GRADX, GRADYM, GRADYP
      COMMON /GRADC/ GRADX(MX),GRADYM(MX,NX),GRADYP(MX,NX)
!
!--   /LGND/: initial. in PARMTR
      REAL SQRHLF, CONSQ, EPSI, REPSI, EMM, ELL
      COMMON /LGND/ SQRHLF,CONSQ(MXP),EPSI(MXP,NXP),REPSI(MXP,NXP),     &
     &              EMM(MXP),ELL(MXP,NXP)
!
!--   /POLS/: initial. in PARMTR
      REAL POLY
      COMMON /POLS/ POLY(MX,NX)
!
!--   /POLS1/: initial. in PARMTR
      REAL CPOL
      COMMON /POLS1/ CPOL(MX2,NX,IY)
!
!--   /UV/: initial. in PARMTR
      REAL UVDX, UVDYM, UVDYP
      COMMON /UV/ UVDX(MX,NX),UVDYM(MX,NX),UVDYP(MX,NX)
!
!--   /VD/: initial. in PARMTR
      REAL VDDYM, VDDYP
      COMMON /VD/ VDDYM(MX,NX),VDDYP(MX,NX)
