!  /LFLAG1/: Logical flags to activate processes throughout 
!            the integration (initial. in INDYNS)
!   LPPRES = Flag to post-process upper-air fields
!            on pressure levels (.false. for model level p.p.)
!   LCO2 = Flag for CO2 optical thickness increase

       LOGICAL LPPRES, LCO2
       COMMON /LFLAG1/ LPPRES, LCO2
!
!   /LFLAG2/: Logical flags to activate processes in selected 
!             time steps (initial. in STEPONE, updated in STLOOP)
!    LRADSW = Flag for shortwave radiation routine
!    LRANDF = Flag for random diabatic forcing

       LOGICAL LRADSW, LRANDF
       COMMON /LFLAG2/ LRADSW, LRANDF
