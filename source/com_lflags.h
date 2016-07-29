C--:
C--:  /LFLAG1/: Logical flags to activate processes throughout 
C--:            the integration (initial. in INDYNS)
C--:   LPPRES = Flag to post-process upper-air fields
C--:            on pressure levels (.false. for model level p.p.)
C--:   LCO2 = Flag for CO2 optical thickness increase

       COMMON /LFLAG1/ LPPRES, LCO2
C--
C--   /LFLAG2/: Logical flags to activate processes in selected 
C--             time steps (initial. in STEPONE, updated in STLOOP)
C--    LRADSW = Flag for shortwave radiation routine
C--    LRANDF = Flag for random diabatic forcing

       COMMON /LFLAG2/ LRADSW, LRANDF

               LOGICAL LPPRES, LCO2, 
     &         LRADSW, LRANDF
