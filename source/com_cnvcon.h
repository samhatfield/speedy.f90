C--
C--:  /CNVCON/: Convection constants (initial. in INPHYS)
C--:   PSMIN  = minimum (norm.) sfc. pressure for the occurrence of convection
C--:   TRCNV  = time of relaxation (in hours) towards reference state
C--:   RHBL   = relative hum. threshold in the boundary layer
C--:   RHIL   = rel. hum. threshold in intermed. layers for secondary mass flux	
C--:   ENTMAX = max. entrainment as a fraction of cloud-base mass flux									
C--:   SMF    = ratio between secondary and primary mass flux at cloud-base 
 
      COMMON /CNVCON/ PSMIN, TRCNV, RHBL, RHIL, ENTMAX, SMF
