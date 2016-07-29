C--
C--:  /DYNC0/: Constants for initialization of dynamics (initial. in INDYNS)
C--:   GAMMA  = Ref. temperature lapse rate (-dT/dz in deg/km)
C--:   HSCALE = Ref. scale height for pressure (in km)
C--:   HSHUM  = Ref. scale height for spec. humidity (in km)
C--:   REFRH1 = Ref. relative humidity of near-surface air
C--:   THD    = Max damping time (in hours) for hor. diffusion (del^6)
C--:            of temperature and vorticity
C--:   THDD   = Max damping time (in hours) for hor. diffusion (del^6)
C--:            of divergence
C--:   THDS   = Max damping time (in hours) for extra diffusion (del^2)
C--:            in the stratosphere 
C--:   TDRS   = Damping time (in hours) for drag on zonal-mean wind
C--:            in the stratosphere 
									
      COMMON /DYNC0/ GAMMA, HSCALE, HSHUM, REFRH1, 
     &               THD, THDD, THDS, TDRS
