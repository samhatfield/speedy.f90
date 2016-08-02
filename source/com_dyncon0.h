!
!:  /DYNC0/: Constants for initialization of dynamics (initial. in INDYNS)
!:   GAMMA  = Ref. temperature lapse rate (-dT/dz in deg/km)
!:   HSCALE = Ref. scale height for pressure (in km)
!:   HSHUM  = Ref. scale height for spec. humidity (in km)
!:   REFRH1 = Ref. relative humidity of near-surface air
!:   THD    = Max damping time (in hours) for hor. diffusion (del^6)
!:            of temperature and vorticity
!:   THDD   = Max damping time (in hours) for hor. diffusion (del^6)
!:            of divergence
!:   THDS   = Max damping time (in hours) for extra diffusion (del^2)
!:            in the stratosphere 
!:   TDRS   = Damping time (in hours) for drag on zonal-mean wind
!:            in the stratosphere 
      REAL GAMMA, HSCALE, HSHUM, REFRH1, THD, THDD, THDS, TDRS
      COMMON /DYNC0/ GAMMA, HSCALE, HSHUM, REFRH1, THD, THDD, THDS, TDRS
