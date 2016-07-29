C--
C--   /LSMASK/ land-sea masks (initial. in INBCON)
      common /LSMASK/ fmask, fmask1
      real fmask(ix,il)            ! original (fractional) land-sea mask
      real fmask1(ix,il)           ! model-defined land fraction
C--									
C--   /SFCFIX/ Time invariant surface fields 
C--            (initial. in INBCON, phis0 initial. in INVARS)
      common /SFCFIX/ phi0, phis0, alb0
      real phi0(ix,il)             ! unfiltered surface geopotential
      real phis0(ix,il)            ! spectrally-filtered sfc. geopotential
      real alb0(ix,il)             ! bare-land annual-mean albedo
C--
C--   /SOILMP/ Soil moisture parameters (initial. in INPHYS)
C--:   SWCAP  = Soil wetness at field capacity (volume fraction)
C--:   SWWIL  = Soil wetness at wilting point  (volume fraction)
C--:   SD2SC  = Snow depth (mm water) corresponding to snow cover = 1
      common /SOILMP/ swcap, swwil, sd2sc
      real            swcap, swwil, sd2sc

