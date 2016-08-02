!
!   /LSMASK/ land-sea masks (initial. in INBCON)
      common /LSMASK/ fmask, fmask1
      real fmask(ix,il)            ! original (fractional) land-sea mask
      real fmask1(ix,il)           ! model-defined land fraction
!									
!   /SFCFIX/ Time invariant surface fields 
!            (initial. in INBCON, phis0 initial. in INVARS)
      common /SFCFIX/ phi0, phis0, alb0
      real phi0(ix,il)             ! unfiltered surface geopotential
      real phis0(ix,il)            ! spectrally-filtered sfc. geopotential
      real alb0(ix,il)             ! bare-land annual-mean albedo
!
!   /SOILMP/ Soil moisture parameters (initial. in INPHYS)
!:   SWCAP  = Soil wetness at field capacity (volume fraction)
!:   SWWIL  = Soil wetness at wilting point  (volume fraction)
!:   SD2SC  = Snow depth (mm water) corresponding to snow cover = 1
      common /SOILMP/ swcap, swwil, sd2sc
      real            swcap, swwil, sd2sc

