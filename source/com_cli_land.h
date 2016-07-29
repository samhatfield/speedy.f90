C--
C--   /LMASKS/ Land masks
      common /LMASKS/ fmask_l, bmask_l
      real fmask_l(ix,il)                ! fraction of land
      real bmask_l(ix,il)                ! binary land mask
C--
C--   /LCLIM/ Monthly-mean climatological fields over land
      common /LCLIM/ stl12, snowd12, soilw12
      real   stl12(ix,il,12)             ! land surface temperature
      real snowd12(ix,il,12)             ! snow depth (water equiv.)
      real soilw12(ix,il,12)             ! soil water availability
