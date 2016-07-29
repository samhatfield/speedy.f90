
C--
C--   /SMASKS/ Sea masks
      common /SMASKS/ fmask_s, bmask_s, deglat_s
      real fmask_s(ix,il)                ! fraction of sea
      real bmask_s(ix,il)                ! binary sea mask
      real deglat_s(il)                  ! grid latitudes 
C--
C--   /SCLIM/ Monthly-mean climatological fields over sea
      common /SCLIM/ sst12, sice12
      real  sst12(ix,il,12)             ! sea/ice surface temperature
      real sice12(ix,il,12)             ! sea ice fraction
C--
C--   /SSTANOM/ SST anomaly fields
      common /SSTANOM/ sstan3
      real sstan3(ix,il,3)              ! sst anomaly in 3 consecutive months
C--
C--   /SMCLIM/ Climatological fields from model output
      common /SMCLIM/ hfseacl, sstom12
      real hfseacl(ix,il)               ! annual-mean heat flux into sea sfc.
      real sstom12(ix,il,12)            ! ocean model SST climatology
