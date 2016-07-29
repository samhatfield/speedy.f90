
C--
C--   /LAND_MC/ Constant parameters and fields in land sfc. model
      common /LAND_MC/ rhcapl, cdland

      real  rhcapl(nlon,nlat)           ! 1./heat_capacity (land)
      real  cdland(nlon,nlat)           ! 1./dissip_time (land)


