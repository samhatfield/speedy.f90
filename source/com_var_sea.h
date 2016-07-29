
C--
C--   /SDAYCL/ Daily observed climatological fields over sea
      common /SDAYCL/ sstcl_ob, sicecl_ob, ticecl_ob
      real  sstcl_ob(ngp)               ! observed clim. SST 
      real sicecl_ob(ngp)               ! clim. sea ice fraction
      real ticecl_ob(ngp)               ! clim. sea ice temperature
C--
C--   /SDAYAN/ Daily observed SST anomaly
      common /SDAYAN/ sstan_ob
      real sstan_ob(ngp)                ! observed sst anomaly 
C--
C--   /SDAYCL_OM/ Daily climatological fields from ocean model
      common /SDAYCL_OM/ sstcl_om
      real sstcl_om(ngp)                ! ocean model clim. SST 
C--
C--   /SVAR_AM/ Sea sfc. fields used by atmospheric model
      common /SVAR_AM/ sst_am, sstan_am, sice_am, tice_am
      real   sst_am(ngp)                 ! SST (full field)
      real sstan_am(ngp)                 ! SST anomaly 
      real  sice_am(ngp)                 ! sea ice fraction
      real  tice_am(ngp)                 ! sea ice temperature
C--
C--   /SVAR_OM/ Sea sfc. fields from ocean/sea-ice model
      common /SVAR_OM/ sst_om, sice_om, tice_om, ssti_om
      real   sst_om(ngp)                 ! ocean model SST
      real  sice_om(ngp)                 ! model sea ice fraction
      real  tice_om(ngp)                 ! model sea ice temperature
      real  ssti_om(ngp)                 ! model SST + sea ice temp.
C--
C--   /WOBSST/ Weight for obs. SST anomaly in coupled runs
      common /WOBSST/ wsst_ob
      real wsst_ob(ngp)                  ! weight mask for obs. SST 
