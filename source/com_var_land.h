
C--
C--   /LDAYCL/ Daily observed climatological fields over land
      common /LDAYCL/ stlcl_ob, snowdcl_ob, soilwcl_ob
      real   stlcl_ob(ngp)              ! clim. land sfc. temperature 
      real snowdcl_ob(ngp)              ! clim. snow depth (water equiv)
      real soilwcl_ob(ngp)              ! clim. soil water availability
C--
C--   /LVAR_AM/ Land sfc. fields used by atmospheric model
      common /LVAR_AM/ stl_am, snowd_am, soilw_am
      real   stl_am(ngp)                 ! land sfc. temperature
      real snowd_am(ngp)                 ! snow depth (water equiv)
      real soilw_am(ngp)                 ! soil water availability
C--
C--   /LVAR_LM/ Land sfc. fields from land model
      common /LVAR_LM/ stl_lm
      real   stl_lm(ngp)                 ! land-model sfc. temperature 

