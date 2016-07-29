
C--
C--   /LANDFLUX/ Fluxes at land surface (all downward, except evaporation)
      common /LANDFLUX/ prec_l, snowf_l, evap_l, ustr_l, vstr_l,
     &                  ssr_l, slr_l, shf_l, ehf_l, hflux_l

      real  prec_l(ngp)                 ! precipitation (land)
      real snowf_l(ngp)                 ! snowfall (land)
      real  evap_l(ngp)                 ! evaporation (land)
      real  ustr_l(ngp)                 ! u-wind stress (land)
      real  vstr_l(ngp)                 ! v-wind stress (land)
      real   ssr_l(ngp)                 ! sfc short-wave radiation (land)
      real   slr_l(ngp)                 ! sfc long-wave radiation (land)
      real   shf_l(ngp)                 ! sensible heat flux (land)
      real   ehf_l(ngp)                 ! latent heat flux (land)
      real hflux_l(ngp)                 ! net heat flux into land sfc.
