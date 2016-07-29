
C--
C--   /SEAFLUX/ Fluxes at sea surface (all downward, except evaporation)
      common /SEAFLUX/ prec_s, snowf_s, evap_s, ustr_s, vstr_s,
     &                 ssr_s, slr_s, shf_s, ehf_s, hflux_s, hflux_i

      real  prec_s(ngp)                 ! precipitation (sea)
      real snowf_s(ngp)                 ! snowfall (sea)
      real  evap_s(ngp)                 ! evaporation (sea)
      real  ustr_s(ngp)                 ! u-wind stress (sea)
      real  vstr_s(ngp)                 ! v-wind stress (sea)
      real   ssr_s(ngp)                 ! sfc short-wave radiation (sea)
      real   slr_s(ngp)                 ! sfc long-wave radiation (sea)
      real   shf_s(ngp)                 ! sensible heat flux (sea)
      real   ehf_s(ngp)                 ! latent heat flux (sea)
      real hflux_s(ngp)                 ! net heat flux into sea sfc.
      real hflux_i(ngp)                 ! net heat flux into sea-ice sfc.
