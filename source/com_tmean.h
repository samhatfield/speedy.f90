C--
C--   /TMSAVE/ : Arrays for the computation of time-means
C--              (initial./used by TMOUT, updated in TMINC and DMFLUX)
      common /TMSAVE/ save3d, save2d_1, save2d_2, rnsave

      real save3d  (ngp,nlev,ns3d)     ! 3-D fields saved every post-proc. step
      real save2d_1(ngp,ns2d_1)        ! 2-D fields saved every post-proc. step
      real save2d_2(ngp,ns2d_2)        ! 2-D fields saved every step (fluxes)
      real rnsave                      ! post-processing counter

C--   /DMSAVE/ : Arrays for the computation of daily-means
C--              (initial./used by DMOUT, updated in TMINC and DMFLUX)
      common /DMSAVE/ save2d_d1, save2d_d2
   
      real save2d_d1(ngp,ns2d_d1)     ! Daily output saved every post-proc step
      real save2d_d2(ngp,ns2d_d2)     ! Daily output saved every step (fluxes)
