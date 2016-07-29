
C--
C--   /ICOUPL/ Flags to set coupling options (see doc_instep.txt)
      common /ICOUPL/ icland, icsea, icice, isstan

      integer icland                    ! flag for land coupling
      integer icsea                     ! flag for sea (SST) coupling
      integer icice                     ! flag for sea-ice coupling
      integer isstan                    ! flag for observed SST anomaly
