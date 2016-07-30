!   /ICOUPL/ Flags to set coupling options (see doc_instep.txt)
      COMMON /ICOUPL/ ICLAND, ICSEA, ICICE, ISSTAN

      INTEGER ICLAND                    ! flag for land coupling
      INTEGER ICSEA                     ! flag for sea (SST) coupling
      INTEGER ICICE                     ! flag for sea-ice coupling
      INTEGER ISSTAN                    ! flag for observed SST anomaly
