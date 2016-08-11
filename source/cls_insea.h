
!     ocean mixed layer depth: d + (d0-d)*(cos_lat)^3
      depth_ml = 60.               ! High-latitude depth
      dept0_ml = 40.               ! Minimum depth (tropics)

!     sea-ice depth : d + (d0-d)*(cos_lat)^2
!      depth_ice = 1.8
      depth_ice = 2.5              ! High-latitude depth
      dept0_ice = 1.5              ! Minimum depth 

!     Dissipation time (days) for sea-surface temp. anomalies
      tdsst  = 90.

!     Dissipation time (days) for sea-ice temp. anomalies
!      tdice = 20.
      tdice = 30.

!     Heat flux coefficient at sea/ice interface [(W/m^2)/deg]
      beta = 1.

!     Minimum fraction of sea for the definition of anomalies
      fseamin = 1./3.

!     Geographical domain
!     note : more than one regional domain may be set .true.

      l_globe  =  .true.         ! global domain
      l_northe = .false.         ! Northern hem. oceans (lat > 20N)
      l_natlan = .false.         ! N. Atlantic (lat 20-80N, lon 100W-45E)
      l_npacif = .false.         ! N. Pacific  (lat 20-80N, lon 100E-100W)
      l_tropic = .false.         ! Tropics (lat 30S-30N)
      l_indian = .false.         ! Indian Ocean (lat 30S-30N, lon 30-120E)
