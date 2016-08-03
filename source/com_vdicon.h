!
!  /VDICON/: Constants for vertical diffusion and shallow convection 
!            (initial. in INPHYS)
!   TRSHC  = relaxation time (in hours) for shallow convection
!   TRVDI  = relaxation time (in hours) for moisture diffusion
!   TRVDS  = relaxation time (in hours) for super-adiab. conditions
!   REDSHC = reduction factor of shallow conv. in areas of deep conv.
!   RHGRAD = maximum gradient of relative humidity (d_RH/d_sigma)
!   SEGRAD = minimum gradient of dry static energy (d_DSE/d_phi)

      REAL TRSHC, TRVDI, TRVDS, REDSHC, RHGRAD, SEGRAD
      COMMON /VDICON/ TRSHC, TRVDI, TRVDS, REDSHC, RHGRAD, SEGRAD
