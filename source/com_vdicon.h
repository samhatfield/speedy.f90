C--
C--:  /VDICON/: Constants for vertical diffusion and shallow convection 
C--:            (initial. in INPHYS)
C--:   TRSHC  = relaxation time (in hours) for shallow convection
C--:   TRVDI  = relaxation time (in hours) for moisture diffusion
C--:   TRVDS  = relaxation time (in hours) for super-adiab. conditions
C--:   REDSHC = reduction factor of shallow conv. in areas of deep conv.
C--:   RHGRAD = maximum gradient of relative humidity (d_RH/d_sigma)
C--:   SEGRAD = minimum gradient of dry static energy (d_DSE/d_phi)

      COMMON /VDICON/ TRSHC, TRVDI, TRVDS, REDSHC, RHGRAD, SEGRAD
