
C--
C--   /CLPV_L/ Input and output land variables exchanged by coupler
      common /CPLV_L/ vland_input, vland_output

      real vland_input(ngp,4)            ! land model input variables
      real vland_output(ngp,2)           ! land model output variables
