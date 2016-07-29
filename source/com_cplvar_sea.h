
C--
C--   /CLPV_S/ Input and output sea variables exchanged by coupler
      common /CPLV_S/ vsea_input, vsea_output

      real vsea_input(ngp,8)            ! ocean model input variables
      real vsea_output(ngp,3)           ! ocean model output variables
