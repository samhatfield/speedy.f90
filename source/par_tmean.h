!
!   Post-processing and output parameters (time-mean and variance fields) 
!   NS2D_1 : no. of 2-d time-mean fields, incremented at post-proc. steps 
!   NS2D_2 : no. of 2-d time-mean fields, incremented every step (fluxes)
!   NS2D_D1: no. of 2-d daily-mean fields, incremented at post-proc. steps 
!   NS2D_D2: no. of 2-d daily-mean fields, incremented every step (fluxes)
!   NS3D1  : no. of 3-d time-mean model variables
!   NS3D2  : no. of 3-d time-mean variances and covariances
!   NS3D3  : no. of 3-d time-mean diabatic heating fields

      INTEGER, PARAMETER :: NS2D_1 = 18
      INTEGER, PARAMETER :: NS2D_2 = 12
      INTEGER, PARAMETER :: NS2D_D1 = 8
      INTEGER, PARAMETER :: NS2D_D2 = 7
      INTEGER, PARAMETER :: NS3D1 = 9
      INTEGER, PARAMETER :: NS3D2 = 6
      INTEGER, PARAMETER :: NS3D3 = 5
      INTEGER, PARAMETER :: NS3D = NS3D1 + NS3D2 + NS3D3
