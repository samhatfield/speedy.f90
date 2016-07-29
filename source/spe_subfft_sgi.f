
      SUBROUTINE INIFFT

C     Initialize FFTs

      include "atparam.h"

      COMMON /FFTCOM/ COEFFT(IX+15)

      CALL DZFFTM1DUI (IX,COEFFT)

      RETURN
      END

*********************************************************************

      SUBROUTINE GRIDX(VARM,VORG,KCOS)

C     From Fourier coefficients to grid-point data

      include "atparam.h"

      PARAMETER (IMAX=IX+2 )

      include "com_spectral.h"

      COMMON /FFTCOM/ COEFFT(IX+15)

      REAL VORG(IX,IL), VARM(MX2,IL)
      REAL FVAR(IMAX,IL)

C     Copy Fourier coefficients into working array

      MX3=MX2+1
      DO J=1,IL
        DO I=1,MX2
          FVAR(I,J)=VARM(I,J)
        ENDDO
        DO I=MX3,IMAX
          FVAR(I,J)=0.
        ENDDO
      ENDDO

C     Inverse FFT

      CALL ZDFFTM1DU (+1,IX,IL,FVAR,1,IMAX,COEFFT)

C     Copy output into grid-point field, scaling by cos(lat) if needed

      IF (KCOS.EQ.1) THEN

        DO J=1,IL
          DO I=1,IX
            VORG(I,J)=FVAR(I,J)
          ENDDO
        ENDDO

      ELSE

        DO J=1,IL
          DO I=1,IX
            VORG(I,J)=FVAR(I,J)*COSGR(J)
          ENDDO
        ENDDO

      ENDIF

      RETURN
      END

******************************************************************

      SUBROUTINE SPECX(VORG,VARM)

C     From grid-point data to Fourier coefficients

      include "atparam.h"

      PARAMETER (IMAX=IX+2 )

      COMMON /FFTCOM/ COEFFT(IX+15)

      REAL VORG(IX,IL), VARM(MX2,IL)
      REAL FVAR(IMAX,IL)

C     Copy grid-point data into working array

      IM1=IX+1
      DO J=1,IL
        DO I=1,IX
          FVAR(I,J)=VORG(I,J)
        ENDDO
        FVAR(IM1 ,J)=0.
        FVAR(IMAX,J)=0.
      ENDDO

C     Direct FFT

      CALL DZFFTM1DU (-1,IX,IL,FVAR,1,IMAX,COEFFT)

C     Copy output into spectral field, dividing by no. of long.

      SCALE=1./FLOAT(IX)

      DO J=1,IL 
        DO M=1,MX2
          VARM(M,J)=FVAR(M,J)*SCALE
        ENDDO
      ENDDO

      RETURN
      END
