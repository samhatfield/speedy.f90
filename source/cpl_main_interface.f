
      SUBROUTINE INI_COUPLER (istart)
C--
C--   SUBROUTINE INI_COUPLER (istart)
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_cli_land.h"
      include "com_cli_sea.h"

      include "com_surfcon.h"

C--   1.1 Initialize land model constants

      CALL LAND_MODEL_INIT (fmask_l,alb0)

C--   1.2 Initialize land model variables

      CALL INI_LAND (istart)

C--   2.1 Initialize sea and ice model constants

      CALL SEA_MODEL_INIT (fmask_s,deglat_s)

C--   2.2 Initialize sea and ice model variables

      CALL INI_SEA (istart)

      RETURN
      END


      SUBROUTINE AGCM_TO_COUPLER (jday)
C--
C--   SUBROUTINE AGCM_TO_COUPLER (jday)
C--

C--   1. Send fields to land model

      CALL ATM2LAND (jday)

C--   2. Send fields to sea and ice model

      CALL ATM2SEA (jday)

      RETURN
      END


      SUBROUTINE COUPLER_TO_AGCM (jday)
C--
C--   SUBROUTINE COUPLER_TO_AGCM (jday)
C--

C--   1. Get updated fields from land model

      CALL LAND2ATM (jday)

C--   2. Get updated fields from sea and ice model

      CALL SEA2ATM (jday)

      RETURN
      END
