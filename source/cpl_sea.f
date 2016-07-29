
      SUBROUTINE INI_SEA (istart)
C--
C--   SUBROUTINE INI_SEA (istart)
C-- 
C--   Input : istart = restart flag ( 0 = no, 1 = yes)

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

      include "com_cpl_flags.h"

      include "com_cli_sea.h"
      include "com_var_sea.h"

C--   1. Compute climatological fields for initial date

      CALL ATM2SEA (0)

C--   2. Initialize prognostic variables of ocean/ice model
C--      in case of no restart or no coupling

      if (istart.le.0) then

        sst_om(:)  = sstcl_ob(:)      ! SST 
        tice_om(:) = ticecl_ob(:)     ! sea ice temperature
        sice_om(:) = sicecl_ob(:)     ! sea ice fraction

      endif

      if (icsea.le.0) sst_om(:) = 0.

C--   3. Compute additional sea/ice variables

      wsst_ob(:) = 0.
      if (icsea.ge.4) call SEA_DOMAIN ('elnino',deglat_s,wsst_ob)

      CALL SEA2ATM (0)

      RETURN
      END


      SUBROUTINE ATM2SEA (jday)
C--
C--   SUBROUTINE ATM2SEA (jday)
C-- 
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

      include "com_date.h"
      include "com_cpl_flags.h"

      include "com_cli_sea.h" 
      include "com_var_sea.h"
      include "com_flx_sea.h"

      include "com_cplvar_sea.h"

      real fmasks(ngp)                  ! sea fraction
      equivalence (fmasks,fmask_s)

      real hfyearm(ngp)                 ! annual mean heat flux into the ocean
      equivalence (hfyearm,hfseacl)


C--   1. Interpolate climatological fields and obs. SST anomaly
C--      to actual date

C     Climatological SST
      call FORIN5 (ngp,imont1,tmonth,sst12,sstcl_ob)

C     Climatological sea ice fraction
      call FORINT (ngp,imont1,tmonth,sice12,sicecl_ob)

C     SST anomaly
      if (isstan.gt.0) then 
         if (iday.eq.1.and.jday.gt.0) call OBS_SSTA
         call FORINT (ngp,2,tmonth,sstan3,sstan_ob)
      endif

C     Ocean model climatological SST
      if (icsea.ge.3) then
         call FORIN5 (ngp,imont1,tmonth,sstom12,sstcl_om)
      endif

C     Adjust climatological fields over sea ice

c     SST at freezing point
      sstfr = 273.2-1.8

      do j=1,ngp

         sstcl0 = sstcl_ob(j)

         if (sstcl_ob(j).gt.sstfr) then

           sicecl_ob(j) = min(0.5,sicecl_ob(j))
           ticecl_ob(j) = sstfr
           if (sicecl_ob(j).gt.0.)
     &        sstcl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/(1.-sicecl_ob(j))

         else

           sicecl_ob(j) = max(0.5,sicecl_ob(j))
           ticecl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/sicecl_ob(j)
c           ticecl_ob(j) = sstcl_ob(j)
           sstcl_ob(j)  = sstfr

         endif

         if (icsea.ge.3) sstcl_om(j) = sstcl_om(j)+(sstcl_ob(j)-sstcl0)

      enddo

      if (jday.le.0) RETURN

C--   2. Set input variables for mixed-layer/ocean model

      if (icsea.gt.0.or.icice.gt.0) then

        VSEA_INPUT(:,1) = sst_om(:)
        VSEA_INPUT(:,2) = tice_om(:)
        VSEA_INPUT(:,3) = sicecl_ob(:)
c        VSEA_INPUT(:,4) = hflux_s(:)*fmasks(:)
c        VSEA_INPUT(:,5) = hflux_i(:)*fmasks(:)
        VSEA_INPUT(:,4) = hflux_s(:)
        VSEA_INPUT(:,5) = hflux_i(:)
        VSEA_INPUT(:,6) = sstcl_ob(:)
        VSEA_INPUT(:,7) = ticecl_ob(:)
c        VSEA_INPUT(:,8) = hfyearm(:)*fmasks(:)
        VSEA_INPUT(:,8) = hfyearm(:)

      endif

C--   3. Call message-passing routines to send data (if needed)

      RETURN
      END


      SUBROUTINE SEA2ATM (jday)
C--
C--   SUBROUTINE SEA2ATM (jday)
C-- 
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

      include "com_cpl_flags.h"

      include "com_var_sea.h"

      include "com_cplvar_sea.h"

      if (jday.gt.0.and.(icsea.gt.0.or.icice.gt.0)) then

C--   1. Run ocean mixed layer or 
C--      call message-passing routines to receive data from ocean model

         CALL SEA_MODEL 

C--   2. Get updated variables for mixed-layer/ocean model

         sst_om(:)   = VSEA_OUTPUT(:,1)      ! SST
         tice_om(:)  = VSEA_OUTPUT(:,2)      ! sea ice temperature 
         sice_om(:)  = VSEA_OUTPUT(:,3)      ! sea ice fraction

c        sice_om(:)  = sicecl_ob(:)

      endif

C--   3. Compute sea-sfc. anomalies and full fields for atm. model

C     3.1 SST

      sstan_am(:) = 0.

      if (icsea.le.1) then

         if (isstan.gt.0) sstan_am(:) = sstan_ob(:)

C        Use observed SST (climatological or full field)
         sst_am(:) = sstcl_ob(:) + sstan_am(:)
         
      else if (icsea.eq.2) then

C        Use full ocean model SST
         sst_am(:) = sst_om(:)

      else if (icsea.ge.3) then

C        Define SST anomaly from ocean model ouput and climatology 
         sstan_am(:) = sst_om(:) - sstcl_om(:)

C        Merge with observed SST anomaly in selected area
         if (icsea.ge.4) sstan_am(:) = sstan_am(:) +
     &                   wsst_ob(:)*(sstan_ob(:)-sstan_am(:))

C        Add observed SST climatology to model SST anomaly 
         sst_am(:) = sstcl_ob(:) + sstan_am(:)

      endif

C     3.2 Sea ice fraction and temperature

      if ( icice.gt.0 ) then

         sice_am(:) = sice_om(:)
         tice_am(:) = tice_om(:)

      else

         sice_am(:) = sicecl_ob(:)
         tice_am(:) = ticecl_ob(:)

      endif

      sst_am(:)  = sst_am(:)+
     &             sice_am(:)*(tice_am(:)-sst_am(:))
      ssti_om(:) = sst_om(:)+
     &             sice_am(:)*(tice_am(:)-sst_om(:))


      RETURN
      END


      SUBROUTINE REST_SEA (imode)
C--
C--   SUBROUTINE REST_SEA (imode)
C--
C--   Purpose : read/write sea variables from/to a restart file
C--   Input :   IMODE = 0 : read model variables from a restart file
C--                   = 1 : write model variables  to a restart file

C-- 
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

      include "com_cpl_flags.h"

      include "com_var_sea.h"

      real sst_c(ngp)              ! SST corrected for sea-ice values

      if (imode.eq.0) then

         read (3)  sst_om(:)       ! SST 
         read (3) tice_om(:)       ! sea ice temperature
         read (3) sice_om(:)       ! sea ice fraction

      else

C        write sea/ice model variables from coupled runs,
C        otherwise write fields used by atmospheric model

         sstfr = 273.2-1.8

         if (icsea.gt.0) then
            write (10) sst_om(:) 
         else
            sst_c(:) = max(sst_am(:),sstfr)
            write (10) sst_c(:)
         endif

         if (icice.gt.0) then
            write (10) tice_om(:) 
            write (10) sice_om(:) 
         else
            write (10) tice_am(:)
            write (10) sice_am(:) 
         endif

      endif

      RETURN
      END


      SUBROUTINE OBS_SSTA 
C--
C--   SUBROUTINE OBS_SSTA 
C--
C--   Purpose : update observed SST anomaly array
 
      include "atparam.h"

      include "com_cli_sea.h"

      real*4 r4inp(ix,il)
   
        do j = 1,il
          do i = 1,ix
             sstan3(i,j,1) = sstan3(i,j,2)
             sstan3(i,j,2) = sstan3(i,j,3)
          enddo
        enddo

        read(30,end=100) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
             sstan3(i,j,3) = r4inp(i,j)
          enddo
        enddo 

      CALL FORCHK (bmask_s,sstan3(1,1,3),ix*il,1,-50.,50.,0.)

      RETURN

 100  continue

      print *, ' WARNING: end-of-file reached on SST anomaly file'
      print *, ' SST anomaly will be kept constant'

C--
      RETURN
      END  


