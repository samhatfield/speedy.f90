
      SUBROUTINE DMFLUX (IADD)

C--   SUBROUTINE DMFLUX (IADD)
C--
C--   Purpose: Add up fluxes to provide daily averages 
C--            used in sea/land models and daily/time-mean output
C--   Input: IADD = 0 to initialize storage arrays to 0
C--               > 0 to increment arrays with current flux values  

      USE mod_tsteps, only: nsteps
      USE mod_atparam

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Parameters for post-processing arrays
      include "par_tmean.h"

      include "com_surfcon.h"

      include "com_physcon.h"
      include "com_physvar.h"

      include "com_radcon.h"

      include "com_flx_land.h"
      include "com_flx_sea.h"

      include "com_var_sea.h"

      include "com_tmean.h"

      real prec(ngp), difice(ngp)

      real fland(ngp)
      equivalence (fland,fmask1)

C--   1. Initialization

      if (iadd.le.0) then

C--     Set all daily-mean arrays to zero
C       NB storage arrays for time-mean output are re-initialized 
C       by subroutines TMOUT and DMOUT after write-up

        prec_l(:)  = 0.
        snowf_l(:) = 0.
        evap_l(:)  = 0.
        hflux_l(:) = 0.

        prec_s(:)  = 0.
        snowf_s(:) = 0.
        evap_s(:)  = 0.
        ustr_s(:)  = 0.
        vstr_s(:)  = 0.
        ssr_s(:)   = 0.
        slr_s(:)   = 0.
        shf_s(:)   = 0.
        ehf_s(:)   = 0.
        hflux_s(:) = 0.
        hflux_i(:) = 0.

        return

      endif

      rsteps = 1./real(nsteps)
      rstep1 = rsteps*0.001
      rstep2 = rsteps*alhc

C     SST at freezing point
      sstfr  = 273.2-1.8

      sstfr4 = sstfr**4
      esbc   = emisfc*sbc

c     total precipitation 
      prec(:) = precls(:)+precnv(:)

C--   2. Store fluxes over land (SI units, all heat fluxes downw.)

      prec_l(:) = prec_l(:) + prec(:)  *rstep1
      evap_l(:) = evap_l(:) + evap(:,1)*rstep1

      hflux_l(:) = hflux_l(:) + hfluxn(:,1)*rsteps

C--   3. Store fluxes over sea (SI units, all heat fluxes downw.)

      prec_s(:) = prec_s(:) + prec(:)  *rstep1
      evap_s(:) = evap_s(:) + evap(:,2)*rstep1

      ustr_s(:) = ustr_s(:) - ustr(:,2)*rsteps
      vstr_s(:) = vstr_s(:) - vstr(:,2)*rsteps

      ssr_s(:) = ssr_s(:) + ssr(:)   *rsteps
      slr_s(:) = slr_s(:) - slr(:)   *rsteps
      shf_s(:) = shf_s(:) - shf(:,2) *rsteps
      ehf_s(:) = ehf_s(:) - evap(:,2)*rstep2

c     difference in net (downw.) heat flux between ice and sea surface
      difice(:) = (albsea-albice)*ssrd(:)
     &            + esbc*(sstfr4-tice_am(:)**4)
     &            + shf(:,2)+evap(:,2)*alhc

      hflux_s(:) = hflux_s(:) + rsteps* hfluxn(:,2)
      hflux_i(:) = hflux_i(:) + rsteps*(hfluxn(:,2)
     &                                 +difice(:)*(1.-sice_am(:)))


C--   4.1 Store fluxes for daily-mean output

c     multiply net heat fluxes by land or sea fractions
      hfluxn(:,1) = hfluxn(:,1)*fland(:)
      hfluxn(:,2) = hfluxn(:,2)*(1.-fland(:))

c     surface water budget (in mm/day)
      save2d_d2(:,1) = save2d_d2(:,1) + prec(:)  *86.400
      save2d_d2(:,2) = save2d_d2(:,2) + evap(:,3)*86.400

c     surface momentum budget 
      save2d_d2(:,3) = save2d_d2(:,3) - ustr(:,3)
      save2d_d2(:,4) = save2d_d2(:,4) - vstr(:,3)

c     OLR
      save2d_d2(:,5) = save2d_d2(:,5) + olr(:)

c     surface energy budget
      save2d_d2(:,6) = save2d_d2(:,6) + hfluxn(:,1)
      save2d_d2(:,7) = save2d_d2(:,7) + hfluxn(:,2)

C--   4.2 Store fluxes for time-mean output

c     surface water budget (in mm/day)
      save2d_2(:,1) = save2d_2(:,1) + precls(:)*86.400
      save2d_2(:,2) = save2d_2(:,2) + precnv(:)*86.400
      save2d_2(:,3) = save2d_2(:,3) + evap(:,3)*86.400

c     surface momentum budget 
      save2d_2(:,4) = save2d_2(:,4) - ustr(:,3)
      save2d_2(:,5) = save2d_2(:,5) - vstr(:,3)

c     top-of-atmosphere energy budget
      save2d_2(:,6) = save2d_2(:,6) + tsr(:)
      save2d_2(:,7) = save2d_2(:,7) + olr(:)

c     surface energy budget
      save2d_2(:,8)  = save2d_2(:,8)  + ssr(:)
      save2d_2(:,9)  = save2d_2(:,9)  + slr(:)
      save2d_2(:,10) = save2d_2(:,10) + shf(:,3)
      save2d_2(:,11) = save2d_2(:,11) + hfluxn(:,1)
      save2d_2(:,12) = save2d_2(:,12) + hfluxn(:,2)

c     end of flux increment
      return

      end
