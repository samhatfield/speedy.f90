
      SUBROUTINE XS_RDF (TT1,TT2,IVM)
C--
c--   SUBROUTINE XS_RDF (TT1,TT2,IVM)
C--
C--   Purpose: compute zonal-mean cross-sec. of random diabatic forcing
C--   Input: TT1, TT2 = diabatic heating fields
C--          IVM      = index of vertical mode (1 or 2)
C--   Modified common block: RANDF
C--
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_physcon.h"

      include "com_randfor.h"

      real tt1(nlon,nlat,nlev), tt2(nlon,nlat,nlev)

      real rand1(0:nlat+1)

      rnlon = 1./float(nlon)
      pigr2 = 4.*asin(1.)

C--   1. Compute cross sections

      do k=1,nlev

         if (ivm.eq.1) then
            rnsig = rnlon
         else
            rnsig = rnlon*sin(pigr2*sig(k))
         endif

         do j=1,nlat
           randfv(j,k,ivm) = 0.
           do i=1,nlon
              randfv(j,k,ivm) = randfv(j,k,ivm)+tt1(i,j,k)+tt2(i,j,k)
           enddo
           randfv(j,k,ivm) = randfv(j,k,ivm)*rnsig 
         enddo

      enddo

C--   2. Perform smoothing in latitude

      do nsmooth=1,2
        do k=1,nlev

          do j=1,nlat
             rand1(j) = randfv(j,k,ivm)
          enddo
          rand1(0) = rand1(2)
          rand1(nlat+1) = rand1(nlat-1)
             
          do j=1,nlat
             randfv(j,k,ivm) = 0.5*rand1(j)+0.25*(rand1(j-1)+rand1(j+1))
          enddo

        enddo
      enddo

      return
      end

      SUBROUTINE SETRDF (TT_RDF)
C--
c--   SUBROUTINE SETRDF (TT_RDF)
C--
C--   Purpose: compute 3-D pattern of random diabatic forcing
C--   Output: TT_RDF = random diabatic forcing
C--
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_randfor.h"

      real tt_rdf(nlon,nlat,nlev)

      do k=1,nlev
        do j=1,nlat
           do i=1,nlon
               tt_rdf(i,j,k) = randfh(i,j,1)*randfv(j,k,1)+
     &                         randfh(i,j,2)*randfv(j,k,2)
           enddo
         enddo
       enddo

       return
       end

