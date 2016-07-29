      SUBROUTINE INIRDF (INDRDF)
C--
C--   SUBROUTINE INIRDF (INDRDF)
C--
C--   Purpose : Initialize random diabatic forcing 
C--   Input :   INIRDF = index of forcing perturbation
C--   Initialized common blocks: RANDF
C--
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_physcon.h"

      include "com_randfor.h"

      real    redgrd(0:36,0:18), randf2(nlon,nlat),
     &        rnlon(0:18), colat(nlat)

      integer nlonrg(0:18)

      data nlonrg / 1,  6, 12, 18, 24, 28, 32, 34, 36, 36, 
     &             36, 34, 32, 28, 24, 18, 12,  6,  1 /


c     rms aplitude of (non-null) horizontal perturbation
      ampl = 0.5

c     frequency of grid points with null perturbation
      freq0 = 0.

c     ntrfor = spectral truncation of random forcing
      ntrfor = 18

C--   1. Initialization 

      iseed = -abs(indrdf)
      if (indrdf.lt.0.) ampl = -ampl

      do i=1,20
         fran = RAN1(iseed)
      enddo

      do jlat=0,18
         rnlon(jlat) = float(nlonrg(jlat))/float(nlon)
      enddo

      rdeg = 9./asin(1.)
      do j=1,nlat
         colat(j)=rdeg*asin(slat(j))+9.
      enddo

      do nf=1,2

C--     2. Fill reduced grid with normally-distributed random numbers

        do jlat=0,18

           call GAUSTS (nlonrg(jlat),0.,ampl,0.,0,iseed,
     &                  redgrd(1,jlat))

           if (freq0.gt.0.) then
             do jlon=1,nlonrg(jlat)
                fran = RAN1(iseed)
                if (fran.lt.freq0) redgrd(jlon,jlat) = 0.
             enddo
           endif

           redgrd(0,jlat) = redgrd(nlonrg(jlat),jlat)

        enddo

C--     3. Interpolate random field to gaussian grid

        do j=1,nlat

           jlat1 = int(colat(j))
           jlat2 = jlat1+1

           do i=1,nlon

              rlon  = (i-1)*rnlon(jlat1)
              jlon  = int(rlon)
              flat1 = redgrd(jlon,jlat1)+(rlon-jlon)*
     &                (redgrd(jlon+1,jlat1)-redgrd(jlon,jlat1))

              rlon  = (i-1)*rnlon(jlat2)
              jlon1 = int(rlon)
              flat2 = redgrd(jlon,jlat2)+(rlon-jlon)*
     &                (redgrd(jlon+1,jlat2)-redgrd(jlon,jlat2))

              randf2(i,j) = flat1+(colat(j)-jlat1)*(flat2-flat1)

           enddo

        enddo

C--     4. Spectral filter of gaussian-grid field

        call TRUNCG (ntrfor,randf2,randfh(1,1,nf))

      enddo

      return
      end


      SUBROUTINE GAUSTS (nt,av,sd,ac,ndis,iseed,
     &                   ts)

C--   SUBROUTINE GAUSTS (nt,av,sd,ac,ndis,iseed,
C--  &                   ts)
C--   Computes a gaussian-dist. time series (ts) of (nt) random values, with
C--   assigned average (av), stand. dev. (sd), and lag-1 autocorrelation (ac). 
C--   Autocor. may be discontinued at the limits between (ndis) sub-series. 
C--   Uses function RAN1 to generate uniform deviates from seed (iseed)
C--   Adapted from Numerical Recipes, Chapter 7.2 

      real ts(nt)

C--   1. Generate a time series of (nt) gaussian deviates

      do j=2,nt,2

 1      continue
        v1=2.*RAN1(iseed)-1.
        v2=2.*RAN1(iseed)-1.
        rsq=v1*v1+v2*v2
        if (rsq.gt.1.or.rsq.eq.0.) go to 1

        fact=sqrt(-2.*log(rsq)/rsq)
        ts(j-1)=v1*fact
        ts(j)  =v2*fact

      enddo

C--   2. Introduce autocorrelation (if requested)

      if (ac.ne.0.) then

        nt2=nt/max(1,ndis)
        sd2=sqrt(1.-ac*ac)

        j=0
        do jd=1,ndis
          j=j+1
          do j2=2,nt2
            j=j+1
            ts(j)=ac*ts(j-1)+sd2*ts(j)
          enddo
        enddo

      endif

C--   3. Set assigned average and standard deviation

      do j=1,nt
        ts(j)=sd*ts(j)+av
      enddo

      return
      end

      FUNCTION RAN1 (IDUM)
C--
C--   FUNCTION RAN1 (IDUM)
C--   Returns a uniform random deviate between 0.0 and 1.0
C--   Set IDUM to any negative value to (re)initialize the sequence
C--   From Numerical Recipes, Chapter 7.1 
C--
      PARAMETER ( IM=714025, IA=1366, IC=150889 )
      PARAMETER ( RM=1./IM )

      INTEGER IR(97)

      DATA IY/ -1/, IR/ 97*0/

      IF (IDUM.LT.0.OR.IY.LT.0) THEN

C       Initialize the shuffle array

        IDUM=MOD(IC+ABS(IDUM),IM)

        DO J=1,97
          IDUM=MOD(IA*IDUM+IC,IM)
          IR(J)=IDUM
        ENDDO

        IDUM=MOD(IA*IDUM+IC,IM)
        IY=IDUM

      ENDIF

C     Get one integer number from the shuffle table 
      J=1+(97*IY)/IM
c      IF (J.GT.97.OR.J.LT.1) STOP ' Error in random no. generator'
      IY=IR(J)

C     Turn the selected integer into a real no. between 0 and 1
      RAN1=IY*RM

C     Replace the selected integer with another random integer
      IDUM=MOD(IA*IDUM+IC,IM)
      IR(J)=IDUM

      RETURN
      END
