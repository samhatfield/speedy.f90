
      SUBROUTINE FORINT (NGP,IMON,FMON,FOR12,FOR1)  

C--   Aux. routine FORINT : linear interpolation of monthly-mean forcing

      REAL FOR12(NGP,*), FOR1(NGP)

      IF (FMON.LE.0.5) THEN
        IMON2 = IMON-1
        IF (IMON.EQ.1) IMON2 = 12
        WMON = 0.5-FMON
      ELSE
        IMON2 = IMON+1
        IF (IMON.EQ.12) IMON2 = 1
        WMON = FMON-0.5
      ENDIF

      DO J=1,NGP
        FOR1(J) = FOR12(J,IMON)+WMON*(FOR12(J,IMON2)-FOR12(J,IMON))
      ENDDO
C--
      RETURN
      END

      subroutine FORIN5 (ngp,imon,fmon,for12,for1)

C--   Aux. routine FORIN5 : non-linear, mean-conserving interpolation 
C--                         of monthly-mean forcing fields

      real for12(ngp,12), for1(ngp)

      im2 = imon-2
      im1 = imon-1
      ip1 = imon+1
      ip2 = imon+2

      if (im2.lt.1)  im2 = im2+12
      if (im1.lt.1)  im1 = im1+12
      if (ip1.gt.12) ip1 = ip1-12
      if (ip2.gt.12) ip2 = ip2-12
 
      c0 = 1./12.
      t0 = c0*fmon
      t1 = c0*(1.-fmon)
      t2 = 0.25*fmon*(1-fmon)

      wm2 =        -t1   +t2
      wm1 =  -c0 +8*t1 -6*t2
      w0  = 7*c0      +10*t2     
      wp1 =  -c0 +8*t0 -6*t2
      wp2 =        -t0   +t2 

      do j=1,ngp
        for1(j) = wm2*for12(j,im2)+wm1*for12(j,im1)
     &           + w0*for12(j,imon)
     &           +wp1*for12(j,ip1)+wp2*for12(j,ip2)
      enddo

      return
      end

