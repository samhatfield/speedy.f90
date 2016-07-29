************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      PARAMETER (NMAX=100,TINY=1.0E-20)

      DIMENSION A(NP,NP),INDX(N),VV(NMAX)

      D=1.

      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J)) 
   11   CONTINUE
        IF(AAMAX.EQ.0.) PAUSE 'singular'
        VV(I)=1./AAMAX
   12 CONTINUE

      DO 19 J=1,N

        IF(J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF(I.GT.1) THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
   13         CONTINUE
              A(I,J)=SUM
            ENDIF
   14     CONTINUE
        ENDIF

        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF(J.GT.1) THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
   15       CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
         IF(DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
   16   CONTINUE

        IF(J.NE.IMAX) THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
   17     CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF

        INDX(J)=IMAX
        IF(J.NE.N) THEN
          IF(A(J,J).EQ.0) A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
   18     CONTINUE
        ENDIF

   19 CONTINUE

      IF(A(N,N).EQ.0.) A(N,N)=TINY

      RETURN
      END
***************************************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      DIMENSION A(NP,NP),INDX(N),B(N)

      II=0

      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF(II.NE.0) THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
   11     CONTINUE
        ELSE IF(SUM.NE.0) THEN
          II=I
        ENDIF
        B(I)=SUM
   12 CONTINUE

      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N) THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
   13     CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
   14 CONTINUE

      RETURN
      END
****************************************************************
      SUBROUTINE INV(A,Y,INDX,N)

      DIMENSION A(N,N),Y(N,N),INDX(N) 

      DO 12 I=1,N
        DO 11 J=1,N
          Y(I,J)=0.
   11   CONTINUE
        Y(I,I)=1.
   12 CONTINUE

      CALL LUDCMP(A,N,N,INDX,D)

      DO 13 J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
   13 CONTINUE

      RETURN
      END
****************************************************************
