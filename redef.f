      SUBROUTINE REDEF(K,CD,ICON0)
C
C   Aquesta subrutina serveix per redefinir els coeficients que surten de la
C   subrutine NEWDER amb la opcio ICON=0, per tal de que quedin definits
C   segons les diferents opcions (ICON) de NEWDER diferents de 1 i de 8
C
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL*1 TORNA
      DIMENSION CD(*)
      IF(ICON0.EQ.0)RETURN
      KK=K*K
      KKP1=KK+1
      KO2=K/2
      KM1=K-1
      KP1=K+1
      KO2P1=KO2+1
      KO2P2=KO2+2
      ICON=ICON0
10    CONTINUE
      IF(ICON.EQ.4)IP=1
      IF(ICON.EQ.5)IP=-1
      IF(ICON.EQ.4.OR.ICON.EQ.5)THEN
        TORNA=.FALSE.
        I0=KO2P1
        I0M1K=(I0-1)*K
        DO I=1,KO2
          IM1K=(I-1)*K
          J=1
          IJ=J+IM1K
          J0=KO2P1-I+1
          J1=J0
          I0J0=J0+I0M1K
          CD(IJ)=CD(I0J0)
C          WRITE(6,1000)I,J,I0,J0
          DO J=2,K
            IJ=J+IM1K
            J0=J0+1
            I0J0=J0+I0M1K
            J1=J1-1
            I0J1=J1+I0M1K
            IF(J0.LE.K)THEN
              IF(J1.GE.1)THEN
                CD(IJ)=CD(I0J0)+IP*CD(I0J1)
C                WRITE(6,2000)I,J,I0,J0,I0,J1
              ELSE
                CD(IJ)=CD(I0J0)
C                WRITE(6,1000)I,J,I0,J0
              ENDIF
            ELSE
              CD(IJ)=0.D0
C              WRITE(6,3000)I,J
            ENDIF
          ENDDO
        ENDDO
      ELSE IF(ICON.EQ.2)THEN
        TORNA=.FALSE.
        DO I=K,KO2P2,-1
          I0=I-1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=K,2,-1
            J0=J-1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=1
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
      ELSE IF(ICON.EQ.3)THEN
C  Fa el mateix que ICON=2 per els primers punts i ICON=2
        ICON=2
        TORNA=.TRUE.
        DO I=1,KO2
          I0=I+1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=1,KM1
            J0=J+1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=K
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
      ELSE IF(ICON.EQ.6)THEN
        TORNA=.TRUE.
C  Fa el mateix que ICON=4 i ICON=2
        DO I=K,KO2P2,-1
          I0=I-1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=K,2,-1
            J0=J-1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=1
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
        ICON=4
      ELSE IF(ICON.EQ.7)THEN
        TORNA=.TRUE.
C  Fa el mateix que ICON=5 i ICON=2
        DO I=K,KO2P2,-1
          I0=I-1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=K,2,-1
            J0=J-1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=1
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
        ICON=5
      ELSE IF(ICON.EQ.12)THEN
        TORNA=.FALSE.
        DO I=KO2P2,K
          I0=I-1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=K,2,-1
            J0=J-1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=1
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
      ELSE IF(ICON.EQ.13)THEN
C  Fa el mateix que ICON=12 per els primers punts i ICON=12
        ICON=12
        TORNA=.TRUE.
        DO I=KO2,1,-1
          I0=I+1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=1,KM1
            J0=J+1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=K
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
      ELSE IF(ICON.EQ.16)THEN
C Fa el mateix que ICON=4 i ICON=12
        TORNA=.TRUE.
        DO I=KO2P2,K
          I0=I-1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=K,2,-1
            J0=J-1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=1
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
        ICON=4
      ELSE IF(ICON.EQ.17)THEN
C Fa el mateix que ICON=5 i ICON=12
        TORNA=.TRUE.
        DO I=KO2P2,K
          I0=I-1
          IM1K=(I-1)*K
          I0M1K=(I0-1)*K
          DO J=K,2,-1
            J0=J-1
            IJ=J+IM1K
            I0J0=J0+I0M1K
            CD(IJ)=CD(I0J0)
C            WRITE(6,1000)I,J,I0,J0
          ENDDO
          J=1
          IJ=J+IM1K
          CD(IJ)=0.D0
C          WRITE(6,3000)I,J
        ENDDO
        ICON=5
      ELSE IF(ICON.EQ.20)THEN
C     Condicions periodiques de Wigner_Seitz ( N+i <--> N-i )
        TORNA=.FALSE.
        I0=KO2P1
        I0M1K=(I0-1)*K
        DO I=K,KO2P2,-1
          IM1K=(I-1)*K
          J=K
          IJ=J+IM1K
          J0=KO2P1-I+K
          J1=J0
          I0J0=J0+I0M1K
          CD(IJ)=CD(I0J0)
C          WRITE(6,1000)I,J,I0,J0
          DO J=K-1,1,-1
            IJ=J+IM1K
            J0=J0-1
            I0J0=J0+I0M1K
            J1=J1+1
            I0J1=J1+I0M1K
            IF(J0.GE.1)THEN
              IF(J1.LE.K)THEN
                CD(IJ)=CD(I0J0)+CD(I0J1)
C                WRITE(6,2000)I,J,I0,J0,I0,J1
              ELSE
                CD(IJ)=CD(I0J0)
C                WRITE(6,1000)I,J,I0,J0
              ENDIF
            ELSE
              CD(IJ)=0.D0
C              WRITE(6,3000)I,J
            ENDIF
          ENDDO
        ENDDO
      ELSE IF(ICON.EQ.24)THEN
C     El mateix que ICON=20 i ICON=4
        I0=KO2P1
        I0M1K=(I0-1)*K
        DO I=K,KO2P2,-1
          IM1K=(I-1)*K
          J=K
          IJ=J+IM1K
          J0=KO2P1-I+K
          J1=J0
          I0J0=J0+I0M1K
          CD(IJ)=CD(I0J0)
C          WRITE(6,1000)I,J,I0,J0
          DO J=K-1,1,-1
            IJ=J+IM1K
            J0=J0-1
            I0J0=J0+I0M1K
            J1=J1+1
            I0J1=J1+I0M1K
            IF(J0.GE.1)THEN
              IF(J1.LE.K)THEN
                CD(IJ)=CD(I0J0)+CD(I0J1)
C                WRITE(6,2000)I,J,I0,J0,I0,J1
              ELSE
                CD(IJ)=CD(I0J0)
C                WRITE(6,1000)I,J,I0,J0
              ENDIF
            ELSE
              CD(IJ)=0.D0
C              WRITE(6,3000)I,J
            ENDIF
          ENDDO
        ENDDO
        TORNA=.TRUE.
        ICON=4
      ELSE IF(ICON.EQ.25)THEN
C     El mateix que ICON=20 i ICON=5
        I0=KO2P1
        I0M1K=(I0-1)*K
        DO I=K,KO2P2,-1
          IM1K=(I-1)*K
          J=K
          IJ=J+IM1K
          J0=KO2P1-I+K
          J1=J0
          I0J0=J0+I0M1K
          CD(IJ)=CD(I0J0)
C          WRITE(6,1000)I,J,I0,J0
          DO J=K-1,1,-1
            IJ=J+IM1K
            J0=J0-1
            I0J0=J0+I0M1K
            J1=J1+1
            I0J1=J1+I0M1K
            IF(J0.GE.1)THEN
              IF(J1.LE.K)THEN
                CD(IJ)=CD(I0J0)+CD(I0J1)
C                WRITE(6,2000)I,J,I0,J0,I0,J1
              ELSE
                CD(IJ)=CD(I0J0)
C                WRITE(6,1000)I,J,I0,J0
              ENDIF
            ELSE
              CD(IJ)=0.D0
C              WRITE(6,3000)I,J
            ENDIF
          ENDDO
        ENDDO
        TORNA=.TRUE.
        ICON=5
      ELSE 
        WRITE(6,*)' Aquest ICON encara no esta definit...',ICON
        TORNA=.FALSE.
      ENDIF
      IF(TORNA)GO TO 10
      RETURN
1000  FORMAT(' (',I2,',',I2,')<-->(',I2,',',I2,')')
2000  FORMAT(' (',I2,',',I2,')<-->(',I2,',',I2,')','(',I2,',',I2,')')
3000  FORMAT(' (',I2,',',I2,')<-->0.D0')
      END
