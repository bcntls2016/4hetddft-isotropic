      SUBROUTINE  NEWDER(K,M,N,H,F,DKF,ICON)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL COND
      DIMENSION F(*),DKF(*)
      COMMON/SDERA/ A(25792)
C
C     W. G. Bickley, Formulae for numerical differentiation, Math. Gaz.
C     25, 19-27, (1941).
C
C     Entradas:
C         K-----> Orden de la derivada
C         M-----> Numero de puntos con los que queremos calcular
C                 dicha derivada
C         N-----> Numero de puntos en donde conocemos la funcion
C         H-----> Intervalo que se ha utilizado para el calculo de
C                 la funcion
C         F-----> Vector de entrada de la funcion
C         DKF---> Vector de salida de la derivada (para ICON=0)
C         ICON--> Parametro de control
C
C     Funciones especiales:
C      ICON:
C         0 Calcula la derivada deseada
C
C         1 Devuelve los parametros de la matriz que define
C           la derivada en el vector DKF(1<-->M*M), en
C           DKF(M*M+1) devuelve el inverso del valor por el
C           que hay que multiplicar para obtener la derivada
C
C         2 Es com si afagisim un punt mes de la funcio a la
C           dreta del ultim punt del grid, amb valor 0.
C
C         3 El mateix que per ICON=2 i a fa la mateixa hipotesi
C           a l'esquerra del primer punt del grid.
C
C         4 Fa una reflexio parell de la funcio per els primers
C           punts, per tal de tenir els punts de la derivada
C           centrats, per aquests primers punts
C
C         5 El mateix que per ICON=4 pero amb reflexio sanar
C
C         6 El mateix que ICON=4 i ICON=2
C
C         7 El mateix que ICON=5 i ICON=2
C
C         8 Condicions periodicas: ( I<1 <--> I+N, I>N <--> I-N )
C
C        18 Condicions antiperiodicas: ( I<1 <--> -(I+N), I>N <--> -(I-N) )
C
C        12 Soposem la funcio nulla a la dreta del ultim punt
C           del grid per tal que sempre fem servir els coeficients
C           centrats
C
C        13 El mateix que per ICON=12 i fa la mateixa hipotesi
C           a l'esquerra del primer punt del grid.
C
C        16 El mateix que ICON=4 i ICON=12
C
C        17 El mateix que ICON=5 i ICON=12
C
C        20 Condicions periodiques de Wigner_Seitz ( N+i <--> N-i )
C
C        24 El mateix que ICON= 4 i ICON=20
C
C        25 El mateix que ICON= 5 i ICON=20
C
C
C     Restricciones:
C         M Debe ser mayor o igual que 3 y mas pequenyo o igual que 25
C           tambien tiene que ser un numero impar
C         K Puede valer desde 1 hasta M-1
C         N Tiene que ser mayor o igual que M
C
C
      IF(M.LT.3.OR.M.GT.25.OR.M.GT.N.OR.MOD(M,2).EQ.0
     &  .OR.K.LT.0.OR.K.GT.M-1)THEN
        WRITE(6,1100)N,M,K
        STOP
      ENDIF
      MM1=M-1
      COND=K.LT.MM1
      JI=1
      M2=M*M
      MM1O2=MM1/2
      M21O2=(M2+1)/2
      FMOFK=1.0D0
      DO I=3,M-2,2
        JI=JI+(I*I+1)*(I-2)/2+I
      ENDDO
      JI=JI+M21O2*(K-1)
      IF(COND)THEN
        JF=JI+M21O2-1
      ELSE
        JF=JI+MM1
      ENDIF
      N0=MM1O2+1
      NF=N-MM1O2
      SIG=1.0D0-2.0D0*MOD(K,2)
      IF(ICON.EQ.1)THEN
C
C       En aquest cas tornem els coeficients de les derivades
C
        IF(COND)THEN
          L=JI-1
          J1=M21O2-1
          DO J=1,J1
            L=L+1
            DKF(J)=A(L)
            JC=M2-J+1
            DKF(JC)=SIG*A(L)
          ENDDO
          L=L+1
          DKF(M21O2)=A(L)
          DKF(M2+1)=FMOFK
          RETURN
        ELSE
          DO I=1,M
            J0=1+(I-1)*M
            J1=J0+M
            L=JI-1
            DO J=J0,J1
              L=L+1
              DKF(J)=A(L)
            ENDDO
            DKF(M2+1)=FMOFK
          ENDDO
          RETURN
        ENDIF
      ENDIF
      HK=H**K
      FACT=1.0D0/(HK*FMOFK)
      IF(ICON.EQ.0)THEN
C
C          Calcula la derivada desitjada
C
        IF(COND)THEN
          DO I=1,MM1O2
            L0=JI+(I-1)*M
            L1=L0+MM1
            SUMA=0.0D0
            L=L0-1
            DO J=1,M
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO II=MM1O2,1,-1
            L1=JI+(II-1)*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO II=1,MM1O2
            I=N0-II
            DKF(I)=DKF(N0)
            I=NF+II
            DKF(I)=DKF(NF)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.2)THEN
C
C          Es com si afagisim un punt mes de la funcio a la
C          dreta del ultim punt del grid, amb valor 0.
C
        IF(COND)THEN
          DO I=1,MM1O2
            L0=JI+(I-1)*M
            L1=L0+MM1
            SUMA=0.0D0
            L=L0-1
            DO J=1,M
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I-1
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          DO II=MM1O2-1,1,-1
            L1=JI+II*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2+1,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO II=1,MM1O2
            I=N0-II
            DKF(I)=DKF(N0)
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          DO I=NF+2,N
            DKF(I)=DKF(NF+1)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.3)THEN
C
C          El mateix que per ICON=2 i a fa la mateixa hipotesi
C          a l'esquerra del primer punt del grid.
C
        IF(COND)THEN
          DO I=1,MM1O2-1
            L0=JI+I*M
            L1=L0+MM1
            SUMA=0.0D0
            L=L0
            DO J=1,MM1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          L0=JI+MM1O2*M
          L1=L0+MM1O2
            I=MM1O2
            SUMA=0.0D0
            L=L0
            DO J=I-MM1O2+1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I-1
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          DO II=MM1O2-1,1,-1
            L1=JI+II*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2+1,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
            I=1
            SUMA=0.0D0
            L=JI
            DO J=1,MM1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          DO I=2,MM1O2
            DKF(I)=DKF(1)
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          DO I=NF+2,N
            DKF(I)=DKF(NF+1)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.4)THEN
C
C           Fa una reflexio parell de la funcio per els primers
C           punts, per tal de tenir els punts de la derivada
C           centrats, per aquests primers punts
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO II=MM1O2,1,-1
            L1=JI+(II-1)*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO II=1,MM1O2
            I=NF+II
            DKF(I)=DKF(NF)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.5)THEN
C
C          Fa una reflexio sanar de la funcio per els primers
C          punts, per tal de tenir els punts de la derivada
C          sempre centrats.
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO II=MM1O2,1,-1
            L1=JI+(II-1)*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO II=1,MM1O2
            I=NF+II
            DKF(I)=DKF(NF)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.6)THEN
C
C          El mateix que per ICON=4 i ICON=2
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I-1
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          DO II=MM1O2-1,1,-1
            L1=JI+II*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2+1,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          DO I=NF+2,N
            DKF(I)=DKF(NF+1)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.7)THEN
C
C          El mateix que per ICON=5 i ICON=2
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I-1
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          DO II=MM1O2-1,1,-1
            L1=JI+II*M
            L0=L1+MM1
            SUMA=0.0D0
            L=L0+1
            I=NF+MM1O2-II+1
            DO J=NF-MM1O2+1,N
              L=L-1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*SIG*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
            I=NF+1
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          DO I=NF+2,N
            DKF(I)=DKF(NF+1)
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.8)THEN
C
C         8 Condicions periodicas: ( I<1 <--> I+N, I>N <--> I-N )
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=I-MM1O2,0
              L=L+1
              J=JJ+N
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DO JJ=N+1,MM1O2+I
              L=L-1
              J=JJ-N
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            SUMA=0.0D0
            L=JI-1
            DO JJ=I-MM1O2,0
              L=L+1
              J=JJ+N
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO JJ=N+1,I+MM1O2
              L=L+1
              J=JJ-N
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.18)THEN
C
C        18 Condicions antiperiodicas: ( I<1 <--> -(I+N), I>N <--> -(I-N) )
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=I-MM1O2,0
              L=L+1
              J=JJ+N
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DO JJ=N+1,MM1O2+I
              L=L-1
              J=JJ-N
              SUMB=SUMB-A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            SUMA=0.0D0
            L=JI-1
            DO JJ=I-MM1O2,0
              L=L+1
              J=JJ+N
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO JJ=N+1,I+MM1O2
              L=L+1
              J=JJ-N
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.12)THEN
C
C           Soposem la funcio nulla a la dreta del ultim punt
C           del grid per tal que sempre fem servir els coeficients
C           centrats
C
        IF(COND)THEN
          DO I=1,MM1O2
            L0=JI+(I-1)*M
            L1=L0+MM1
            SUMA=0.0D0
            L=L0-1
            DO J=1,M
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO     
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO II=1,MM1O2
            I=N0-II
            DKF(I)=DKF(N0)
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.13)THEN
C
C         El mateix que per ICON=12 i fa la mateixa hipotesi
C         a l'esquerra del primer punt del grid.
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-I+MM1O2
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            SUMA=0.0D0
            L=JI-I+MM1O2
            DO J=1,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.16)THEN
C
C          El mateix que ICON= 4 i ICON=12
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.17)THEN
C
C          El mateix que ICON= 5 i ICON=12
C
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO     
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.20)THEN
C
C        Condicions periodiques de Wigner_Seitz ( N+i <--> N-i )
C
        IF(COND)THEN
          DO I=1,MM1O2
            L0=JI+(I-1)*M
            L1=L0+MM1
            SUMA=0.0D0
            L=L0-1
            DO J=1,M
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DO JJ=N+1,MM1O2+I
              L=L-1
              J=2*N-JJ
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO II=1,MM1O2
            I=N0-II
            DKF(I)=DKF(N0)
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO JJ=N+1,I+MM1O2
              L=L+1
              J=2*N-JJ
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.24)THEN
C
C        El mateix que ICON=4 i ICON=20
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DO JJ=N+1,MM1O2+I
              L=L-1
              J=2*N-JJ
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO JJ=N+1,I+MM1O2
              L=L+1
              J=2*N-JJ
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE IF(ICON.EQ.25)THEN
C
C        El mateix que ICON=5 i ICON=20
C
        IF(COND)THEN
          L0=JI+MM1O2*M
          L1=L0+MM1O2
          DO I=1,MM1O2
            SUMA=0.0D0
            L=L0-1
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            SUMB=0.0D0
            L=L1
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,MM1O2+I
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=L0-1
            DO J=I-MM1O2,I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            L=L1
            SUMB=0.0D0
            DO J=I+1,N
              L=L-1
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DO JJ=N+1,MM1O2+I
              L=L-1
              J=2*N-JJ
              SUMB=SUMB+A(L)*F(J)
            ENDDO
            DKF(I)=(SUMA+SUMB*SIG)*FACT
          ENDDO
        ELSE
          DO I=1,MM1O2
            L=JI-1
            SUMA=0.0D0
            DO JJ=MM1O2+1-I,1,-1
              J=JJ+1
              L=L+1
              SUMA=SUMA-A(L)*F(J)
            ENDDO
            DO J=1,MM1O2+I
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=N0,NF
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,I+MM1O2
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
          DO I=NF+1,N
            SUMA=0.0D0
            L=JI-1
            DO J=I-MM1O2,N
              L=L+1
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DO JJ=N+1,I+MM1O2
              L=L+1
              J=2*N-JJ
              SUMA=SUMA+A(L)*F(J)
            ENDDO
            DKF(I)=SUMA*FACT
          ENDDO
        ENDIF
      ELSE
        WRITE(*,*) ' Subroutine newder mal utilitzada: ICON-->',ICON
        RETURN
      ENDIF
      RETURN
1000  FORMAT(11F11.0)
1100  FORMAT(' Subroutine Newder mal utilizada: N,M,K------->',3I5)
      END
