      SUBROUTINE CUDUDV4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVP1,
     *     KTOPM2,   KTYPE,    KCBOT,    PAPHP1,   LDCUM,
     *     PUEN,     PVEN,     PVOM,     PVOL,
     *     PUU,      PUD,      PVU,      PVD,
     *     PMFU,     PMFD,     PSDISS)
C
C
C**** *CUDUDV* - UPDATES U AND V TENDENCIES,
C                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C**   INTERFACE.
C     ----------
C
C          *CUDUDV* IS CALLED FROM *CUMASTR*
C
CJHC*CALL COMCON
      INCLUDE 'comcon.h'
C
C
      REAL     PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PVOL(K2LP2,KLEV),       PVOM(K2LP2,KLEV),
     *         PAPHP1(KLP2,KLEVP1)
      REAL     PUU(KLP2,KLEV),         PUD(KLP2,KLEV),
     *         PVU(KLP2,KLEV),         PVD(KLP2,KLEV),
     *         PMFU(KLP2,KLEV),        PMFD(KLP2,KLEV)
      INTEGER  KTYPE(KLP2),            KCBOT(KLP2)
      LOGICAL  LDCUM(KLP2)
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      REAL     ZMFUU(JPHR,MLEV),       ZMFDU(JPHR,MLEV),
     *         ZMFUV(JPHR,MLEV),       ZMFDV(JPHR,MLEV),
     *         ZDISS(JPHR)
C
C
C----------------------------------------------------------------------
C
C*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
C                  ----------------------------------------------
C
  100 CONTINUE
      DO 120 JK=KTOPM2,KLEV
      IK=JK-1
      DO 110 JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
         ZMFUU(JL,JK)=PMFU(JL,JK)*(PUU(JL,JK)-PUEN(JL,IK))
         ZMFUV(JL,JK)=PMFU(JL,JK)*(PVU(JL,JK)-PVEN(JL,IK))
         ZMFDU(JL,JK)=PMFD(JL,JK)*(PUD(JL,JK)-PUEN(JL,IK))
         ZMFDV(JL,JK)=PMFD(JL,JK)*(PVD(JL,JK)-PVEN(JL,IK))
      END IF
  110 CONTINUE
  120 CONTINUE
      DO 140 JK=KTOPM2,KLEV
CDIR$ IVDEP
      DO 130 JL=KIDIA,KFDIA
      IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
         IKB=KCBOT(JL)
         ZZP=((PAPHP1(JL,KLEVP1)-PAPHP1(JL,JK))/
     1        (PAPHP1(JL,KLEVP1)-PAPHP1(JL,IKB)))
         ZZP=CVMGT(ZZP**2,ZZP,KTYPE(JL).EQ.3)
         ZMFUU(JL,JK)=ZMFUU(JL,IKB)*ZZP
         ZMFUV(JL,JK)=ZMFUV(JL,IKB)*ZZP
         ZMFDU(JL,JK)=ZMFDU(JL,IKB)*ZZP
         ZMFDV(JL,JK)=ZMFDV(JL,IKB)*ZZP
      END IF
  130 CONTINUE
  140 CONTINUE
C
      DO 150 JL=KIDIA,KFDIA
      ZDISS(JL)=0.
  150 CONTINUE
C
      DO 190 JK=KTOPM2,KLEV
C
      IF(JK.LT.KLEV) THEN
         DO 160 JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               ZDUDT=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     1               (ZMFUU(JL,JK+1)-ZMFUU(JL,JK)+
     1                ZMFDU(JL,JK+1)-ZMFDU(JL,JK))
               ZDVDT=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     1               (ZMFUV(JL,JK+1)-ZMFUV(JL,JK)+
     1                ZMFDV(JL,JK+1)-ZMFDV(JL,JK))
               ZDISS(JL)=ZDISS(JL)+
     1                   PUEN(JL,JK)*(ZMFUU(JL,JK+1)-ZMFUU(JL,JK)+
     1                                ZMFDU(JL,JK+1)-ZMFDU(JL,JK))+
     1                   PVEN(JL,JK)*(ZMFUV(JL,JK+1)-ZMFUV(JL,JK)+
     1                                ZMFDV(JL,JK+1)-ZMFDV(JL,JK))
               PVOM(JL,JK)=PVOM(JL,JK)+ZDUDT
               PVOL(JL,JK)=PVOL(JL,JK)+ZDVDT
            END IF
  160    CONTINUE
C
      ELSE
         DO 170 JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               ZDUDT=-(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     1                   (ZMFUU(JL,JK)+ZMFDU(JL,JK))
               ZDVDT=-(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     1                   (ZMFUV(JL,JK)+ZMFDV(JL,JK))
               ZDISS(JL)=ZDISS(JL)-
     1                   (PUEN(JL,JK)*(ZMFUU(JL,JK)+ZMFDU(JL,JK))+
     1                    PVEN(JL,JK)*(ZMFUV(JL,JK)+ZMFDV(JL,JK)))
               PVOM(JL,JK)=PVOM(JL,JK)+ZDUDT
               PVOL(JL,JK)=PVOL(JL,JK)+ZDVDT
            END IF
  170    CONTINUE
       END IF
C
  190 CONTINUE
C
      ZSUM=SSUM(KLON,ZDISS(1),1)
      PSDISS=PSDISS+ZSUM
C
      RETURN
      END
