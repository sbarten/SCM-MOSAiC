      SUBROUTINE CUDUDV
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     KTOPM2,   KTYPE,    KCBOT,    PAPH,     LDCUM,
     *     PUEN,     PVEN,     PVOM,     PVOL,
     *     PUU,      PUD,      PVU,      PVD,
     *     PMFU,     PMFD,     PSDISS
     *   )
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
      INCLUDE 'comcon.h'
C
C
      REAL     PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PVOL(KHOR2,KLEV),       PVOM(KHOR2,KLEV),
     *         PAPH(KHOR,KLEVP1)
      REAL     PUU(KHOR,KLEV),         PUD(KHOR,KLEV),
     *         PVU(KHOR,KLEV),         PVD(KHOR,KLEV),
     *         PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV)
      INTEGER  KTYPE(KHOR),            KCBOT(KHOR)
      LOGICAL  LDCUM(KHOR)
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      REAL
     *       ZDISS (JPHR),
     *       ZMFDU (JPHR,MLEV),
     *       ZMFDV (JPHR,MLEV),
     *       ZMFUU (JPHR,MLEV),
     *       ZMFUV (JPHR,MLEV)
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
      DO 110 JL=KSTART,KSTOP
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
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 130 JL=KSTART,KSTOP
      IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
         IKB=KCBOT(JL)
         ZZP=((PAPH(JL,KLEVP1)-PAPH(JL,JK))/
     1        (PAPH(JL,KLEVP1)-PAPH(JL,IKB)))
         ZZP=CVMGT(ZZP**2,ZZP,KTYPE(JL).EQ.3)
         ZMFUU(JL,JK)=ZMFUU(JL,IKB)*ZZP
         ZMFUV(JL,JK)=ZMFUV(JL,IKB)*ZZP
      END IF
  130 CONTINUE
  140 CONTINUE
C
      DO 150 JL=KSTART,KSTOP
      ZDISS(JL)=0.
  150 CONTINUE
C
      DO 190 JK=KTOPM2,KLEV
C
      IF(JK.LT.KLEV) THEN
         DO 160 JL=KSTART,KSTOP
            IF(LDCUM(JL)) THEN
               ZDUDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1               (ZMFUU(JL,JK+1)-ZMFUU(JL,JK)+
     1                ZMFDU(JL,JK+1)-ZMFDU(JL,JK))
               ZDVDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
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
         DO 170 JL=KSTART,KSTOP
            IF(LDCUM(JL)) THEN
               ZDUDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1                   (ZMFUU(JL,JK)+ZMFDU(JL,JK))
               ZDVDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
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
      ZSUM=SSUM(KLEN,ZDISS(KSTART),1)
      PSDISS=PSDISS+ZSUM
C
      RETURN
      END
