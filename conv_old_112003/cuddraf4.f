      SUBROUTINE CUDDRAF4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVP1,
     *     PTENH,    PQENH,    PUEN,     PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   PXTENH,   PXTD,     PMFDXT,

C LG- end

     *     PGEOH,    PAPHP1,   PRFL,     PLU,
     *     KTYPE,    KCBOT,
     *     PTD,      PQD,      PUD,      PVD,
     *     PMFD,     PMFDS,    PMFDQ,    PDMFDP,
     *     KDTOP,    LDDRAF)
C
C          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
C
C          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
C          (I.E. T,Q,U AND V AND FLUXES)
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
C          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
C          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
C
C          METHOD.
C          --------
C          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
C          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
C          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
C
C          EXTERNALS
C          ---------
C          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
C          SATURATED DESCENT
C
C          REFERENCE
C          ---------
C          (TIEDTKE,1989)
C
CJHC*CALL COMCON
CJHC*CALL COMCUMF
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KLP2,KLEV),       PQENH(KLP2,KLEV),
     *         PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PGEOH(KLP2,KLEV),       PAPHP1(KLP2,KLEVP1)
C
      REAL     PTD(KLP2,KLEV),         PQD(KLP2,KLEV),
     *         PUD(KLP2,KLEV),         PVD(KLP2,KLEV),
     *         PMFD(KLP2,KLEV),        PMFDS(KLP2,KLEV),
     *         PMFDQ(KLP2,KLEV),       PDMFDP(KLP2,KLEV),
     *         PRFL(KLP2)
      INTEGER  KTYPE(KLP2),            KCBOT(KLP2),
     *         KDTOP(KLP2)
      LOGICAL  LDDRAF(KLP2)
C
      INCLUDE 'paramh.h'
      REAL     ZDMFEN(JPHR),           ZDMFDE(JPHR),
     *         ZCOND(JPHR)
      REAL     ZPH(JPHR)
      LOGICAL  LLO2(JPHR)

C LG- adding the tracers

      REAL   PXTENH(KLON,KLEV,KTRAC),  PXTD(KLON,KLEV,KTRAC),
     *       PMFDXT(KLON,KLEV,KTRAC)

C LG- end

      LOGICAL  LLO1
C
C
C----------------------------------------------------------------------
C
C     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
C                     (A) CALCULATING ENTRAINMENT RATES, ASSUMING
C                         LINEAR DECREASE OF MASSFLUX IN PBL
C                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
C                         AND MOISTENING IS CALCULATED IN *CUADJTQ*
C                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
C                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
C                    -------------------------------------------------
C
  100 CONTINUE
      DO 180 JK=3,KLEV
      IS=0
      DO 110 JL=KIDIA,KFDIA
      ZPH(JL)=PAPHP1(JL,JK)
      LLO2(JL)=LDDRAF(JL).AND.PMFD(JL,JK-1).LT.0.
      IS=IS+ICVMGT(1,0,LLO2(JL))
  110 CONTINUE
      IF(IS.EQ.0) GO TO 180
      DO 122 JL=KIDIA,KFDIA
      IF(LLO2(JL)) THEN
         ZENTR=ENTRDD*PMFD(JL,JK-1)*RD*PTENH(JL,JK-1)/
     1         (G*PAPHP1(JL,JK-1))*(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))
         ZDMFEN(JL)=ZENTR
         ZDMFDE(JL)=ZENTR
      END IF
  122 CONTINUE
      ITOPDE=KLEV-2
         IF(JK.GT.ITOPDE) THEN
            DO 124 JL=KIDIA,KFDIA
            IF(LLO2(JL)) THEN
               ZDMFEN(JL)=0.
               ZDMFDE(JL)=PMFD(JL,ITOPDE)*
     1         (PAPHP1(JL,JK)-PAPHP1(JL,JK-1))/
     1         (PAPHP1(JL,KLEVP1)-PAPHP1(JL,ITOPDE))
            END IF
  124       CONTINUE
         END IF
C
      DO 126 JL=KIDIA,KFDIA
         IF(LLO2(JL)) THEN
            PMFD(JL,JK)=PMFD(JL,JK-1)+ZDMFEN(JL)-ZDMFDE(JL)
            ZSEEN=(CPD*PTENH(JL,JK-1)+PGEOH(JL,JK-1))*ZDMFEN(JL)
            ZQEEN=PQENH(JL,JK-1)*ZDMFEN(JL)
            ZSDDE=(CPD*PTD(JL,JK-1)+PGEOH(JL,JK-1))*ZDMFDE(JL)
            ZQDDE=PQD(JL,JK-1)*ZDMFDE(JL)
            ZMFDSK=PMFDS(JL,JK-1)+ZSEEN-ZSDDE
            ZMFDQK=PMFDQ(JL,JK-1)+ZQEEN-ZQDDE
            PQD(JL,JK)=ZMFDQK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
            PTD(JL,JK)=(ZMFDSK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))-
     1                  PGEOH(JL,JK))*RCPD
            PTD(JL,JK)=MIN(400.,PTD(JL,JK))
            PTD(JL,JK)=MAX(100.,PTD(JL,JK))
            ZCOND(JL)=PQD(JL,JK)
         END IF
  126 CONTINUE
C

C LG- adding the tracers

      DO 1264 JT=1,KTRAC
      DO 1262 JL=KIDIA,KFDIA
      IF(LLO2(JL)) THEN
       ZXTEEN=PXTENH(JL,JK-1,JT)*ZDMFEN(JL)
       ZXTDDE=PXTD(JL,JK-1,JT)*ZDMFDE(JL)
       ZMFDXTK=PMFDXT(JL,JK-1,JT)+ZXTEEN-ZXTDDE
       PXTD(JL,JK,JT)=ZMFDXTK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
      ENDIF
 1262 CONTINUE
 1264 CONTINUE
 
C LG- end

C
C
      IK=JK
      ICALL=2
      CALL CUADJTQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     IK,
     *     ZPH,      PTD,      PQD,      LLO2,     ICALL)
C
      DO 150 JL=KIDIA,KFDIA
         IF(LLO2(JL)) THEN
            ZCOND(JL)=ZCOND(JL)-PQD(JL,JK)
            ZBUO=PTD(JL,JK)*(1.+VTMPC1*PQD(JL,JK))-
     1      PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))
            LLO1=ZBUO.LT.0..AND.(PRFL(JL)-PMFD(JL,JK)*ZCOND(JL).GT.0.)
            PMFD(JL,JK)=CVMGT(PMFD(JL,JK),0.,LLO1)
            PMFDS(JL,JK)=(CPD*PTD(JL,JK)+PGEOH(JL,JK))*PMFD(JL,JK)
            PMFDQ(JL,JK)=PQD(JL,JK)*PMFD(JL,JK)
            ZDMFDP=-PMFD(JL,JK)*ZCOND(JL)
            PDMFDP(JL,JK-1)=ZDMFDP
            PRFL(JL)=PRFL(JL)+ZDMFDP
         END IF
  150 CONTINUE
C

C LG- adding the tracers

      DO 1504 JT=1,KTRAC
      DO 1502 JL=KIDIA,KFDIA
      IF(LLO2(JL)) THEN
       PMFDXT(JL,JK,JT)=PXTD(JL,JK,JT)*PMFD(JL,JK)
      ENDIF
 1502 CONTINUE
 1504 CONTINUE
 
C LG- end

        IF(LMFDUDV) THEN
          DO 160 JL=KIDIA,KFDIA
             IF(LLO2(JL).AND.PMFD(JL,JK).LT.0.) THEN
                ZMFDUK=PMFD(JL,JK-1)*PUD(JL,JK-1)+
     1          ZDMFEN(JL)*PUEN(JL,JK-1)-ZDMFDE(JL)*PUD(JL,JK-1)
                ZMFDVK=PMFD(JL,JK-1)*PVD(JL,JK-1)+
     1          ZDMFEN(JL)*PVEN(JL,JK-1)-ZDMFDE(JL)*PVD(JL,JK-1)
                PUD(JL,JK)=ZMFDUK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
                PVD(JL,JK)=ZMFDVK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
             END IF
  160     CONTINUE
        END IF
C
  180 CONTINUE
C
      RETURN
      END
