      SUBROUTINE CUDTDQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     CONACC,
     *     KTOPM2,   PAPH,     PGEO,     PTS,      LDLAND,
     *     LDCUM,    PTEN,     PTTE,     PQTE,
     *     PMFUS,    PMFDS,    PMFUQ,    PMFDQ,
     *     PMFUL,    PDMFUP,   PDMFDP,   PLUDE,
     *     PRAIN,    PRFL,     PSRAIN,   PSEVAP,   PSHEAT,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,
C
     *     TWODT
     *   )
C
C
C**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
C                DOES GLOBAL DIAGNOSTICS
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C**   INTERFACE.
C     ----------
C
C          *CUDTDQ* IS CALLED FROM *CUMASTR*
C
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      LOGICAL  LLO1,LLO2,LLO3
C
      REAL     PTTE(KHOR,KLEV),        PQTE(KHOR,KLEV),
     *         PTEN(KHOR,KLEV),        PTS(KHOR),
     *         PGEO(KHOR,KLEV),        PAPH(KHOR,KLEVP1),
     *         PAPRC(KHOR),            PAPRS(KHOR),
     *         PRSFC(KHOR),            PSSFC(KHOR)
      REAL     PMFUS(KHOR,KLEV),       PMFDS(KHOR,KLEV),
     *         PMFUQ(KHOR,KLEV),       PMFDQ(KHOR,KLEV),
     *         PMFUL(KHOR,KLEV),       PLUDE(KHOR,KLEV),
     *         PDMFUP(KHOR,KLEV),      PDMFDP(KHOR,KLEV),
     *         PRFL(KHOR),             PRAIN(KHOR)
      LOGICAL  LDLAND(KHOR),           LDCUM(KHOR)
C
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING SPACE
C
      INCLUDE 'paramh.h'
C
      INTEGER
     *        ITMELT(JPHR)
      REAL
     *        ZSHEAT(JPHR)
C
C
C----------------------------------------------------------------------
C
C*    1.0          SPECIFY PARAMETERS
C                  ------------------
C
  100 CONTINUE
      ZDIAGT=CONACC*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
C
C
C----------------------------------------------------------------------
C
C*    2.0          INCREMENTATION OF T AND Q TENDENCIES
C                  ------------------------------------
C
  200 CONTINUE
      DO 210 JL=KSTART,KSTOP
      ITMELT(JL)=KLEV
      ZSHEAT(JL)=0.
  210 CONTINUE
C
      DO 250 JK=KTOPM2,KLEV
C
      IF(JK.LT.KLEV) THEN
         DO 220 JL=KSTART,KSTOP
         IF(LDCUM(JL)) THEN
            LLO1=(PTEN(JL,JK)-TMELT).GT.0.
            ZALV=CVMGT(ALV,ALS,LLO1)
            ZDTDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     1               (PMFUS(JL,JK+1)-PMFUS(JL,JK)+
     1                PMFDS(JL,JK+1)-PMFDS(JL,JK)
     1         -ZALV*(PMFUL(JL,JK+1)-PMFUL(JL,JK)-
     1               (PDMFUP(JL,JK)+PDMFDP(JL,JK))))
            PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
            ZDQDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1               (PMFUQ(JL,JK+1)-PMFUQ(JL,JK)+
     1                PMFDQ(JL,JK+1)-PMFDQ(JL,JK)+
     1                PMFUL(JL,JK+1)-PMFUL(JL,JK)-
     1           (PDMFUP(JL,JK)+PDMFDP(JL,JK)))
            PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
            LLO2=PGEO(JL,JK).LT.3000.AND.ITMELT(JL).EQ.KLEV
            ITMELT(JL)=ICVMGT(JK,ITMELT(JL),LLO2)
            ZSHEAT(JL)=ZSHEAT(JL)+ZALV*(PDMFUP(JL,JK)+PDMFDP(JL,JK))
         END IF
  220 CONTINUE
C
      ELSE
         DO 230 JL=KSTART,KSTOP
         IF(LDCUM(JL)) THEN
            LLO1=(PTEN(JL,JK)-TMELT).GT.0.
            ZALV=CVMGT(ALV,ALS,LLO1)
            ZDTDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     1                (PMFUS(JL,JK)+PMFDS(JL,JK)-ZALV*
     1                (PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))
            PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
            ZDQDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1                (PMFUQ(JL,JK)+PMFDQ(JL,JK)+
     1                (PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))
            PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
            LLO2=PGEO(JL,JK).LT.3000.AND.ITMELT(JL).EQ.KLEV
            ITMELT(JL)=ICVMGT(JK,ITMELT(JL),LLO2)
            ZSHEAT(JL)=ZSHEAT(JL)+ZALV*(PDMFUP(JL,JK)+PDMFDP(JL,JK))
         END IF
  230    CONTINUE
      END IF
C
  250 CONTINUE
C
C
C---------------------------------------------------------------------
C
C      3.          UPDATE SURFACE FIELDS AND DO GLOBAL BUDGETS
C                  -------------------------------------------
C
  300 CONTINUE
      DO 310 JL=KSTART,KSTOP
      LLO1=PTS(JL).GE.TMELT
      LLO2=PTEN(JL,ITMELT(JL)).GT.(TMELT-3.)
      LLO3=LCVMGT(LLO1.OR.LLO2,LLO2,LDLAND(JL))
      PRSFC(JL)=CVMGT(PRFL(JL),0.,LLO3)
      PSSFC(JL)=CVMGT(0.,PRFL(JL),LLO3)
      PAPRC(JL)=PAPRC(JL)+ZDIAGW*PRFL(JL)
      PAPRS(JL)=PAPRS(JL)+CVMGT(0.,ZDIAGW*PRFL(JL),LLO3)
      PSHEAT=PSHEAT+ZSHEAT(JL)
      PSRAIN=PSRAIN+PRAIN(JL)
      PSEVAP=PSEVAP-PRFL(JL)
  310 CONTINUE
C
      PSEVAP=PSEVAP+PSRAIN
C
C
      RETURN
      END
