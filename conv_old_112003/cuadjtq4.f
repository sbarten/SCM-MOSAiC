      SUBROUTINE CUADJTQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KK,
     *     PP,       PT,       PQ,       LDFLAG,   KCALL)
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C          D.SALMOND         CRAY(UK))      12/8/91
C
C          PURPOSE.
C          --------
C          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM SUBROUTINES:
C              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
C              *CUASC*    (T AND Q AT CLOUD LEVELS)
C              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
C          INPUT ARE UNADJUSTED T AND Q VALUES,
C          IT RETURNS ADJUSTED VALUES OF T AND Q
C          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
C               KCALL=0    ENV. T AND QS IN*CUINI*
C               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
C               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
C
C          EXTERNALS
C          ---------
C          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
C          FOR CONDENSATION CALCULATIONS.
C          THE TABLES ARE INITIALISED IN *SETPHYS*.
C
CJHC*CALL COMCON
CJHC*CALL YOTLUC
      INCLUDE 'comcon.h'
      INCLUDE 'yotluc.h'
C
      REAL     PT(KLP2,KLEV),          PQ(KLP2,KLEV),
     *         PP(KLP2)
      LOGICAL  LDFLAG(KLP2)
C
CDIR$ VFUNCTION EXP
C
      INCLUDE 'paramh.h'
      REAL     ZCOND(JPHR),            ZQP(JPHR)
C
C----------------------------------------------------------------------
C
C     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
C                  -----------------------------------------------------
C
  200 CONTINUE
      IF (KCALL.EQ.1 ) THEN
C
         ISUM=0
CDIR$    IVDEP
         DO 210 JL=KIDIA,KFDIA
         IF(LDFLAG(JL)) THEN
            ZQP(JL)=1./PP(JL)
            IT=PT(JL,KK)*1000.
C
            ZQSAT=TLUCUA(IT)*ZQP(JL)
C
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            ZCOND(JL)=MAX(ZCOND(JL),0.)
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
            PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
            IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
         ELSE
            ZCOND(JL)=0.0
         END IF
  210    CONTINUE
C
         IF(ISUM.EQ.0) GO TO 230
C
CDIR$    IVDEP
         DO 220 JL=KIDIA,KFDIA
         IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
            IT=PT(JL,KK)*1000.
            ZQSAT=TLUCUA(IT)*ZQP(JL)
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
            PQ(JL,KK)=PQ(JL,KK)-ZCOND1
         END IF
  220    CONTINUE
C
  230    CONTINUE
C
      END IF

      IF(KCALL.EQ.2) THEN
C
         ISUM=0
CDIR$    IVDEP
         DO 310 JL=KIDIA,KFDIA
         IF(LDFLAG(JL)) THEN
            IT=PT(JL,KK)*1000.
            ZQP(JL)=1./PP(JL)
            ZQSAT=TLUCUA(IT)*ZQP(JL)
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            ZCOND(JL)=MIN(ZCOND(JL),0.)
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
            PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
            IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
         ELSE
            ZCOND(JL)=0.0
         END IF
  310    CONTINUE
C
         IF(ISUM.EQ.0) GO TO 330
C
CDIR$    IVDEP
         DO 320 JL=KIDIA,KFDIA
         IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
            IT=PT(JL,KK)*1000.
            ZQSAT=TLUCUA(IT)*ZQP(JL)
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
            PQ(JL,KK)=PQ(JL,KK)-ZCOND1
         END IF
  320    CONTINUE
C
  330    CONTINUE
C
      END IF

      IF(KCALL.EQ.0) THEN
C
         ISUM=0
CDIR$    IVDEP
         DO 410 JL=KIDIA,KFDIA
         IT=PT(JL,KK)*1000.
            ZQP(JL)=1./PP(JL)
         ZQSAT=TLUCUA(IT)*ZQP(JL)
         ZQSAT=MIN(0.5,ZQSAT)
         ZCOR=1./(1.-VTMPC1*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
         PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
         PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
         IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
  410    CONTINUE
C
         IF(ISUM.EQ.0) GO TO 430
C
CDIR$    IVDEP
         DO 420 JL=KIDIA,KFDIA
         IT=PT(JL,KK)*1000.
         ZQSAT=TLUCUA(IT)*ZQP(JL)
         ZQSAT=MIN(0.5,ZQSAT)
         ZCOR=1./(1.-VTMPC1*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
         PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
         PQ(JL,KK)=PQ(JL,KK)-ZCOND1
  420    CONTINUE
C
  430    CONTINUE
C
      END IF

      IF(KCALL.EQ.4) THEN

CDIR$    IVDEP
         DO 510 JL=KIDIA,KFDIA
         IT=PT(JL,KK)*1000.
         ZQP(JL)=1./PP(JL)
C
         ZQSAT=TLUCUA(IT)*ZQP(JL)
C
         ZQSAT=MIN(0.5,ZQSAT)
         ZCOR=1./(1.-VTMPC1*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
         PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
         PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
  510    CONTINUE

CDIR$    IVDEP
         DO 520 JL=KIDIA,KFDIA
         IT=PT(JL,KK)*1000.
C
         ZQSAT=TLUCUA(IT)*ZQP(JL)
C
         ZQSAT=MIN(0.5,ZQSAT)
         ZCOR=1./(1.-VTMPC1*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
         PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
         PQ(JL,KK)=PQ(JL,KK)-ZCOND1
  520    CONTINUE

      END IF

C
      RETURN
      END
