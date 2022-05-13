      SUBROUTINE CUADJTQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   KK,
     *     KSTART,   KSTOP,    KLEN,
     *     PAPH,     PT,       PQ,       LDFLAG,   KCALL
     *  )
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
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
C          NONE
C
      INCLUDE 'comcon.h'
C
      REAL     PT(KHOR,KLEV),          PQ(KHOR,KLEV),
     *         PAPH(KHOR,KLEVP1)
      LOGICAL  LDFLAG(KHOR)
C
CDIR$ VFUNCTION EXP
C
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
      REAL
     *        ZCOND  (JPHR),
     *        ZQP    (JPHR)
C
C
C----------------------------------------------------------------------
C
C     1.           DEFINE CONSTANTS
C                  ----------------
C
  100 CONTINUE
      Z5ALVCP=C5LES*ALV/CPD
      Z5ALSCP=C5IES*ALS/CPD
      ZALVDCP=ALV/CPD
      ZALSDCP=ALS/CPD
C
C
C----------------------------------------------------------------------
C
C     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
C                  -----------------------------------------------------
C
  200 CONTINUE
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 210 JL=KSTART,KSTOP
      IF(LDFLAG(JL)) THEN
         IF(PT(JL,KK)-TMELT.GT.0.) THEN
            ZCVM3=C3LES
            ZCVM4=C4LES
            ZCVM5=Z5ALVCP
            ZLDCP=ZALVDCP
         ELSE
            ZCVM3=C3IES
            ZCVM4=C4IES
            ZCVM5=Z5ALSCP
            ZLDCP=ZALSDCP
         END IF
         ZQP(JL)=1./PAPH(JL,KK)
         ZQSAT=C2ES*EXP(ZCVM3*(PT(JL,KK)-TMELT)*
     1              (1./(PT(JL,KK)-ZCVM4)))*ZQP(JL)
         ZQSAT=MIN(0.5,ZQSAT)
         ZCOR=1./(1.-VTMPC1*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZCVM5*ZQSAT*ZCOR*
     1              (1./(PT(JL,KK)-ZCVM4))**2)
         IF(KCALL.EQ.1) ZCOND(JL)=MAX(ZCOND(JL),0.)
         IF(KCALL.EQ.2) ZCOND(JL)=MIN(ZCOND(JL),0.)
         PT(JL,KK)=PT(JL,KK)+ZLDCP*ZCOND(JL)
         PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
      END IF
  210 CONTINUE
C
      DO 215 JL=KSTART,KSTOP
      ISUM=ISUM+ICVMGT(1,0,ZCOND(JL).NE.0.)
  215 CONTINUE
      IF(ISUM.EQ.0) GO TO 230
C
      DO 220 JL=KSTART,KSTOP
      IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
         IF(PT(JL,KK)-TMELT.GT.0.) THEN
            ZCVM3=C3LES
            ZCVM4=C4LES
            ZCVM5=Z5ALVCP
            ZLDCP=ZALVDCP
         ELSE
            ZCVM3=C3IES
            ZCVM4=C4IES
            ZCVM5=Z5ALSCP
            ZLDCP=ZALSDCP
         END IF
         ZQSAT=C2ES*EXP(ZCVM3*(PT(JL,KK)-TMELT)*
     1              (1./(PT(JL,KK)-ZCVM4)))*ZQP(JL)
         ZQSAT=MIN(0.5,ZQSAT)
         ZCOR=1./(1.-VTMPC1*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZCVM5*ZQSAT*ZCOR*
     1              (1./(PT(JL,KK)-ZCVM4))**2)
         PT(JL,KK)=PT(JL,KK)+ZLDCP*ZCOND1
         PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      END IF
  220 CONTINUE
C
  230 CONTINUE
C
      RETURN
      END
