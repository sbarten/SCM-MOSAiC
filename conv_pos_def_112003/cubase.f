      SUBROUTINE CUBASE
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   KLEVM1,
     *     KSTART,   KSTOP,    KLEN,
     *     PTENH,    PQENH,    PGEOH,    PAPH,
     *     PTU,      PQU,      PLU,
     *     LDCUM,    KCBOT,    KLAB
     *    )
C
C          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
C          FOR CUMULUS PARAMETERIZATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
C          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
C                 KLAB=1 FOR SUBCLOUD LEVELS
C                 KLAB=2 FOR CONDENSATION LEVEL
C
C          METHOD.
C          --------
C          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
C          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
C
C          EXTERNALS
C          ---------
C          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
C
      INCLUDE 'comcon.h'
C
      REAL     PTENH(KHOR,KLEV),       PQENH(KHOR,KLEV),
     *         PGEOH(KHOR,KLEV),       PAPH(KHOR,KLEVP1)
C
      REAL     PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PLU(KHOR,KLEV)
      INTEGER  KLAB(KHOR,KLEV),        KCBOT(KHOR)
      LOGICAL  LDCUM(KHOR)
C
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
C
      LOGICAL
     *        LOFLAG (JPHR)
      REAL
     *        ZQOLD  (JPHR)
C
C----------------------------------------------------------------------
C
C     1.           INITIALIZE VALUES AT LIFTING LEVEL
C                  ----------------------------------
C

  100 CONTINUE
      DO 110 JL=KSTART,KSTOP
      KLAB(JL,KLEV)=1
      KCBOT(JL)=KLEVM1
      LDCUM(JL)=.FALSE.
  110 CONTINUE
C
C
C----------------------------------------------------------------------
C
C     2.0          DO ASCENT IN SUBCLOUD LAYER,
C                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
C                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
C                  CHECK FOR BUOYANCY AND SET FLAGS
C                  -------------------------------------
C
  200 CONTINUE
      DO 290 JK=KLEVM1,2,-1
      IS=0
      DO 210 JL=KSTART,KSTOP
      IS=IS+ICVMGT(1,0,KLAB(JL,JK+1).EQ.1)
      LOFLAG(JL)=LCVMGT(.TRUE.,.FALSE.,KLAB(JL,JK+1).EQ.1)
  210 CONTINUE
      IF(IS.EQ.0) GO TO 290
      DO 220 JL=KSTART,KSTOP
      IF(LOFLAG(JL)) THEN
         PQU(JL,JK)=PQU(JL,JK+1)
         PTU(JL,JK)=(CPD*PTU(JL,JK+1)+PGEOH(JL,JK+1)
     1              -PGEOH(JL,JK))*RCPD
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
         IF(ZBUO.GT.0.) KLAB(JL,JK)=1
         ZQOLD(JL)=PQU(JL,JK)
      END IF
  220 CONTINUE
C
      IK=JK
      ICALL=1
      CALL CUADJTQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *     KSTART,   KSTOP,    KLEN,
     *     PAPH,     PTU,      PQU,      LOFLAG,   ICALL
     *  )
C
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 240 JL=KSTART,KSTOP
      IF(LOFLAG(JL).AND.PQU(JL,JK).NE.ZQOLD(JL)) THEN
         KLAB(JL,JK)=2
         PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
         IF(ZBUO.GT.0.) THEN
            KCBOT(JL)=JK
      LDCUM(JL)=.TRUE.
         END IF
      END IF
  240 CONTINUE
C
  290 CONTINUE
C
      RETURN
      END
