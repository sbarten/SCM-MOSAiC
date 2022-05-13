      SUBROUTINE CUBASE4
     *    (KCUCALL,  KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,   KLEVM1,
     *     K2LP2,
     *     PTENH,    PQENH,    PGEOH,    PAPH,

C LG- adding the tracers

     *   KTRAC,
     *   PXTU,

C LG- end     

     *     PTU,      PQU,      PLU,
     *     PUEN,    PVEN,     PUU,      PVU,
     *     LDCUM,    KCBOT,    KLAB 
CTCL ...
     *   , KBOTSC  , LDSC )
C
C          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
C          FOR CUMULUS PARAMETERIZATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C          E.VAN MEIJGAARD   KNMI                MODIF. 01/97
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
CJHC*CALL COMCON
CJHC*CALL COMCUMF
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KLP2,KLEV),       PQENH(KLP2,KLEV),
     *         PGEOH(KLP2,KLEV),       PAPH(KLP2,KLEVP1)
C
      REAL     PTU(KLP2,KLEV),         PQU(KLP2,KLEV),
     *         PLU(KLP2,KLEV)
      REAL     PUEN(K2LP2,KLEV),        PVEN(K2LP2,KLEV),
     *         PUU(KLP2,KLEV),          PVU(KLP2,KLEV)
      INTEGER  KLAB(KLP2,KLEV),        KCBOT(KLP2)
      LOGICAL  LDCUM(KLP2)
CTCL ...
      INTEGER KBOTSC(KLP2)
      LOGICAL LDSC(KLP2)
CTCL ...
C
      INCLUDE 'paramh.h'
      REAL     ZQOLD(JPHR)
      REAL     ZPH(JPHR)
      LOGICAL  LOFLAG(JPHR)
COBC
     *        , LCVMGT

C LG- adding the tracers

      REAL   PXTU(KLON,KLEV,KTRAC)

C LG- end

CTCL ...
      LOGICAL LLO1
CTCL ...
C
C
C----------------------------------------------------------------------
C
C     1.           INITIALIZE VALUES AT LIFTING LEVEL
C                  ----------------------------------
C
  100 CONTINUE
      DO 110 JL=KIDIA,KFDIA
      KLAB(JL,KLEV)=1
      KCBOT(JL)=KLEVM1
      LDCUM(JL)=.FALSE.
CTCL ...
      KBOTSC(JL)=KLEVP1
      LDSC(JL)=.FALSE.
CTCL ...
      PUU(JL,KLEV)=PUEN(JL,KLEV)*(PAPH(JL,KLEVP1)-PAPH(JL,KLEV))
      PVU(JL,KLEV)=PVEN(JL,KLEV)*(PAPH(JL,KLEVP1)-PAPH(JL,KLEV))
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
      DO 210 JL=KIDIA,KFDIA
      IS=IS+ICVMGT(1,0,KLAB(JL,JK+1).EQ.1)
CTCL ...
      if (kcucall.eq.4) then
        LOFLAG(JL)=LCVMGT(.TRUE.,.FALSE.,KLAB(JL,JK+1).EQ.1)
      else if (kcucall.eq.6) then
        LOFLAG(JL)=LCVMGT(.TRUE.,.FALSE.,
     S     (KLAB(JL,JK+1).EQ.1.OR.(LDCUM(JL).AND.KCBOT(JL).EQ.JK+1)))
      endif
CTCL ...
      ZPH(JL)=PAPH(JL,JK)
  210 CONTINUE
      IF(IS.EQ.0) GO TO 290
      DO 220 JL=KIDIA,KFDIA
      IF(LOFLAG(JL)) THEN
         PQU(JL,JK)=PQU(JL,JK+1)
         PTU(JL,JK)=(CPD*PTU(JL,JK+1)+PGEOH(JL,JK+1)
     1              -PGEOH(JL,JK))*RCPD
        if (kcucall.eq.4) then
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
         IF(ZBUO.GT.0.) KLAB(JL,JK)=1
        endif
         ZQOLD(JL)=PQU(JL,JK)
      END IF
  220 CONTINUE
C

C LG- adding the tracers

      DO 2204 JT=1,KTRAC
      DO 2202 JL=1,KLON
      IF(LOFLAG(JL)) THEN
       PXTU(JL,JK,JT)=PXTU(JL,JK+1,JT)
      ENDIF
 2202 CONTINUE
 2204 CONTINUE
 
C LG- end

      IK=JK
      ICALL=1
      CALL CUADJTQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     IK,
     *     ZPH,      PTU,      PQU,      LOFLAG,   ICALL)
C
CDIR$ IVDEP
      DO 240 JL=KIDIA,KFDIA
CTCL ...
      if (kcucall.eq.4) then
CTCL ...
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
CTCL ...
      else if (kcucall.eq.6) then
       IF(LOFLAG(JL)) THEN
       IF(PQU(JL,JK).EQ.ZQOLD(JL)) THEN
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
	 IF(ZBUO.GT.0) KLAB(JL,JK)=1
       ELSE
	 KLAB(JL,JK)=2
	 PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
	 LLO1=ZBUO.GT.0..AND.KLAB(JL,JK+1).EQ.1
	 IF(LLO1) THEN
	   LDSC(JL)=.TRUE.
	   KBOTSC(JL)=JK
	   KCBOT(JL)=JK
	   LDCUM(JL)=.TRUE.
	 ENDIF
       ENDIF
       ENDIF
      endif
CTCL ...
  240 CONTINUE
C
C             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
C             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
C
      IF(LMFDUDV) THEN
         DO 250 JL=KIDIA,KFDIA
         IF(JK.GE.KCBOT(JL)) THEN
            PUU(JL,KLEV)=PUU(JL,KLEV)+
     *                   PUEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
            PVU(JL,KLEV)=PVU(JL,KLEV)+
     *                   PVEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
         END IF
 250     CONTINUE
      END IF
C
  290 CONTINUE
C
C
      IF(LMFDUDV) THEN
         DO 310 JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            IKB=KCBOT(JL)
            ZZ=1./(PAPH(JL,KLEVP1)-PAPH(JL,IKB))
            PUU(JL,KLEV)=PUU(JL,KLEV)*ZZ
            PVU(JL,KLEV)=PVU(JL,KLEV)*ZZ
         ELSE
            PUU(JL,KLEV)=PUEN(JL,KLEVM1)
            PVU(JL,KLEV)=PVEN(JL,KLEVM1)
         END IF
 310     CONTINUE
      END IF
C
      RETURN
      END
