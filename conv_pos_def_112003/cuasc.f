      SUBROUTINE CUASC
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     PTENH,    PQENH,    PUEN,     PVEN,
     *     PTEN,     PQEN,     PQSEN,
     *     PGEO,     PGEOH,    PAP,      PAPH,
     *     PQTE,     PVERV,    KLWMIN,   PDQPBL,
     *     LDLAND,   LDCUM,    KTYPE,    KLAB,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     PMFU,     PMFD,     PMFUB,    PENTR,
     *     PMFUS,    PMFUQ,
     *     PMFUL,    PLUDE,    PDMFUP,
     *     KCBOT,    KCTOP,    KCTOP0,   KCUM
     *    )
C
C          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
C          FOR CUMULUS PARAMETERIZATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
C          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
C           FLUXES AS WELL AS PRECIPITATION RATES)
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C
C          METHOD.
C          --------
C          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
C          AND THEN CALCULATE MOIST ASCENT FOR
C          ENTRAINING/DETRAINING PLUME.
C          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
C          SHALLOW AND DEEP CUMULUS CONVECTION.
C          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
C          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
C          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
C
C          EXTERNALS
C          ---------
C          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
C          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
C          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
C
C          REFERENCE
C          ---------
C          (TIEDTKE,1989)
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KHOR,KLEV),       PQENH(KHOR,KLEV),
     *         PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PTEN(KHOR,KLEV),        PQEN(KHOR,KLEV),
     *         PGEO(KHOR,KLEV),        PGEOH(KHOR,KLEV),
     *         PAP(KHOR,KLEV),         PAPH(KHOR,KLEVP1),
     *         PQSEN(KHOR,KLEV),       PQTE(KHOR,KLEV),
     *         PVERV(KHOR,KLEV),       PDQPBL(KHOR)
C
      REAL     PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PUU(KHOR,KLEV),         PVU(KHOR,KLEV),
     *         PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV),
     *         PMFUB(KHOR),            PENTR(KHOR),
     *         PMFUS(KHOR,KLEV),       PMFUQ(KHOR,KLEV),
     *         PLU(KHOR,KLEV),         PLUDE(KHOR,KLEV),
     *         PMFUL(KHOR,KLEV),       PDMFUP(KHOR,KLEV)
      INTEGER  KLWMIN(KHOR),           KTYPE(KHOR),
     *         KLAB(KHOR,KLEV),        KCBOT(KHOR),
     *         KCTOP(KHOR),            KCTOP0(KHOR)
      LOGICAL  LDLAND(KHOR),           LDCUM(KHOR)
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
     *        ZDMFDE (JPHR),
     *        ZDMFEN (JPHR),
     *        ZMFUU  (JPHR),
     *        ZMFUV  (JPHR),
     *        ZPBASE (JPHR),
     *        ZQOLD  (JPHR)
C
C----------------------------------------------------------------------
C
C*    1.           SPECIFY PARAMETERS
C                  ------------------
C
  100 CONTINUE
      ZDMFMAX=CMFCMAX/(KLEV/4)
C
C
C----------------------------------------------------------------------
C
C     2.           SET DEFAULT VALUES
C                  ------------------
C
  200 CONTINUE
      DO 210 JL=KSTART,KSTOP
      IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
  210 CONTINUE
      DO 230 JK=1,KLEV
      DO 220 JL=KSTART,KSTOP
      PLU(JL,JK)=0.
      PMFU(JL,JK)=0.
      PMFUS(JL,JK)=0.
      PMFUQ(JL,JK)=0.
      PMFUL(JL,JK)=0.
      PLUDE(JL,JK)=0.
      PDMFUP(JL,JK)=0.
      IF(.NOT.LDCUM(JL).OR.KTYPE(JL).EQ.3) KLAB(JL,JK)=0
      IF(.NOT.LDCUM(JL).AND.PAPH(JL,JK).LT.4.E4) KCTOP0(JL)=JK
  220 CONTINUE
  230 CONTINUE
C
C
C----------------------------------------------------------------------
C
C     3.0          INITIALIZE VALUES AT LIFTING LEVEL
C                  ----------------------------------
C
  300 CONTINUE
      DO 310 JL=KSTART,KSTOP
      KCTOP(JL)=KLEVM1
         IF(.NOT.LDCUM(JL)) THEN
            KCBOT(JL)=KLEVM1
            PMFUB(JL)=0.
            PQU(JL,KLEV)=0.
         END IF
      PMFU(JL,KLEV)=PMFUB(JL)
      PMFUS(JL,KLEV)=PMFUB(JL)*(CPD*PTU(JL,KLEV)+PGEOH(JL,KLEV))
      PMFUQ(JL,KLEV)=PMFUB(JL)*PQU(JL,KLEV)
         IF(LMFDUDV) THEN
            ZMFUU(JL)=PMFUB(JL)*PUU(JL,KLEV)
            ZMFUV(JL)=PMFUB(JL)*PVU(JL,KLEV)
         END IF
  310 CONTINUE
      DO 320 JL=KSTART,KSTOP
      LDCUM(JL)=.FALSE.
  320 CONTINUE
C
C
C----------------------------------------------------------------------
C
C     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
C                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
C                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
C                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
C                  -------------------------------------------------
C
  400 CONTINUE
      DO 480 JK=KLEVM1,2,-1
C
C                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
C                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
C                  ----------------------------------------------------
C
      IK=JK
      CALL CUBASMC
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVM1,   IK,
     *     KSTART,   KSTOP,    KLEN,
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     *     PVERV,    PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,
     *     PMFU,     PMFUB,    PENTR,    KCBOT,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   ZMFUU,    ZMFUV)
C
      IS=0
      DO 410 JL=KSTART,KSTOP
      IS=IS+KLAB(JL,JK+1)
      IF(KLAB(JL,JK+1).EQ.0) KLAB(JL,JK)=0
      LOFLAG(JL)=LCVMGT(.TRUE.,.FALSE.,KLAB(JL,JK+1).GT.0)
  410 CONTINUE
      IF(IS.EQ.0) GO TO 480
C
C
C*                  SPECIFY ENTRAINMENT RATES IN *CUENTR*
C                   -------------------------------------
C
      IK=JK
      CALL CUENTR
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *     KSTART,   KSTOP,    KLEN,     KLAB,
     *     PTENH,    PQENH,    PQTE,     PAPH,     PAP,
     *     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     *     ZPBASE,   PMFU,     PENTR,
     *     ZDMFEN,   ZDMFDE
     *  )
C
C
C
C                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
C                  ---------------------------------------------------
C
      DO 420 JL=KSTART,KSTOP
      IF(LOFLAG(JL)) THEN
         ZDMFEN(JL)=MIN(ZDMFEN(JL),ZDMFMAX)
         IF(LMFNEW)THEN
	   ZDMFDE(JL)=MIN(ZDMFDE(JL),1.0 *PMFU(JL,JK+1))
         ELSE
	   ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75*PMFU(JL,JK+1))
	 ENDIF
         PMFU(JL,JK)=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
         ZQEEN=PQENH(JL,JK+1)*ZDMFEN(JL)
         ZSEEN=(CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFEN(JL)
         ZSCDE=(CPD*PTU(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFDE(JL)
         ZQUDE=PQU(JL,JK+1)*ZDMFDE(JL)
         PLUDE(JL,JK)=PLU(JL,JK+1)*ZDMFDE(JL)
         ZMFUSK=PMFUS(JL,JK+1)+ZSEEN-ZSCDE
         ZMFUQK=PMFUQ(JL,JK+1)+ZQEEN-ZQUDE
         ZMFULK=PMFUL(JL,JK+1)    -PLUDE(JL,JK)
         PLU(JL,JK)=ZMFULK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
         PQU(JL,JK)=ZMFUQK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
         PTU(JL,JK)=(ZMFUSK*(1./MAX(CMFCMIN,PMFU(JL,JK)))-
     1               PGEOH(JL,JK))*RCPD
         PTU(JL,JK)=MAX(100.,PTU(JL,JK))
         PTU(JL,JK)=MIN(400.,PTU(JL,JK))
         ZQOLD(JL)=PQU(JL,JK)
      END IF
  420 CONTINUE
C
C
C                  DO CORRECTIONS FOR MOIST ASCENT
C                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
C                  -----------------------------------
C
      IK=JK
      ICALL=1
      CALL CUADJTQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *     KSTART,   KSTOP,    KLEN,
     *     PAPH,     PTU,      PQU,      LOFLAG,  ICALL
     *  )
C
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 440 JL=KSTART,KSTOP
      IF(LOFLAG(JL).AND.PQU(JL,JK).NE.ZQOLD(JL)) THEN
         KLAB(JL,JK)=2
         PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     1   PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))
         IF(KLAB(JL,JK+1).EQ.1) ZBUO=ZBUO+0.5
         IF(ZBUO.GT.0..AND.PMFU(JL,JK).GE.0.1*PMFUB(JL)) THEN
            KCTOP(JL)=JK
            LDCUM(JL)=.TRUE.
            ZDNOPRC=CVMGT(3.0E4,1.5E4,LDLAND(JL))
            ZPRCON=CVMGT(0.,CPRCON,ZPBASE(JL)-PAPH(JL,JK).LT.ZDNOPRC)
            ZLNEW=PLU(JL,JK)/(1.+ZPRCON*(PGEOH(JL,JK)-PGEOH(JL,JK+1)))
            IF(KTYPE(JL).EQ.3) ZLNEW=0
            PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
            PLU(JL,JK)=ZLNEW
         ELSE
            KLAB(JL,JK)=0
            PMFU(JL,JK)=0.
         END IF
      END IF
  440 CONTINUE
      DO 455 JL=KSTART,KSTOP
      IF(LOFLAG(JL)) THEN
         PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
         PMFUS(JL,JK)=(CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
         PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
      END IF
  455 CONTINUE
        IF(LMFDUDV) THEN
           DO 460 JL=KSTART,KSTOP
           IF(LOFLAG(JL)) THEN
              ZMFUU(JL)=ZMFUU(JL)+
     1                  ZDMFEN(JL)*PUEN(JL,JK)-ZDMFDE(JL)*PUU(JL,JK+1)
              ZMFUV(JL)=ZMFUV(JL)+
     1                  ZDMFEN(JL)*PVEN(JL,JK)-ZDMFDE(JL)*PVU(JL,JK+1)
              IF(PMFU(JL,JK).GT.0.) THEN
                 PUU(JL,JK)=ZMFUU(JL)*(1./PMFU(JL,JK))
                 PVU(JL,JK)=ZMFUV(JL)*(1./PMFU(JL,JK))
              END IF
           END IF
  460      CONTINUE
        END IF
C
  480 CONTINUE
C
C
C----------------------------------------------------------------------
C
C     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
C                  ----------------------------------------------------
		  
  500 CONTINUE
c
      DO 510 JL=KSTART,KSTOP
      IF(KCTOP(JL).EQ.KLEVM1) LDCUM(JL)=.FALSE.
      KCBOT(JL)=MAX(KCBOT(JL),KCTOP(JL))
  510 CONTINUE
      IS=0
      DO 520 JL=KSTART,KSTOP
      IS=IS+ICVMGT(1,0,LDCUM(JL))
  520 CONTINUE
      KCUM=IS
      IF(IS.EQ.0) GO TO 800
C
C-------------------------------------------------------------------
C								   |
C     NEW RESOLUTION-INDEPENDENT TREATMENT OF THE OVERSHOOT-REGIME |
C     ONLY ACTIVE IN HIGH-RESOLUTION MODE: LMFHRES = .TRUE.	   |
C								   |
C-------------------------------------------------------------------
C
      JL = 1
      IF (LMFHRES.AND.(KTYPE(JL).LT.3)) THEN
C
C-------------------------------------------------------------------
C								   |
C     DETERMINE LEVELS WHERE OVERSHOOT NEEDS TO BE CALCULATED	   |
C	           LABEL WITH KLAB = 3				   |
C								   |
C								   |
C     ICTOP: HIGHEST LEVEL WITH OVERSHOOT .LT. DP_SHOOT		   |
C            HIGHEST LEVEL WITH MASS_FLUX .GT. 0.0		   |
C								   |
C-------------------------------------------------------------------
C
      JK = KCTOP(JL)
525   CONTINUE
      ICTOP         = JK
      JK            = JK-1
      GUARD = PAPH(JL,KCTOP(JL)) - PAPH(JL,JK)
      IF (GUARD.LT.DP_SHOOT ) GOTO 525
C
      DZTOT   = (PGEOH(JL,ICTOP-1) - PGEOH(JL,KCTOP(JL)))/G
      DECAY   = 0.5*DZTOT/
     &          LOG((1+SQRT(1-4*CMFCTOP*(1-CMFCTOP)))/(2*CMFCTOP))
C
C-------------------------------------------------------------------
C								   |
C     DO MOIST ADIABATIC ASCENT IN OVERSHOOT LAYER		   |
C     WITH TURBULENT ENTRAINMENT = TURBULENT DETRAINMENT	   |
C     PLUS MASSIVE DETRAINMENT CAUSING				   | 
C     AN EXPONENTIAL DEACY OF THE MASS FLUX			   |
C-------------------------------------------------------------------
C
      DO 530 JK=KCTOP(JL)-1,ICTOP,-1
C
C-------------------------------------------------------------------
C								   |
C     LABEL OVERSHOOT LAYER: KLAB = 3				   |
C     DO TURBULENT ENTRAINMENT + DETRAINMENT			   | 
C								   |
C-------------------------------------------------------------------
C
c      KLAB  (JL,JK) = 3
      IK=JK
      CALL CUENTR
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *     KSTART,   KSTOP,    KLEN,     KLAB,
     *     PTENH,    PQENH,    PQTE,     PAPH,     PAP,
     *     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     *     ZPBASE,   PMFU,     PENTR,
     *     ZDMFEN,   ZDMFDE
     *  )
C
C-------------------------------------------------------------------
C								   |
C      CALCULATE A EXPONENTIAL DECAY OF THE MASS FLUX		   |
C      DUE TO MASSIVE DETRAINMENT:						   |
C								   |
C								   |
C	  M(Z(KCTOP))                           = M_TOP		   |
C         M(Z(KCTOP)+0.5*(Z(ICTOP-1)-(Z(KCTOP)) = CMFCTOP*M_TOP	   |
C         M(Z(ICTOP-1))                         = 0.0		   |
C								   |
C								   |
C     ICTOP: HIGHEST LEVEL WITH OVERSHOOT .LT. DP_SHOOT		   |
C            HIGHEST LEVEL WITH MASS_FLUX .GT. 0.0		   |
C								   |
C-------------------------------------------------------------------
C
       PMFU(JL,JK)   = PMFU(JL,KCTOP(JL))
     &  *(EXP((PGEOH(JL,ICTOP-1)-PGEOH(JL,JK))       /(DECAY*G)) - 1.0)
     &  /(EXP((PGEOH(JL,ICTOP-1)-PGEOH(JL,KCTOP(JL)))/(DECAY*G)) - 1.0) 
                      
C
C-------------------------------------------------------------------
C								   |
C      DO DRY-ADIABTIC ASCENT WITH ENTRAINMENT AND DETRAINMENT	   |
C      FOR Q, T AND L 						   |
C								   |
C-------------------------------------------------------------------
C
       ZDMFDE(JL)    = PMFU  (JL,JK+1) - PMFU(JL,JK) + ZDMFDE(JL)
       ZDMFEN(JL)=MIN(ZDMFEN(JL),ZDMFMAX)
       IF(LMFNEW)THEN
           IF ( ZDMFDE(JL).GT.PMFU(JL,JK+1) ) THEN
	    ZDMFDE(JL)  = PMFU(JL,JK+1)
            PMFU(JL,JK) = ZDMFEN(JL)
           ENDIF
       ELSE
	   ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75*PMFU(JL,JK+1))
       ENDIF
       ZQEEN       = PQENH(JL,JK+1)*ZDMFEN(JL)
       ZSEEN       = (CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFEN(JL)
       ZSCDE       = (CPD*PTU  (JL,JK+1)+PGEOH(JL,JK+1))*ZDMFDE(JL)
       ZQUDE       = PQU(JL,JK+1)*ZDMFDE(JL)
       PLUDE(JL,JK)= PLU(JL,JK+1)*ZDMFDE(JL)
       ZMFUSK      = PMFUS(JL,JK+1)+ZSEEN-ZSCDE
       ZMFUQK      = PMFUQ(JL,JK+1)+ZQEEN-ZQUDE
       ZMFULK      = PMFUL(JL,JK+1)-PLUDE(JL,JK)
       PLU(JL,JK)  = ZMFULK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
       PQU(JL,JK)  = ZMFUQK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
       PTU(JL,JK)  = (ZMFUSK*(1./MAX(CMFCMIN,PMFU(JL,JK)))-
     1                  PGEOH(JL,JK))*RCPD
       PTU(JL,JK)  = MAX(100.,PTU(JL,JK))
       PTU(JL,JK)  = MIN(400.,PTU(JL,JK))
       ZQOLD(JL)   = PQU(JL,JK)
       IF(LMFDUDV) THEN
            PUU(JL,JK) = PUU(JL,JK+1)
            PVU(JL,JK) = PVU(JL,JK+1)
       ENDIF
C
C-------------------------------------------------------------------
C								   |
C          DO CORRECTIONS FOR MOIST ASCENT			   |
C          BY ADJUSTING T,Q AND L IN *CUADJTQ*			   |
C								   |
C-------------------------------------------------------------------
C
           ZQOLD(JL) = PQU(JL,JK)
           IK        = JK
           ICALL     = 1
           LOFLAG(JL) = .TRUE.
c
           CALL CUADJTQ
     *      (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *       KSTART,   KSTOP,    KLEN,
     *       PAPH,     PTU,      PQU,      LOFLAG,  ICALL
     *  )
C
           PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
           ZDNOPRC=CVMGT(3.0E4,1.5E4,LDLAND(JL))
           ZPRCON=CVMGT(0.,CPRCON,ZPBASE(JL)-PAPH(JL,JK).LT.ZDNOPRC)
           ZLNEW=PLU(JL,JK)/(1.+ZPRCON*(PGEOH(JL,JK)-PGEOH(JL,JK+1)))
           PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
           PLU(JL,JK)=ZLNEW
C
           PMFUL(JL,JK) =  PLU(JL,JK)                  *PMFU(JL,JK)
           PMFUS(JL,JK) = (CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
           PMFUQ(JL,JK) =  PQU(JL,JK)                  *PMFU(JL,JK)
C
 530    CONTINUE
C
C
C-------------------------------------------------------------------
C								   |
C          DETRAIN REMAINING LIQUID WATER AT THE FIRST LEVEL	   |
C	   WHERE MASS_FLUX = 0.0				   |
C								   |
C-------------------------------------------------------------------
C
          PLUDE(JL,ICTOP-1) = PMFUL(JL,ICTOP)
C
        ENDIF
C
C**********  END OF RESOLUTION_INDEPENDENT TREATMENT ***********


C       DO ORIGINAL TIEDTKE OVERSHOOT-TREATMENT IF
C       MODEL IS NOT IN HIGH_RES MODE
C       -------------------------------------------
C
      IF (.NOT.LMFHRES.OR.KTYPE(JL).EQ.3) THEN

CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 540 JL=KSTART,KSTOP
      IF(LDCUM(JL)) THEN
         JK=KCTOP(JL)-1
         ZZDMF=CMFCTOP
         ZDMFDE(JL)=(1.-ZZDMF)*PMFU(JL,JK+1)
         PLUDE(JL,JK)=ZDMFDE(JL)*PLU(JL,JK+1)
         PMFU(JL,JK)=PMFU(JL,JK+1)-ZDMFDE(JL)
         ZLNEW=PLU(JL,JK)
         IF(KTYPE(JL).EQ.3) ZLNEW=0
         PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
         PLU(JL,JK)=ZLNEW
         PMFUS(JL,JK)=(CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
         PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
         PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
         PLUDE(JL,JK-1)=PMFUL(JL,JK)
      END IF
  540 CONTINUE
        IF(LMFDUDV) THEN
CDIR$      IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
           DO 550 JL=KSTART,KSTOP
           IF(LDCUM(JL)) THEN
              JK=KCTOP(JL)-1
              PUU(JL,JK)=PUU(JL,JK+1)
              PVU(JL,JK)=PVU(JL,JK+1)
           END IF
  550      CONTINUE
        END IF
C
      END IF
  800 CONTINUE
C
      RETURN
      END
