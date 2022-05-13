      SUBROUTINE CUMASTR
     *    (KHOR,     KHOR2,    KLON,     KLEV,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEVP1,   KLEVM1,   KROW,     ILAB,
     *     CONACC,
     *     PTEN,     PQEN,     PUEN,     PVEN,
     *     PVERV,    PQSEN,    PTS,      PQHFL,
     *     PAP,      PAPH,     PGEO,     LDLAND,
     *     PTTE,     PQTE,     PVOM,     PVOL,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,
     *     PRFLCK,
     *     LDCUM,    KTYPE,    KCBOT,    KCTOP,
     *     PTU,      PQU,      PLU,      PLUDE,
     *     PMFU,     PMFD,     PRAIN,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,
CHL --- lines added
     *     NSTEP , NSTART,   NSTOP,    NPRINT,  TWODT ,
CEVM---    CONDENSATION RATES
     *     TEMFCD ,   QEMFCD , XEMFCD
     *     )
C
C**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
C
C     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
C
C
C     PURPOSE
C     -------
C
C          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
C     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
C     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
C     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
C     SATURATED CUMULUS DOWNDRAFTS.
C
C**   INTERFACE.
C     ----------
C
C          *CUMASTR* IS CALLED FROM *CUCALL*
C     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
C     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
C     IT RETURNS ITS OUTPUT TO THE SAME SPACE
C      1.MODIFIED TENDENCIES OF MODEL VARIABLES
C      2.RATES OF CONVECTIVE PRECIPITATION
C        (USED IN SUBROUTINE SURF)
C      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
C        (USED IN SUBROUTINE CLOUD)
C
C     METHOD
C     -------
C
C     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
C        (1) DEFINE CONSTANTS AND PARAMETERS
C        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
C            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
C        (3) CALCULATE CLOUD BASE IN 'CUBASE'
C            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
C        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
C        (5) DO DOWNDRAFT CALCULATIONS:
C              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
C              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
C              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
C                  EFFECT OF CU-DOWNDRAFTS
C        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
C        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
C            DO EVAPORATION IN SUBCLOUD LAYER
C        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
C        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
C
C     EXTERNALS.
C     ----------
C
C       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
C       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
C       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
C       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
C       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
C       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
C       CUDQDT: UPDATES TENDENCIES FOR T AND Q
C       CUDUDV: UPDATES TENDENCIES FOR U AND V
C
C     SWITCHES.
C     --------
C
C          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
C          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
C          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
C          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
C          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
C
C
C     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
C     ------------------------------------------------
C     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
C     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
C     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
C     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
C     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVE
C     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
C     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
C     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
C     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
C
C     REFERENCE.
C     ----------
C
C          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
C
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
      INCLUDE 'comph2.h'
C
      LOGICAL LLO1
C
      REAL     PTEN(KHOR,KLEV),        PQEN(KHOR,KLEV),
     *         PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PTTE(KHOR,KLEV),        PQTE(KHOR,KLEV),
     *         PVOM(KHOR2,KLEV),       PVOL(KHOR2,KLEV),
     *         PQSEN(KHOR,KLEV),       PGEO(KHOR,KLEV),
     *         PAP(KHOR,KLEV),         PAPH(KHOR,KLEVP1),
     *         PVERV(KHOR,KLEV),       PTS(KHOR),
     *         PQHFL(KHOR)
      REAL     PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PLU(KHOR,KLEV),         PLUDE(KHOR,KLEV),
     *         PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV),
     *         PAPRC(KHOR),            PAPRS(KHOR),
     *         PRSFC(KHOR),            PSSFC(KHOR),
     *         PRAIN(KHOR)
      INTEGER  KCBOT(KHOR),            KCTOP(KHOR),
     *         KTYPE(KHOR)
      INTEGER  ILAB(KHOR,KLEV)
      LOGICAL  LDLAND(KHOR),           LDCUM(KHOR)
C
      REAL PRFLCK(KHOR,KLEV)
C  
      REAL TEMFCD(KLEV),QEMFCD(KLEV),XEMFCD(KLEV)
C
C
C     -------------------------------------------------------
C     DECLARATION OF
C        LOGICAL AND INTEGER ONE-DIMENSIONAL WORKING ARRAYS
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
C
      LOGICAL
     *       LODDRAF(JPHR)
      INTEGER
     *       IDTOP  (JPHR),
     *       ICTOP0 (JPHR),
     *       ILWMIN (JPHR),
     *       ITMELT (JPHR)
C
C     -------------------------------------------------------
C     DECLARATION OF REAL ONE-DIMENSIONAL WORKING ARRAYS
C
      REAL
     *       ZDQCV  (JPHR),
     *       ZDQPBL (JPHR),
     *       ZENTR  (JPHR),
     *       ZHCBASE(JPHR),
     *       ZMFUB  (JPHR),
     *       ZMFUB1 (JPHR),
     *       ZRFL   (JPHR),

C
C     -------------------------------------------------------
C     DECLARATION OF REAL TWO-DIMENSIONAL WORKING ARRAYS
C
     *       ZDMFDP (JPHR,MLEV),
     *       ZDMFUP (JPHR,MLEV),
     *       ZGEOH  (JPHR,MLEV),
     *       ZMFDQ  (JPHR,MLEV),
     *       ZMFDS  (JPHR,MLEV),
     *       ZMFUL  (JPHR,MLEV),
     *       ZMFUQ  (JPHR,MLEV),
     *       ZMFUS  (JPHR,MLEV),
     *       ZQD    (JPHR,MLEV),
     *       ZQENH  (JPHR,MLEV),
     *       ZQSENH (JPHR,MLEV),
     *       ZTD    (JPHR,MLEV),
     *       ZTENH  (JPHR,MLEV),
     *       ZUD    (JPHR,MLEV),
     *       ZUU    (JPHR,MLEV),
     *       ZVD    (JPHR,MLEV),
     *       ZVU    (JPHR,MLEV)
C
C
C
C
C---------------------------------------------------------------------
C
C     1.           SPECIFY CONSTANTS AND PARAMETERS
C                  --------------------------------
C
  100 CONTINUE
C
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZCONS2=1./(G*ZTMST)
C
C
C----------------------------------------------------------------------
C
C*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
C                  ---------------------------------------------------
C
  200 CONTINUE
      CALL CUINI
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     *     PVERV,    PGEO,     PAPH,     ZGEOH,
     *     ZTENH,    ZQENH,    ZQSENH,   ILWMIN,
     *     PTU,      PQU,      ZTD,      ZQD,
     *     ZUU,      ZVU,      ZUD,      ZVD,
     *     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     *     ZMFUQ,   ZMFDQ,   ZDMFUP,   ZDMFDP,   PRFLCK,
     *     PLU,      PLUDE,    ILAB
     *   )
C
C
C---------------------------------------------------------------------
C
C*    3.0          CLOUD BASE CALCULATIONS
C                  -----------------------
C
  300 CONTINUE
C
C*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
C                  ---------------------------------------
C
      CALL CUBASE
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   KLEVM1,
     *     KSTART,   KSTOP,    KLEN,
     *     ZTENH,    ZQENH,    ZGEOH,    PAPH,
     *     PTU,      PQU,      PLU,
     *     LDCUM,    KCBOT,    ILAB
     *   )
C
C*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
C*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
C                  -----------------------------------------
C
       JK=1
       DO 310 JL=KSTART,KSTOP
       ZDQCV(JL) =PQTE(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
       ZDQPBL(JL)=0.0
       IDTOP(JL)=0
  310  CONTINUE
       DO 320 JK=2,KLEV
       DO 315 JL=KSTART,KSTOP
       ZDQCV(JL)=ZDQCV(JL)+PQTE(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
       IF(JK.GE.KCBOT(JL)) ZDQPBL(JL)=ZDQPBL(JL)+PQTE(JL,JK)
     1                               *(PAPH(JL,JK+1)-PAPH(JL,JK))
  315 CONTINUE
  320 CONTINUE
      DO 330 JL=KSTART,KSTOP
      KTYPE(JL)=ICVMGT(1,2,ZDQCV(JL).GT.MAX(0.,-1.1*PQHFL(JL)*G))
  330 CONTINUE
C
C     IF (LMFQHFL = TRUE) INVOKE ALTERNATIVE CLOSURE
C
      IF (LMFQHFL) THEN
	DO 321 JL=KSTART,KSTOP
	  ZDQPBL(JL)=-PQHFL(JL)*G
 321    CONTINUE
      ENDIF
C
C*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
C*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
C*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
C                  ------------------------------------------
C
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 340 JL=KSTART,KSTOP
      IKB=KCBOT(JL)
      ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-ZQENH(JL,IKB)
      ZDQMIN=MAX(0.01*ZQENH(JL,IKB),1.E-10)
      LLO1=ZDQPBL(JL).GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM(JL)
      ZMFUB(JL)=CVMGT(ZDQPBL(JL)/(G*MAX(ZQUMQE,ZDQMIN)),0.01,LLO1)
      ZMFUB(JL)=MIN(ZMFUB(JL),CMFCMAX)
      IF(.NOT.LLO1) LDCUM(JL)=.FALSE.
      ZENTR(JL)=CVMGT(ENTRPEN,ENTRSCV,KTYPE(JL).EQ.1)
  340 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
C                  -------------------------------------------
C
C
  400 CONTINUE
C
C*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
C*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
C*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
C                  -------------------------------------------------
C
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 410 JL=KSTART,KSTOP
      IKB=KCBOT(JL)
      ZHCBASE(JL)=CPD*PTU(JL,IKB)+ZGEOH(JL,IKB)+ALV*PQU(JL,IKB)
      ICTOP0(JL)=KCBOT(JL)-1
  410 CONTINUE
      ZALVDCP=ALV/CPD
      ZQALV=1./ALV
      DO 430 JK=KLEVM1,3,-1
c$dir no_recurrence, force_vector, force_parallel_ext
CDIR$ IVDEP
      DO 420 JL=KSTART,KSTOP
      ZHSAT=CPD*ZTENH(JL,JK)+ZGEOH(JL,JK)+ALV*ZQSENH(JL,JK)
      ZGAM=C5LES*ZALVDCP*ZQSENH(JL,JK)/
     1     ((1.-VTMPC1*ZQSENH(JL,JK))*(ZTENH(JL,JK)-C4LES)**2)
      ZZZ=CPD*ZTENH(JL,JK)*0.608
      ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)*
     1            MAX(ZQSENH(JL,JK)-ZQENH(JL,JK),0.)
      IF(JK.LT.ICTOP0(JL).AND.ZHCBASE(JL).GT.ZHHAT) ICTOP0(JL)=JK
  420 CONTINUE
  430 CONTINUE
C
C*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
C                  --------------------------------------------
C
      CALL CUASC
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     ZTENH,    ZQENH,    PUEN,     PVEN,
     *     PTEN,     PQEN,     PQSEN,
     *     PGEO,     ZGEOH,    PAP,      PAPH,
     *     PQTE,     PVERV,    ILWMIN,   ZDQPBL,
     *     LDLAND,   LDCUM,    KTYPE,    ILAB,
     *     PTU,      PQU,      PLU,      ZUU,      ZVU,
     *     PMFU,     PMFD,     ZMFUB,    ZENTR,
     *     ZMFUS,    ZMFUQ,
     *     ZMFUL,    PLUDE,    ZDMFUP,
     *     KCBOT,    KCTOP,    ICTOP0,   ICUM
     *    )
C
      IF(ICUM.EQ.0) THEN
	ITOPM2=1
        GO TO 1000
      ENDIF
C
C*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
C              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
C              -----------------------------------------------------
C
C
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 440 JL=KSTART,KSTOP
      ZPBMPT=PAPH(JL,KCBOT(JL))-PAPH(JL,KCTOP(JL))
      IF(LDCUM(JL).AND.KTYPE(JL).EQ.1.AND.ZPBMPT.LT.2.E4) KTYPE(JL)=2
      IF(LDCUM(JL)) ICTOP0(JL)=KCTOP(JL)
      IF(KTYPE(JL).EQ.2) ZENTR(JL)=ENTRSCV
      ZRFL(JL)=ZDMFUP(JL,1)
  440 CONTINUE
      DO 460 JK=2,KLEV
      DO 450 JL=KSTART,KSTOP
      ZRFL(JL)=ZRFL(JL)+ZDMFUP(JL,JK)
  450 CONTINUE
  460 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
C                  ------------------------------
C
C
  500 CONTINUE
C
      IF(LMFDD) THEN
C
C*             (A) DETERMINE LFS IN 'CUDLFS'
C                  -------------------------
C
      CALL CUDLFS
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     ZTENH,    ZQENH,    PUEN,     PVEN,
     *     ZGEOH,    PAPH,     LDLAND,   ZDQPBL,
     *     PTU,      PQU,      PLU,      ZUU,      ZVU,
     *     LDCUM,    KCBOT,    KCTOP,    ZMFUB,    ZRFL,
     *     ZTD,      ZQD,      ZUD,      ZVD,
     *     PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,
     *     IDTOP,    LODDRAF
     *   )
C
C*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
C                  -----------------------------------------------
C
      CALL CUDDRAF
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     ZTENH,    ZQENH,    PUEN,     PVEN,
     *     ZGEOH,    PAPH,     ZRFL,     PLU,
     *     KTYPE,    KCBOT,
     *     ZTD,      ZQD,      ZUD,      ZVD,
     *     PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,
     *     IDTOP,    LODDRAF
     *   )
C
C*            (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
C                  DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
C                  --------------------------------------------
C
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
         DO 520 JL=KSTART,KSTOP
         IF(LODDRAF(JL)) THEN
            IKB=KCBOT(JL)
            LLO1=PMFD(JL,IKB).LT.0.
            ZEPS=CVMGT(CMFDEPS,0.,LLO1)
            ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-
     1      ZEPS*ZQD(JL,IKB)-(1.-ZEPS)*ZQENH(JL,IKB)
            ZDQMIN=MAX(0.01*ZQENH(JL,IKB),1.E-10)
            LLO1=ZDQPBL(JL).GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM(JL)
            ZMFUB1(JL)=CVMGT(ZDQPBL(JL)/(G*MAX(ZQUMQE,ZDQMIN)),
     1                 ZMFUB(JL),LLO1)
            ZMFUB1(JL)=CVMGT(ZMFUB1(JL),ZMFUB(JL),
     1                 (KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.2).AND.
     1                 ABS(ZMFUB1(JL)-ZMFUB(JL)).LT.0.2*ZMFUB(JL))
         END IF
  520    CONTINUE
         DO 540 JK=1,KLEV
c$dir no_recurrence, force_vector, force_parallel_ext
         DO 530 JL=KSTART,KSTOP
         IF(LODDRAF(JL)) THEN
            ZFAC=ZMFUB1(JL)/MAX(ZMFUB(JL),1.E-10)
            PMFD(JL,JK)=PMFD(JL,JK)*ZFAC
            ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZFAC
            ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZFAC
            ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZFAC
         END IF
  530    CONTINUE
  540    CONTINUE
         DO 550 JL=KSTART,KSTOP
         IF(LODDRAF(JL)) ZMFUB(JL)=ZMFUB1(JL)
  550    CONTINUE
C
      END IF
C
C
C-----------------------------------------------------------------------
C
C*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
C*                 FOR PENETRATIVE CONVECTION (TYPE=1),
C*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
C*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
C                  -------------------------------------------------
C
  600 CONTINUE
      CALL CUASC
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     ZTENH,    ZQENH,    PUEN,     PVEN,
     *     PTEN,     PQEN,     PQSEN,
     *     PGEO,     ZGEOH,    PAP,      PAPH,
     *     PQTE,     PVERV,    ILWMIN,   ZDQPBL,
     *     LDLAND,   LDCUM,    KTYPE,    ILAB,
     *     PTU,      PQU,      PLU,      ZUU,      ZVU,
     *     PMFU,     PMFD,     ZMFUB,    ZENTR,
     *     ZMFUS,    ZMFUQ,
     *     ZMFUL,    PLUDE,    ZDMFUP,
     *     KCBOT,    KCTOP,    ICTOP0,   ICUM
     *    )
C
C
C-----------------------------------------------------------------------
C
C*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
C*                 AND EVAPORATION OF RAIN IN SUBCLOUD LAYER
C                  ------------------------------------------
C
  700 CONTINUE
      CALL CUFLX
     *    (KHOR,     KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     PQEN,     PQSEN,    ZTENH,    ZQENH,
     *     PAPH,     LDLAND,   ZGEOH,
     *     KCBOT,    KCTOP,    IDTOP,
     *     KTYPE,    LODDRAF,  LDCUM,
     *     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     *     ZMFUQ,    ZMFDQ,    ZMFUL,    PLUDE,
     *     ZDMFUP,  ZDMFDP,  ZRFL,    PRAIN,   PRFLCK,
     *     ITOPM2)
C
      IBOT=KLEV
      DO 765 JL=KSTART,KSTOP
      ZRFL(JL)=MAX(0.,ZRFL(JL))
      IBOT1=ICVMGT(KCBOT(JL),KLEV,LDCUM(JL).AND.ZRFL(JL).GT.0.)
      IBOT=MIN(IBOT,IBOT1)
  765 CONTINUE
      IF(IBOT.EQ.KLEV) GO TO 800
      ZCUCOV=0.05
      DO 780 JK=IBOT,KLEV
      DO 775 JL=KSTART,KSTOP
      IF(LDCUM(JL).AND.JK.GE.KCBOT(JL).AND.ZRFL(JL).GT.0.) THEN
         ZRNEW=(MAX(0.,SQRT(ZRFL(JL)/ZCUCOV)-
     1         CEVAPCU(JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
     1         *MAX(0.,PQSEN(JL,JK)-PQEN(JL,JK))))**2*ZCUCOV
         ZRMIN=ZRFL(JL)-ZCUCOV*MAX(0.,0.8*PQSEN(JL,JK)-PQEN(JL,JK))
     1         *ZCONS2*(PAPH(JL,JK+1)-PAPH(JL,JK))
         ZRNEW=MAX(ZRNEW,ZRMIN)
         ZRFLN=MAX(ZRNEW,0.)
      ZRFLN=MIN(ZRNEW,ZRFL(JL))
      IF((ZDMFUP(JL,JK)+ZRFLN-ZRFL(JL)).LT.0.) THEN
       ZRFLN=ZRFL(JL)-ZDMFUP(JL,JK)
       ZDMFUP(JL,JK)=0.
      ELSE
         ZDMFUP(JL,JK)=ZDMFUP(JL,JK)+ZRFLN-ZRFL(JL)
      END IF
         ZRFL(JL)=ZRFLN
      PRFLCK(JL,JK)=ZRFL(JL)
C
      END IF
  775 CONTINUE
  780 CONTINUE
C
C
C----------------------------------------------------------------------
C
C*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
C                  --------------------------------------------------
C
  800 CONTINUE
      CALL CUDTDQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     CONACC,
     *     ITOPM2,   PAPH,     PGEO,     PTS,      LDLAND,
     *     LDCUM,    PTEN,     PTTE,     PQTE,
     *     ZMFUS,    ZMFDS,    ZMFUQ,    ZMFDQ,
     *     ZMFUL,    ZDMFUP,   ZDMFDP,   PLUDE,
     *     PRAIN,    ZRFL,     PSRAIN,   PSEVAP,   PSHEAT,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,
C
     *     TWODT
     *   )
C
C
C----------------------------------------------------------------------
C
C*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
C                  --------------------------------------------------
C
  900 CONTINUE
      IF(LMFDUDV) THEN
      CALL CUDUDV
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     ITOPM2,   KTYPE,    KCBOT,    PAPH,     LDCUM,
     *     PUEN,     PVEN,     PVOM,     PVOL,
     *     ZUU,      ZUD,      ZVU,      ZVD,
     *     PMFU,     PMFD,     PSDISS
     *   )
C
      END IF
C
 1000 CONTINUE
C
C
C
C----------------------------------------------------------------------
C
C		PRINT HALF LEVEL FIELDS
C		-----------------------
C
      IF (LMFPRINT) THEN
C
      CALL CUPR
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,   NSTEP,    NPRINT,
     *     KCBOT,    KCTOP,    KTYPE,    ITOPM2, 
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     *     PVERV,    PGEO,     PAP,	 PAPH,     ZGEOH,
     *     ZTENH,    ZQENH,    ZQSENH,   ILWMIN,
     *     PTU,      PQU,      ZTD,      ZQD,
     *     ZUU,      ZVU,      ZUD,      ZVD,
     *     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     *     ZMFUQ,    ZMFDQ,    ZDMFUP,   ZDMFDP,   PRFLCK,
     *     ZMFUL,    PLU,      PLUDE,    ILAB,
CEVM---    CONDENSATION RATES
     *     TEMFCD ,   QEMFCD , XEMFCD
     *   )
C
      ENDIF
C
      RETURN
      END
