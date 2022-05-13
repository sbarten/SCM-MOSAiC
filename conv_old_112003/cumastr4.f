      SUBROUTINE CUMASTR4
     *     (KCUCALL, KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,
     *     KLEVP1,   KLEVM1,   KROW,     ILAB,
     *     NSTEP ,   NSTART, NPRINT,    TWODT,   CONACC,
     *     PTEN,   PQEN,   PXEN,   PUEN,   PVEN,

C LG- adding the tracers

     *   KTRAC,    PDTIME,
     *   PXTEN,    PXTU,     PXTTE,

C LG- end

     *     PVERV,    PQSEN,    PTS,      PQHFL,
     *     PAPP1,    PAPHP1,   PGEO,     LDLAND,
     *     PTTE,     PQTE,     PVOM,     PVOL,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,  PXTEC,
     *     PRFLCK,
     *     LDCUM,    KTYPE,    KCBOT,    KCTOP,
     *     PTU,      PQU,      PLU,      PLUDE,
     *     PMFU,     PMFD,     PRAIN,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,   PSMELT,
CEVM---    CONDENSATION RATES
     *     TEMFCD ,   QEMFCD , XEMFCD 
CTCL---    ADDITIONAL ARRAYS FOR PROGNOSTIC CLOUD SCHEME
     *    ,KBOTSC,   LDSC)
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
cjhc*CALL PARAM
cjhc*CALL COMSDS
cjhc*CALL COMCTL
cjhc*CALL COMCON
cjhc*CALL COMCUMF
cjhc*CALL COMPH2
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
      INCLUDE 'comph2.h'
C
C
      REAL     PTEN(KLP2,KLEV),        PQEN(KLP2,KLEV),
     *     PXEN(KLP2,KLEV),
     *         PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PTTE(KLP2,KLEV),        PQTE(KLP2,KLEV),
     *         PVOM(K2LP2,KLEV),       PVOL(K2LP2,KLEV),
     *         PQSEN(KLP2,KLEV),       PGEO(KLP2,KLEV),
     *         PAPP1(KLP2,KLEV),       PAPHP1(KLP2,KLEVP1),
     *         PVERV(KLP2,KLEV),       PTS(KLP2),
     *         PQHFL(KLP2)
      REAL     PTU(KLP2,KLEV),         PQU(KLP2,KLEV),
     *         PLU(KLP2,KLEV),         PLUDE(KLP2,KLEV),
     *         PMFU(KLP2,KLEV),        PMFD(KLP2,KLEV),
     *         PAPRC(KLP2),            PAPRS(KLP2),
     *         PRSFC(KLP2),            PSSFC(KLP2),
     *         PRAIN(KLP2)
      INTEGER  KCBOT(KLP2),            KCTOP(KLP2),
     *         KTYPE(KLP2)
      LOGICAL  LDLAND(KLP2),           LDCUM(KLP2)
C
      REAL PXTEC(KLP2,KLEV)
C
      REAL PRFLCK(KLP2,KLEV)
C
      REAL TEMFCD(KLEV),QEMFCD(KLEV),XEMFCD(KLEV)
CTCL
      INTEGER KBOTSC(KLP2)
      LOGICAL LDSC(KLP2)
CTCL
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'

C LG- 

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

C LG- end

      REAL     ZTENH(JPHR,MLEV),       ZQENH(JPHR,MLEV),
     *      ZXENH(JPHR,MLEV),
     *         ZGEOH(JPHR,MLEV),       ZQSENH(JPHR,MLEV),
     *         ZTD(JPHR,MLEV),         ZQD(JPHR,MLEV),
     *         ZMFUS(JPHR,MLEV),       ZMFDS(JPHR,MLEV),
     *         ZMFUQ(JPHR,MLEV),       ZMFDQ(JPHR,MLEV),
     *         ZDMFUP(JPHR,MLEV),      ZDMFDP(JPHR,MLEV),
     *         ZMFUL(JPHR,MLEV),       ZRFL(JPHR),
     *         ZUU(JPHR,MLEV),         ZVU(JPHR,MLEV),
     *         ZUD(JPHR,MLEV),         ZVD(JPHR,MLEV)
      REAL     ZENTR(JPHR),            ZHCBASE(JPHR),
     *         ZMFUB(JPHR),            ZMFUB1(JPHR),
     *         ZDQPBL(JPHR),           ZDQCV(JPHR)
      REAL     ZSFL(JPHR),               ZDPMEL(JPHR,MLEV)
      REAL     ZCAPE(JPHR),             ZHEAT(JPHR)
      REAL     ZHMIN(JPHR)
      REAL     ZHHATT(JPHR,MLEV)
      INTEGER  IHMIN(JPHR)
      INTEGER  ILAB(JPHR,MLEV),        IDTOP(JPHR),
     *         ICTOP0(JPHR),
     *          ILWMIN(JPHR)
      LOGICAL  LODDRAF(JPHR)

C LG- adding the tracers

      REAL PXTEN(KLON,KLEV,KTRAC),     PXTTE(KLON,NLEVT,KTRAC),
     *     PXTU(KLON,KLEV,KTRAC),
     *     ZXTENH(KLON,KLEV,KTRAC),    ZXTD(KLON,KLEV,KTRAC),
     *     ZMFUXT(KLON,KLEV,KTRAC),    ZMFDXT(KLON,KLEV,KTRAC)

C LG- end

      LOGICAL  LLO1

C LG-

*/ ----- prepare NO-emission from flash ---------
*I CUMASTR.162
C -----------------------------------------------
C Calculate flash intensity (Price & Rind)
C -----------------------------------------------

C LG- This call is to routine which only contains common block,
C     which is now in common block parameter file 'parchem.h'
 
*CALL GJFLASH

C LG- LOLAND has the parameter name LDLAND in this routine
c      LOGICAL LOLAND(KLON)

C LG- end

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

C LG-

*I CUMASTR.210

      CALL RESETR (EMFLASH,KLON,-1.)
      
C LG- end      

C
C
C----------------------------------------------------------------------
C
C*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
C                  ---------------------------------------------------
C
  200 CONTINUE

      CALL CUINI4
     *     (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,
     *     KLEV,     KLEVP1,   KLEVM1,
     *      PTEN,   PQEN,   PQSEN,   PXEN,   PUEN,   PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   PXTEN,    ZXTENH,   PXTU,     ZXTD,     ZMFUXT,   ZMFDXT,

C LG- end

     *     PVERV,    PGEO,     PAPHP1,   ZGEOH,
     *      ZTENH,  ZQENH,  ZQSENH,  ZXENH,  ILWMIN,
     *     PTU,      PQU,      ZTD,      ZQD,
     *     ZUU,      ZVU,      ZUD,      ZVD,
     *     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     *     ZMFUQ,    ZMFDQ,    ZDMFUP,   ZDMFDP, PRFLCK,
     *     ZDPMEL,   PLU,      PLUDE,    ILAB)

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

      CALL CUBASE4
     *     (KCUCALL, KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,   KLEVM1,
     *     K2LP2,
     *     ZTENH,    ZQENH,    ZGEOH,    PAPHP1,

C LG- adding the tracers

     *   KTRAC,
     *   PXTU,

C LG- end

     *     PTU,      PQU,      PLU,
     *     PUEN,    PVEN,     ZUU,      ZVU,
     *     LDCUM,    KCBOT,    ILAB 
CTCL       
     *   , KBOTSC, LDSC)

C
C*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
C*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
C                  -----------------------------------------
C
       JK=1
       DO 310 JL=KIDIA,KFDIA
       ZDQCV(JL) =PQTE(JL,JK)*(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))
       ZDQPBL(JL)=0.0
       IDTOP(JL)=0
  310  CONTINUE
       DO 320 JK=2,KLEV
       DO 315 JL=KIDIA,KFDIA
       ZDQCV(JL)=ZDQCV(JL)+PQTE(JL,JK)*(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))
       IF(JK.GE.KCBOT(JL)) ZDQPBL(JL)=ZDQPBL(JL)+PQTE(JL,JK)
     1                               *(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))
  315 CONTINUE
  320 CONTINUE
      DO 330 JL=KIDIA,KFDIA
      KTYPE(JL)=ICVMGT(1,2,ZDQCV(JL).GT.MAX(0.,-1.1*PQHFL(JL)*G))
  330 CONTINUE
C
CEVM
C     IF (LMFQHFL = TRUE) INVOKE ALTERNATIVE CLOSURE
C
      IF (LMFQHFL) THEN
        DO 321 JL=KIDIA,KFDIA
	  ZDQPBL(JL)=-PQHFL(JL)*G
 321    CONTINUE
      ENDIF
CEVM
C
C*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
C*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
C*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
C                  ------------------------------------------
C
CDIR$ IVDEP
      DO 340 JL=KIDIA,KFDIA
      IKB=KCBOT(JL)
      ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-ZQENH(JL,IKB)
      ZDQMIN=MAX(0.01*ZQENH(JL,IKB),1.E-10)
      LLO1=ZDQPBL(JL).GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM(JL)
      ZMFUB(JL)=CVMGT(ZDQPBL(JL)/(G*MAX(ZQUMQE,ZDQMIN)),0.01,LLO1)
      ZMFMAX=(PAPHP1(JL,IKB)-PAPHP1(JL,IKB-1))*ZCONS2
      ZMFUB(JL)=MIN(ZMFUB(JL),ZMFMAX)
      IF(.NOT.LLO1) LDCUM(JL)=.FALSE.
      ZCAPE(JL)=0.
      ZHEAT(JL)=0.
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
      DO 410 JL=KIDIA,KFDIA
      IKB=KCBOT(JL)
      ZHCBASE(JL)=CPD*PTU(JL,IKB)+ZGEOH(JL,IKB)+ALV*PQU(JL,IKB)
      ICTOP0(JL)=KCBOT(JL)-1
  410 CONTINUE
      ZALVDCP=ALV/CPD
      ZQALV=1./ALV
      DO 430 JK=KLEVM1,3,-1
CDIR$ IVDEP
      DO 420 JL=KIDIA,KFDIA
      ZHSAT=CPD*ZTENH(JL,JK)+ZGEOH(JL,JK)+ALV*ZQSENH(JL,JK)
      ZGAM=C5LES*ZALVDCP*ZQSENH(JL,JK)/
     1     ((1.-VTMPC1*ZQSENH(JL,JK))*(ZTENH(JL,JK)-C4LES)**2)
      ZZZ=CPD*ZTENH(JL,JK)*0.608
      ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)*
     1            MAX(ZQSENH(JL,JK)-ZQENH(JL,JK),0.)
      ZHHATT(JL,JK)=ZHHAT
      IF(JK.LT.ICTOP0(JL).AND.ZHCBASE(JL).GT.ZHHAT) ICTOP0(JL)=JK
  420 CONTINUE
  430 CONTINUE
C
      DO JL=KIDIA,KFDIA
       JK=KCBOT(JL)
       ZHSAT=CPD*ZTENH(JL,JK)+ZGEOH(JL,JK)+ALV*ZQSENH(JL,JK)
       ZGAM=C5LES*ZALVDCP*ZQSENH(JL,JK)/
     1     ((1.-VTMPC1*ZQSENH(JL,JK))*(ZTENH(JL,JK)-C4LES)**2)
       ZZZ=CPD*ZTENH(JL,JK)*0.608
       ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)*
     1       MAX(ZQSENH(JL,JK)-ZQENH(JL,JK),0.)
       ZHHATT(JL,JK)=ZHHAT
      ENDDO
C
C
C                  FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
C                  -------------------------------------------
C
      DO JL=KIDIA,KFDIA
      LLO1=LDCUM(JL).AND.KTYPE(JL).EQ.1
      IF(LLO1) THEN
       IKB=KCBOT(JL)
       ZHMIN(JL)=0.
       IHMIN(JL)=IKB
      ENDIF
      ENDDO
C
      ZB=25.
      ZBI=1./(ZB*G)
      DO JK=KLEV,1,-1
      DO JL=KIDIA,KFDIA
      LLO1=LDCUM(JL).AND.KTYPE(JL).EQ.1.AND.IHMIN(JL).EQ.KCBOT(JL)
      IF(LLO1.AND.JK.LT.KCBOT(JL).AND.JK.GE.ICTOP0(JL)) THEN
         IKB=KCBOT(JL)
         ZRO=PAPHP1(JL,JK)/(RD*ZTENH(JL,JK))
         ZDZ=(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))/(G*ZRO)
         ZDHDZ=( CPD*(PTEN(JL,JK-1)-PTEN(JL,JK))+
     *           ALV*(PQEN(JL,JK-1)-PQEN(JL,JK))+
     *               (PGEO(JL,JK-1)-PGEO(JL,JK)) )*G/
     *               (PGEO(JL,JK-1)-PGEO(JL,JK))
         ZDEPTH=ZGEOH(JL,JK)-ZGEOH(JL,IKB)
         ZFAC=SQRT(1.+ZDEPTH*ZBI)
         ZHMIN(JL)=ZHMIN(JL) + ZDHDZ*ZFAC*ZDZ
         ZRH=-ALV*(ZQSENH(JL,JK)-ZQENH(JL,JK))*ZFAC
         IF(ZHMIN(JL).GT.ZRH) IHMIN(JL)=JK
      ENDIF
      ENDDO
      ENDDO
C
      DO JL=KIDIA,KFDIA
      IF(LDCUM(JL).AND.KTYPE(JL).EQ.1) THEN
       IF(IHMIN(JL).LT.ICTOP0(JL)) IHMIN(JL)=ICTOP0(JL)
      ENDIF
      ENDDO
C
C
C
C*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
C                  --------------------------------------------
C

      CALL CUASC4
     *     (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     NSTEP,    NSTART,   TWODT,
     *      ZTENH,  ZQENH,  ZXENH,   PUEN,   PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   ZXTENH,   PXTEN,    PXTU,     ZMFUXT,

C LG- end

     *     PTEN,     PQEN,     PQSEN,
     *     PGEO,     ZGEOH,    PAPP1,    PAPHP1,
     *     PQTE,     PVERV,    ILWMIN,   ZDQPBL,
     *     LDLAND,   LDCUM,    KTYPE,    ILAB,
     *     PTU,      PQU,      PLU,      ZUU,      ZVU,
     *     PMFU,     PMFD,     ZMFUB,    ZENTR,
     *     ZMFUS,    ZMFUQ,
     *     ZMFUL,    PLUDE,    ZDMFUP,
     *     IHMIN,    ZHHATT,   ZHCBASE,   ZQSENH,
     *     KCBOT,    KCTOP,    ICTOP0,   ICUM)

C
C
C*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
C              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
C              -----------------------------------------------------
C
C
CDIR$ IVDEP
      DO 440 JL=KIDIA,KFDIA
      ZPBMPT=PAPHP1(JL,KCBOT(JL))-PAPHP1(JL,KCTOP(JL))
      IF(LDCUM(JL).AND.KTYPE(JL).EQ.1.AND.ZPBMPT.LT.2.E4) KTYPE(JL)=2
C
      IF(KTYPE(JL).EQ.2) ZENTR(JL)=ENTRSCV
      ZRFL(JL)=ZDMFUP(JL,1)
  440 CONTINUE
      DO 460 JK=2,KLEV
      DO 450 JL=KIDIA,KFDIA
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

         CALL CUDLFS4
     *     (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVP1,
     *     ZTENH,    ZQENH,    PUEN,     PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   ZXTENH,   PXTU,     ZXTD,     ZMFDXT,

C LG- end

     *     ZGEOH,    PAPHP1,   LDLAND,   ZDQPBL,
     *     PTU,      PQU,      PLU,      ZUU,      ZVU,
     *     LDCUM,    KCBOT,    KCTOP,    ZMFUB,    ZRFL,
     *     ZTD,      ZQD,      ZUD,      ZVD,
     *     PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,
     *     IDTOP,    LODDRAF)
C
C*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
C                  -----------------------------------------------
C

         CALL CUDDRAF4
     *     (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVP1,
     *     ZTENH,    ZQENH,    PUEN,     PVEN,

C LG- adding hte tracers

     *   KTRAC,
     *   ZXTENH,   ZXTD,     ZMFDXT,

C LG- end

     *     ZGEOH,    PAPHP1,   ZRFL,     PLU,
     *     KTYPE,    KCBOT,
     *     ZTD,      ZQD,      ZUD,      ZVD,
     *     PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,
     *     IDTOP,    LODDRAF)
C
C
      END IF
C
C
C*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
C*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
C*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
C*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
C                  -------------------------------------------
C
C
      DO JL=KIDIA,KFDIA
        ZHEAT(JL)=0.
        ZCAPE(JL)=0.
        ZMFUB1(JL)=ZMFUB(JL)
      ENDDO
C
CEVM  971023
      IF (.NOT.LCLOS_VAPCVG) THEN
C
      DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
      LLO1=LDCUM(JL).AND.KTYPE(JL).EQ.1
      IF(LLO1.AND.JK.LE.KCBOT(JL).AND.JK.GT.KCTOP(JL)) THEN
         IKB=KCBOT(JL)
         ZRO=PAPHP1(JL,JK)/(RD*ZTENH(JL,JK))
         ZDZ=(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))/(G*ZRO)
         ZHEAT(JL)=ZHEAT(JL) +
     I    (  (PTEN(JL,JK-1)-PTEN(JL,JK) + G*ZDZ/CPD)/ZTENH(JL,JK)
     I        +  0.608*(PQEN(JL,JK-1)-PQEN(JL,JK))  ) *
     I       (G*(PMFU(JL,JK)+PMFD(JL,JK)))/ZRO
         ZCAPE(JL)=ZCAPE(JL) +
     I        (G*(PTU(JL,JK)-ZTENH(JL,JK))/ZTENH(JL,JK)
     I             +G*0.608*(PQU(JL,JK)-ZQENH(JL,JK))
     I             -G*PLU(JL,JK) ) * ZDZ
      ENDIF
      ENDDO
      ENDDO
C
CJHC      IF(NN.EQ.21) THEN
CJHC        ZTAU=2.*3600.
CJHC      ELSE
CJHC        ZTAU=1.*3600.
CJHC      ENDIF
CJHC from comph2:      ZTAU=CTAUMF
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL).AND.KTYPE(JL).EQ.1) THEN
           IKB=KCBOT(JL)
           ZMFUB1(JL)=(ZCAPE(JL)*ZMFUB(JL))/(ZHEAT(JL)*CTAUMF)
           ZMFUB1(JL)=MAX(ZMFUB1(JL),0.001)
           ZMFMAX=(PAPHP1(JL,IKB)-PAPHP1(JL,IKB-1))*ZCONS2
           ZMFUB1(JL)=MIN(ZMFUB1(JL),ZMFMAX)
        ENDIF
      ENDDO
CEVM  971023
C     ENDIF (.NOT.LCLOS_VAPCVG) THEN
      ENDIF
C
C
C
C
C
C*                 RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
C*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
C*                 FOR SHALLOW CONVECTION (KTYPE=2)
C                  --------------------------------------------
C
C
C
CDIR$ IVDEP
         DO 520 JL=KIDIA,KFDIA
CEVM  971023
C         IF(KTYPE(JL).NE.1) THEN
          IF(KTYPE(JL).NE.1 .OR. LCLOS_VAPCVG) THEN
            IKB=KCBOT(JL)
            LLO1=PMFD(JL,IKB).LT.0..AND.LODDRAF(JL)
            ZEPS=CVMGT(CMFDEPS,0.,LLO1)
            ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-
     1      ZEPS*ZQD(JL,IKB)-(1.-ZEPS)*ZQENH(JL,IKB)
            ZDQMIN=MAX(0.01*ZQENH(JL,IKB),1.E-10)
      ZMFMAX=(PAPHP1(JL,IKB)-PAPHP1(JL,IKB-1))*ZCONS2
            LLO1=ZDQPBL(JL).GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM(JL)
     1      .AND.ZMFUB(JL).LT.ZMFMAX
            ZMFUB1(JL)=CVMGT(ZDQPBL(JL)/(G*MAX(ZQUMQE,ZDQMIN)),
     1                 ZMFUB(JL),LLO1)
            ZMFUB1(JL)=CVMGT(ZMFUB1(JL),ZMFUB(JL),
     1                 (KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.2).AND.
     1                 ABS(ZMFUB1(JL)-ZMFUB(JL)).LT.0.2*ZMFUB(JL))
         END IF
  520    CONTINUE
         DO 540 JK=1,KLEV
         DO 530 JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            ZFAC=ZMFUB1(JL)/MAX(ZMFUB(JL),1.E-10)
            PMFD(JL,JK)=PMFD(JL,JK)*ZFAC
            ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZFAC
            ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZFAC
            ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZFAC
         END IF
  530    CONTINUE
C

C LG- adding the tracers

      DO 5304 JT=1,KTRAC
      DO 5302 JL=KIDIA,KFDIA
       IF(LDCUM(JL)) THEN
        ZFAC=ZMFUB1(JL)/MAX(ZMFUB(JL),1.E-10)
        ZMFDXT(JL,JK,JT)=ZMFDXT(JL,JK,JT)*ZFAC
       ENDIF
 5302 CONTINUE
 5304 CONTINUE

C LG- end

C
  540    CONTINUE
C
C
C*                 NEW VALUES OF CLOUD BASE MASS FLUX
C                  ----------------------------------
C

         DO 550 JL=KIDIA,KFDIA
         IF(LDCUM(JL)) ZMFUB(JL)=ZMFUB1(JL)
  550    CONTINUE
C
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

      CALL CUASC4
     *     (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     NSTEP,    NSTART,   TWODT,
     *      ZTENH,  ZQENH,  ZXENH,   PUEN,   PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   ZXTENH,   PXTEN,    PXTU,     ZMFUXT,

C LG- end

     *     PTEN,     PQEN,     PQSEN,
     *     PGEO,     ZGEOH,    PAPP1,    PAPHP1,
     *     PQTE,     PVERV,    ILWMIN,   ZDQPBL,
     *     LDLAND,   LDCUM,    KTYPE,    ILAB,
     *     PTU,      PQU,      PLU,      ZUU,      ZVU,
     *     PMFU,     PMFD,     ZMFUB,    ZENTR,
     *     ZMFUS,    ZMFUQ,
     *     ZMFUL,    PLUDE,    ZDMFUP,
     *     IHMIN,    ZHHATT,   ZHCBASE,   ZQSENH,
     *     KCBOT,    KCTOP,    ICTOP0,   ICUM)

C LG- lightning calculations

*I CUMASTR.455
c - - -
c ZNEMFAC in 5.14e26 molec/flash/(1x1degree) * XxX degree
c this value is tuned to yield appr. 4 Tg NO yr-1
c Array EMFLASH has unit: molec/s/gridbox
c
c 22/11/95: For T30, ZNEMFAC=multiplied by 2
c - - -
      ZNEMFAC=5.14E26*(360./REAL(NLONT30))**2.
      ZNEMFAC=ZNEMFAC*2.
      DO 604 JL=1,KLON

C LG- KTYPE indicates the convection type

      IF (KTYPE(JL).EQ.1.AND.KCTOP(JL).LT.KLEV) THEN
       JKT(JL)=KCTOP(JL)
       JKB(JL)=KCBOT(JL)
       TOPHEIGHT=MAX(0.,PGEO(JL,JKT(JL))/(G*1000.))
       GJFF=0.
       IF (LDLAND(JL)) THEN
        GJFF=3.44e-5*TOPHEIGHT**4.9
        JKB(JL)=KLEV
       ELSE
        GJFF=6.4e-4*TOPHEIGHT**1.73
       ENDIF
       EMFLASH(JL)=GJFF/60.*ZNEMFAC
      ENDIF
  604 CONTINUE
*/ ----- end flsh emission -------

C LG- end

C
C
C-----------------------------------------------------------------------
C
C*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
C                  ------------------------------------------
C
  700 CONTINUE

      CALL CUFLX4
     *     (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,
     *     NSTEP ,   NSTART,  TWODT,   
     *     PQEN,     PQSEN,    ZTENH,    ZQENH,

C LG- adding the tracers

     *   KTRAC,   PDTIME,
     *   PXTEN,    ZXTENH,   ZMFUXT,   ZMFDXT,

C LG-
C*I US310195.30
     * PXTU, PLU,PXTD,

C LG- end

     *     PAPHP1,   LDLAND,   ZGEOH,
     *     KCBOT,    KCTOP,    IDTOP,
     *     KTYPE,    LODDRAF,  LDCUM,
     *     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     *     ZMFUQ,    ZMFDQ,    ZMFUL,    PLUDE,
     *     ZDMFUP,   ZDMFDP,   ZRFL,     PRAIN,  PRFLCK, 
     *     PTEN,     ZSFL,     ZDPMEL,   ITOPM2)

C
C
C----------------------------------------------------------------------
C
C*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
C                  --------------------------------------------------
C

  800 CONTINUE

      CALL CUDTDQ4
     *     (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,
     *     NSTEP ,   NSTART,  TWODT,   CONACC,
     *     ITOPM2,   PAPHP1,   PGEO,     PTS,      LDLAND,
     *     LDCUM,    PTEN,     PTTE,     PQTE,

C LG- adding the tracers

     *     KTRAC,   PDTIME,
     *     PXTTE,    ZMFUXT,   ZMFDXT,

C LG-
*I HF240591.153
     *   PXTEN,

C LG- end

     *     PXTEC,
     *     ZMFUS,    ZMFDS,    ZMFUQ,    ZMFDQ,
     *     ZMFUL,    ZDMFUP,   ZDMFDP,   PLUDE,
     *     ZDPMEL,   PRAIN,    ZRFL,     ZSFL,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSMELT,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS)

C
C
C----------------------------------------------------------------------
C
C*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
C                  --------------------------------------------------
C
  900 CONTINUE
      IF(LMFDUDV) THEN
      CALL CUDUDV4
     *     (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVP1,
     *     ITOPM2,   KTYPE,    KCBOT,    PAPHP1,   LDCUM,
     *     PUEN,     PVEN,     PVOM,     PVOL,
     *     ZUU,      ZUD,      ZVU,      ZVD,
     *     PMFU,     PMFD,     PSDISS)
C
      END IF
C
 1000 CONTINUE
C
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
     *    (KLP2,     K2LP2,    KLON,
     *     KIDIA,    KFDIA,    KFDIA-KIDIA+1,
     *     KLEV,     KLEVP1,   KLEVM1,   NSTEP,   NPRINT,
     *     KCBOT,    KCTOP,    KTYPE,    ITOPM2, 
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     *     PVERV,    PGEO,     PAPP1,    PAPHP1,   ZGEOH,
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
