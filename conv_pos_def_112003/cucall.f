      SUBROUTINE CUCALL
     *    (KHOR,     KHOR2,    KLON,     KLEV,
     *     KSTART,   KSTOP,    KLEN,
     *     CONACC,
     *     KLEVP1,   KLEVM1,   KROW,     ILAB,
     *      PTM1,   PQM1,   PUM1,   PVM1,    PXM1,
     *      PTTE,   PQTE,   PVOM,   PVOL,    PXTE,
     *     PVERV,    PTS,      PQHFL,
     *       PACLCM,
     *     PAP,      PAPH,     PGEO,     LDLAND,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,
     *     KTYPE,    KTOPC,    KBASEC,   PARPRC,
     *     PTOPMAX,  PTOPMAXM,  PGEOSPM,  PRFLCK,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,
CHL ---    NEXT LINE REPLACES *COMSDS*   ---
     *     NSTEP ,  NSTART ,   NSTOP ,   NPRINT,
CHL ---    NEXT LINE REPLACES *COMCTL*   ---
     *     TWODT,
CHL ---    NEXT LINE REPLACES *COMRSW*   ---
     *     LRAD  ,   NRADFR,
CEVM---    CONDENSATION RATES 
     *     TEMFCD,   QEMFCD, XEMFCD
     *     )
C
C          *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
C                     *CUMASTR* (CUMULUS PARAMETERIZATION)
C                     *CUCCDIA* (CUMULUS CLOUDS FOR RADIATION)
C                     *STRATCU* (PBL_STRATOCUMULUS)
C
C           M.TIEDTKE      E.C.M.W.F.     12/1989
C
C**   PURPOSE.
C     --------
C
C          *CUCALL* - INTERFACE FOR *CUMASTR*,*CUCCDIA* AND *STRATCU*:
C                     PROVIDES INPUT FOR CUMASTR, CUCCDIA AND STRATCU.
C                     RECEIVES UPDATED TENDENCIES, PRECIPITATION
C                     AND CONVECTIVE CLOUD PARAMETERS FOR RADIATION.
C
C**   INTERFACE.
C     ----------
C
C          *CUCALL* IS CALLED FROM *PHYSC*
C
C     EXTERNALS.
C     ----------
C
C          CUMASTR
C          CUCCDIA
C          STRATCU
C
      INCLUDE 'comcon.h'
C
CDIR$ VFUNCTION EXP
C
      REAL     PTM1(KHOR2,KLEV),       PQM1(KHOR2,KLEV),
     *         PUM1(KHOR2,KLEV),       PVM1(KHOR2,KLEV),
     *         PTTE(KHOR,KLEV),        PQTE(KHOR,KLEV),
     *         PVOM(KHOR2,KLEV),       PVOL(KHOR2,KLEV),
     *         PVERV(KHOR,KLEV),       PGEO(KHOR,KLEV),
     *         PAP(KHOR,KLEV),         PAPH(KHOR,KLEVP1),
     *         PTS(KHOR),              PQHFL(KHOR)
      REAL     PAPRC(KHOR),            PAPRS(KHOR),
     *         PRSFC(KHOR),            PSSFC(KHOR),
     *         PARPRC(KHOR)
      INTEGER  KTOPC(KHOR),            KBASEC(KHOR),
     *         KTYPE(KHOR)
      REAL    PTOPMAX(KHOR),   PTOPMAXM(KHOR),
     *     PGEOSPM(KHOR), PRFLCK(KHOR,KLEV)
      INTEGER ILAB(KHOR,KLEV)
      LOGICAL  LDLAND(KHOR)
      REAL PACLCM(KHOR,KLEV)
      REAL PXM1(KHOR,KLEV), PXTE(KHOR,KLEV)
C
      REAL  TEMFCD(KLEV),QEMFCD(KLEV),XEMFCD(KLEV)
C
      LOGICAL  LRAD
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
     *       LOCUM  (JPHR)
      INTEGER
     *       ICBOT  (JPHR),
     *       ICTOP  (JPHR),
     *       ITOPEC2(JPHR)
C
C     -------------------------------------------------------
C     DECLARATION OF REAL ONE-DIMENSIONAL WORKING ARRAYS
C
      REAL
     *       ZRAIN  (JPHR),
C
C     -------------------------------------------------------
C     DECLARATION OF REAL TWO-DIMENSIONAL WORKING ARRAYS
C
     *       ZLU    (JPHR,MLEV),
     *       ZLUDE  (JPHR,MLEV),
     *       ZMFD   (JPHR,MLEV),
     *       ZMFU   (JPHR,MLEV),
     *       ZQP1   (JPHR,MLEV),
     *       ZQSAT  (JPHR,MLEV),
     *       ZQU    (JPHR,MLEV),
     *       ZTP1   (JPHR,MLEV),
     *       ZTU    (JPHR,MLEV),
     *       ZTVP1  (JPHR,MLEV),
C
     *       ZGEOP  (JPHR,MLEVP1)
C
C-----------------------------------------------------------------------
C
C*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
C*                 -----------------------------------
C
C
  100 CONTINUE
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      DO 120 JK=1,KLEV
      DO 110 JL=KSTART,KSTOP
      ZTP1(JL,JK)=PTM1(JL,JK)+PTTE(JL,JK)*ZTMST
      ZQP1(JL,JK)=PQM1(JL,JK)+PQTE(JL,JK)*ZTMST
         IF(ZTP1(JL,JK)-TMELT.GT.0.) THEN
            ZCVM3=C3LES
            ZCVM4=C4LES
         ELSE
            ZCVM3=C3IES
            ZCVM4=C4IES
         END IF
      ZQSAT(JL,JK)=C2ES*EXP(ZCVM3*(ZTP1(JL,JK)-TMELT)/
     1                        (ZTP1(JL,JK)-ZCVM4))/PAP(JL,JK)
      ZQSAT(JL,JK)=MIN(0.5,ZQSAT(JL,JK))
      ZQSAT(JL,JK)=ZQSAT(JL,JK)/(1.-VTMPC1*ZQSAT(JL,JK))
  110 CONTINUE
  120 CONTINUE
      DO 130 JL=KSTART,KSTOP
      ZRAIN(JL)=0.
      LOCUM(JL)=.FALSE.
  130 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
C*           -----------------------------------------------------------
C
C
  200 CONTINUE
      CALL CUMASTR
     *    (KHOR,     KHOR2,    KLON,     KLEV,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEVP1,   KLEVM1,   KROW,     ILAB,
     *     CONACC,
     *     ZTP1,     ZQP1,     PUM1,     PVM1,
     *     PVERV,    ZQSAT,    PTS,      PQHFL,
     *     PAP,      PAPH,     PGEO,     LDLAND,
     *     PTTE,     PQTE,     PVOM,     PVOL,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,
     *     PRFLCK,
     *     LOCUM,    KTYPE,    ICBOT,    ICTOP,
     *     ZTU,      ZQU,      ZLU,      ZLUDE,
     *     ZMFU,     ZMFD,     ZRAIN,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,
CHL --- lines added
     *     NSTEP ,  NSTART,   NSTOP,    NPRINT,   TWODT ,
CEVM---    CONDENSATION RATES 
     *     TEMFCD ,   QEMFCD , XEMFCD
     *     )
C
C
C ---------------------------------------------------------------------
C
C*    2.1 GEOPOTENTIAL OF CONVECTIVE CLOUD TOPS (KUO0)
C
C
  210 CONTINUE
C
      DO 2101 JL=KSTART,KSTOP
      ZGEOP(JL,KLEVP1)=PGEOSPM(JL)
 2101 CONTINUE
C
      DO 2113 JK=KLEV,1,-1
      DO 2112 JL=KSTART,KSTOP
      ZTVP1(JL,JK)=ZTP1(JL,JK)*(1.+VTMPC1*PQM1(JL,JK))
 2112 CONTINUE
 2113 CONTINUE
C
      DO 2103 JK=KLEV,2,-1
      DO 2102 JL=KSTART,KSTOP
      ZGEOP(JL,JK)=ZGEOP(JL,JK+1)
     >            +RD*ZTVP1(JL,JK)*LOG(PAPH(JL,JK+1)/PAPH(JL,JK))
 2102 CONTINUE
 2103 CONTINUE
C
      DO 2104 JL=KSTART,KSTOP
      ITOPEC2(JL)=KLEVP1
 2104 CONTINUE
C
      DO 2106 JK=1,KLEV
      DO 2105 JL=KSTART,KSTOP
      IF(ILAB(JL,JK).EQ.2 .AND. ITOPEC2(JL).EQ.KLEVP1) THEN
         ITOPEC2(JL)=JK
      END IF
 2105 CONTINUE
 2106 CONTINUE
C
      DO 2107 JL=KSTART,KSTOP
      IF(ITOPEC2(JL).EQ.1) THEN
         PTOPMAX(JL)=35000.
      ELSE IF(ITOPEC2(JL).NE.KLEVP1) THEN
         PTOPMAX(JL)=ZGEOP(JL,ITOPEC2(JL))/G
      ELSE
         PTOPMAX(JL)=-99.
      END IF
      PTOPMAX(JL)=MAX(PTOPMAX(JL),PTOPMAXM(JL))
 2107 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C*    3.0       CALL 'CUCCDIA' TO UPDATE CLOUD PARAMETERS FOR RADIATION
C               -------------------------------------------------------
C
  300 CONTINUE
      CALL CUCCDIA
     *    (KHOR,     KLON,     KLEV,
     *     KSTART,   KSTOP,    KLEN,
     *     LOCUM,    ZQU,      ZLU,      ZMFU,
     *     ZRAIN,    ICBOT,    ICTOP,
     *     PARPRC,   KTOPC,    KBASEC,
CHL
     *     NSTEP ,   LRAD ,    NRADFR)
C
C
C---------------------------------------------------------------------
C
C*    4.0          CALL 'STRATCU' FOR PARAMETERIZATION OF PBL-CLOUDS
C                  -------------------------------------------------
C
  400 CONTINUE
      CALL STRATCU
     *    (KHOR,     KLON,     KLEV,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEVP1,   KLEVM1,   KROW,
     *     PAP,      PAPH,     PGEO,
     *     ZTP1,     ZQP1,     ZQSAT,    KTYPE,
     *     PTTE,     PQTE,     PXTE,     PXM1,
CHL --- lines added
     *     NSTEP , NSTART,   NSTOP,    TWODT
     *   )
C
      RETURN
      END
