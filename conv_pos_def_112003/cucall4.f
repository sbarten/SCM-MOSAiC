      SUBROUTINE CUCALL4
     *    (KCUCALL,  KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,
     *     KLEVP1,   KLEVM1,   KROW,     ILAB,
     *     NSTEP ,   NSTART, NPRINT,    TWODT, CONACC,

C LG- adding the tracers

     *   KTRAC,   PDTIME,
     *   PXTM1,    PXTTE,

C LG- end

     *      PTM1,   PQM1,   PUM1,   PVM1,    PXM1,
     *      PTTE,   PQTE,   PVOM,   PVOL,    PXTE,
     *     PVERV,    PTS,      PQHFL,
     *     PACLCM,    PXTEC,
     *     PAPP1,    PAPHP1,   PGEO,     LDLAND,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,
     *     PRFLCK,
     *     KTYPE,
     *      PTOPMAX,     PTOPMAXM,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,   PSMELT,
CEVM---    CONDENSATION RATES
     *     TEMFCD,   QEMFCD, XEMFCD 
CTCL---    VARIABLES FROM MASS FLUX TO ECMWF CLOUD SCHEME
     *   , PLU, PLUDE, PMFU, PMFD
     *   , KCTOP , KCBOT , KBOTSC
     *   , LDCUM , LDSC
CTCL---    
     *     )
C
C          *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
C                     *CUMASTR* (CUMULUS PARAMETERIZATION)
C
C           M.TIEDTKE      E.C.M.W.F.     12/1989
C           E.VAN MEIJGAARD   KNMI        01/1997
C
C**   PURPOSE.
C     --------
C
C          *CUCALL* - INTERFACE FOR *CUMASTR*:
C                     PROVIDES INPUT FOR CUMASTR
C                     RECEIVES UPDATED TENDENCIES, PRECIPITATION.
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
C
CJHC*CALL PARAM
CJHC*CALL COMSDS
CJHC*CALL COMCTL
CJHC*CALL COMCON
C
CJHC*CALL YOTLUC
      INCLUDE 'comcon.h'
      INCLUDE 'yotluc.h'
C
CDIR$ VFUNCTION EXP
C
      REAL     PTM1(K2LP2,KLEV),       PQM1(K2LP2,KLEV),
     *         PUM1(K2LP2,KLEV),       PVM1(K2LP2,KLEV),
     *         PTTE(KLP2,KLEV),        PQTE(KLP2,KLEV),
     *         PVOM(K2LP2,KLEV),       PVOL(K2LP2,KLEV),
     *         PVERV(KLP2,KLEV),       PGEO(KLP2,KLEV),
     *         PAPP1(KLP2,KLEV),       PAPHP1(KLP2,KLEVP1),
     *         PTS(KLP2),              PQHFL(KLP2)
      REAL     PAPRC(KLP2),            PAPRS(KLP2),
     *         PRSFC(KLP2),            PSSFC(KLP2)
      REAL     PRFLCK(KLP2,KLEV)
CEVM  INTEGER  KTOPC(KLP2),            KBASEC(KLP2),
      INTEGER  KTYPE(KLP2)
      REAL    PTOPMAX(KLP2),   PTOPMAXM(KLP2)
      INTEGER ILAB(KLP2,KLEV)
      LOGICAL  LDLAND(KLP2)
      REAL PACLCM(KLP2,KLEV)
      REAL PXTEC(KLP2,KLEV)
      REAL PXM1(KLP2,KLEV), PXTE(KLP2,KLEV)
C
      REAL  TEMFCD(KLEV),QEMFCD(KLEV),XEMFCD(KLEV)
CTCL---    VARIABLES FROM MASS FLUX TO ECMWF CLOUD SCHEME
      REAL
     *     PLU(KLP2,KLEV), PLUDE(KLP2,KLEV)
     *   , PMFU(KLP2,KLEV), PMFD(KLP2,KLEV)
      INTEGER
     *     KCTOP(KLP2) , KCBOT(KLP2) , KBOTSC(KLP2)
      LOGICAL
     *     LDCUM(KLP2) , LDSC(KLP2)
CTCL---    
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'

C LG- 

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'
      
C LG- end

      REAL     ZTP1(JPHR,MLEV),        ZQP1(JPHR,MLEV),
     *         ZXP1(JPHR,MLEV),
     *         ZUP1(JPHR,MLEV),        ZVP1(JPHR,MLEV),
     *         ZTU(JPHR,MLEV),         ZQU(JPHR,MLEV),
     *         ZQSAT(JPHR,MLEV),       ZRAIN(JPHR)
      INTEGER ITOPEC2(JPHR)

C LG- adding the tracers

      REAL   ZXTP1(KLON,KLEV,KTRAC),   ZXTU(KLON,KLEV,KTRAC),
     *       PXTM1(KLON,NLEVT,KTRAC),   PXTTE(KLON,NLEVT,KTRAC)

C LG- end

C
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
      DO 110 JL=KIDIA,KFDIA
      ZTP1(JL,JK)=PTM1(JL,JK)+PTTE(JL,JK)*ZTMST
      ZQP1(JL,JK)=PQM1(JL,JK)+PQTE(JL,JK)*ZTMST
      ZXP1(JL,JK)=PXM1(JL,JK)+PXTE(JL,JK)*ZTMST
      ZUP1(JL,JK)=PUM1(JL,JK)+PVOM(JL,JK)*ZTMST
      ZVP1(JL,JK)=PVM1(JL,JK)+PVOL(JL,JK)*ZTMST
      IT=ZTP1(JL,JK)*1000.
      ZQSAT(JL,JK)=TLUCUA(IT)/PAPP1(JL,JK)
      ZQSAT(JL,JK)=MIN(0.5,ZQSAT(JL,JK))
      ZQSAT(JL,JK)=ZQSAT(JL,JK)/(1.-VTMPC1*ZQSAT(JL,JK))
  110 CONTINUE
C

C LG- adding the tracers

      DO 1104 JT=1,KTRAC
      DO 1102 JL=KIDIA,KFDIA
       ZXTP1(JL,JK,JT)=PXTM1(JL,JK,JT)+PXTTE(JL,JK,JT)*ZTMST
 1102 CONTINUE
 1104 CONTINUE

C LG- end

C
  120 CONTINUE
      DO 130 JL=KIDIA,KFDIA
      ZRAIN(JL)=0.
      LDCUM(JL)=.FALSE.
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

      CALL CUMASTR4
     *    (KCUCALL,  KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,
     *     KLEVP1,   KLEVM1,   KROW,     ILAB,
     *     NSTEP ,   NSTART, NPRINT,    TWODT, CONACC,
     *     ZTP1,   ZQP1,   ZXP1,   ZUP1,   ZVP1,

C LG- adding the tracers

     *   KTRAC,    PDTIME,
     *   ZXTP1,    ZXTU,      PXTTE,

C LG- end

     *     PVERV,    ZQSAT,    PTS,      PQHFL,
     *     PAPP1,    PAPHP1,   PGEO,     LDLAND,
     *     PTTE,     PQTE,     PVOM,     PVOL,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS,  PXTEC,
     *     PRFLCK,
     *     LDCUM,    KTYPE,    KCBOT,    KCTOP,
     *     ZTU,      ZQU,      PLU,      PLUDE,
     *     PMFU,     PMFD,     ZRAIN,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,   PSMELT,
CEVM---    CONDENSATION RATES
     *     TEMFCD,   QEMFCD, XEMFCD 
CTCL---    ADDITIONAL ARRAYS FOR PROGNOSTIC CLOUD SCHEME
     *    ,KBOTSC,   LDSC)

C
C
C ------------------------------------------------------------------
C
C*     3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
C             -------- -------- -- ---------- ----- -----
C

  300 CONTINUE
C
      ILEVMIN=KLEV-4
C
      DO 301 JL=KIDIA,KFDIA
      ITOPEC2(JL)=KLEVP1
  301 CONTINUE
C
      DO 303 JK=1,ILEVMIN
      DO 302 JL=KIDIA,KFDIA
      IF(ILAB(JL,JK).EQ.2 .AND. ITOPEC2(JL).EQ.KLEVP1) THEN
         ITOPEC2(JL)=JK
      END IF
  302 CONTINUE
  303 CONTINUE
C
      DO 304 JL=KIDIA,KFDIA
      IF(ITOPEC2(JL).EQ.1) THEN
         PTOPMAX(JL)=PAPP1(JL,1)
      ELSE IF(ITOPEC2(JL).NE.KLEVP1) THEN
         PTOPMAX(JL)=PAPHP1(JL,ITOPEC2(JL))
      ELSE
         PTOPMAX(JL)=99999.
      END IF
      PTOPMAX(JL)=MIN(PTOPMAX(JL),PTOPMAXM(JL))
  304 CONTINUE
C
C
C---------------------------------------------------------------------
C
      RETURN
      END
