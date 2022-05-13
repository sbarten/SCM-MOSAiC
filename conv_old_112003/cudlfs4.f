      SUBROUTINE CUDLFS4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVP1,
     *     PTENH,    PQENH,    PUEN,     PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   PXTENH,   PXTU,     PXTD,     PMFDXT,

C LG- end

     *     PGEOH,    PAPH,     LDLAND,   PDQPBL,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     LDCUM,    KCBOT,    KCTOP,    PMFUB,    PRFL,
     *     PTD,      PQD,      PUD,      PVD,
     *     PMFD,     PMFDS,    PMFDQ,    PDMFDP,
     *     KDTOP,    LDDRAF)
C
C          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
C          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
C
C          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
C          FOR MASSFLUX CUMULUS PARAMETERIZATION
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
C          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
C          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
C          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
C
C          METHOD.
C          --------
C          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
C          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
C
C          EXTERNALS
C          ---------
C          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
C
cjhc*CALL COMCON
cjhc*CALL COMCUMF
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KLP2,KLEV),       PQENH(KLP2,KLEV),
     *         PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PGEOH(KLP2,KLEV),       PAPH(KLP2,KLEVP1),
     *         PTU(KLP2,KLEV),         PQU(KLP2,KLEV),
     *         PUU(KLP2,KLEV),         PVU(KLP2,KLEV),
     *         PLU(KLP2,KLEV),         PDQPBL(KLP2),
     *         PMFUB(KLP2),            PRFL(KLP2)
C
      REAL     PTD(KLP2,KLEV),         PQD(KLP2,KLEV),
     *         PUD(KLP2,KLEV),         PVD(KLP2,KLEV),
     *         PMFD(KLP2,KLEV),        PMFDS(KLP2,KLEV),
     *         PMFDQ(KLP2,KLEV),       PDMFDP(KLP2,KLEV)
      INTEGER  KCBOT(KLP2),            KCTOP(KLP2),
     *         KDTOP(KLP2)
      LOGICAL  LDLAND(KLP2),           LDCUM(KLP2),
     *         LDDRAF(KLP2)
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      REAL     ZTENWB(JPHR,MLEV),      ZQENWB(JPHR,MLEV),
     *         ZCOND(JPHR)
      REAL     ZPH(JPHR)
      LOGICAL  LLO2(JPHR)

C LG- adding the tracers

      REAL   PXTENH(KLON,KLEV,KTRAC),  PXTU(KLON,KLEV,KTRAC),
     *       PXTD(KLON,KLEV,KTRAC),    PMFDXT(KLON,KLEV,KTRAC)

C LG- end

      LOGICAL   LLO3(JPHR)
C
C
C----------------------------------------------------------------------
C
C     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
C                  ---------------------------------
C
  100 CONTINUE
      DO 110 JL=KIDIA,KFDIA
      LDDRAF(JL)=.FALSE.
      KDTOP(JL)=KLEVP1
  110 CONTINUE
C
      IF(.NOT.LMFDD) GO TO 300
C
C
C----------------------------------------------------------------------
C
C     2.           DETERMINE LEVEL OF FREE SINKING BY
C                  DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
C
C                  FOR EVERY POINT AND PROCEED AS FOLLOWS:
C
C                    (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
C                    (2) DO MIXING WITH CUMULUS CLOUD AIR
C                    (3) CHECK FOR NEGATIVE BUOYANCY
C
C                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
C                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
C                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
C                  EVAPORATION OF RAIN AND CLOUD WATER)
C                  ----------------------------------------------------
C
  200 CONTINUE
C
      KE=KLEV-3
      DO 290 JK=3,KE
C
C
C     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
C                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
C                  -------------------------------------------
C
  210 CONTINUE
      IS=0
      DO 212 JL=KIDIA,KFDIA
      ZTENWB(JL,JK)=PTENH(JL,JK)
      ZQENWB(JL,JK)=PQENH(JL,JK)
      ZPH(JL)=PAPH(JL,JK)
      LLO2(JL)=LDCUM(JL).AND.PRFL(JL).GT.0..AND..NOT.LDDRAF(JL).AND.
     1         (JK.LT.KCBOT(JL).AND.JK.GT.KCTOP(JL))
      IS=IS+ICVMGT(1,0,LLO2(JL))
  212 CONTINUE
      IF(IS.EQ.0) GO TO 290
C
      IK=JK
      ICALL=2
      CALL CUADJTQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     IK,
     *     ZPH,      ZTENWB,   ZQENWB,   LLO2,     ICALL)
C
C
C     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
C                  AND CHECK FOR NEGATIVE BUOYANCY.
C                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
C                  ----------------------------------------
C
  220 CONTINUE
CDIR$ IVDEP
      DO 222 JL=KIDIA,KFDIA
      LLO3(JL)=.FALSE.
      IF(LLO2(JL)) THEN
         ZTTEST=0.5*(PTU(JL,JK)+ZTENWB(JL,JK))
         ZQTEST=0.5*(PQU(JL,JK)+ZQENWB(JL,JK))
         ZBUO=ZTTEST*(1.+VTMPC1*ZQTEST)-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))
         ZCOND(JL)=PQENH(JL,JK)-ZQENWB(JL,JK)
         ZMFTOP=-CMFDEPS*PMFUB(JL)
         IF(ZBUO.LT.0..AND.PRFL(JL).GT.10.*ZMFTOP*ZCOND(JL)) THEN
      LLO3(JL)=.TRUE.
            KDTOP(JL)=JK
            LDDRAF(JL)=.TRUE.
            PTD(JL,JK)=ZTTEST
            PQD(JL,JK)=ZQTEST
            PMFD(JL,JK)=ZMFTOP
            PMFDS(JL,JK)=PMFD(JL,JK)*(CPD*PTD(JL,JK)+PGEOH(JL,JK))
            PMFDQ(JL,JK)=PMFD(JL,JK)*PQD(JL,JK)
            PDMFDP(JL,JK-1)=-0.5*PMFD(JL,JK)*ZCOND(JL)
            PRFL(JL)=PRFL(JL)+PDMFDP(JL,JK-1)
         END IF
      END IF
  222 CONTINUE
C

C LG- adding the tracers

      DO 2224 JT=1,KTRAC
      DO 2222 JL=KIDIA,KFDIA
      IF(LLO3(JL)) THEN
       PXTD(JL,JK,JT)=0.5*(PXTU(JL,JK,JT)+PXTENH(JL,JK,JT))
       PMFDXT(JL,JK,JT)=PMFD(JL,JK)*PXTD(JL,JK,JT)
      ENDIF
 2222 CONTINUE
 2224 CONTINUE

C LG- end

C
         IF(LMFDUDV) THEN
            DO 224 JL=KIDIA,KFDIA
            IF(PMFD(JL,JK).LT.0.) THEN
               PUD(JL,JK)=0.5*(PUU(JL,JK)+PUEN(JL,JK-1))
               PVD(JL,JK)=0.5*(PVU(JL,JK)+PVEN(JL,JK-1))
            END IF
  224       CONTINUE
         END IF
C
  290 CONTINUE
C
 300  CONTINUE
      RETURN
      END
