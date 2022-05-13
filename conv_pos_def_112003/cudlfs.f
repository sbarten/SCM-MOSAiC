      SUBROUTINE CUDLFS
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     PTENH,    PQENH,    PUEN,     PVEN,
     *     PGEOH,    PAPH,     LDLAND,   PDQPBL,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     LDCUM,    KCBOT,    KCTOP,    PMFUB,    PRFL,
     *     PTD,      PQD,      PUD,      PVD,
     *     PMFD,     PMFDS,    PMFDQ,    PDMFDP,
     *     KDTOP,    LDDRAF
     *   )
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
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KHOR,KLEV),       PQENH(KHOR,KLEV),
     *         PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PGEOH(KHOR,KLEV),       PAPH(KHOR,KLEVP1),
     *         PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PUU(KHOR,KLEV),         PVU(KHOR,KLEV),
     *         PLU(KHOR,KLEV),         PDQPBL(KHOR),
     *         PMFUB(KHOR),            PRFL(KHOR)
C
      REAL     PTD(KHOR,KLEV),         PQD(KHOR,KLEV),
     *         PUD(KHOR,KLEV),         PVD(KHOR,KLEV),
     *         PMFD(KHOR,KLEV),        PMFDS(KHOR,KLEV),
     *         PMFDQ(KHOR,KLEV),       PDMFDP(KHOR,KLEV)
      INTEGER  KCBOT(KHOR),            KCTOP(KHOR),
     *         KDTOP(KHOR)
      LOGICAL  LDLAND(KHOR),           LDCUM(KHOR),
     *         LDDRAF(KHOR)
C
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      LOGICAL
     *        LLO2   (JPHR)
      REAL
     *        ZCOND  (JPHR),
     *        ZQENWB (JPHR,MLEV),
     *        ZTENWB (JPHR,MLEV)
C
C
C----------------------------------------------------------------------
C
C     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
C                  ---------------------------------
C
  100 CONTINUE
      DO 110 JL=KSTART,KSTOP
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
C                  DOING A SCAN FROM BASE TO TOP OF CUMULUS CLOUDS
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
      DO 212 JL=KSTART,KSTOP
      ZTENWB(JL,JK)=PTENH(JL,JK)
      ZQENWB(JL,JK)=PQENH(JL,JK)
      LLO2(JL)=LDCUM(JL).AND.PRFL(JL).GT.0..AND..NOT.LDDRAF(JL).AND.
     1         (JK.LT.KCBOT(JL).AND.JK.GT.KCTOP(JL))
      IS=IS+ICVMGT(1,0,LLO2(JL))
  212 CONTINUE
      IF(IS.EQ.0) GO TO 290
C
      IK=JK
      ICALL=2
      CALL CUADJTQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *     KSTART,   KSTOP,    KLEN,
     *     PAPH,     ZTENWB,   ZQENWB,   LLO2,     ICALL
     *  )
C
C
C     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
C                  AND CHECK FOR NEGATIVE BUOYANCY.
C                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
C                  ----------------------------------------
C
  220 CONTINUE
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 222 JL=KSTART,KSTOP
      IF(LLO2(JL)) THEN
         ZTTEST=0.5*(PTU(JL,JK)+ZTENWB(JL,JK))
         ZQTEST=0.5*(PQU(JL,JK)+ZQENWB(JL,JK))
         ZBUO=ZTTEST*(1.+VTMPC1*ZQTEST)-
     1        PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))
         ZCOND(JL)=PQENH(JL,JK)-ZQENWB(JL,JK)
         ZMFTOP=-CMFDEPS*PMFUB(JL)
         IF(ZBUO.LT.0..AND.PRFL(JL).GT.10.*ZMFTOP*ZCOND(JL)) THEN
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
         IF(LMFDUDV) THEN
            DO 224 JL=KSTART,KSTOP
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
