      SUBROUTINE CUBASMC4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVM1,   KK,
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   PXTEN,    PXTU,     PMFUXT,

C LG- end

     *     PVERV,    PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,
     *     PMFU,     PMFUB,    PENTR,    KCBOT,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   PMFUU,    PMFUV)
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE.
C          --------
C          THIS ROUTINE CALCULATES CLOUD BASE VALUES
C          FOR MIDLEVEL CONVECTION
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CUASC*.
C          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
C          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
C
C          METHOD.
C          --------
C          S. TIEDTKE (1989)
C
C          EXTERNALS
C          ---------
C          NONE
C
cjhc*CALL COMCON
cjhc*CALL COMCUMF
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTEN(KLP2,KLEV),        PQEN(KLP2,KLEV),
     *         PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PQSEN(KLP2,KLEV),       PVERV(KLP2,KLEV),
     *         PGEO(KLP2,KLEV),        PGEOH(KLP2,KLEV)
C
      REAL     PTU(KLP2,KLEV),         PQU(KLP2,KLEV),
     *         PUU(KLP2,KLEV),         PVU(KLP2,KLEV),
     *         PLU(KLP2,KLEV),         PMFU(KLP2,KLEV),
     *         PMFUB(KLP2),            PENTR(KLP2),
     *         PMFUS(KLP2,KLEV),       PMFUQ(KLP2,KLEV),
     *         PMFUL(KLP2,KLEV),       PDMFUP(KLP2,KLEV),
     *         PMFUU(KLP2),            PMFUV(KLP2)
      INTEGER  KTYPE(KLP2),            KCBOT(KLP2),
     *         KLAB(KLP2,KLEV)
      LOGICAL  LDCUM(KLP2)
C
      LOGICAL  LLO2

C LG- adding the tracers

      REAL PXTEN(KLON,KLEV,KTRAC),     PXTU(KLON,KLEV,KTRAC),
     *     PMFUXT(KLON,KLEV,KTRAC)
     
C LG- end     
    
      INCLUDE 'paramh.h'
      LOGICAL LLO3(JPHR)
C
C
C----------------------------------------------------------------------
C
C*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
C                  -------------------------------------------
C
  100 CONTINUE
      IF(LMFMID.AND.KK.LT.KLEVM1.AND.KK.GT.KLEV/2) THEN
CDIR$ IVDEP
         DO 150 JL=KIDIA,KFDIA
      LLO3(JL)=.FALSE.
         IF(.NOT.LDCUM(JL).AND.KLAB(JL,KK+1).EQ.0.
     1   .AND.PQEN(JL,KK).GT.0.90*PQSEN(JL,KK)) THEN
      LLO3(JL)=.TRUE.
            PTU(JL,KK+1)=(CPD*PTEN(JL,KK)+PGEO(JL,KK)-PGEOH(JL,KK+1))
     1                          *RCPD
            PQU(JL,KK+1)=PQEN(JL,KK)
            PLU(JL,KK+1)=0.
            ZZZMB=MAX(CMFCMIN,-PVERV(JL,KK)/G)
            ZZZMB=MIN(ZZZMB,CMFCMAX)
            PMFUB(JL)=ZZZMB
            PMFU(JL,KK+1)=PMFUB(JL)
            PMFUS(JL,KK+1)=PMFUB(JL)*(CPD*PTU(JL,KK+1)+PGEOH(JL,KK+1))
            PMFUQ(JL,KK+1)=PMFUB(JL)*PQU(JL,KK+1)
            PMFUL(JL,KK+1)=0.
            PDMFUP(JL,KK+1)=0.
            KCBOT(JL)=KK
            KLAB(JL,KK+1)=1
            KTYPE(JL)=3
            PENTR(JL)=ENTRMID
               IF(LMFDUDV) THEN
                  PUU(JL,KK+1)=PUEN(JL,KK)
                  PVU(JL,KK+1)=PVEN(JL,KK)
                  PMFUU(JL)=PMFUB(JL)*PUU(JL,KK+1)
                  PMFUV(JL)=PMFUB(JL)*PVU(JL,KK+1)
               END IF
         END IF
  150    CONTINUE
cjhcCDIR$ IVDEP

C LG- adding the tracers

      DO 1504 JT=1,KTRAC
      DO 1502 JL=1,KLON
      IF (LLO3(JL)) THEN
       PXTU(JL,KK+1,JT)=PXTEN(JL,KK,JT)
       PMFUXT(JL,KK+1,JT)=PMFUB(JL)*PXTU(JL,KK+1,JT)
      ENDIF
 1502 CONTINUE
 1504 CONTINUE

C LG- end

C
      END IF
C
      RETURN
      END
