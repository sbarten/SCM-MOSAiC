      SUBROUTINE CUBASMC
     *    (KHOR,     KHOR2,    KLON,     KLEV,     KLEVM1,   KK,
     *     KSTART,   KSTOP,    KLEN,
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
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
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTEN(KHOR,KLEV),        PQEN(KHOR,KLEV),
     *         PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PQSEN(KHOR,KLEV),       PVERV(KHOR,KLEV),
     *         PGEO(KHOR,KLEV),        PGEOH(KHOR,KLEV)
C
      REAL     PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PUU(KHOR,KLEV),         PVU(KHOR,KLEV),
     *         PLU(KHOR,KLEV),         PMFU(KHOR,KLEV),
     *         PMFUB(KHOR),            PENTR(KHOR),
     *         PMFUS(KHOR,KLEV),       PMFUQ(KHOR,KLEV),
     *         PMFUL(KHOR,KLEV),       PDMFUP(KHOR,KLEV),
     *         PMFUU(KHOR),            PMFUV(KHOR)
      INTEGER  KTYPE(KHOR),            KCBOT(KHOR),
     *         KLAB(KHOR,KLEV)
      LOGICAL  LDCUM(KHOR)
C
      LOGICAL  LLO2
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
c$dir no_recurrence, force_vector, force_parallel_ext
         DO 150 JL=KSTART,KSTOP
         IF(.NOT.LDCUM(JL).AND.KLAB(JL,KK+1).EQ.0.
     1   .AND.PQEN(JL,KK).GT.0.90*PQSEN(JL,KK)) THEN
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
      END IF
C
      RETURN
      END
