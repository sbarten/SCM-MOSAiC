      SUBROUTINE CUCCDIA
     *    (KHOR,     KLON,     KLEV,
     *     KSTART,   KSTOP,    KLEN,
     *     LDCUM,    PQU,      PLU,      PMFU,
     *     PRAIN,    KCBOT,    KCTOP,
     *     PARPRC,   KTOPC,    KBASEC,
CHL
     *     NSTEP ,   LRAD ,    NRADFR)
C
C
C**** *CUCCDIA*- UPDATES PRECIPITAION, CLOUD BASE AND CLOUD TOP
C                FOR DIAGNOSTIC SCHEME FOR CONVECTIVE CLOUDS
C
C          M.TIEDTKE         E.C.M.W.F.    12/89
C
C**   INTERFACE.
C     ----------
C
C          *CUCCDIA* IS CALLED FROM *CUCALL*
C
      INCLUDE 'comcumf.h'
C
C
      REAL     PQU(KHOR,KLEV),         PLU(KHOR,KLEV),
     *         PMFU(KHOR,KLEV)
      REAL     PRAIN(KHOR)
      INTEGER  KCBOT(KHOR),            KCTOP(KHOR)
      LOGICAL  LDCUM(KHOR)
C
      REAL     PARPRC(KHOR)
      INTEGER  KTOPC(KHOR),            KBASEC(KHOR)
      LOGICAL  LRAD
C
C
C
C---------------------------------------------------------------------
C
C*    1.0          STORE CLOUD PARAMETERS FOR RADIATION CALCULATION
C                  -----------------------------------------------
C
  100 CONTINUE
C
      IF(LMFPEN.AND.LRAD.AND.MOD(NSTEP+1,NRADFR).NE.0) THEN
         ZNORMR=1./FLOAT(NRADFR-1)
CDIR$ IVDEP
c$dir no_recurrence, force_vector, force_parallel_ext
         DO 120 JL=KSTART,KSTOP
         IF(LDCUM(JL)) THEN
            KBASEC(JL)=MAX(KBASEC(JL),KCBOT(JL))
            KTOPC(JL)=MIN(KTOPC(JL),KCTOP(JL))
            IF(KTOPC(JL).EQ.1) KTOPC(JL)=KCTOP(JL)
            IKB=KCBOT(JL)
            IKT=KCTOP(JL)
            ZDMFQ=PMFU(JL,IKB)*MAX(PQU(JL,IKB)+PLU(JL,IKB)-PQU(JL,IKT),
     1                           PQU(JL,IKB)*0.05)
            PARPRC(JL)=PARPRC(JL)+MAX(ZDMFQ,PRAIN(JL))*ZNORMR
         END IF
  120    CONTINUE
      END IF
C
      RETURN
      END
