      SUBROUTINE CUFLX
     *    (KHOR,     KLON,     KLEV,     KLEVP1,
     *     KSTART,   KSTOP,    KLEN,
     *     PQEN,     PQSEN,    PTENH,    PQENH,
     *     PAPH,     LDLAND,   PGEOH,
     *     KCBOT,    KCTOP,    KDTOP,
     *     KTYPE,    LDDRAF,   LDCUM,
     *     PMFU,     PMFD,     PMFUS,    PMFDS,
     *     PMFUQ,    PMFDQ,    PMFUL,    PLUDE,
     *     PDMFUP,  PDMFDP,  PRFL,    PRAIN,   PRFLCK,
     *     KTOPM2)
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE
C          -------
C
C          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
C          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C
C          EXTERNALS
C          ---------
C          NONE
C
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PQEN(KHOR,KLEV),        PQSEN(KHOR,KLEV),
     *         PTENH(KHOR,KLEV),       PQENH(KHOR,KLEV),
     *         PAPH(KHOR,KLEVP1),      PGEOH(KHOR,KLEV)
C
      REAL     PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV),
     *         PMFUS(KHOR,KLEV),       PMFDS(KHOR,KLEV),
     *         PMFUQ(KHOR,KLEV),       PMFDQ(KHOR,KLEV),
     *         PDMFUP(KHOR,KLEV),      PDMFDP(KHOR,KLEV),
     *         PMFUL(KHOR,KLEV),       PLUDE(KHOR,KLEV),
     *         PRFL(KHOR),             PRAIN(KHOR)
      REAL PRFLCK(KHOR,KLEV)
      INTEGER  KCBOT(KHOR),            KCTOP(KHOR),
     *         KDTOP(KHOR),            KTYPE(KHOR)
      LOGICAL  LDDRAF(KHOR),           LDLAND(KHOR),
     *         LDCUM(KHOR)
C
C
C
C*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
C                  ---------------------------------
C
  100 CONTINUE
      ITOP=KLEV
      DO 110 JL=KSTART,KSTOP
      ITOP=MIN(ITOP,KCTOP(JL))
      IF(.NOT.LDCUM(JL).OR.KDTOP(JL).LT.KCTOP(JL)) LDDRAF(JL)=.FALSE.
      IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
  110 CONTINUE
      KTOPM2=ITOP-2
C
C------------------------------------------------------------------
C								  |
C     DETERMINATION OF KTOPM2 IN CASE OF HIGH RESOLUTION MODE	  |
C								  |
C------------------------------------------------------------------
C
      JL   = 1
      IF (LMFHRES.AND.(KTYPE(JL).LT.3)) THEN
C
C
C-----------------------------------------------
C     KTOPM2: LOWEST LEVEL WITH ZERO MASS FLUX  |
C						|
C     ORIGINAL TIEDTKE: KTOPM2 = ITOP-2		|
C------------------------------------------------
C
C
      ICTOP = KLEV
      DO 102 JK = 2,KCTOP(JL)
       IF ((PMFU(JL,JK).EQ.0.0).AND.(PMFU(JL,JK+1).GT.0.0)) THEN
           ICTOP = JK+1
       ENDIF
  102 CONTINUE
      KTOPM2 = ICTOP-1
C
C
      ENDIF
C
C--------------- END OF LMFHIRES  -------------
C
C
      DO 120 JK=KTOPM2,KLEV
c$dir no_recurrence, force_vector, force_parallel_ext
CDIR$ IVDEP
      DO 115 JL=KSTART,KSTOP
      IF(LDCUM(JL).AND.JK.GE.(KTOPM2+1)) THEN
         PMFUS(JL,JK)=PMFUS(JL,JK)-PMFU(JL,JK)*
     1                (CPD*PTENH(JL,JK)+PGEOH(JL,JK))
         PMFUQ(JL,JK)=PMFUQ(JL,JK)-PMFU(JL,JK)*PQENH(JL,JK)
         ZDP=CVMGT(3.E4,1.5E4,LDLAND(JL))
         IF(PAPH(JL,KCBOT(JL))-PAPH(JL,KCTOP(JL)).GE.ZDP.AND.
     1      PQEN(JL,JK-1).GT.0.8*PQSEN(JL,JK-1))
     1      PDMFUP(JL,JK-1)=PDMFUP(JL,JK-1)+PLUDE(JL,JK-1)
         IF(LDDRAF(JL).AND.JK.GE.KDTOP(JL)) THEN
            PMFDS(JL,JK)=PMFDS(JL,JK)-PMFD(JL,JK)*
     1                   (CPD*PTENH(JL,JK)+PGEOH(JL,JK))
            PMFDQ(JL,JK)=PMFDQ(JL,JK)-PMFD(JL,JK)*PQENH(JL,JK)
         ELSE
            PMFD(JL,JK)     = 0.
            PMFDS(JL,JK)    = 0.
            PMFDQ(JL,JK)    = 0.
            PDMFDP(JL,JK-1) = 0.
         END IF
      ELSE
         PMFU(JL,JK)     = 0.
         PMFD(JL,JK)     = 0.
         PMFUS(JL,JK)    = 0.
         PMFDS(JL,JK)    = 0.
         PMFUQ(JL,JK)    = 0.
         PMFDQ(JL,JK)    = 0.
         PMFUL(JL,JK)    = 0.
         PDMFUP(JL,JK-1) = 0.
         PDMFDP(JL,JK-1) = 0.
         PLUDE(JL,JK-1)  = 0.
      END IF
  115 CONTINUE
  120 CONTINUE
      DO 130 JK=KTOPM2,KLEV
c$dir no_recurrence, force_vector, force_parallel_ext
CDIR$ IVDEP
      DO 125 JL=KSTART,KSTOP
      IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
         IKB=KCBOT(JL)
         ZZP=((PAPH(JL,KLEVP1)-PAPH(JL,JK))/
     1        (PAPH(JL,KLEVP1)-PAPH(JL,IKB)))
         ZZP=CVMGT(ZZP**2,ZZP,KTYPE(JL).EQ.3)
         PMFUS(JL,JK)=PMFUS(JL,IKB)*ZZP
         PMFUQ(JL,JK)=PMFUQ(JL,IKB)*ZZP
         PMFUL(JL,JK)=PMFUL(JL,IKB)*ZZP
      END IF
  125 CONTINUE
  130 CONTINUE
C
      DO 140 JL=KSTART,KSTOP
      PRFL(JL)=0.
      PRAIN(JL)=0.
  140 CONTINUE
      DO 160 JK=KTOPM2,KLEV
      DO 150 JL=KSTART,KSTOP
      IF(LDCUM(JL)) THEN
         PRFL(JL)=PRFL(JL)+PDMFUP(JL,JK)+PDMFDP(JL,JK)
      PRFLCK(JL,JK)=MAX(0.,PRFL(JL))
         PRAIN(JL)=PRAIN(JL)+PDMFUP(JL,JK)
      END IF
  150 CONTINUE
  160 CONTINUE
      RETURN
      END
