      SUBROUTINE CUINI
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     *     PVERV,    PGEO,     PAPH,     PGEOH,
     *     PTENH,    PQENH,    PQSENH,   KLWMIN,
     *     PTU,      PQU,      PTD,      PQD,
     *     PUU,      PVU,      PUD,      PVD,
     *     PMFU,     PMFD,     PMFUS,    PMFDS,
     *     PMFUQ,   PMFDQ,   PDMFUP,   PDMFDP,   PRFLCK,
     *     PLU,      PLUDE,    KLAB
     *   )
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE
C          -------
C
C          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
C          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
C          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
C          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C
C          METHOD.
C          --------
C          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
C
C          EXTERNALS
C          ---------
C          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
C
      INCLUDE 'comcon.h'
C
      REAL     PTEN(KHOR,KLEV),        PQEN(KHOR,KLEV),
     *         PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PQSEN(KHOR,KLEV),       PVERV(KHOR,KLEV),
     *         PGEO(KHOR,KLEV),        PGEOH(KHOR,KLEV),
     *         PAPH(KHOR,KLEVP1),      PTENH(KHOR,KLEV),
     *         PQENH(KHOR,KLEV),       PQSENH(KHOR,KLEV)
C
      REAL     PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PTD(KHOR,KLEV),         PQD(KHOR,KLEV),
     *         PUU(KHOR,KLEV),         PUD(KHOR,KLEV),
     *         PVU(KHOR,KLEV),         PVD(KHOR,KLEV),
     *         PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV),
     *         PMFUS(KHOR,KLEV),       PMFDS(KHOR,KLEV),
     *         PMFUQ(KHOR,KLEV),       PMFDQ(KHOR,KLEV),
     *         PDMFUP(KHOR,KLEV),      PDMFDP(KHOR,KLEV),
     *         PLU(KHOR,KLEV),         PLUDE(KHOR,KLEV)
      INTEGER  KLAB(KHOR,KLEV),        KLWMIN(KHOR)
C
      REAL PRFLCK(KHOR,KLEV)
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
      LOGICAL
     *        LOFLAG (JPHR)
      REAL
     *        ZWMAX  (JPHR)
C
C
C----------------------------------------------------------------------
C
C*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
C*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
C*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
C                  ----------------------------------------------
C
  100 CONTINUE
      ZDP=0.5
      DO 130 JK=2,KLEV
      DO 110 JL=KSTART,KSTOP
      PGEOH(JL,JK)=PGEO(JL,JK)+(PGEO(JL,JK-1)-PGEO(JL,JK))*ZDP
      SKM1=CPD*PTEN(JL,JK-1) + PGEO(JL,JK-1)
      SK  =CPD*PTEN(JL,JK)   + PGEO(JL,JK)
      PTENH(JL,JK)=(MAX(SK,SKM1)-PGEOH(JL,JK))*RCPD
CPS    IF (SKM1.GT.SK) THEN
         PQSENH(JL,JK)=PQSEN(JL,JK-1)
CPS    ELSE
CPS      PQSENH(JL,JK)=PQSEN(JL,JK)
CPS      WRITE (6,*) JK, 'modified'
CPS    ENDIF
      LOFLAG(JL)=.TRUE.
  110 CONTINUE
C
      IK=JK
      ICALL=0
      CALL CUADJTQ
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   IK,
     *     KSTART,   KSTOP,    KLEN,
     *     PAPH,     PTENH,    PQSENH,   LOFLAG,   ICALL
     *  )
C
      DO 120 JL=KSTART,KSTOP
      PQENH(JL,JK)=MIN(PQEN(JL,JK-1),PQSEN(JL,JK-1))
     1            +(PQSENH(JL,JK)-PQSEN(JL,JK-1))
      PQENH(JL,JK)=MAX(PQENH(JL,JK),0.)
  120 CONTINUE
  130 CONTINUE
C
      DO 140 JL=KSTART,KSTOP
      PTENH(JL,KLEV)=(CPD*PTEN(JL,KLEV)+PGEO(JL,KLEV)-
     1                PGEOH(JL,KLEV))*RCPD
      PQENH(JL,KLEV)=PQEN(JL,KLEV)
      PTENH(JL,1)=PTEN(JL,1)
      PQENH(JL,1)=PQEN(JL,1)
      PGEOH(JL,1)=PGEO(JL,1)
      KLWMIN(JL)=KLEV
      ZWMAX(JL)=0.
  140 CONTINUE
C
      DO 160 JK=KLEVM1,2,-1
      DO 150 JL=KSTART,KSTOP
      ZZS=MAX(CPD*PTENH(JL,JK)+PGEOH(JL,JK),
     1        CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))
      PTENH(JL,JK)=(ZZS-PGEOH(JL,JK))*RCPD
  150 CONTINUE
  160 CONTINUE
C
      DO 190 JK=KLEV,1,-1
c$dir no_recurrence, force_vector, force_parallel_ext
CDIR$ IVDEP
      DO 180 JL=KSTART,KSTOP
      IF(PVERV(JL,JK).LT.ZWMAX(JL)) THEN
         ZWMAX(JL)=PVERV(JL,JK)
         KLWMIN(JL)=JK
      END IF
  180 CONTINUE
  190 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
C*                 ---------------------------------------------
C
  200 CONTINUE
      DO 230 JK=1,KLEV
      IK=JK-1
      IF(JK.EQ.1) IK=1
c$dir no_recurrence, force_vector, force_parallel_ext
      DO 220 JL=KSTART,KSTOP
      PTU(JL,JK)=PTENH(JL,JK)
      PTD(JL,JK)=PTENH(JL,JK)
      PQU(JL,JK)=PQENH(JL,JK)
      PQD(JL,JK)=PQENH(JL,JK)
      PLU(JL,JK)=0.
      PUU(JL,JK)=PUEN(JL,IK)
      PUD(JL,JK)=PUEN(JL,IK)
      PVU(JL,JK)=PVEN(JL,IK)
      PVD(JL,JK)=PVEN(JL,IK)
      PMFU(JL,JK)=0.
      PMFD(JL,JK)=0.
      PMFUS(JL,JK)=0.
      PMFDS(JL,JK)=0.
      PMFUQ(JL,JK)=0.
      PMFDQ(JL,JK)=0.
      PDMFUP(JL,JK)=0.
      PDMFDP(JL,JK)=0.
      PLUDE(JL,JK)=0.
      KLAB(JL,JK)=0
      PRFLCK(JL,JK)=0.
  220 CONTINUE
  230 CONTINUE
C
      RETURN
      END
