      SUBROUTINE CUINI4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     PTEN,   PQEN,   PQSEN,   PXEN,   PUEN,   PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   PXTEN,    PXTENH,   PXTU,     PXTD,     PMFUXT,   PMFDXT,

C LG- end

     *     zpmfun,                        ! op_ck_20031001

     *     PVERV,    PGEO,     PAPHP1,   PGEOH,
     *     PTENH,    PQENH,    PQSENH,   PXENH,  KLWMIN,
     *     PTU,      PQU,      PTD,      PQD,
     *     PUU,      PVU,      PUD,      PVD,
     *     PMFU,     PMFD,     PMFUS,    PMFDS,
     *     PMFUQ,    PMFDQ,    PDMFUP,   PDMFDP, PRFLCK,

     *     pcpen,    pcpcu,               ! mz_lg_20031117+

     *     PDPMEL,   PLU,      PLUDE,    KLAB)
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
CJHC*CALL COMCON
      INCLUDE 'comcon.h'
C
      REAL     PTEN(KLP2,KLEV),        PQEN(KLP2,KLEV),
     *         PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PQSEN(KLP2,KLEV),       PVERV(KLP2,KLEV),
     *         PGEO(KLP2,KLEV),        PGEOH(KLP2,KLEV),
     *         PAPHP1(KLP2,KLEVP1),    PTENH(KLP2,KLEV),
     *         PXENH(KLP2,KLEV),   PXEN(KLP2,KLEV),
     *         PQENH(KLP2,KLEV),       PQSENH(KLP2,KLEV)
C
      REAL     pcpen(KLP2,klev),       pcpcu(KLP2,klev) ! mz_lg-20031117+

      REAL     PTU(KLP2,KLEV),         PQU(KLP2,KLEV),
     *         PTD(KLP2,KLEV),         PQD(KLP2,KLEV),
     *         PUU(KLP2,KLEV),         PUD(KLP2,KLEV),
     *         PVU(KLP2,KLEV),         PVD(KLP2,KLEV),
     *         PMFU(KLP2,KLEV),        PMFD(KLP2,KLEV),
     *         PMFUS(KLP2,KLEV),       PMFDS(KLP2,KLEV),
     *         PMFUQ(KLP2,KLEV),       PMFDQ(KLP2,KLEV),
     *         PDMFUP(KLP2,KLEV),      PDMFDP(KLP2,KLEV),
     *         PLU(KLP2,KLEV),         PLUDE(KLP2,KLEV)
      REAL     PDPMEL(KLP2,KLEV)
      REAL     PRFLCK(KLP2,KLEV)
      INTEGER  KLAB(KLP2,KLEV),        KLWMIN(KLP2)
C
      INCLUDE 'paramh.h'
      REAL     ZWMAX(JPHR)

      REAL     zpmfun(JPHR,klev) ! op_ck_20031001

      REAL     ZPH(JPHR)
      LOGICAL  LOFLAG(JPHR)

C LG- adding the tracers

      REAL   PXTEN(KLON,KLEV,KTRAC),   PXTENH(KLON,KLEV,KTRAC),
     *       PXTU(KLON,KLEV,KTRAC),    PXTD(KLON,KLEV,KTRAC),
     *       PMFUXT(KLON,KLEV,KTRAC),  PMFDXT(KLON,KLEV,KTRAC)

C LG- end

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

! mz_lg_20031117+
      DO 101 jk=1,klev
        DO 102 jl=1,KIDIA,KFDIA
          pcpen(jl,jk)=cpd
102     END DO
101   END DO
! mz_lg_20031117-

      ZDP=0.5
      DO 130 JK=2,KLEV
      DO 110 JL=KIDIA,KFDIA
      PGEOH(JL,JK)=PGEO(JL,JK)+(PGEO(JL,JK-1)-PGEO(JL,JK))*ZDP
      PTENH(JL,JK)=(MAX(CPD*PTEN(JL,JK-1)+PGEO(JL,JK-1),
     1             CPD*PTEN(JL,JK)+PGEO(JL,JK))-PGEOH(JL,JK))*RCPD
      PQSENH(JL,JK)=PQSEN(JL,JK-1)
      ZPH(JL)=PAPHP1(JL,JK)
      LOFLAG(JL)=.TRUE.
  110 CONTINUE
C

C LG- adding the tracers

      DO 1104 JT=1,KTRAC
      DO 1102 JL=KIDIA,KFDIA
       PXTENH(JL,JK,JT)=(PXTEN(JL,JK,JT)+PXTEN(JL,JK-1,JT))
     *                    *ZDP
 1102 CONTINUE
 1104 CONTINUE

C LG- end

C
C
      IK=JK
      ICALL=0
      CALL CUADJTQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     IK,
     *     ZPH,      PTENH,    PQSENH,   LOFLAG,   ICALL)
C
      DO 120 JL=KIDIA,KFDIA
      PXENH(JL,JK)=(PXEN(JL,JK)+PXEN(JL,JK-1))*ZDP
      PQENH(JL,JK)=MIN(PQEN(JL,JK-1),PQSEN(JL,JK-1))
     1            +(PQSENH(JL,JK)-PQSEN(JL,JK-1))
      PQENH(JL,JK)=MAX(PQENH(JL,JK),0.)

      pcpcu(jl,jk)=cpd  ! mz_lg_20031117+

  120 CONTINUE
  130 CONTINUE
C
      DO 140 JL=KIDIA,KFDIA
      PTENH(JL,KLEV)=(CPD*PTEN(JL,KLEV)+PGEO(JL,KLEV)-
     1                PGEOH(JL,KLEV))*RCPD
      PXENH(JL,KLEV)=PXEN(JL,KLEV)
      PQENH(JL,KLEV)=PQEN(JL,KLEV)
      PTENH(JL,1)=PTEN(JL,1)
      PXENH(JL,1)=PXEN(JL,1)
      PQENH(JL,1)=PQEN(JL,1)
      PGEOH(JL,1)=PGEO(JL,1)
      KLWMIN(JL)=KLEV
      ZWMAX(JL)=0.
  140 CONTINUE
C

C LG- adding the tracers

      DO 1404 JT=1,KTRAC
      DO 1402 JL=KIDIA,KFDIA
       PXTENH(JL,KLEV,JT)=PXTEN(JL,KLEV,JT)
       PXTENH(JL,1,JT)=PXTEN(JL,1,JT)
 1402 CONTINUE
 1404 CONTINUE

C LG- end

C
C
      DO 160 JK=KLEVM1,2,-1
      DO 150 JL=KIDIA,KFDIA
      ZZS=MAX(CPD*PTENH(JL,JK)+PGEOH(JL,JK),
     1        CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))
      PTENH(JL,JK)=(ZZS-PGEOH(JL,JK))*RCPD
  150 CONTINUE
  160 CONTINUE
C
      DO 190 JK=KLEV,3,-1
CDIR$ IVDEP
      DO 180 JL=KIDIA,KFDIA
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
      DO 220 JL=KIDIA,KFDIA
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

      zpmfun(jl,jk)=0. ! op_ck_20031001

      PMFUS(JL,JK)=0.
      PMFDS(JL,JK)=0.
      PMFUQ(JL,JK)=0.
      PMFDQ(JL,JK)=0.
      PDMFUP(JL,JK)=0.
      PDMFDP(JL,JK)=0.
      PDPMEL(JL,JK)=0.
      PLUDE(JL,JK)=0.
      KLAB(JL,JK)=0
      PRFLCK(JL,JK)=0.
  220 CONTINUE
C

C LG- adding the tracers

      DO 2204 JT=1,KTRAC
      DO 2202 JL=KIDIA,KFDIA
       PXTU(JL,JK,JT)=PXTENH(JL,JK,JT)
       PXTD(JL,JK,JT)=PXTENH(JL,JK,JT)
       PMFUXT(JL,JK,JT)=0.
       PMFDXT(JL,JK,JT)=0.
 2202 CONTINUE
 2204 CONTINUE

C LG- end

C
  230 CONTINUE
C
      RETURN
      END
