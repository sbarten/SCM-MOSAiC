      SUBROUTINE CUASC4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,
     *     KLEV,     KLEVP1,   KLEVM1,
     *     NSTEP,    NSTART,   TWODT,
     *     PTENH,    PQENH,    PXENH,    PUEN,     PVEN,

C LG- adding the tracers

     *   KTRAC,
     *   PXTENH,   PXTEN,    PXTU,     PMFUXT,

C LG- end

     *     zpmfun,   ! op_ck_20031001

     *     PTEN,     PQEN,     PQSEN,
     *     PGEO,     PGEOH,    PAPP1,    PAPHP1,
     *     PQTE,     PVERV,    KLWMIN,   PDQPBL,
     *     LDLAND,   LDCUM,    KTYPE,    KLAB,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     PMFU,     PMFD,     PMFUB,    PENTR,
     *     PMFUS,    PMFUQ,
     *     PMFUL,    PLUDE,    PDMFUP,
     *     KHMIN,    PHHATT,   PHCBASE,  PQSENH,

     *     pcpen,    pcpcu,               ! mz_lg_20031117+

     *     KCBOT,    KCTOP,    KCTOP0,   KCUM)
C
C          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
C          FOR CUMULUS PARAMETERIZATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
C          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
C           FLUXES AS WELL AS PRECIPITATION RATES)
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C
C          METHOD.
C          --------
C          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
C          AND THEN CALCULATE MOIST ASCENT FOR
C          ENTRAINING/DETRAINING PLUME.
C          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
C          SHALLOW AND DEEP CUMULUS CONVECTION.
C          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
C          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
C          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
C
C          EXTERNALS
C          ---------
C          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
C          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
C          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
C
C          REFERENCE
C          ---------
C          (TIEDTKE,1989)
C
CJHC*CALL PARAM
CJHC*CALL COMSDS
CJHC*CALL COMCTL
CJHC*CALL COMCON
CJHC*CALL COMCUMF
CJHC*CALL COMTRC
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KLP2,KLEV),       PQENH(KLP2,KLEV),
     *         PXENH(KLP2,KLEV),
     *         PUEN(K2LP2,KLEV),       PVEN(K2LP2,KLEV),
     *         PTEN(KLP2,KLEV),        PQEN(KLP2,KLEV),
     *         PGEO(KLP2,KLEV),        PGEOH(KLP2,KLEV),
     *         PAPP1(KLP2,KLEV),       PAPHP1(KLP2,KLEVP1),
     *         PQSEN(KLP2,KLEV),       PQTE(KLP2,KLEV),
     *         PVERV(KLP2,KLEV),       PDQPBL(KLP2)
C
      REAL     PTU(KLP2,KLEV),         PQU(KLP2,KLEV),
     *         PUU(KLP2,KLEV),         PVU(KLP2,KLEV),
     *         PMFU(KLP2,KLEV),        PMFD(KLP2,KLEV),
     *         PMFUB(KLP2),            PENTR(KLP2),
     *         PMFUS(KLP2,KLEV),       PMFUQ(KLP2,KLEV),
     *         PLU(KLP2,KLEV),         PLUDE(KLP2,KLEV),
     *         PMFUL(KLP2,KLEV),       PDMFUP(KLP2,KLEV)

      REAL     zpmfun(KLP2,klev) ! op_ck_20031001

      REAL     pcpen(KLP2,klev),       pcpcu(KLP2,klev) ! mz_lg_20031117+

      INTEGER  KLWMIN(KLP2),           KTYPE(KLP2),
     *         KLAB(KLP2,KLEV),        KCBOT(KLP2),
     *         KCTOP(KLP2),            KCTOP0(KLP2)
      INTEGER  KHMIN(KLP2)
      REAL     PHHATT(KLP2,KLEV)
      REAL     PHCBASE(KLP2)
      REAL     PQSENH(KLP2,KLEV)
      LOGICAL  LDLAND(KLP2),           LDCUM(KLP2)
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      REAL     ZDMFEN(JPHR),           ZDMFDE(JPHR),
     *         ZMFUU(JPHR),            ZMFUV(JPHR),
     *         ZPBASE(JPHR),           ZQOLD(JPHR)
      REAL      ZDLAND(JPHR)
      REAL     ZPH(JPHR)
      REAL     ZODETR(JPHR,MLEV)
      REAL     ZOENTR(JPHR,MLEV)
      REAL     ZBUOY(JPHR)
      LOGICAL  LOFLAG(JPHR)
COBC
     *        , LCVMGT, LOTEST

C LG- adding the tracers

      REAL   PXTENH(KLON,KLEV,KTRAC),  PXTEN(KLON,KLEV,KTRAC),
     *       PXTU(KLON,KLEV,KTRAC),    PMFUXT(KLON,KLEV,KTRAC)

C LG- end

      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
C
C
C----------------------------------------------------------------------
C
C*    1.           SPECIFY PARAMETERS
C                  ------------------
C
  100 CONTINUE
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZCONS2=1./(G*ZTMST)
      ZTGLACE=TMELT-13.
C
C
C----------------------------------------------------------------------
C
C     2.           SET DEFAULT VALUES
C                  ------------------
C
  200 CONTINUE
      DO 210 JL=KIDIA,KFDIA
       ZMFUU(JL)=0.
       ZMFUV(JL)=0.
      IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
  210 CONTINUE
      DO 230 JK=1,KLEV
      DO 220 JL=KIDIA,KFDIA
      PLU(JL,JK)=0.
      PMFU(JL,JK)=0.

      zpmfun(jl,jk)=0. ! op_ck_20031001

      PMFUS(JL,JK)=0.
      PMFUQ(JL,JK)=0.
      PMFUL(JL,JK)=0.
      PLUDE(JL,JK)=0.
      PDMFUP(JL,JK)=0.
      IF(.NOT.LDCUM(JL).OR.KTYPE(JL).EQ.3) KLAB(JL,JK)=0
      IF(.NOT.LDCUM(JL).AND.PAPHP1(JL,JK).LT.4.E4) KCTOP0(JL)=JK
  220 CONTINUE

C LG- adding the tracers

      DO 2204 JT=1,KTRAC
      DO 2202 JL=KIDIA,KFDIA
       PMFUXT(JL,JK,JT)=0.
 2202 CONTINUE
 2204 CONTINUE

C LG- end

C
  230 CONTINUE
CDIR$ IVDEP
      DO 240 JL=KIDIA,KFDIA
      IF(LDLAND(JL)) THEN
         ZDLAND(JL)=3.0E4
         ZDPHI=PGEOH(JL,KCTOP0(JL))-PGEOH(JL,KCBOT(JL))
         IF(PTU(JL,KCTOP0(JL)).GE.ZTGLACE) ZDLAND(JL)=ZDPHI
         ZDLAND(JL)=MAX(3.0E4,ZDLAND(JL))
         ZDLAND(JL)=MIN(5.0E4,ZDLAND(JL))
      ZDLAND(JL)=1.5E4
      END IF
  240 CONTINUE
      DO JK=1,KLEV
       DO JL=KIDIA,KFDIA
        ZOENTR(JL,JK)=0.
        ZODETR(JL,JK)=0.
       ENDDO
      ENDDO
C
C
C----------------------------------------------------------------------
C
C     3.0          INITIALIZE VALUES AT LIFTING LEVEL
C                  ----------------------------------
C
  300 CONTINUE
      DO 310 JL=KIDIA,KFDIA
      KCTOP(JL)=KLEVM1
         IF(.NOT.LDCUM(JL)) THEN
            KCBOT(JL)=KLEVM1
            PMFUB(JL)=0.
            PQU(JL,KLEV)=0.
         END IF
      PMFU(JL,KLEV)=PMFUB(JL)
      PMFUS(JL,KLEV)=PMFUB(JL)*(CPD*PTU(JL,KLEV)+PGEOH(JL,KLEV))
      PMFUQ(JL,KLEV)=PMFUB(JL)*PQU(JL,KLEV)
         IF(LMFDUDV) THEN
            ZMFUU(JL)=PMFUB(JL)*PUU(JL,KLEV)
            ZMFUV(JL)=PMFUB(JL)*PVU(JL,KLEV)
         END IF
  310 CONTINUE
C
! op_ck_20031001+
      DO 3101 jk=1,klev
!DIR$ IVDEP
        DO 3102 jl=1,KIDIA,KFDIA
!
! LINEARISATION OF MASS FLUX FROM CLOUD BASE TO THE SURFACE:
! NEW VARIABLE FOR MASS FLUX (TRACER ONLY AT THIS STAGE): ZPMFUN
!
          IF (ktype(jl).ne.3) then
            IF (jk.eq.kcbot(jl)) zpmfun(jl,jk)=pmfub(jl)
            IF (jk.gt.kcbot(jl)) then
              ikb=kcbot(jl)
              zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/                
     *           (paphp1(jl,klevp1)-paphp1(jl,ikb))
              zpmfun(jl,jk)=zpmfun(jl,ikb)*zzp
            END IF
          END IF
3102    END DO
3101  END DO
!
      DO 311 jl=1,KIDIA,KFDIA
        pmfu(jl,klev)=pmfub(jl)
        pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)          
     *                                  +pgeoh(jl,klev))
        pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
        IF(lmfdudv) THEN
          zmfuu(jl)=pmfub(jl)*puu(jl,klev)
          zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
        END IF
311   END DO
!
C LG- adding the tracers

      DO 3112 JT=1,KTRAC
      DO 3110 JL=KIDIA,KFDIA
      IF(.NOT.LDCUM(JL)) THEN
        PXTU(JL,KLEV,JT)=0.
      ENDIF
C        PMFUXT(JL,KLEV,JT)=PMFUB(JL)*PXTU(JL,KLEV,JT)

        pmfuxt(jl,klev,jt)=zpmfun(jl,klev)*pxtu(jl,klev,jt) ! op_ck_20031001

 3110 CONTINUE
 3112 CONTINUE
 
C LG- end 

C
      DO 320 JL=KIDIA,KFDIA
      LDCUM(JL)=.FALSE.
  320 CONTINUE
C
C
C
C----------------------------------------------------------------------
C
C     3.5          FIND ORGANIZED ENTRAINMENT AT CLOUD BASE
C                  ----------------------------------------
C
  350 CONTINUE
      DO JL=KIDIA,KFDIA
       LOTEST=(KTYPE(JL).EQ.1.AND..NOT.LNORDSCV) .OR.
     1        (KTYPE(JL).NE.3.AND.LNORDSCV)
       IF(LOTEST) THEN
        IKB=KCBOT(JL)
        ZBUOY(JL)=G*(PTU(JL,IKB)-PTENH(JL,IKB))/PTENH(JL,IKB) +
     1            G*0.608*(PQU(JL,IKB)-PQENH(JL,IKB))
        IF(ZBUOY(JL).GT.0.) THEN
         ZDZ=(PGEO(JL,IKB-1)-PGEO(JL,IKB))/G
         ZDRODZ=-LOG(PTEN(JL,IKB-1)/PTEN(JL,IKB))/ZDZ
     1         -G/(RD*PTENH(JL,IKB))
         ! nb zoentr is here a fractional value
         ZOENTR(JL,IKB-1)=ZBUOY(JL)*0.5/(1.+ZBUOY(JL)*ZDZ)
     1                    + ZDRODZ
	 IF (.NOT.LNORDSCV) THEN
           IF(ZOENTR(JL,IKB-1).GT.1.E-3) ZOENTR(JL,IKB-1)=1.E-3
         ELSE
           IF(ZOENTR(JL,IKB-1).GT.1.E-2) ZOENTR(JL,IKB-1)=1.E-2
	 ENDIF
         IF(ZOENTR(JL,IKB-1).LT.0.0) ZOENTR(JL,IKB-1)=0.0
        ENDIF
       ENDIF
       ENDDO
C
C
C----------------------------------------------------------------------
C
C     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
C                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
C                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
C                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
C                  -------------------------------------------------
C
  400 CONTINUE
      DO 480 JK=KLEVM1,2,-1
C
C                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
C                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
C                  ----------------------------------------------------
C
      IK=JK
      CALL CUBASMC4
     *    (KIDIA, KFDIA,
     *     KLP2,     K2LP2,    KLON,     KLEV,     KLEVM1,   IK,
     *     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,

C LG- adding the tracers

     *   KTRAC,

     *   paphp1,   zpmfun,   ! op_ck_20031001

     *   PXTEN,    PXTU,     PMFUXT,

C LG- end

     *     PVERV,    PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,
     *     PMFU,     PMFUB,    PENTR,    KCBOT,
     *     PTU,      PQU,      PLU,      PUU,      PVU,
     *     PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   ZMFUU,    ZMFUV)
C
      IS=0
      DO 410 JL=KIDIA,KFDIA
      IS=IS+KLAB(JL,JK+1)
      IF(KLAB(JL,JK+1).EQ.0) KLAB(JL,JK)=0
      LOFLAG(JL)=LCVMGT(.TRUE.,.FALSE.,KLAB(JL,JK+1).GT.0)
      ZPH(JL)=PAPHP1(JL,JK)
      IF(KTYPE(JL).EQ.3.AND.JK.EQ.KCBOT(JL)) THEN
         ZMFMAX=(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))*ZCONS2
         IF(PMFUB(JL).GT.ZMFMAX) THEN
            ZFAC=PMFUB(JL)/ZMFMAX
            PMFU(JL,JK+1)=PMFU(JL,JK+1)*ZFAC
            PMFUS(JL,JK+1)=PMFUS(JL,JK+1)*ZFAC
            PMFUQ(JL,JK+1)=PMFUQ(JL,JK+1)*ZFAC
            ZMFUU(JL)=ZMFUU(JL)*ZFAC
            ZMFUV(JL)=ZMFUV(JL)*ZFAC
         END IF
      END IF
  410 CONTINUE

C LG- adding the tracers

      DO 4102 JT=1,KTRAC
      DO 4101 JL=KIDIA,KFDIA
      IF(KTYPE(JL).EQ.3.AND.JK.EQ.KCBOT(JL)) THEN
         ZMFMAX=(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))*ZCONS2
         IF(PMFUB(JL).GT.ZMFMAX) THEN
C             ZFAC=PMFUB(JL)/ZMFMAX

            zfac=zmfmax/pmfub(jl)  !!!!!!

            PMFUXT(JL,JK+1,JT)=PMFUXT(JL,JK+1,JT)*ZFAC
         END IF
      END IF
 4101 CONTINUE
 4102 CONTINUE

C LG- end

C
C RESET PMFUB IF NECESSARY
C
      DO 4103 JL=KIDIA,KFDIA
      IF(KTYPE(JL).EQ.3.AND.JK.EQ.KCBOT(JL)) THEN
         ZMFMAX=(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))*ZCONS2
         IF(PMFUB(JL).GT.ZMFMAX) THEN
            PMFUB(JL)=ZMFMAX
         END IF
      END IF
 4103 CONTINUE

      IF(IS.EQ.0) GO TO 480
C
C
C*                 SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
C                  RATES PLUS ORGANIZED DETRAINMENT RATES IN *CUENTR*
C                   -------------------------------------
C
      IK=JK
      CALL CUENTR4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,   IK,
     *     KLAB,
     *     PTENH,    PQENH,    PQTE,     PAPHP1,   PAPP1,
     *     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     *     ZPBASE,   PMFU,     PENTR,    ZODETR,  ZOENTR,
     *     KHMIN,    ZBUOY,    PGEOH,
     *     ZDMFEN,   ZDMFDE)
C
C
C
C                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
C                  THE CLOUD ENSEMBLE ENTRAINS ENVIRONMENTAL VALUES
C                  IN TURBULENT DETRAINMENT CLOUD ENSEMBLE VALUES
C                  ARE DETRAINED
C                  IN ORGANIZED DETRAINMENT THE DRY STATIC ENERGY AND
C                  MOISTURE THAT ARE NEUTRAL COMPARED TO THE
C                  ENVIRONMENTAL AIR ARE DETRAINED
C                  ---------------------------------------------------
C
      DO 420 JL=KIDIA,KFDIA
      IF(LOFLAG(JL)) THEN
      IF(JK.LT.KCBOT(JL)) THEN
	 IF (LNORDSCV) THEN
	   ZDMFEN(JL)=0.
	   ZDMFDE(JL)=0.
	 ENDIF
         ZMFTEST=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
         ZMFMAX=MIN(ZMFTEST,(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))*ZCONS2)
         ZDMFEN(JL)=MAX(ZDMFEN(JL)-MAX(ZMFTEST-ZMFMAX,0.),0.)
      END IF
CEVM     ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75*PMFU(JL,JK+1))
CEVM     A.P. SIEBESMA modification on shallow convection
         IF(LMFNEW)THEN
	   ZDMFDE(JL)=MIN(ZDMFDE(JL),1.0 *PMFU(JL,JK+1))
         ELSE
	   ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75*PMFU(JL,JK+1))
	 ENDIF
CEVM    
         PMFU(JL,JK)=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
      IF(JK.LT.KCBOT(JL)) THEN
         ZDPRHO=(PGEOH(JL,JK)-PGEOH(JL,JK+1))/G
         ZOENTR(JL,JK)=ZOENTR(JL,JK)*ZDPRHO*PMFU(JL,JK+1)
         ZMFTEST=PMFU(JL,JK)+ZOENTR(JL,JK)-ZODETR(JL,JK)
         ZMFMAX=MIN(ZMFTEST,(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))*ZCONS2)
         ZOENTR(JL,JK)=MAX(ZOENTR(JL,JK)-MAX(ZMFTEST-ZMFMAX,0.),0.)
      ENDIF
         LOTEST=(KTYPE(JL).EQ.1.AND..NOT.LNORDSCV) .OR.
     1          (KTYPE(JL).NE.3.AND.LNORDSCV)
         IF(LOTEST.AND.JK.LT.KCBOT(JL).
     I         AND.JK.LE.KHMIN(JL)) THEN
         ! limit organized detrainment to not allowing for too
         ! deep clouds
         ZMSE =CPD*PTU(JL,JK+1)+ALV*PQU(JL,JK+1)+PGEOH(JL,JK+1)
         IKT=KCTOP0(JL)
         ZNEVN=(PGEOH(JL,IKT)-PGEOH(JL,JK+1))*
     I         (ZMSE-PHHATT(JL,JK+1))/G
               IF(ZNEVN.LE.0.) ZNEVN=1.
         ZDPRHO=(PGEOH(JL,JK)-PGEOH(JL,JK+1))/G
         ZODMAX=((PHCBASE(JL)-ZMSE)/ZNEVN)*ZDPRHO*PMFU(JL,JK+1)
         ZODMAX=MAX(ZODMAX,0.)
         ZODETR(JL,JK)=MIN(ZODETR(JL,JK),ZODMAX)
         ENDIF
         ZODETR(JL,JK)=MIN(ZODETR(JL,JK),0.75*PMFU(JL,JK))
         PMFU(JL,JK)=PMFU(JL,JK)+ZOENTR(JL,JK)-ZODETR(JL,JK)
         ZQEEN=PQENH(JL,JK+1)*ZDMFEN(JL)
         ZQEEN=ZQEEN+
     I         PQENH(JL,JK+1)*ZOENTR(JL,JK)
         ZSEEN=(CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFEN(JL)
         ZSEEN=ZSEEN+
     I         (CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZOENTR(JL,JK)
         ZSCDE=(CPD*PTU(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFDE(JL)
         ! find moist static energy that give nonbuoyant air
         ZGA=ALV*PQSENH(JL,JK+1)/(RV*(PTENH(JL,JK+1)**2))
         ZDT=(PLU(JL,JK+1)-0.608*(PQSENH(JL,JK+1)-PQENH(JL,JK+1)))/
     I         (1./PTENH(JL,JK+1) + 0.608*ZGA)
         ZSCOD=CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1)+CPD*ZDT
         ZSCDE=ZSCDE+ZODETR(JL,JK)*ZSCOD
         ZQUDE=PQU(JL,JK+1)*ZDMFDE(JL)
         ZQCOD=PQSENH(JL,JK+1)+ZGA*ZDT
         ZQUDE=ZQUDE+ZODETR(JL,JK)*ZQCOD
         PLUDE(JL,JK)=PLU(JL,JK+1)*ZDMFDE(JL)
         PLUDE(JL,JK)=PLUDE(JL,JK)+PLU(JL,JK+1)*ZODETR(JL,JK)
         ZMFUSK=PMFUS(JL,JK+1)+ZSEEN-ZSCDE
         ZMFUQK=PMFUQ(JL,JK+1)+ZQEEN-ZQUDE
         ZMFULK=PMFUL(JL,JK+1)    -PLUDE(JL,JK)
         PLU(JL,JK)=ZMFULK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
         PQU(JL,JK)=ZMFUQK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
         PTU(JL,JK)=(ZMFUSK*(1./MAX(CMFCMIN,PMFU(JL,JK)))-
     1               PGEOH(JL,JK))*RCPD
         PTU(JL,JK)=MAX(100.,PTU(JL,JK))
         PTU(JL,JK)=MIN(400.,PTU(JL,JK))
         ZQOLD(JL)=PQU(JL,JK)
      END IF
  420 CONTINUE
C
C

C LG- adding the tracers

      DO 4204 JT=1,KTRAC
      DO 4202 JL=KIDIA,KFDIA
      IF(LOFLAG(JL)) THEN

! op_ck_20031001+
!        ZXTEEN=PXTENH(JL,JK+1,JT)*(ZDMFEN(JL)+ZOENTR(JL,JK))
!        ZXTUDE=PXTU(JL,JK+1,JT)*(ZDMFDE(JL)+ZODETR(JL,JK))
!        ZMFUXTK=PMFUXT(JL,JK+1,JT)+ZXTEEN-ZXTUDE
!        PXTU(JL,JK,JT)=ZMFUXTK*(1./MAX(CMFCMIN,PMFU(JL,JK)))

             IF (jk.ge.kcbot(jl)) then
               IF (ktype(jl).ne.3) then
                 pmfuxt(jl,jk,jt)=(pxtu(jl,jk+1,jt)*zpmfun(jl,jk+1)    
     *              + (zpmfun(jl,jk)-zpmfun(jl,jk+1))*pxten(jl,jk,jt))
                 pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)                       
     *                          /max(cmfcmin,zpmfun(jl,jk))
               ELSE
                 zpmfun(jl,jk)=pmfu(jl,jk)
                 pmfuxt(jl,jk,jt)=pmfuxt(jl,jk+1,jt)
                 pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)                       
     *                          /max(cmfcmin,pmfu(jl,jk))
               ENDIF
             ELSE
               zpmfun(jl,jk)=pmfu(jl,jk)
               zxteen=pxten(jl,jk,jt)*(zdmfen(jl)+zoentr(jl,jk))
               zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen
               pmfuxt(jl,jk,jt)=zmfuxtk/(1.+(zdmfde(jl)+zodetr(jl,jk)) 
     *                          /max(cmfcmin,pmfu(jl,jk)))
               pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)/max(cmfcmin,pmfu(jl,jk))
             ENDIF
! op_ck_20031001-

      ENDIF
 4202 CONTINUE
 4204 CONTINUE
 
C LG- end 

C
C
C                  DO CORRECTIONS FOR MOIST ASCENT
C                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
C                  -----------------------------------
C
      IK=JK
      ICALL=1
      CALL CUADJTQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     IK,
     *     ZPH,      PTU,      PQU,      LOFLAG,   ICALL)
C
CDIR$ IVDEP
      DO 440 JL=KIDIA,KFDIA
      IF(LOFLAG(JL).AND.PQU(JL,JK).NE.ZQOLD(JL)) THEN
         KLAB(JL,JK)=2
         PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
         ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK)-PLU(JL,JK))-
     1   PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK)-PXENH(JL,JK))
         IF(KLAB(JL,JK+1).EQ.1) ZBUO=ZBUO+0.5
         IF(ZBUO.GT.0..AND.PMFU(JL,JK).GE.0.01*PMFUB(JL).AND.
     *      JK.GE.KCTOP0(JL)) THEN
            KCTOP(JL)=JK
            LDCUM(JL)=.TRUE.
      ZDNOPRC=CVMGT(ZDLAND(JL),1.5E4,LDLAND(JL))
            ZPRCON=CVMGT(0.,CPRCON,ZPBASE(JL)-PAPHP1(JL,JK).LT.ZDNOPRC)
            ZLNEW=PLU(JL,JK)/(1.+ZPRCON*(PGEOH(JL,JK)-PGEOH(JL,JK+1)))
            PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
            PLU(JL,JK)=ZLNEW
         ELSE

            zpmfun(jl,jk)=0. ! op_ck_20031001

            KLAB(JL,JK)=0
            PMFU(JL,JK)=0.
         END IF
      END IF
  440 CONTINUE
      DO 455 JL=KIDIA,KFDIA
      IF(LOFLAG(JL)) THEN
         PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
         PMFUS(JL,JK)=(CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
         PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
      END IF
  455 CONTINUE

C LG- adding the tracers

      DO 4554 JT=1,KTRAC
      DO 4552 JL=KIDIA,KFDIA
      IF(LOFLAG(JL)) THEN
C         PMFUXT(JL,JK,JT)=PXTU(JL,JK,JT)*PMFU(JL,JK)

        pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*zpmfun(jl,jk) ! op_ck_20031001

      ENDIF
 4552 CONTINUE
 4554 CONTINUE

C LG- end

C
        IF(LMFDUDV) THEN
         DO JL=KIDIA,KFDIA
           ZDMFEN(JL)=ZDMFEN(JL)+ZOENTR(JL,JK)
           ZDMFDE(JL)=ZDMFDE(JL)+ZODETR(JL,JK)
         ENDDO
           DO 460 JL=KIDIA,KFDIA
           IF(LOFLAG(JL)) THEN
              IF(KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.3) THEN
                 ZZ=CVMGT(3.,2.,ZDMFEN(JL).EQ.0.)
              ELSE
                 ZZ=CVMGT(1.,0.,ZDMFEN(JL).EQ.0.)
              END IF
                 ZDMFEU=ZDMFEN(JL)+ZZ*ZDMFDE(JL)
                 ZDMFDU=ZDMFDE(JL)+ZZ*ZDMFDE(JL)
                 ZDMFDU=MIN(ZDMFDU,0.75*PMFU(JL,JK+1))
              ZMFUU(JL)=ZMFUU(JL)+
     1         ZDMFEU*PUEN(JL,JK)-ZDMFDU*PUU(JL,JK+1)
              ZMFUV(JL)=ZMFUV(JL)+
     1         ZDMFEU*PVEN(JL,JK)-ZDMFDU*PVU(JL,JK+1)
              IF(PMFU(JL,JK).GT.0.) THEN
                 PUU(JL,JK)=ZMFUU(JL)*(1./PMFU(JL,JK))
                 PVU(JL,JK)=ZMFUV(JL)*(1./PMFU(JL,JK))
              END IF
           END IF
  460      CONTINUE
        END IF
C
C
C
C                  COMPUTE ORGANIZED ENTRAINMENT
C                  FOR USE AT NEXT LEVEL
C                  ------------------------------
C
      DO JL=KIDIA,KFDIA
      LOTEST=(KTYPE(JL).EQ.1.AND..NOT.LNORDSCV) .OR.
     1       (KTYPE(JL).NE.3.AND.LNORDSCV)
      IF(LOFLAG(JL).AND.LOTEST) THEN
        ZBUOYZ=G*(PTU(JL,JK)-PTENH(JL,JK))/PTENH(JL,JK) +
     1         G*0.608*(PQU(JL,JK)-PQENH(JL,JK))-G*PLU(JL,JK)
        ZBUOYZ=MAX(ZBUOYZ,0.0)
        ZDZ=(PGEO(JL,JK-1)-PGEO(JL,JK))/G
        ZDRODZ=-LOG(PTEN(JL,JK-1)/PTEN(JL,JK))/ZDZ
     1         -G/(RD*PTENH(JL,JK))
        ZBUOY(JL)=ZBUOY(JL)+ZBUOYZ*ZDZ
        ZOENTR(JL,JK-1)=ZBUOYZ*0.5/(1.+ZBUOY(JL))
     1                 + ZDRODZ
	IF (.NOT.LNORDSCV) THEN
          IF(ZOENTR(JL,JK-1).GT.1.E-3) ZOENTR(JL,JK-1)=1.E-3
        ELSE
          IF(ZOENTR(JL,JK-1).GT.1.E-2) ZOENTR(JL,JK-1)=1.E-2
	ENDIF
        IF(ZOENTR(JL,JK-1).LT.0.0)   ZOENTR(JL,JK-1)=0.
C
      ENDIF
      ENDDO
C
C
  480 CONTINUE
C
C
C----------------------------------------------------------------------
C
C     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
C                  ----------------------------------------------------
C                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
C                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
C                         FROM PREVIOUS CALCULATIONS ABOVE)
C
  500 CONTINUE
      DO 510 JL=KIDIA,KFDIA
      IF(KCTOP(JL).EQ.KLEVM1) LDCUM(JL)=.FALSE.
      KCBOT(JL)=MAX(KCBOT(JL),KCTOP(JL))
  510 CONTINUE
      IS=0
      DO 520 JL=KIDIA,KFDIA
      IS=IS+ICVMGT(1,0,LDCUM(JL))
  520 CONTINUE
      KCUM=IS

      IF(IS.EQ.0) GO TO 800
C=================================================================
C
C-------------------------------------------------------------------
C								   |
C     NEW RESOLUTION-INDEPENDENT TREATMENT OF THE OVERSHOOT-REGIME |
C     ONLY ACTIVE IN HIGH-RESOLUTION MODE: LMFHRES = .TRUE.	   |
C								   |
C-------------------------------------------------------------------
C
      JL = 1
      IF (LMFHRES.AND.(KTYPE(JL).LT.3)) THEN
C
C-------------------------------------------------------------------
C								   |
C     DETERMINE LEVELS WHERE OVERSHOOT NEEDS TO BE CALCULATED	   |
C	           LABEL WITH KLAB = 3				   |
C								   |
C								   |
C     ICTOP: HIGHEST LEVEL WITH OVERSHOOT .LT. DP_SHOOT		   |
C            HIGHEST LEVEL WITH MASS_FLUX .GT. 0.0		   |
C								   |
C-------------------------------------------------------------------
C
      JK = KCTOP(JL)
525   CONTINUE
      ICTOP         = JK
      JK            = JK-1
      GUARD = PAPHP1(JL,KCTOP(JL)) - PAPHP1(JL,JK)
      IF (GUARD.LT.DP_SHOOT ) GOTO 525
C
      DZTOT   = (PGEOH(JL,ICTOP-1) - PGEOH(JL,KCTOP(JL)))/G
      DECAY   = 0.5*DZTOT/
     &          LOG((1+SQRT(1-4*CMFCTOP*(1-CMFCTOP)))/(2*CMFCTOP))
C
C-------------------------------------------------------------------
C								   |
C     DO MOIST ADIABATIC ASCENT IN OVERSHOOT LAYER		   |
C     WITH TURBULENT ENTRAINMENT = TURBULENT DETRAINMENT	   |
C     PLUS MASSIVE DETRAINMENT CAUSING				   | 
C     AN EXPONENTIAL DEACY OF THE MASS FLUX			   |
C-------------------------------------------------------------------
C
      DO 528 JK=KCTOP(JL)-1,ICTOP,-1
C
C-------------------------------------------------------------------
C								   |
C     LABEL OVERSHOOT LAYER: KLAB = 3				   |
C     DO TURBULENT ENTRAINMENT + DETRAINMENT			   | 
C								   |
C-------------------------------------------------------------------
C
c      KLAB  (JL,JK) = 3
      IK=JK
      CALL CUENTR4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,   IK,
     *     KLAB,
     *     PTENH,    PQENH,    PQTE,     PAPHP1,   PAPP1,
     *     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     *     ZPBASE,   PMFU,     PENTR,    ZODETR,  ZOENTR,
     *     KHMIN,    ZBUOY,    PGEOH,
     *     ZDMFEN,   ZDMFDE)
C
C-------------------------------------------------------------------
C								   |
C      CALCULATE A EXPONENTIAL DECAY OF THE MASS FLUX		   |
C      DUE TO MASSIVE DETRAINMENT:						   |
C								   |
C								   |
C	  M(Z(KCTOP))                           = M_TOP		   |
C         M(Z(KCTOP)+0.5*(Z(ICTOP-1)-(Z(KCTOP)) = CMFCTOP*M_TOP	   |
C         M(Z(ICTOP-1))                         = 0.0		   |
C								   |
C								   |
C     ICTOP: HIGHEST LEVEL WITH OVERSHOOT .LT. DP_SHOOT		   |
C            HIGHEST LEVEL WITH MASS_FLUX .GT. 0.0		   |
C								   |
C-------------------------------------------------------------------
C
       PMFU(JL,JK)   = PMFU(JL,KCTOP(JL))
     &  *(EXP((PGEOH(JL,ICTOP-1)-PGEOH(JL,JK))       /(DECAY*G)) - 1.0)
     &  /(EXP((PGEOH(JL,ICTOP-1)-PGEOH(JL,KCTOP(JL)))/(DECAY*G)) - 1.0) 
                      
C
C-------------------------------------------------------------------
C								   |
C      DO DRY-ADIABTIC ASCENT WITH ENTRAINMENT AND DETRAINMENT	   |
C      FOR Q, T AND L 						   |
C								   |
C-------------------------------------------------------------------
C
CEVM
CEVM   960410
CEVM   leave out test on size of entrainment (recommended by APS)
CEVM
CEVM   IF(JK.LT.KCBOT(JL)) THEN
CEVM     ZMFTEST=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
CEVM     ZMFMAX=MIN(ZMFTEST,(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))*ZCONS2)
CEVM     ZDMFEN(JL)=MAX(ZDMFEN(JL)-MAX(ZMFTEST-ZMFMAX,0.),0.)
CEVM   END IF
       ZDMFDE(JL)    = PMFU  (JL,JK+1) - PMFU(JL,JK) + ZDMFDE(JL)
CEVM     A.P. SIEBESMA modification on shallow convection
       IF(LMFNEW)THEN
           IF ( ZDMFDE(JL).GT.PMFU(JL,JK+1) ) THEN
	    ZDMFDE(JL)  = PMFU(JL,JK+1)
            PMFU(JL,JK) = ZDMFEN(JL)
           ENDIF
       ELSE
	   ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75*PMFU(JL,JK+1))
       ENDIF
       ZQEEN       = PQENH(JL,JK+1)*ZDMFEN(JL)
       ZSEEN       = (CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFEN(JL)
       ZSCDE       = (CPD*PTU  (JL,JK+1)+PGEOH(JL,JK+1))*ZDMFDE(JL)
       ZQUDE       = PQU(JL,JK+1)*ZDMFDE(JL)
       PLUDE(JL,JK)= PLU(JL,JK+1)*ZDMFDE(JL)
       ZMFUSK      = PMFUS(JL,JK+1)+ZSEEN-ZSCDE
       ZMFUQK      = PMFUQ(JL,JK+1)+ZQEEN-ZQUDE
       ZMFULK      = PMFUL(JL,JK+1)-PLUDE(JL,JK)
       PLU(JL,JK)  = ZMFULK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
       PQU(JL,JK)  = ZMFUQK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
       PTU(JL,JK)  = (ZMFUSK*(1./MAX(CMFCMIN,PMFU(JL,JK)))-
     1                  PGEOH(JL,JK))*RCPD
       PTU(JL,JK)  = MAX(100.,PTU(JL,JK))
       PTU(JL,JK)  = MIN(400.,PTU(JL,JK))
       ZQOLD(JL)   = PQU(JL,JK)
       IF(LMFDUDV) THEN
            PUU(JL,JK) = PUU(JL,JK+1)
            PVU(JL,JK) = PVU(JL,JK+1)
       ENDIF
C
C-------------------------------------------------------------------
C								   |
C          DO CORRECTIONS FOR MOIST ASCENT			   |
C          BY ADJUSTING T,Q AND L IN *CUADJTQ*			   |
C								   |
C-------------------------------------------------------------------
C
           ZQOLD(JL) = PQU(JL,JK)
           IK        = JK
           ICALL     = 1
           LOFLAG(JL) = .TRUE.
c
           CALL CUADJTQ4
     *      (KIDIA, KFDIA,
     *       KLP2,     KLON,     KLEV,     IK,
     *       ZPH,      PTU,      PQU,      LOFLAG,   ICALL)
C
           PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
           ZDNOPRC=CVMGT(3.0E4,1.5E4,LDLAND(JL))
           ZPRCON=CVMGT(0.,CPRCON,ZPBASE(JL)-PAPHP1(JL,JK).LT.ZDNOPRC)
           ZLNEW=PLU(JL,JK)/(1.+ZPRCON*(PGEOH(JL,JK)-PGEOH(JL,JK+1)))
           PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
           PLU(JL,JK)=ZLNEW
C
           PMFUL(JL,JK) =  PLU(JL,JK)                  *PMFU(JL,JK)
           PMFUS(JL,JK) = (CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
           PMFUQ(JL,JK) =  PQU(JL,JK)                  *PMFU(JL,JK)
C
 528    CONTINUE
C
C
C-------------------------------------------------------------------
C								   |
C          DETRAIN REMAINING LIQUID WATER AT THE FIRST LEVEL	   |
C	   WHERE MASS_FLUX = 0.0				   |
C								   |
C-------------------------------------------------------------------
C
          PLUDE(JL,ICTOP-1) = PMFUL(JL,ICTOP)
C
        ENDIF
C
C**********  END OF RESOLUTION_INDEPENDENT TREATMENT ***********


C       DO ORIGINAL TIEDTKE OVERSHOOT-TREATMENT IF
C       MODEL IS NOT IN HIGH_RES MODE
C       -------------------------------------------
C
C=================================================================
      IF (.NOT.LMFHRES.OR.KTYPE(JL).EQ.3) THEN

CDIR$ IVDEP
      DO 530 JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
         JK=KCTOP(JL)-1
         ZZDMF=CMFCTOP
         ZDMFDE(JL)=(1.-ZZDMF)*PMFU(JL,JK+1)
         PLUDE(JL,JK)=ZDMFDE(JL)*PLU(JL,JK+1)
         PMFU(JL,JK)=PMFU(JL,JK+1)-ZDMFDE(JL)
         ZLNEW=PLU(JL,JK)
         PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
         PLU(JL,JK)=ZLNEW
         PMFUS(JL,JK)=(CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
         PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
         PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
         PLUDE(JL,JK-1)=PMFUL(JL,JK)

         zpmfun(jl,jk)=pmfu(jl,jk) ! op_ck_20031001

      END IF
  530 CONTINUE

C LG- adding the tracers

      DO 5312 JT=1,KTRAC
      DO 5310 JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        JK=KCTOP(JL)-1
C         PMFUXT(JL,JK,JT)=PXTU(JL,JK,JT)*PMFU(JL,JK)

        pmfuxt(jl,jk,jt)=zpmfun(jl,jk)*pxtu(jl,jk,jt) ! op_ck_20031001

      ENDIF
 5310 CONTINUE
 5312 CONTINUE

C LG- end

C
        IF(LMFDUDV) THEN
CDIR$      IVDEP
           DO 540 JL=KIDIA,KFDIA
           IF(LDCUM(JL)) THEN
              JK=KCTOP(JL)-1
              PUU(JL,JK)=PUU(JL,JK+1)
              PVU(JL,JK)=PVU(JL,JK+1)
           END IF
  540      CONTINUE
        END IF
C
      END IF
  800 CONTINUE
C
      RETURN
      END
