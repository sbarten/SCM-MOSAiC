      SUBROUTINE CLOUDSC_TCL
C---input
     S    (KIDIA,    KFDIA,    KHOR,     KLON,     KLEV,
     S     KSTEP,    PTSPHY,
     S     PT,       PQ,       PTENT,    PTENQ,
     S     PDTDIAB,  PHRSW,    PHRLW,
     S     PU,       PV,       PTENU,    PTENV,
     S     PVERVEL,  PAP,      PAPH,     PGEO,
     S     PGEOH,    LDCUM,    KCBOT,    KCTOP,   KTYPE,
     S     PLU,      PLUDE,    PMFU,     PMFD,
     S     LDSC,     KBOTSC,   PQFLDIA,  PSFLDIA,
C---prognostoc fields
     S     PL,       PI,       PA,
     S     PTENL,    PTENI,    PTENA     ,
C---resulting fluxes
     S     PFPLSL,   PFPLSN,
C---constants
     S     RD, RCPD, RG, RLMLT, RALVDCP, RALSDCP, RALFDCP, RETV,
     S     RTT, RTWAT, RTICE, RTBER )
C
CEVMC-------------------------------------------------------
C---970120                                               ---
C-- original dummy list replaced by minimum list         ---
C-- modifications:                                       ---
C-- 1) resulting (or diagnostic) fluxes taken out        ---
C--    (except precipitative fluxes)                     ---
C-- 2) LDLAND taken out since it is not used             ---
C-- 3) ECMWF-constants initialized from ECHAM constants  ---
C--                                                      ---
C-- further changes: saturation humidity obtained from   ---
C--                  TLUCUA-functions in YOTLUC          ---
C---input
C    S    (KIDIA,    KFDIA,    KHOR,     KLON,     KLEV,
C    S     PTSPHY,
C    S     PT,       PQ,       PTENT,    PTENQ,
C    S     PDTDIAB,  PHRSW,    PHRLW,
C    S     PU,       PV,       PTENU,    PTENV,
C    S     PVERVEL,  PAP,      PAPH,     PGEO,
C    S     PGEOH,    LDCUM,    KCBOT,    KCTOP,   KTYPE,
C    S     PLU,      PLUDE,    PMFU,     PMFD,
C    S     LDLAND,       
C    S     LDSC,     KBOTSC,   PQFLDIA,  PSFLDIA,
CEVMC---prognostoc fields
C    S     PL,       PI,       PA,
C    S     PTENL,    PTENI,    PTENA,
CEVMC---resulting fluxes
C    S     PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG,
C    S     PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN )
CEVMC-------------------------------------------------------
C
C**** *CLOUDSC* -  ROUTINE FOR PARAMATERIZATION OF CLOUD PROCESSES
C                  FOR PROGNOSTIC CLOUD SCHEME
C
C     M.TIEDTKE      E.C.M.W.F.     8/1988, 2/1990
C     CH. JAKOB      E.C.M.W.F.     2/1994 IMPLEMENTATION INTO IFS
C
C     PURPOSE
C     -------
C
C          THIS ROUTINE UPDATES THE CONV/STRAT CLOUD FIELDS.
C          THE FOLLOWING PROCESSES ARE CONSIDERED:
C
C      (1) DETRAINMENT OF CLOUD WATER FROM CONVECTIVE UPDRAFTS
C      (2) FORMATION OF CLOUDS AT TOP OF CONVECTIVE BOUNDARY LAYERS
C      (3) EVAPORATION/CONDENSATION OF CLOUD WATER IN CONNECTION
C          WITH HEATING/COOLING SUCH AS BY SUBSIDENCE/ASCENT
C      (4) EROSION OF CLOUDS BY TURBULENT MIXING OF CLOUD AIR
C          WITH UNSATURATED ENVIRONMENTAL AIR
C      (5) TURBULENT TRANSPORTS OF S,Q,U,V AT CLOUD TOPS DUE TO
C          BUOYANCY FLUXES AND LW RADIATIVE COOLING
C      (6) CONVERSION OF CLOUD WATER INTO PRECIPITATION
C      (7) EVAPORATION OF RAIN
C      (8) MELTING OF SNOW
C
C**   INTERFACE.
C     ----------
C
C          *CLOUDSC* IS CALLED FROM *CALLPAR*
C     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
C     T,Q,L,PHI AND DETRAINMENT OF CLOUD WATER FROM THE
C     CONVECTIVE CLOUDS (MASSFLUX CONVECTION SCHEME), BOUNDARY
C     LAYER TURBULENT FLUXES OF HEAT AND MOISTURE, RADIATIVE FLUXES,
C     OMEGA.
C     IT RETURNS ITS OUTPUT TO:
C      1.MODIFIED TENDENCIES OF MODEL VARIABLES T AND Q
C        AS WELL AS CLOUD VARIABLES L AND C
C      2.GENERATES PRECIPITATION FLUXES FROM STRATIFORM CLOUDS
C
C
C     EXTERNALS.
C     ----------
C
C          NONE
C
C     SWITCHES.
C     --------
C
C          LCLOUD=.TRUE.   PREDICTION OF CONV/STRAT CLOUD-FIELDS
C
C
C     MODEL PARAMETERS
C     ----------------
C
C     RCLDIFF:    PARAMETER FOR EROSION OF CLOUDS
C     RCLCRIT:    THRESHOLD VALUE FOR AUTOCONVERSION
C     RKCONV:    PARAMETER FOR AUTOCONVERSION OF CLOUDS (KESSLER)
C     RCLDMAX:    MAXIMUM POSSIBLE CLW CONTENT (MASON,1971)
C
C     REFERENCE.
C     ----------
C
C          TIEDTKE (1993)
C
CEVM#include "yomcst.h"
CEVM#include "yoethf.h"
CEVM#include "yoecldp.h"
      INCLUDE 'yoecldp.h'
      INCLUDE 'yotluc.h'
C
      LOGICAL LLO1,LLO2,LLO3,LLO4
      LOGICAL LOH,LOM,LOL
C
C
C---input variables
      REAL     PT(KHOR,KLEV),          PQ(KHOR,KLEV),
     S         PAP(KHOR,KLEV),         PAPH(KHOR,KLEV+1),
     S         PGEO(KHOR,KLEV),        PGEOH(KHOR,KLEV+1),
     S         PTENT(KHOR,KLEV),       PTENQ(KHOR,KLEV),
     S         PU(KHOR,KLEV),          PV(KHOR,KLEV),
     S         PTENU(KHOR,KLEV),       PTENV(KHOR,KLEV),
     S         PVERVEL(KHOR,KLEV)
      REAL     PDTDIAB(KHOR,KLEV),
     S         PHRSW(KHOR,KLEV),       PHRLW(KHOR,KLEV)
      REAL     PLU(KHOR,KLEV),         PLUDE(KHOR,KLEV),
     S         PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV)
      REAL     PQFLDIA(KHOR,KLEV+1),   PSFLDIA(KHOR,KLEV+1)
      INTEGER  KCBOT(KHOR),            KCTOP(KHOR)
      INTEGER  KBOTSC(KHOR),           KTYPE(KHOR)
      LOGICAL  LDCUM(KHOR),            LDSC(KHOR)
C
C---prognostic variables 
      REAL     PL(KHOR,KLEV),        PA(KHOR,KLEV),
     S         PI(KHOR,KLEV), 
     S         PTENL(KHOR,KLEV),     PTENA(KHOR,KLEV),
     S         PTENI(KHOR,KLEV)
C---fluxes for budget
      REAL
CEVM S         PFHPSL(KHOR,KLEV+1),    PFHPSN(KHOR,KLEV+1),
     S         PFPLSL(KHOR,KLEV+1),    PFPLSN(KHOR,KLEV+1) 
CEVM S         PFSQLF(KHOR,KLEV+1),    PFSQIF(KHOR,KLEV+1),
CEVM S         PFCQLNG(KHOR,KLEV+1),   PFCQNNG(KHOR,KLEV+1)
C 
C---local fields
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
C---source terms for budget
      REAL     PLCOND1(JPHR,MLEV),     PLCOND2(JPHR,MLEV),
     S         PLEVAP(JPHR,MLEV),      PLEROS(JPHR,MLEV),
     S         PLADJ(JPHR,MLEV),       PLSCGE(JPHR,MLEV),
     S         PDLFLEN(JPHR,MLEV),     PLCVPR(JPHR,MLEV),
     S         PLCVPS(JPHR,MLEV),      PLSEDIM(JPHR,MLEV)
C
      REAL     ZQS(JPHR,MLEV),         ZLCLD(JPHR),
     S         ZAMEAN(JPHR),           ZACOND(JPHR),
     S         ZAEVAP(JPHR),
     S         ZAEROS(JPHR),
     S         ZADETR(JPHR),
     S         ZASCGE(JPHR),
     S         ZREVAP(JPHR),           ZSEVAP(JPHR),
     S         ZDPREC(JPHR),           ZSMELT(JPHR),
     S         ZPREC(JPHR),            ZPRCADJ(JPHR),
     S         ZCOVPC(JPHR),           ZDQSDT(JPHR),
     S         ZTOLD(JPHR),            ZQOLD(JPHR),
     S         ZQMDL(JPHR),            ZPP(JPHR),
     S         ZDTGDP(JPHR),           ZTRPAUS(JPHR)
      REAL     ZDTFORC(JPHR),          ZDTDIAB(JPHR)
      REAL     ZTP1(JPHR,MLEV),        ZQP1(JPHR,MLEV)
      REAL     ZBUOINT(JPHR),          ZWENTR(JPHR),
     S         ZSFLEN(JPHR),
     S         ZLFLEN(JPHR),           
     S         ZLUSC(JPHR),            ZSUMDP(JPHR)
      REAL     ZALUDEM(JPHR),          ZACONDM(JPHR),
     S         ZASCGEM(JPHR),          ZAEVAPM(JPHR),
     S         ZAEROSM(JPHR),          ZLDIFDT(JPHR,MLEV)
C---for flux calculation
      REAL     ZQTORIG(JPHR,MLEV),      ZTTORIG(JPHR,MLEV),
     S         ZUTORIG(JPHR,MLEV),      ZVTORIG(JPHR,MLEV)
      REAL     ZDL(JPHR,MLEV),          ZDA(JPHR,MLEV)
      REAL     ZL(JPHR,MLEV),           ZA(JPHR,MLEV)
      REAL     ZLNEG(JPHR,MLEV)
      LOGICAL  LLFLAG(JPHR)
C
CDIR$ VFUNCTION EXPHF
CEVM#include "fcttre.h"
#ifdef DOC
C
C     *****************************************************************
C
C           CONSIDERATION OF MIXED PHASES
C
C     *****************************************************************
C 
C     FOEALFA is calculated to distinguish the three cases:
C
C			FOEALFA=1            water phase
C			FOEALFA=0            ice phase
C			0 < FOEALFA < 1      mixed phase
C
C		INPUT : PTARG = TEMPERATURE
#endif
C
      FOEALFA (PTARG) = MIN(1.,((MAX(RTICE,MIN(RTWAT,PTARG))-RTICE)
     S  /(RTWAT-RTICE))**2) 
C----------------------------------------------------------------------
C
C
C        0.    SAVE TENDENCIES FOR FLUX CALCULATIONS,
C              INITIALIZE SOURCE TERMS AND CHECK A AND L
C              FOR MINIMUM LIMITS
C
C
   10 CONTINUE
      DO 12 JK=1,KLEV
      DO 11 JL=KIDIA,KFDIA
      ZL(JL,JK)=PL(JL,JK)+PI(JL,JK)+PTSPHY*(PTENL(JL,JK)+PTENI(JL,JK))
      ZA(JL,JK)=PA(JL,JK)+PTSPHY*PTENA(JL,JK)
      ZTP1(JL,JK)=PT(JL,JK)+PTSPHY*PTENT(JL,JK)
      ZQP1(JL,JK)=PQ(JL,JK)+PTSPHY*PTENQ(JL,JK)
      ZQP1(JL,JK)=MAX(ZQP1(JL,JK),0.)
      ZTTORIG(JL,JK)=PTENT(JL,JK)
      ZQTORIG(JL,JK)=PTENQ(JL,JK)
      ZUTORIG(JL,JK)=PTENU(JL,JK)
      ZVTORIG(JL,JK)=PTENV(JL,JK)
      ZDL(JL,JK)=0.
      ZDA(JL,JK)=0.
      PLCOND1(JL,JK)=0.
      PLCOND2(JL,JK)=0.
      PLEVAP(JL,JK)=0.
      PLEROS(JL,JK)=0.
      PLADJ(JL,JK)=0.
      PLSCGE(JL,JK)=0.
      PDLFLEN(JL,JK)=0.
      PLCVPR(JL,JK)=0.
      PLCVPS(JL,JK)=0.
      PLSEDIM(JL,JK)=0.
      IF(.NOT.LDCUM(JL)) THEN
         PLUDE(JL,JK)=0.
         PMFU(JL,JK)=0.
         PMFD(JL,JK)=0.
      ENDIF
      IF(ZA(JL,JK).LT.RAMIN.OR.ZL(JL,JK).LT.RLMIN) THEN
         ZDL(JL,JK)=-ZL(JL,JK)
         ZDA(JL,JK)=-ZA(JL,JK)
         ZA(JL,JK)=0.
         ZL(JL,JK)=0.
      ENDIF
      IF(ZA(JL,JK).GT.1.) THEN
         ZDA(JL,JK)=1.-ZA(JL,JK)
         ZA(JL,JK)=1.
      ENDIF
      ZLNEG(JL,JK)=0.0
   11 CONTINUE
   12 CONTINUE
C
C
C----------------------------------------------------------------------
C
C
C
C        1.       CONSTANTS AND PARAMETERS
C                 CALCULATE L IN UPDRAFTS OF BL-CLOUDS
C                 SPECIFY QS, P/PS FOR TROPOPAUSE (FOR C2)
C                 AND INITIALIZE VARIABLES
C                 -----------------------------------
C
C
C
  100 CONTINUE
C
      ZQTMST=1./PTSPHY
      ZCPDG=RCPD/RG
      ZQCP=1./RCPD
      ZQCPDT=1./(RCPD*PTSPHY)
      ZGDCP=RG/RCPD
      ZRDCP=RD/RCPD
      ZCONS1A=RCPD/(RLMLT*RG*RTAUMEL)
      ZCONS2=1./(RG*PTSPHY)
      ZCONS3=0.5*RCPD*RETV
      ZEPSEC=1.E-10
      ZCKCODT=RKCONV*PTSPHY
C
      DO 120 JL=KIDIA,KFDIA
      IF(LDCUM(JL)) LDSC(JL)=.FALSE.
      IF(LDSC(JL)) THEN
         IK=KBOTSC(JL)-1
         ZQPH=1./PAP(JL,IK)
         ZQADIAB=ZQP1(JL,KLEV)
         ZTADIAB=(RCPD*ZTP1(JL,KLEV)+PGEO(JL,KLEV)-PGEO(JL,IK))*ZQCP
         IT=NINT(ZTADIAB*1000.)
	 ZQSAT=TLUCUA(IT)*ZQPH
CEVM     ZQSAT=FOEEWM(ZTADIAB)*ZQPH
         ZCOR=1./(1.-RETV*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZCOND=(ZQADIAB-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
CEVM     ZCOND=(ZQADIAB-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(ZTADIAB))
         ZTNEW=ZTADIAB+TLUCUC(IT)*ZCOND
CEVM     ZTNEW=ZTADIAB+FOELDCPM(ZTADIAB)*ZCOND
         IT=NINT(ZTNEW*1000.)
	 ZQSAT=TLUCUA(IT)*ZQPH
CEVM     ZQSAT=FOEEWM(ZTNEW)*ZQPH
         ZCOR=1./(1.-RETV*ZQSAT)
         ZQSAT=ZQSAT*ZCOR
         ZLUSC(JL)=ZQADIAB-ZQSAT
         IF(ZLUSC(JL).LE.1.E-10) LDSC(JL)=.FALSE.
      END IF
  120 CONTINUE
C
      DO 125 JL=KIDIA,KFDIA
      ZTRPAUS(JL)=0.1
  125 CONTINUE
      DO 130 JK=1,KLEV
      DO 130 JL=KIDIA,KFDIA
      ZLDIFDT(JL,JK)=RCLDIFF*PTSPHY
      IF(KTYPE(JL).EQ.2..AND.PLUDE(JL,JK).GT.0.)
     S       ZLDIFDT(JL,JK)=5.*ZLDIFDT(JL,JK)
      IT=NINT(ZTP1(JL,JK)*1000.)
      ZQS(JL,JK)=TLUCUA(IT)/PAP(JL,JK)
CEVM  ZQS(JL,JK)=FOEEWM(ZTP1(JL,JK))/PAP(JL,JK)
      ZQS(JL,JK)=MIN(0.5,ZQS(JL,JK))
      ZQS(JL,JK)=ZQS(JL,JK)/(1.-RETV*ZQS(JL,JK))
  130 CONTINUE
      DO 140 JK=1,KLEV-1
      DO 140 JL=KIDIA,KFDIA
      ZSIG=PAP(JL,JK)/PAPH(JL,KLEV+1)
      LLO1=ZSIG.GT.0.1.AND.ZSIG.LT.0.4.AND.
     S     ZTP1(JL,JK).GT.ZTP1(JL,JK+1)
      IF(LLO1) THEN
         ZTRPAUS(JL)=ZSIG
      END IF
  140 CONTINUE
C
      DO 160 JL=KIDIA,KFDIA
      ZPREC(JL)=0.
      ZCOVPC(JL)=0.
      ZALUDEM(JL)=0.
      ZACONDM(JL)=0.
      ZASCGEM(JL)=0.
      ZAEVAPM(JL)=0.
      ZAEROSM(JL)=0.
      IF(.NOT.LDSC(JL)) ZLUSC(JL)=0.
      PFPLSL(JL,1)=0.
      PFPLSN(JL,1)=0.
  160 CONTINUE
C
C----------------------------------------------------------------------
C
C
C
C        2.       CALCULATE VERTICAL INTEGRATED BUOYANCY FLUXES IN
C                 CONVECTIVE BOUNDARY LAYER FOR TOP ENTRAINMENT IN SC
C                 ---------------------------------------------------
C
C
C
  200 CONTINUE
      DO 210 JL=KIDIA,KFDIA
      ZBUOINT(JL)=-(PAPH(JL,KLEV+1)-PAP(JL,KLEV))*
     S  (PSFLDIA(JL,KLEV+1)+RCPD*ZTP1(JL,KLEV)*RETV*PQFLDIA(JL,KLEV+1))
      ZSUMDP(JL)=PAPH(JL,KLEV+1)-PAP(JL,KLEV)
  210 CONTINUE
      DO 230 JK=KLEV,2,-1
      DO 220 JL=KIDIA,KFDIA
      LLO1=LDSC(JL).AND.JK.GE.KBOTSC(JL)
      IF(LLO1) THEN
         ZBUOINT(JL)=ZBUOINT(JL)-(PAP(JL,JK)-PAP(JL,JK-1))*
     S     (PSFLDIA(JL,JK)+RCPD*0.5*(ZTP1(JL,JK)+ZTP1(JL,JK-1))*
     S                    RETV*PQFLDIA(JL,JK))
         ZSUMDP(JL)=ZSUMDP(JL)+PAP(JL,JK)-PAP(JL,JK-1)
      END IF
  220 CONTINUE
  230 CONTINUE
      DO 240 JL=KIDIA,KFDIA
      ZBUOINT(JL)=ZBUOINT(JL)/ZSUMDP(JL)
  240 CONTINUE
C----------------------------------------------------------------------
C
C
C
C        3.       START CLOUD CALCULATIONS WITH OUTER LOOP OVER LEVELS
C                 ----------------------------------------------------
C
C
C
  300 CONTINUE
C
      DO 590 JK=1,KLEV
C
C
C       INITIALIZE VARIABLES
C       --------------------
C
C
      DO 304 JL=KIDIA,KFDIA
      PFPLSL(JL,JK+1)=PFPLSL(JL,JK)
      PFPLSN(JL,JK+1)=PFPLSN(JL,JK)
      ZLCLD(JL)=0.
      ZPRCADJ(JL)=0.
      ZREVAP(JL)=0.
      ZSEVAP(JL)=0.
      ZACOND(JL)=0.
      ZAEROS(JL)=0.
      ZAEVAP(JL)=0.
      ZADETR(JL)=0.
      ZASCGE(JL)=0.
      ZDPREC(JL)=0.
      ZSMELT(JL)=0.
      ZSMELT(JL)=0.
      ZAMEAN(JL)=0.
      ZWENTR(JL)=0.
      ZSFLEN(JL)=0.
      ZLFLEN(JL)=0.
      ZDTFORC(JL)=0.
      ZDTDIAB(JL)=0.
      ZDTGDP(JL)=PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZLCLD(JL)=ZL(JL,JK)/MAX(ZA(JL,JK),ZEPSEC)
  304 CONTINUE
C----------------------------------------------------------------------
C
C
C        3.1      DETRAINMENT OF CLOUD AIR FROM CONVECTIVE UPDRAFTS
C                 -------------------------------------------------
C
C
  310 CONTINUE
      IF(JK.LT.KLEV) THEN
         DO 315 JL=KIDIA,KFDIA
         PLUDE(JL,JK)=PLUDE(JL,JK)*ZDTGDP(JL)
         LLO1=LDCUM(JL).AND.PLUDE(JL,JK).GT.RLMIN.AND.PLU(JL,JK+1)
     1         .GT.ZEPSEC
         IF(LLO1) THEN
            ZADETR(JL)=PLUDE(JL,JK)/PLU(JL,JK+1)
         ELSE
            PLUDE(JL,JK)=0.
         END IF
  315    CONTINUE
      END IF
C----------------------------------------------------------------------
C
C
C        3.2      ENTRAINMENT VELOCITY DUE TO LW-RADIATIVE COOLING
C                 ------------------------------------------------
C
C
  320 CONTINUE
      IF(JK.GT.1) THEN
         DO 322 JL=KIDIA,KFDIA
         LLO1=ZA(JL,JK).GT.0..AND.ZA(JL,JK-1).EQ.0.
         IF(LLO1) THEN
           ZWENTR(JL)=
     S      -RENTRRA*PHRLW(JL,JK)*(PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZCPDG
            ZWENTR(JL)=MAX(0.,ZWENTR(JL))
         END IF
  322    CONTINUE
      END IF
C----------------------------------------------------------------------
C
C
C        3.3    GENERATION OF CLOUDS AT TOP OF CONVECTIVE BOUNDARY LAYER
C               --------------------------------------------------------
C
C
  330 CONTINUE
      IF(JK.LT.KLEV) THEN
         DO 332 JL=KIDIA,KFDIA
         IK=KBOTSC(JL)-1
         LLFLAG(JL)=LDSC(JL).AND.JK.EQ.IK
         IF(LLFLAG(JL)) THEN
            ZQK=MIN(ZQP1(JL,JK),ZQS(JL,JK))
            ZQKLEV=MIN(ZQP1(JL,KLEV),ZQS(JL,KLEV))
            ZDQBASE=ZQKLEV-
     S      ((1.-ZA(JL,JK))*ZQK+ZA(JL,JK)*(ZQS(JL,JK)+ZLCLD(JL)))
            LLO2=ZDQBASE.GT.1.E-6.AND.PQFLDIA(JL,JK+1).LT.0.
     S      .AND.ZQK.GT.0.2*ZQS(JL,JK)
            IF(LLO2) THEN
               ZMFLBAS=-PQFLDIA(JL,JK+1)/ZDQBASE
               PLSCGE(JL,JK)=ZMFLBAS*ZDTGDP(JL)
               ZASCGE(JL)=ZMFLBAS*ZDTGDP(JL)
               PLSCGE(JL,JK)=MIN(PLSCGE(JL,JK),
     S          0.5*ZQP1(JL,IK)/MAX(ZLUSC(JL),ZEPSEC))
               ZASCGE(JL)=MIN(ZASCGE(JL),1.)
            END IF
            LLO1=ZA(JL,JK).GT.0..AND.ZA(JL,JK-1).EQ.0.
            IF(LLO1) THEN
               ZWENTR(JL)=ZWENTR(JL)+RENTRTU*ZBUOINT(JL)
               ZWENTR(JL)=MAX(0.,ZWENTR(JL))
            END IF
         END IF
  332    CONTINUE
      END IF
         DO 334 JL=KIDIA,KFDIA
         LLO1=ZWENTR(JL).GT.0.
         IF(LLO1) THEN
            ZSD=RCPD*ZTP1(JL,JK)+PGEO(JL,JK)
            ZSU=RCPD*ZTP1(JL,JK-1)+PGEO(JL,JK-1)
            ZDSV=ZSU-ZSD+ZCONS3*(ZTP1(JL,JK-1)+ZTP1(JL,JK))*
     S                           (ZQP1(JL,JK-1)-ZQP1(JL,JK))
            LLO2=ZDSV.GT.ZEPSEC.AND.ZQP1(JL,JK-1).LT.0.9*ZQS(JL,JK-1)
            IF(LLO2) THEN
               ZWENTR(JL)=ZWENTR(JL)/ZDSV
               ZWENTR(JL)=MIN(1.,ZWENTR(JL))
               ZRHO=PAP(JL,JK)/(RD*ZTP1(JL,JK))
               ZWENTR(JL)=ZRHO*ZWENTR(JL)
               ZSFLEN(JL)=-ZWENTR(JL)*MAX(0.,ZSU-ZSD)
               ZLFLEN(JL)= ZWENTR(JL)*ZL(JL,JK)
               PDLFLEN(JL,JK)=ZLFLEN(JL)*ZDTGDP(JL)
            ELSE
               ZWENTR(JL)=0.
            END IF
         END IF
  334 CONTINUE
C----------------------------------------------------------------------
C
C
C        3.4      EROSION OF CLOUDS BY TURBULENT MIXING:
C                 NOTE: THIS PROCESS DECREASES THE CLOUD AREA BUT
C                       LEAVES THE SPECIFIC CLOUD WATER CONTENT
C                       WITHIN CLOUDS UNCHANGED
C                 ---------------------------------------------
C
C
  340 CONTINUE
         DO 345 JL=KIDIA,KFDIA
         IF(ZLCLD(JL).GT.ZEPSEC) THEN
            ZE=ZLDIFDT(JL,JK)*MAX(ZQS(JL,JK)-ZQP1(JL,JK),0.)
            PLEROS(JL,JK)=MIN(ZE,ZLCLD(JL))
            ZAEROS(JL)=ZE/ZLCLD(JL)
         END IF
  345    CONTINUE
C----------------------------------------------------------------------
C
C
C        3.5      CONDENSATION/EVAPORATION DUE TO DQSAT/DT
C                 ----------------------------------------
C
C
  350 CONTINUE
C
C                 CALCULATE FIRST DQS/DT
C                 ----------------------
C
      DO 351 JL=KIDIA,KFDIA
      ZDTDP=ZRDCP*ZTP1(JL,JK)/PAP(JL,JK)
      ZDPMXDT=(PAPH(JL,JK+1)-PAPH(JL,JK))*ZQTMST
      ZMFDN=0.
      IF(JK.LT.KLEV) ZMFDN=PMFU(JL,JK+1)+PMFD(JL,JK+1)
      ZWTOT=PVERVEL(JL,JK)+0.5*RG*(PMFU(JL,JK)+PMFD(JL,JK)+ZMFDN)
      IF(PMFU(JL,JK).GT.0..AND.ZWTOT.LT.0.) ZWTOT=0.
      ZWTOT=MIN(ZDPMXDT,MAX(-ZDPMXDT,ZWTOT))
      ZDTFLEN=-ZSFLEN(JL)*ZDTGDP(JL)
      ZZZDT=PDTDIAB(JL,JK)+PHRSW(JL,JK)+PHRLW(JL,JK)
      ZDTDIAB(JL)=MIN(ZDPMXDT*ZDTDP,MAX(-ZDPMXDT*ZDTDP,ZZZDT))
     S            *PTSPHY+ZDTFLEN*ZQCP
      ZDTFORC(JL)=ZDTDP*ZWTOT*PTSPHY+ZDTDIAB(JL)
      ZQOLD(JL)=ZQS(JL,JK)
      ZTOLD(JL)=ZTP1(JL,JK)
      ZTP1(JL,JK)=ZTP1(JL,JK)+ZDTFORC(JL)
      ZPP(JL)=PAP(JL,JK)
      LLFLAG(JL)=.TRUE.
  351 CONTINUE
C
      IK=JK
      ICALL=4
      CALL CUADJTQ4
     S    (KIDIA,    KFDIA,    KHOR,      KLON,     KLEV,     
     S     IK,
     S     ZPP,      ZTP1,       ZQS,      LLFLAG,   ICALL)
C
      DO 352 JL=KIDIA,KFDIA
      ZDQSDT(JL)=ZQS(JL,JK)-ZQOLD(JL)
      ZQS(JL,JK)=ZQOLD(JL)
      ZTP1(JL,JK)=ZTOLD(JL)
  352 CONTINUE
C
C                 ZDQSDT(JL) > 0:  EVAPORATION OF CLOUDS
C                 NOTE: THIS PROCESS DECREASES CLW CONTENT
C                       BUT LEAVES CLOUD AREA UNCHANGED
C                 ----------------------------------------
C
         DO 355 JL=KIDIA,KFDIA
         PLEVAP(JL,JK)=MAX(ZDQSDT(JL),0.)
         PLEVAP(JL,JK)=MIN(PLEVAP(JL,JK),ZLCLD(JL))
  355    CONTINUE
C
C                 ZDQSDT(JL) < 0: FORMATION OF CLOUDS
C                 (1) INCREASE OF LWC IN EXISTING CLOUDS
C                 --------------------------------------
C
         DO 356 JL=KIDIA,KFDIA
         LLO1=ZA(JL,JK).GT.ZEPSEC.AND.ZDQSDT(JL).LE.-RLMIN
         IF(LLO1) THEN
            PLCOND1(JL,JK)=-ZDQSDT(JL)
            IF(ZA(JL,JK).GT.0.99) THEN
               ZCOR=1./(1.-RETV*ZQS(JL,JK))
	       IT=NINT(ZTP1(JL,JK)*1000.)
               ZCDMAX=(ZQP1(JL,JK)-ZQS(JL,JK))/
     S                (1.+ZCOR*ZQS(JL,JK)*TLUCUB(IT))
CEVM S                (1.+ZCOR*ZQS(JL,JK)*FOEDEM(ZTP1(JL,JK)))
            ELSE
               ZCDMAX=(ZQP1(JL,JK)-ZA(JL,JK)*ZQS(JL,JK))/ZA(JL,JK)
            END IF
            PLCOND1(JL,JK)=MAX(MIN(PLCOND1(JL,JK),ZCDMAX),0.)
         END IF
  356    CONTINUE
C
C                 ZDQSDT(JL) < 0: FORMATION OF CLOUDS
C                 (2) GENERATION OF NEW CLOUDS (DA/DT>0)
C                 --------------------------------------
C
         DO 357 JL=KIDIA,KFDIA
         LLO2=ZDQSDT(JL).LE.-RLMIN.AND.JK.LE.KBOTSC(JL)-1
         IF(LLO2) THEN
            ZRHC=RAMID
            ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
            IF(ZSIGK.GT.0.8) THEN
               ZRHC=RAMID+(1.-RAMID)*((ZSIGK-0.8)/0.2)**2
            END IF
            ZBOTT=ZTRPAUS(JL)+0.2
            IF(ZSIGK.LT.ZBOTT) THEN
               ZRHC=RAMID+(1.-RAMID)*MIN(((ZBOTT-ZSIGK)/0.2)**2,1.)
            END IF
            LLO3=ZQP1(JL,JK).GE.ZRHC*ZQS(JL,JK)
            IF(LLO3) THEN
                ZACOND(JL)=-ZDQSDT(JL)/
     S           MAX(ZQS(JL,JK)-ZQP1(JL,JK),ZEPSEC)
                ZACOND(JL)=MIN(ZACOND(JL),1.)
                PLCOND2(JL,JK)=-ZDQSDT(JL)
                PLCOND2(JL,JK)=MIN(PLCOND2(JL,JK),
     S           (ZQP1(JL,JK)-ZA(JL,JK)*ZQS(JL,JK))/
     S           MAX(ZACOND(JL),ZEPSEC))
                PLCOND2(JL,JK)=MAX(PLCOND2(JL,JK),0.)
                IF(PLCOND2(JL,JK).EQ.0.) ZACOND(JL)=0.
            END IF
         END IF
  357 CONTINUE
C----------------------------------------------------------------------
C
C
C        3.6      INTEGRATE ANALYTICALLY TENDENCY EQUATION FOR
C                 CLOUD COVER A AND SPECIFY TERMS IN TENDENCY
C                 EQUATION OF CLOUD WATER L ACCORDINGLY.
C                 ------------------------------------------------
C
C
  360 CONTINUE
      DO 362 JL=KIDIA,KFDIA
      ZAA=ZADETR(JL)+ZASCGE(JL)+ZACOND(JL)
      ZB=ZAEROS(JL)
      ZINT=EXP(-(ZAA+ZB))
      ZANEW=MIN(1.,ZAA/MAX(ZAA+ZB,ZEPSEC))*(1.-ZINT)+ZA(JL,JK)*ZINT
      ZANEW=MIN(1.,MAX(0.,ZANEW))
      ZDA(JL,JK)=ZANEW-ZA(JL,JK)
      ZAMEAN(JL)=0.5*(ZANEW+ZA(JL,JK))
  362 CONTINUE
      DO 363 JL=KIDIA,KFDIA
      ZADETR(JL)=(1.-ZAMEAN(JL))*ZADETR(JL)
      ZASCGE(JL)=(1.-ZAMEAN(JL))*ZASCGE(JL)
      ZACOND(JL)=(1.-ZAMEAN(JL))*ZACOND(JL)
      ZAEROS(JL)=ZAMEAN(JL)*ZAEROS(JL)
      ZALUDEM(JL)=MAX(ZALUDEM(JL),ZADETR(JL))
      ZASCGEM(JL)=MAX(ZASCGEM(JL),ZASCGE(JL))
      ZACONDM(JL)=MAX(ZACONDM(JL),ZACOND(JL))
      ZAEROSM(JL)=MAX(ZAEROSM(JL),ZAEROS(JL))
  363 CONTINUE
      DO 364 JL=KIDIA,KFDIA
      PLSCGE(JL,JK)=PLSCGE(JL,JK)*(ZLUSC(JL)-ZA(JL,JK)*ZLCLD(JL))
      IF(PLSCGE(JL,JK).LT.RLMIN) PLSCGE(JL,JK)=0.
      PLCOND1(JL,JK)=ZA(JL,JK)*PLCOND1(JL,JK)
      IF(PLCOND1(JL,JK).LT.RLMIN) PLCOND1(JL,JK)=0.
      ZDACOND=(1.-(1.-ZA(JL,JK))*EXP(-ZACOND(JL)))-ZA(JL,JK)
      PLCOND2(JL,JK)=ZDACOND*PLCOND2(JL,JK)
      IF(PLCOND2(JL,JK).LT.RLMIN) PLCOND2(JL,JK)=0.
      PDLFLEN(JL,JK)=MIN(PDLFLEN(JL,JK),ZL(JL,JK))
      PLEROS(JL,JK)=MIN(ZA(JL,JK)*PLEROS(JL,JK),ZL(JL,JK)
     S              -PDLFLEN(JL,JK))
      PLEVAP(JL,JK)=MIN(ZA(JL,JK)*PLEVAP(JL,JK),ZL(JL,JK)
     S              -PDLFLEN(JL,JK)-PLEROS(JL,JK))
      ZDL(JL,JK)=ZDL(JL,JK)+PLUDE(JL,JK)+PLSCGE(JL,JK)
     S           +PLCOND1(JL,JK)+PLCOND2(JL,JK)-PLEVAP(JL,JK)
     S           -PLEROS(JL,JK)-PDLFLEN(JL,JK)
      LLO1=ZL(JL,JK)+ZDL(JL,JK).LE.RLMIN
      IF(LLO1) THEN
         ZDL(JL,JK)=-ZL(JL,JK)
         ZAEVAP(JL)=MAX(0.,ZA(JL,JK)+ZDA(JL,JK))
         ZAEVAPM(JL)=MAX(ZAEVAPM(JL),ZAEVAP(JL))
         ZDA(JL,JK)=-ZA(JL,JK)
      END IF
  364 CONTINUE
C----------------------------------------------------------------------
C
C
C        3.7      FORMATION OF GRID SCALE CLOUD (A=1) IF Q > QS
C                 CHECK FOR SUPERSATURATION AND ADJUST T AND Q.
C                 NOTE THAT THIS ADJUSTMENT IS ONLY TO AVOID
C                 SUPERSATURATED STATES WHICH ARE NOT DUE TO
C                 ADIABATIC AND DIABATIC PROCESSES BUT DUE TO
C                 NUMERICAL PROBLEMS.
C                 ---------------------------------------------
C
C
  370 CONTINUE
      DO 372 JL=KIDIA,KFDIA
      ZZDL=(PLSCGE(JL,JK)+PLCOND1(JL,JK)+PLCOND2(JL,JK)
     S    -PLEVAP(JL,JK)-PLEROS(JL,JK))
      ZQOLD(JL)=ZQP1(JL,JK)
      ZTOLD(JL)=ZTP1(JL,JK)
      ZQMDL(JL)=ZQP1(JL,JK)-ZZDL
      ZQP1(JL,JK)=ZQP1(JL,JK)-ZZDL
      IT=NINT(ZTP1(JL,JK)*1000.)
      ZTP1(JL,JK)=ZTP1(JL,JK)+TLUCUC(IT)*ZZDL
CEVM  ZTP1(JL,JK)=ZTP1(JL,JK)+FOELDCPM(ZTP1(JL,JK))*ZZDL
      ZPP(JL)=PAP(JL,JK)
      LLFLAG(JL)=.TRUE.
  372 CONTINUE
C
      IK=JK
      ICALL=4
      CALL CUADJTQ4
     S    (KIDIA,    KFDIA,    KHOR,      KLON,     KLEV,     
     S     IK,
     S     ZPP,      ZTP1,       ZQP1,       LLFLAG,   ICALL)
C
      DO 374 JL=KIDIA,KFDIA
      ZZDQAD=MAX(ZQMDL(JL)-ZQP1(JL,JK),0.)
      ZQP1(JL,JK)=ZQOLD(JL)
      ZTP1(JL,JK)=ZTOLD(JL)
      ZPRCADJ(JL)=ZZDQAD
      PLCOND1(JL,JK)=PLCOND1(JL,JK)+ZPRCADJ(JL)
  374 CONTINUE
C----------------------------------------------------------------------
C
C
C        3.8      CONVERSION OF CLOUD WATER INTO PRECIPITATION.
C                 CLOUD PHYSICS FOLLOWING H.SUNDQUIST (1988)
C                 FOR WARM CLOUDS AND MIXED WATER/ICE CLOUDS
C                 AND HEYMSFIELD AND DONNER(1990) FOR ICE CLOUDS.
C                 NOTE THAT WE DISTINGUISH ACTUAL CLOUD COVER
C                 AND FRACTIONAL COVER OF PRECIPITATION AREA.
C                 WE ASSUME THAT AT LOWEST MODEL LEVEL ALL CONDENSED
C                 CLOUD WATER IS CONVERTED INTO PRECIPITATION
C                 --------------------------------------------
C
C
  380 CONTINUE
C
C                 SET FRACTIONAL COVER FOR PRECIPITATING CLOUDS
C                 ---------------------------------------------
C
      DO 381 JL=KIDIA,KFDIA
      LLFLAG(JL)=ZL(JL,JK)+ZDL(JL,JK).GT.0..AND.ZAMEAN(JL).GT.ZEPSEC
      LLO1=LLFLAG(JL).AND.ZCOVPC(JL).EQ.0.
      IF(LLO1) THEN
         ZCOVPC(JL)=ZAMEAN(JL)
      END IF
  381 CONTINUE
C
C
C     1)         ICE CLOUDS: SEDIMENTATION OF ICE
C                --------------------------------
C
C
      DO 382 JL=KIDIA,KFDIA
      IF(ZTP1(JL,JK).LE.RTICE) THEN
         ZPRCGEN=0.
         IF(LLFLAG(JL)) THEN
            ZCLDOLD=ZL(JL,JK)/MAX(ZAMEAN(JL),ZEPSEC)
            ZC=ZDL(JL,JK)/MAX(ZAMEAN(JL),ZEPSEC)
            ZRHO=PAP(JL,JK)/(RD*ZTP1(JL,JK))
            ZVICE=1.7*(ZRHO*(ZCLDOLD+ZC))**0.07
     S         +851.7*(ZRHO*(ZCLDOLD+ZC))
            ZD=ZRHO*ZVICE*ZDTGDP(JL)
            ZINT=EXP(-ZD)
            ZCLDNEW=ZCLDOLD*ZINT+ZC/ZD*(1.-ZINT)
            ZCLDNEW=MAX(0.,ZCLDNEW)
            ZCLDNEW=MIN(RCLDMAX,ZCLDNEW)
            ZLSEDIM=MAX(0.,ZCLDOLD+ZC-ZCLDNEW)
            ZDL(JL,JK)=ZDL(JL,JK)-ZLSEDIM*ZAMEAN(JL)
            PLSEDIM(JL,JK)=-ZLSEDIM*ZAMEAN(JL)
            ZPRCGEN=ZLSEDIM*ZAMEAN(JL)
         END IF
         ZDPREC(JL)=(ZPRCGEN+ZPRCADJ(JL))/ZDTGDP(JL)
         ZPREC(JL)=ZPREC(JL)+ZDPREC(JL)
         PLCVPS(JL,JK)=ZPRCGEN+ZPRCADJ(JL)
         PFPLSN(JL,JK+1)=PFPLSN(JL,JK+1)+ZDPREC(JL)
      END IF
  382 CONTINUE
C
C
C
C
C     2)          WARM CLOUDS AND MIXED WATER/ICE CLOUDS
C                 --------------------------------------
C
      DO 384 JL=KIDIA,KFDIA
      IF(ZTP1(JL,JK).GT.RTICE) THEN
         ZPRCGEN=0.
         LLO1=LLFLAG(JL).AND.ZCOVPC(JL).GT.ZEPSEC
         IF(LLO1) THEN
C
C           A) PARAMETERS FOR BERGERON-FINDEISEN PROCESS ( T < -5C )
C
            ZDT=MIN(RTBER-RTICE,MAX(RTBER-ZTP1(JL,JK),0.))
            ZCBF=1.+RPRC2*SQRT(ZDT)
            ZZCO=ZCKCODT*ZCBF
            ZLCRIT=RCLCRIT/ZCBF
C
C           B) PARAMETERS FOR CLOUD COLLECTION BY RAIN DROPS
C
            ZPRECIP=(PFPLSL(JL,JK+1)+PFPLSN(JL,JK+1))/ZCOVPC(JL)
            ZCFPR=1.+RPRC1*SQRT(MAX(ZPRECIP,0.))
            ZZCO=ZZCO*ZCFPR
            ZLCRIT=ZLCRIT/ZCFPR
C
C           C) ANALYTIC INTEGRATION OF EQUATION FOR L
C
            ZCLDOLD=ZL(JL,JK)/MAX(ZAMEAN(JL),ZEPSEC)
            ZC=ZDL(JL,JK)/MAX(ZAMEAN(JL),ZEPSEC)
            ZD=ZZCO*(1.-EXP(-((ZCLDOLD+ZC)/ZLCRIT)**2))
            ZINT=EXP(-ZD)
            ZCLDNEW=ZCLDOLD*ZINT+ZC/ZD*(1.-ZINT)
            ZCLDNEW=MAX(0.,ZCLDNEW)
            ZCLDNEW=MIN(RCLDMAX,ZCLDNEW)
C
C           D) CALCULATE PRECIPITATION RATES
C
            ZPR=(ZCLDOLD+ZC-ZCLDNEW)*ZAMEAN(JL)
            ZPR=MIN(ZPR,ZL(JL,JK)+ZDL(JL,JK))
            ZPR=MAX(ZPR,0.)
            ZDL(JL,JK)=ZDL(JL,JK)-ZPR
            ZPRCGEN=ZPR
         END IF
         ZDPREC(JL)=(ZPRCGEN+ZPRCADJ(JL))/ZDTGDP(JL)
         ZPREC(JL)=ZPREC(JL)+ZDPREC(JL)
         ZALFAW=FOEALFA(ZTP1(JL,JK))
         ZALFAI=1.-ZALFAW
         PLCVPR(JL,JK)=ZALFAW*(ZPRCGEN+ZPRCADJ(JL))
         PLCVPS(JL,JK)=(1.-ZALFAW)*(ZPRCGEN+ZPRCADJ(JL))
         ZDPR=ZALFAW*ZDPREC(JL)
         ZDPS=(1.-ZALFAW)*ZDPREC(JL)
         PFPLSL(JL,JK+1)=PFPLSL(JL,JK+1)+ZDPR
         PFPLSN(JL,JK+1)=PFPLSN(JL,JK+1)+ZDPS
      END IF
  384 CONTINUE
C
C                 UPDATE FRACTIONAL COVER FOR PRECIPITATING CLOUDS
C                 ------------------------------------------------
C
      DO 386 JL=KIDIA,KFDIA
      IF(ZPREC(JL).GT.ZEPSEC) THEN
         ZDPRADJ=ZPRCADJ(JL)/ZDTGDP(JL)
         ZCOVPC(JL)=ZCOVPC(JL)+
     1   MAX(0.,(ZAMEAN(JL)*ZDPREC(JL)+ZDPRADJ*(1.-ZAMEAN(JL))-
     1           ZCOVPC(JL)*ZDPREC(JL))/ZPREC(JL))
         ZCOVPC(JL)=MIN(1.,ZCOVPC(JL))
      END IF
  386 CONTINUE
C----------------------------------------------------------------------
C
C
C        3.9      MELTING OF SNOW AND
C                 EVAPORATION OF RAIN/SNOW
C                 NOTE THAT QSAT CHANGES DUE TO MELTING
C                 -------------------------------------
C
C
  390 CONTINUE
         DO 392 JL=KIDIA,KFDIA
         LLO1=PFPLSN(JL,JK+1).GT.0..AND.ZTP1(JL,JK).GT.RTT
         IF(LLO1) THEN
            ZCONS1=ZCONS1A*(1.+0.5*(ZTP1(JL,JK)-RTT))
            ZFAC=ZCONS1*(PAPH(JL,JK+1)-PAPH(JL,JK))
            ZSNMLT=MIN(PFPLSN(JL,JK+1),ZFAC*(ZTP1(JL,JK)-RTT))
            ZSMELT(JL)=ZSNMLT*ZDTGDP(JL)
            PFPLSN(JL,JK+1)=PFPLSN(JL,JK+1)-ZSNMLT
            PFPLSL(JL,JK+1)=PFPLSL(JL,JK+1)+ZSNMLT
            PFPLSL(JL,JK+1)=MAX(0.,PFPLSL(JL,JK+1))
            PFPLSN(JL,JK+1)=MAX(0.,PFPLSN(JL,JK+1))
            ZTP1(JL,JK)=ZTP1(JL,JK)-ZSNMLT/ZFAC
            IT=NINT(ZTP1(JL,JK)*1000.)
            ZQS(JL,JK)=TLUCUA(IT)/PAP(JL,JK)
CEVM        ZQS(JL,JK)=FOEEWM(ZTP1(JL,JK))/PAP(JL,JK)
         END IF
  392 CONTINUE
      DO 394 JL=KIDIA,KFDIA
      ZZRH=0.8
      IF(LDCUM(JL)) ZZRH=0.7
      LLO1=ZPREC(JL).GT.ZEPSEC.AND.ZCOVPC(JL).GT.ZAMEAN(JL).AND.
     S     ZQP1(JL,JK).LT.ZZRH*ZQS(JL,JK)
      IF(LLO1) THEN
         ZPROLD=ZPREC(JL)
         ZDPR=RPECONS*MAX(0.,ZQS(JL,JK)-ZQP1(JL,JK))*ZCOVPC(JL)*
     S        (SQRT(PAP(JL,JK)/PAPH(JL,KLEV+1))/5.09E-3*
     S         ZPROLD/MAX(ZCOVPC(JL),ZEPSEC))**0.5777*
     S        (PAPH(JL,JK+1)-PAPH(JL,JK))
         ZDPR=MIN(ZDPR,MAX(0.,ZZRH*ZQS(JL,JK)-ZQP1(JL,JK))
     S         *ZCOVPC(JL)*ZCONS2*(PAPH(JL,JK+1)-PAPH(JL,JK)))
         ZDPR=MIN(ZDPR,ZPROLD)
         ZDPEVAP=MAX(0.,ZDPR*
     S         MAX((ZCOVPC(JL)-ZAMEAN(JL))/MAX(ZCOVPC(JL),ZEPSEC),0.))
         ZPRNEW=ZPREC(JL)-ZDPEVAP
         ZDPRSFL=ZDPEVAP*PFPLSL(JL,JK+1)*(1./ZPROLD)
         PFPLSL(JL,JK+1)=PFPLSL(JL,JK+1)-ZDPRSFL
         ZREVAP(JL)=ZDPRSFL*ZDTGDP(JL)
         ZDPSSFL=ZDPEVAP*PFPLSN(JL,JK+1)*(1./ZPROLD)
         PFPLSN(JL,JK+1)=PFPLSN(JL,JK+1)-ZDPSSFL
         ZSEVAP(JL)=ZDPSSFL*ZDTGDP(JL)
         ZPREC(JL)=ZPRNEW
      END IF
  394 CONTINUE
      DO 396 JL=KIDIA,KFDIA
         IF(ZPREC(JL).EQ.0.) ZCOVPC(JL)=0.
  396 CONTINUE
C----------------------------------------------------------------------
C
C
C
C        4.0      UPDATE FIELDS OF L,A AND TENDENCIES OF T AND Q
C                 1. DUE TO CLOUD PROCESSES AND
C                 2. ENTRAINMENT PROCESSES (BUOYANCY + RADIATION)
C                        (NOTE IMPLICIT CALCUCLATION)
C                 -------------------------------------------------
C
C
C
  400 CONTINUE
      DO 410 JL=KIDIA,KFDIA
      ZLNEW=(PL(JL,JK)+PI(JL,JK))+ZDL(JL,JK)+(PTENL(JL,JK)+
     S PTENI(JL,JK))*PTSPHY
      ZANEW=PA(JL,JK)+ZDA(JL,JK)+PTENA(JL,JK)*PTSPHY
      LLO1=ZLNEW.LE.RLMIN.OR.ZANEW.LE.RAMIN.OR.
     S     (ZANEW.GT.0..AND.ZLNEW/MAX(ZANEW,RAMIN).GT.RCLDMAX)
      IF(LLO1) THEN
         ZLNEG(JL,JK)=-ZLNEW
         ZDA(JL,JK)=-PA(JL,JK)-PTENA(JL,JK)*PTSPHY
         ZDL(JL,JK)=-(PL(JL,JK)+PI(JL,JK))-(PTENL(JL,JK)+
     S      PTENI(JL,JK))*PTSPHY
         PLADJ(JL,JK)=ZLNEW
      END IF
      IF(ZANEW.GT.1.) THEN
         ZDA(JL,JK)=1.-PA(JL,JK)-PTENA(JL,JK)*PTSPHY
      END IF
      ZDLSUM=PLSCGE(JL,JK)+PLCOND1(JL,JK)+PLCOND2(JL,JK)
     S      -PLEVAP(JL,JK)-PLEROS(JL,JK)
      ZDQ=ZDLSUM-PLADJ(JL,JK)-ZREVAP(JL)-ZSEVAP(JL)
      PTENQ(JL,JK)=PTENQ(JL,JK)-ZDQ*ZQTMST
      IT=NINT(ZTP1(JL,JK)*1000.)
      ZDT=TLUCUC(IT)*(ZDLSUM-PLADJ(JL,JK))
CEVM  ZDT=FOELDCPM(ZTP1(JL,JK))*(ZDLSUM-PLADJ(JL,JK))
     1   -RALVDCP*ZREVAP(JL)-RALSDCP*ZSEVAP(JL)-RALFDCP*ZSMELT(JL)
      PTENT(JL,JK)=PTENT(JL,JK)+ZDT*ZQTMST
      PTENA(JL,JK)=PTENA(JL,JK)+ZDA(JL,JK)*ZQTMST
      ZALFAW=FOEALFA(ZTP1(JL,JK))
C     IF((PL(JL,JK)+PI(JL,JK)).GT.ZEPSEC) THEN
C        ZALFAOLD=PL(JL,JK)/(PL(JL,JK)+PI(JL,JK))
C     ELSE
C        ZALFAOLD=0.0
C     ENDIF
C     ZPHFLUX=(ZALFAW-ZALFAOLD)*(PL(JL,JK)+PI(JL,JK))
C     PTENL(JL,JK)=PTENL(JL,JK)+ZALFAW     *ZDL(JL,JK)*ZQTMST+
C    S             ZPHFLUX*ZQTMST
C     PTENI(JL,JK)=PTENI(JL,JK)+(1.-ZALFAW)*ZDL(JL,JK)*ZQTMST-
C    S             ZPHFLUX*ZQTMST
      ZLINEW=(PL(JL,JK)+PI(JL,JK))+ZDL(JL,JK)+(PTENL(JL,JK)+
     S PTENI(JL,JK))*PTSPHY
      ZLNEW=    ZALFAW *ZLINEW
      ZINEW=(1.-ZALFAW)*ZLINEW      
      PTENL(JL,JK)=(ZLNEW-PL(JL,JK))*ZQTMST
      PTENI(JL,JK)=(ZINEW-PI(JL,JK))*ZQTMST
  410 CONTINUE
C
      DO 420 JL=KIDIA,KFDIA
      IF(ZWENTR(JL).GT.0.) THEN
      ZSD=RCPD*ZTP1(JL,JK)+PGEO(JL,JK)
      ZSU=RCPD*ZTP1(JL,JK-1)+PGEO(JL,JK-1)
      LLO1=ZSU.GT.ZSD+1.E3.AND.ZQP1(JL,JK-1).LT.ZQP1(JL,JK)
      IF(LLO1) THEN
C
C        TRANSPORT OF CLOUD WATER
C
         ZGDPUP=RG/(PAPH(JL,JK)-PAPH(JL,JK-1))
         ZGDPDN=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
         ZDLTOP=ZLFLEN(JL)*ZGDPUP
         PTENQ(JL,JK-1)=PTENQ(JL,JK-1)+ZDLTOP
	 IT=NINT(ZTP1(JL,JK-1)*1000.)
         PTENT(JL,JK-1)=PTENT(JL,JK-1)-TLUCUC(IT)*ZDLTOP
CEVM     PTENT(JL,JK-1)=PTENT(JL,JK-1)-FOELDCPM(ZTP1(JL,JK-1))*ZDLTOP
C
C        TRANSPORT OF S, Q, U AND V
C
         ZALFA=PTSPHY*ZWENTR(JL)*ZGDPUP
         ZBETA=PTSPHY*ZWENTR(JL)*ZGDPDN
         ZALFA1=ZALFA/(1.+ZALFA)
         ZBETA1=ZBETA/(1.+ZBETA)
         ZXAB=1./(1.-ZALFA1*ZBETA1)
         Z1D1PA=ZXAB/(1.+ZALFA)
         Z1D1PB=ZXAB/(1.+ZBETA)
         ZSUNEW=ZSU*Z1D1PA+ZSD*ZALFA1*Z1D1PB
         ZSDNEW=ZSD*Z1D1PB+ZSU*ZBETA1*Z1D1PA
         ZQUNEW=ZQP1(JL,JK-1)*Z1D1PA+ZQS(JL,JK)*ZALFA1*Z1D1PB
         ZQDNEW=ZQS(JL,JK)*Z1D1PB+ZQP1(JL,JK-1)*ZBETA1*Z1D1PA
         ZUUNEW=PU(JL,JK-1)*Z1D1PA+PU(JL,JK)*ZALFA1*Z1D1PB
         ZUDNEW=PU(JL,JK)*Z1D1PB+PU(JL,JK-1)*ZBETA1*Z1D1PA
         ZVUNEW=PV(JL,JK-1)*Z1D1PA+PV(JL,JK)*ZALFA1*Z1D1PB
         ZVDNEW=PV(JL,JK)*Z1D1PB+PV(JL,JK-1)*ZBETA1*Z1D1PA
         ZACPDT=ZAMEAN(JL)*ZQCPDT
         ZADT=ZAMEAN(JL)*ZQTMST
         PTENT(JL,JK)=PTENT(JL,JK)+ZACPDT*(ZSDNEW-ZSD)
         PTENT(JL,JK-1)=PTENT(JL,JK-1)+ZACPDT*(ZSUNEW-ZSU)
         PTENQ(JL,JK)=PTENQ(JL,JK)+ZADT*(ZQDNEW-ZQS(JL,JK))
         PTENQ(JL,JK-1)=PTENQ(JL,JK-1)+ZADT*(ZQUNEW-ZQP1(JL,JK-1))
         PTENU(JL,JK)=PTENU(JL,JK)+ZADT*(ZUDNEW-PU(JL,JK))
         PTENU(JL,JK-1)=PTENU(JL,JK-1)+ZADT*(ZUUNEW-PU(JL,JK-1))
         PTENV(JL,JK)=PTENV(JL,JK)+ZADT*(ZVDNEW-PV(JL,JK))
         PTENV(JL,JK-1)=PTENV(JL,JK-1)+ZADT*(ZVUNEW-PV(JL,JK-1))
      END IF
      END IF
  420 CONTINUE
C****************************************************************
C
C               END OF VERTICAL LOOP
C
C****************************************************************
  590 CONTINUE
C
C----------------------------------------------------------------------
C
C        6.       FLUX COMPUTATIONS
C
C
  600 CONTINUE
CEVMC
CEVM      DO 601 JL=KIDIA,KFDIA
CEVM      PFSQLF(JL,1) =0.
CEVM      PFSQIF(JL,1) =0.
CEVM      PFCQLNG(JL,1)=0.
CEVM      PFCQNNG(JL,1)=0.
CEVM  601 CONTINUE
CEVMC
CEVM      DO 612 JK=1,KLEV
CEVM      DO 611 JL=KIDIA,KFDIA
CEVMC
CEVM      ZGDPH=-RG/(PAPH(JL,JK+1)-PAPH(JL,JK))*PTSPHY
CEVM      ZALFAW=FOEALFA(ZTP1(JL,JK))
CEVM      PFSQLF(JL,JK+1)=ZALFAW*(ZDL(JL,JK)-PLUDE(JL,JK))/ZGDPH+
CEVM     S  PFSQLF(JL,JK)
CEVM      PFSQIF(JL,JK+1)=(1.-ZALFAW)*(ZDL(JL,JK)-PLUDE(JL,JK))/ZGDPH+
CEVM     S  PFSQIF(JL,JK)
CEVM      PFCQLNG(JL,JK+1)=ZALFAW*ZLNEG(JL,JK)/ZGDPH+
CEVM     S  PFCQLNG(JL,JK)
CEVM      PFCQNNG(JL,JK+1)=(1.-ZALFAW)*ZLNEG(JL,JK)/ZGDPH+
CEVM     S  PFCQNNG(JL,JK)
CEVMC     
CEVM  611 CONTINUE
CEVM  612 CONTINUE
CEVMC
CEVMC---enthalpy flux due to precipitation
CEVMC
CEVM      DO 614 JK=1,KLEV+1
CEVM      DO 613 JL=KIDIA,KFDIA
CEVM      PFHPSL(JL,JK)=-RLVTT*PFPLSL(JL,JK)
CEVM      PFHPSN(JL,JK)=-RLSTT*PFPLSN(JL,JK)
CEVM  613 CONTINUE
CEVM  614 CONTINUE
C
C
C----------------------------------------------------------------------
C
      RETURN
      END
