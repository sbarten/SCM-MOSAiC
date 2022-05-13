      SUBROUTINE CUENTR4
     *     (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,   KK,
     *     KLAB,
     *     PTENH,    PQENH,    PQTE,     PAPHP1,   PAPP1,
     *     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     *     PPBASE,   PMFU,     PENTR,    PODETR,   POENTR,
     *     KHMIN,    PBUOY,    PGEOH,
     *     PDMFEN,   PDMFDE)
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE.
C          --------
C          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
C          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CUASC*.
C          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
C          AND UPDRAFT VALUES T,Q ETC
C          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
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
      REAL     PTENH(KLP2,KLEV),       PQENH(KLP2,KLEV),
     *         PAPP1(KLP2,KLEV),       PAPHP1(KLP2,KLEVP1),
     *         PQTE(KLP2,KLEV),        PMFU(KLP2,KLEV),
     *         PENTR(KLP2),            PPBASE(KLP2)
      REAL     PODETR(KLP2,KLEV)
      REAL     POENTR(KLP2,KLEV)
      REAL     PGEOH (KLP2,KLEV)
      REAL     PBUOY (KLP2)
      INTEGER  KHMIN (KLP2)
      INTEGER  KLWMIN(KLP2),           KTYPE(KLP2),
     *         KCBOT(KLP2),            KCTOP0(KLP2)
      INTEGER  KLAB (KLP2,KLEV)
      LOGICAL  LDCUM(KLP2)
C
      REAL     PDMFEN(KLP2),           PDMFDE(KLP2)
C
      LOGICAL  LLO1,LLO2
C
C
C----------------------------------------------------------------------
C
C*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
C                  -------------------------------------------
C
  100 CONTINUE
C
C
C*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
C                  --------------------------------------------
C
  110 CONTINUE
C
C
C*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
C                  -----------------------------------------
C
  120 CONTINUE
      ZRG=1./G
      DO 125 JL=1,KLON
      PPBASE(JL)=PAPHP1(JL,KCBOT(JL))
      ZRRHO=(RD*PTENH(JL,KK+1))/PAPHP1(JL,KK+1)
      ZDPRHO=(PAPHP1(JL,KK+1)-PAPHP1(JL,KK))*ZRG
      ZPMID=0.5*(PPBASE(JL)+PAPHP1(JL,KCTOP0(JL)))
      ZENTR=PENTR(JL)*PMFU(JL,KK+1)*ZDPRHO*ZRRHO
      LLO1=KK.LT.KCBOT(JL).AND.LDCUM(JL)
CEVM  PDMFDE(JL)=CVMGT(ZENTR,0.,LLO1)
CEVM     aps modification on shallow convection
      IF(LLO1) THEN
        IF (LMFNEW.AND.KTYPE(JL).EQ.2) THEN
          PDMFDE(JL)=1.35*ZENTR
          IF (KLAB(JL,KK).EQ.0) THEN
            PDMFDE(JL)=ZENTR
          ENDIF
        ELSE
          PDMFDE(JL)=ZENTR
        ENDIF
      ELSE
        PDMFDE(JL)=0.
      ENDIF
CEVM
      LLO2=LLO1.AND.KTYPE(JL).EQ.2.AND.
     1        (PPBASE(JL)-PAPHP1(JL,KK).LT.0.2E5.OR.
     1         PAPHP1(JL,KK).GT.ZPMID)
      PDMFEN(JL)=CVMGT(ZENTR,0.,LLO2)
      IKLWMIN=MAX(KLWMIN(JL),KCTOP0(JL)+2)
         LLO2=LLO1.AND.KTYPE(JL).EQ.3.AND.
     1        (KK.GE.IKLWMIN.OR.PAPP1(JL,KK).GT.ZPMID)
      IF(LLO2) PDMFEN(JL)=ZENTR
         IF (.NOT.LNORDSCV) THEN
           LLO2=LLO1.AND.KTYPE(JL).EQ.1
         ELSE
           LLO2=LLO1.AND.KTYPE(JL).NE.3
         ENDIF
         IF(LLO2) PDMFEN(JL)=ZENTR !turbulent entrainment
         ! organized detrainment, detrainment starts at khmin
         IKB=KCBOT(JL)
         PODETR(JL,KK)=0.
         IF(LLO2.AND.KK.LE.KHMIN(JL).AND.KK.GE.KCTOP0(JL)) THEN
           IKT=KCTOP0(JL)
           IKH=KHMIN(JL)
           IF(IKH.GT.IKT) THEN
           ZZMZK  =-(PGEOH(JL,IKH)-PGEOH(JL,KK))*ZRG
           ZTMZK  =-(PGEOH(JL,IKH)-PGEOH(JL,IKT))*ZRG
           ARG  =3.1415*(ZZMZK/ZTMZK)*0.5
           ZORGDE=TAN(ARG)*3.1415*0.5/ZTMZK
           ZDPRHO=(PAPHP1(JL,KK+1)-PAPHP1(JL,KK))*(ZRG*ZRRHO)
           IF (.NOT.LNORDSCV) THEN
             PODETR(JL,KK)=MIN(ZORGDE,1.E-3)*PMFU(JL,KK+1)*ZDPRHO
           ELSE
             PODETR(JL,KK)=MIN(ZORGDE,1.E-2)*PMFU(JL,KK+1)*ZDPRHO
           ENDIF
           ENDIF
         ENDIF
  125 CONTINUE
C
      RETURN
      END
