      SUBROUTINE CUENTR
     *    (KHOR,     KLON,     KLEV,     KLEVP1,   KK,
     *     KSTART,   KSTOP,    KLEN,     KLAB,
     *     PTENH,    PQENH,    PQTE,     PAPH,     PAP,
     *     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     *     PPBASE,   PMFU,     PENTR,
     *     PDMFEN,   PDMFDE  
     *   )
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
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
C
      REAL     PTENH(KHOR,KLEV),       PQENH(KHOR,KLEV),
     *         PAP(KHOR,KLEV),         PAPH(KHOR,KLEVP1),
     *         PQTE(KHOR,KLEV),        PMFU(KHOR,KLEV),
     *         PENTR(KHOR),            PPBASE(KHOR)
      INTEGER  KLWMIN(KHOR),           KTYPE(KHOR),
     *         KCBOT(KHOR),            KCTOP0(KHOR),
     *         KLAB (KHOR,KLEV)
      LOGICAL  LDCUM(KHOR)
C
      REAL     PDMFEN(KHOR),           PDMFDE(KHOR)
C
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
C
      REAL
     *        ZPTOP  (JPHR),
     *        ZRHO   (JPHR)
      LOGICAL  LLO1,LLO2
C
C
C----------------------------------------------------------------------
C
C*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
C                  -------------------------------------------
C
  100 CONTINUE
      DO 105 JL=KSTART,KSTOP
      PDMFEN(JL)=0.
      PDMFDE(JL)=0.
      ZRHO(JL)=PAPH(JL,KK+1)/(RD*PTENH(JL,KK+1))
      PPBASE(JL)=PAPH(JL,KCBOT(JL))
      ZPTOP(JL)=PAPH(JL,KCTOP0(JL))
  105 CONTINUE
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
      DO 125 JL=KSTART,KSTOP
      IF(LDCUM(JL)) THEN
         ZDPRHO=(PAPH(JL,KK+1)-PAPH(JL,KK))/(G*ZRHO(JL))
         ZENTR=PENTR(JL)*PMFU(JL,KK+1)*ZDPRHO
         LLO1=KK.LT.KCBOT(JL)
C        IF(LLO1) PDMFDE(JL)=ZENTR
         IF(LLO1) THEN
	   IF (LMFNEW.AND.KTYPE(JL).EQ.2) THEN
	     PDMFDE(JL)=1.35*ZENTR
             IF (KLAB(JL,KK).EQ.0) THEN
	        PDMFDE(JL)=ZENTR
	      ENDIF
	   ELSE
	     PDMFDE(JL)=ZENTR
	   ENDIF
	 ENDIF
         ZPMID=0.5*(PPBASE(JL)+ZPTOP(JL))
         LLO2=LLO1.AND.KTYPE(JL).EQ.2.AND.
     1        (PPBASE(JL)-PAPH(JL,KK).LT.0.2E5.OR.
     1         PAPH(JL,KK).GT.ZPMID)
         IF(LLO2) PDMFEN(JL)=ZENTR
         IKLWMIN=MAX(KLWMIN(JL),KCTOP0(JL)+2)
         LLO2=LLO1.AND.(KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.3).AND.
     1        (KK.GE.IKLWMIN.OR.PAP(JL,KK).GT.ZPMID)
         IF(LLO2) PDMFEN(JL)=ZENTR
         IF(LLO2.AND.PQENH(JL,KK+1).GT.1.E-5)
     1   PDMFEN(JL)=ZENTR+MAX(PQTE(JL,KK),0.)/PQENH(JL,KK+1)*
     1              ZRHO(JL)*ZDPRHO
      END IF
  125    CONTINUE
C
      RETURN
      END
