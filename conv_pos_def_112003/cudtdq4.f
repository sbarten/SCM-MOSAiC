      SUBROUTINE CUDTDQ4
     *    (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,
     *     NSTEP ,   NSTART,   TWODT,    CONACC,
     *     KTOPM2,   PAPHP1,   PGEO,     PTS,      LDLAND,
     *     LDCUM,    PTEN,     PTTE,     PQTE,

C LG- adding the tracers

     *     KTRAC,    PDTIME,
     *     PXTTE,    PMFUXT,   PMFDXT,

C LG-
*I HF240591.222
! pxten also needed for Brinkop tracer update, ! op_ck_20031001
     *     PXTEN,

C LG- end

     *     PXTEC,
     *     PMFUS,    PMFDS,    PMFUQ,    PMFDQ,
     *     PMFUL,    PDMFUP,   PDMFDP,   PLUDE,
     *     PDPMEL,   PRAIN,    PRFL,     PSFL,
     *     PSRAIN,   PSEVAP,   PSHEAT,   PSMELT,
     *     PRSFC,    PSSFC,    PAPRC,    PAPRS)
C
C
C**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
C                DOES GLOBAL DIAGNOSTICS
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C**   INTERFACE.
C     ----------
C
C          *CUDTDQ* IS CALLED FROM *CUMASTR*
C
CJHC*CALL PARAM
CJHC*CALL COMCTL
CJHC*CALL COMCON
CJHC*CALL COMCUMF
C
CJHC*CALL COMTRC
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'

C LG-

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

C LG- end


      LOGICAL  LLO1,LLO2,LLO3
C
      REAL     PTTE(KLP2,KLEV),        PQTE(KLP2,KLEV),
     *         PTEN(KLP2,KLEV),        PTS(KLP2),
     *         PGEO(KLP2,KLEV),        PAPHP1(KLP2,KLEVP1),
     *         PAPRC(KLP2),            PAPRS(KLP2),
     *         PRSFC(KLP2),            PSSFC(KLP2)
      REAL     PMFUS(KLP2,KLEV),       PMFDS(KLP2,KLEV),
     *         PMFUQ(KLP2,KLEV),       PMFDQ(KLP2,KLEV),
     *         PMFUL(KLP2,KLEV),       PLUDE(KLP2,KLEV),
     *         PDMFUP(KLP2,KLEV),      PDMFDP(KLP2,KLEV),
     *         PXTEC(KLP2,KLEV),
     *         PRFL(KLP2),             PRAIN(KLP2)
      REAL     PDPMEL(KLP2,KLEV),      PSFL(KLP2)
      LOGICAL  LDLAND(KLP2),           LDCUM(KLP2)
C
      REAL     zalpha(KLP2,klev,ktrac),zgamma(KLP2,ktrac),      ! op_ck_20031001
     *         ztend(KLP2,klev,ktrac), ztenm(KLP2,klev,ktrac)   ! op_ck_20031001

      INCLUDE 'paramh.h'
      REAL     ZMELT(JPHR)
      REAL     ZSHEAT(JPHR)
      REAL     ztmst   ! mz_lg_20031117+

C LG- adding the tracers

      REAL     PXTTE(KLON,NLEVT,KTRAC),   PMFUXT(KLON,KLEV,KTRAC),
     *         PMFDXT(KLON,KLEV,KTRAC)

C LG-
*I HF240591.224

      REAL PXTEN(KLON,KLEV,KTRAC),DHETCH(KTRAC)

      LOGICAL LO
C
C LG-

*I CUDTDQ.57
C -----------------------------------------------
C Correct Tiedtke's scheme for heterogeneously
C adapted mass fluxes
C -----------------------------------------------
      CALL RESETR(CVDPREC,KLON*KLEV,0.)
      CALL RESETR(DHETCH,KTRAC,0.)

C LG- end

C
C
C----------------------------------------------------------------------
C
C*    1.0          SPECIFY PARAMETERS
C                  ------------------
C
  100 CONTINUE

      zeps1=(1.-1.E-10)   ! op_ck_20031001

      ! mz_lg_20031118+
      IF (NSTEP.EQ.0) THEN
       ztmst=PDTIME
      ELSE
       ztmst=TWODT
      ENDIF
      ! mz_lg_20031118-

      ZDIAGT=CONACC*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
C
C
C----------------------------------------------------------------------
C
C*    2.0          INCREMENTATION OF T AND Q TENDENCIES
C                  ------------------------------------
C
  200 CONTINUE
      DO 210 JL=KIDIA,KFDIA
      ZMELT(JL)=0.
      ZSHEAT(JL)=0.
  210 CONTINUE

! op_ck_20031001+
      DO 2100 jt=1,ktrac
        DO 2101 jl=1,KIDIA,KFDIA
          ztend(jl,1,jt)=(pmfuxt(jl,1,jt)+pmfdxt(jl,1,jt))             
     *                *ztmst*g/(paphp1(jl,2)-paphp1(jl,1))
          ztenm(jl,1,jt)=max(ztend(jl,1,jt),pxten(jl,1,jt)*zeps1)
          ztend(jl,klev,jt)=(pmfuxt(jl,klev,jt)+pmfdxt(jl,klev,jt))    
     *                   *ztmst*g/(paphp1(jl,klev+1)-paphp1(jl,klev))
          ztenm(jl,klev,jt)=max(ztend(jl,klev,jt)                      
     *                       ,pxten(jl,klev,jt)*zeps1)
          zgamma(jl,jt)=999.
2101    ENDDO
2100  ENDDO
!
      DO 2110 jt=1,ktrac
        DO 2111 jk=2,klev-1
           DO 2112 jl=1,KIDIA,KFDIA
             ztend(jl,jk,jt)=(pmfuxt(jl,jk,jt)-pmfuxt(jl,jk+1,jt)     
     *                    +pmfdxt(jl,jk,jt)-pmfdxt(jl,jk+1,jt))     
     *                    *ztmst*g/(paphp1(jl,jk+1)-paphp1(jl,jk))
             ztenm(jl,jk,jt)=max(ztend(jl,jk,jt),pxten(jl,jk,jt)*zeps1)
2112       ENDDO
2111    ENDDO
2110  ENDDO
!
      DO 2120 jt=1,ktrac
        DO 2121 jk=1,klev
          DO 2122 jl=1,KIDIA,KFDIA
            IF (pxten(jl,jk,jt).le.0..and.ztenm(jl,jk,jt).eq.0.) THEN
               zalpha(jl,jk,jt)=1.
            ELSE
               zalpha(jl,jk,jt)=(pxten(jl,jk,jt)*zeps1)/(ztenm(jl,jk,jt))
            ENDIF
            zgamma(jl,jt)=min(zgamma(jl,jt),zalpha(jl,jk,jt))
            zgamma(jl,jt)=min(zgamma(jl,jt),1.)
            zgamma(jl,jt)=max(zgamma(jl,jt),0.)
2122      ENDDO 
2121    ENDDO  
2120  ENDDO
! op_ck_20031001-
C

C LG-

*I CUDTDQ.78
      ZXPERAT=0.
      ZXPERBT=0.
      ZXO3AT=0.
      ZXO3BT=0.

C LG- end

      DO 250 JK=KTOPM2,KLEV

C
      IF(JK.LT.KLEV) THEN
         DO 220 JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            LLO1=(PTEN(JL,JK)-TMELT).GT.0.
            ZALV=CVMGT(ALV,ALS,LLO1)
            ZDTDT=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*RCPD*
     1               (PMFUS(JL,JK+1)-PMFUS(JL,JK)+
     1                PMFDS(JL,JK+1)-PMFDS(JL,JK)
     1      -ALF*PDPMEL(JL,JK)
     1         -ZALV*(PMFUL(JL,JK+1)-PMFUL(JL,JK)-
     1                  PLUDE(JL,JK)-
     1               (PDMFUP(JL,JK)+PDMFDP(JL,JK))))
            PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
            ZDQDT=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     1               (PMFUQ(JL,JK+1)-PMFUQ(JL,JK)+
     1                PMFDQ(JL,JK+1)-PMFDQ(JL,JK)+
     1                PMFUL(JL,JK+1)-PMFUL(JL,JK)-
     1                  PLUDE(JL,JK)-
     1           (PDMFUP(JL,JK)+PDMFDP(JL,JK)))
            PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
       PXTEC(JL,JK)=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*PLUDE(JL,JK)
            ZSHEAT(JL)=ZSHEAT(JL)+ZALV*(PDMFUP(JL,JK)+PDMFDP(JL,JK))
            ZMELT(JL)=ZMELT(JL)+PDPMEL(JL,JK)
         END IF
  220 CONTINUE
C

C LG- the calculations yield strange results over the poles and 
C     therefore the calculations are only performed for IROW >10

      IF (LXTCONV.AND.IROW.GE.10) THEN

      DO 2204 JT=1,KTRAC

C LG-

*I HF240591.229
      DXTHET=0.

      DO 2202 JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN

C LG-
*I HF240591.234
       ZDXDEP=0.

! ======================================================================
! mz_lg_20031118+ modified by removing IF (JT.EQ... etc.) statement.
!      Calculation is being done for all the tracers in order to secure
!      positive definite convective tracer transport scheme!
! ======================================================================
       ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     *       DCONCCV(JL,JK,JT)/1000.

! op_ck_20031001+
C        ZDXTDT=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))
C      *       *(PMFUXT(JL,JK+1,JT)-PMFUXT(JL,JK,JT)
C      *        +PMFDXT(JL,JK+1,JT)-PMFDXT(JL,JK,JT))
       zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        
     *                 *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) 
     *                 +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt) 
     *                 +zdxdep)  ! Geert-Jan's code
     *                 *zgamma(jl,jt)
! op_ck_20031001-

       DXTHET=DXTHET+ZDXDEP*GRMASS(JL,JK)*PDTIME

C LG-  update of the tracer tendency considering both the upward and
C      and downward mass flux and the wet deposition processes (ZDXDEP)

       PXTTE(JL,JK,JT)=PXTTE(JL,JK,JT)+ZDXTDT

      ENDIF

 2202 CONTINUE

C LG-
*I HF240591.237
      DHETCH(JT)=DHETCH(JT)+DXTHET

 2204 CONTINUE

C LG-
*I HF240591.238
      DO 2205 JL=1,NLON
       IF (LDCUM(JL)) THEN
         GPAPH=G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))/1000.
         ZXPERA=-1.*GPAPH*XPERA(JL,JK)
         ZXPERB=-1.*GPAPH*XPERB(JL,JK)
         ZXO3A=-1.*GPAPH*XO3A(JL,JK)
         ZXO3B=-1.*GPAPH*XO3B(JL,JK)
         ZXPERAT=ZXPERAT+ZXPERA*GRMASS(JL,JK)*ZDIAGT
         ZXPERBT=ZXPERBT+ZXPERB*GRMASS(JL,JK)*ZDIAGT
         ZXO3AT=ZXO3AT+ZXO3A*GRMASS(JL,JK)*ZDIAGT
         ZXO3BT=ZXO3BT+ZXO3B*GRMASS(JL,JK)*ZDIAGT
       ENDIF
 2205 CONTINUE

      ENDIF

C LG- end

C
      ELSE
         DO 230 JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            LLO1=(PTEN(JL,JK)-TMELT).GT.0.
            ZALV=CVMGT(ALV,ALS,LLO1)
            ZDTDT=-(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*RCPD*
     1             (PMFUS(JL,JK)+PMFDS(JL,JK)+ALF*PDPMEL(JL,JK)-ZALV*
     1          (PMFUL(JL,JK)+PDMFUP(JL,JK)
     1           +PDMFDP(JL,JK)+PLUDE(JL,JK)))
            PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
            ZDQDT=-(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*
     1                (PMFUQ(JL,JK)+PMFDQ(JL,JK)+
     1           PLUDE(JL,JK)+
     1                (PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))
            PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
         PXTEC(JL,JK)=(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*PLUDE(JL,JK)
            ZSHEAT(JL)=ZSHEAT(JL)+ZALV*(PDMFUP(JL,JK)+PDMFDP(JL,JK))
            ZMELT(JL)=ZMELT(JL)+PDPMEL(JL,JK)
         END IF
  230    CONTINUE
C

C LG- adding the tracers

      IF (LXTCONV.AND.IROW.GT.10) THEN

      DO 2304 JT=1,KTRAC

C LG-

*I HF240591.243
      DXTHET=0.

      DO 2302 JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN

C LG-

*I HF240591.247
       ZDXDEP=0.

! ======================================================================
! mz_lg_20031118+ modified by removing IF (JT.EQ... etc.) statement.
!      Calculation is being done for all the tracers in order to secure
!      positive definite convective tracer transport scheme!
! ======================================================================

       ZDXDEP=-(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))*DCONCCV(JL,JK,JT)
     * /1000.

! op_ck_20031001+
C        ZDXTDT=-(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))
C      *        *(PMFUXT(JL,JK,JT)+PMFDXT(JL,JK,JT))

       zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))       
     *                  *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt)       
     *                  +zdxdep)  ! Geert-Jan's code
     *                  *zgamma(jl,jt)

       DXTHET=DXTHET+ZDXDEP*GRMASS(JL,JK)*PDTIME

C LG-  update of the tracer tendency considering both the upward and
C      and downward mass flux and the wet deposition processes (ZDXDEP)

       PXTTE(JL,JK,JT)=PXTTE(JL,JK,JT)+ZDXTDT

      ENDIF
 2302 CONTINUE

C LG-

*I HF240591.250
      DHETCH(JT)=DHETCH(JT)+DXTHET

 2304 CONTINUE

C LG-

*I HF240591.251
      DO 2305 JL=1,NLON
       IF (LDCUM(JL)) THEN
         GPAPH=G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))/1000.
         ZXPERA=-1.*GPAPH*XPERA(JL,JK)
         ZXPERB=-1.*GPAPH*XPERB(JL,JK)
         ZXO3A=-1.*GPAPH*XO3A(JL,JK)
         ZXO3B=-1.*GPAPH*XO3B(JL,JK)
         ZXPERAT=ZXPERAT+ZXPERA*GRMASS(JL,JK)*ZDIAGT
         ZXPERBT=ZXPERBT+ZXPERB*GRMASS(JL,JK)*ZDIAGT
         ZXO3AT=ZXO3AT+ZXO3A*GRMASS(JL,JK)*ZDIAGT
         ZXO3BT=ZXO3BT+ZXO3B*GRMASS(JL,JK)*ZDIAGT
       ENDIF
 2305 CONTINUE

      ENDIF

C
      END IF
C

C LG-

*I CUDTDQ.123

      DO 2405 JL=1,KLON
        CVDPREC(JL,JK)=PDMFUP(JL,JK)+PDMFDP(JL,JK)
 2405 CONTINUE

  250 CONTINUE

C LG-

*I CUDTDQ.126    

      IW=IWHERE(1,KLEV)
      PEROXT=ZXPERAT+ZXPERBT
      O3OXT=ZXO3AT+ZXO3BT
      TOTOX=PEROXT+O3OXT
      SCAVOX=ZXPERAT+ZXO3AT
      RLISOX=TOTOX-SCAVOX
      BXT(IW,IBWDC,ihno3)=BXT(IW,IBWDC,ihno3)-DHETCH(ihno3)
      BXT(IW,IBWDC,ih2o2)=BXT(IW,IBWDC,ih2o2)-DHETCH(ih2o2)+PEROXT
      BXT(IW,IBWDC,iso2)=BXT(IW,IBWDC,iso2)-DHETCH(iso2)+TOTOX
      BXT(IW,IBWDC,iso4)=BXT(IW,IBWDC,iso4)-DHETCH(iso4)-TOTOX
      BXT(IW,IBWDC,ihplus)=BXT(IW,IBWDC,ihplus)-
     *   DHETCH(ihno3)-DHETCH(iso2)-DHETCH(iso4)

      BXT(IW,IBCHEM,iso2)=BXT(IW,IBCHEM,iso2)+TOTOX
      BXT(IW,IBCHEM,iso4)=BXT(IW,IBCHEM,iso4)-TOTOX
      BXT(IW,IBCHEM,ih2o2)=BXT(IW,IBCHEM,ih2o2)+PEROXT
      BRX(IW,45)=BRX(IW,45)-PEROXT
      BRX(IW,46)=BRX(IW,46)-O3OXT

C LG- added some trace gases of the hydrocarbon chemistry scheme

      DO JT=1,KTRAC
       IF (HENRY(JT).GT.0.) THEN
         BXT(IW,IBWDC,JT)=BXT(IW,IBWDC,JT)-DHETCH(JT)
       ENDIF
      ENDDO

C LG- end



*/ -------

C LG- end


C
C
C---------------------------------------------------------------------
C
C      3.          UPDATE SURFACE FIELDS AND DO GLOBAL BUDGETS
C                  -------------------------------------------
C
  300 CONTINUE
      DO 310 JL=KIDIA,KFDIA
      PRSFC(JL)=PRFL(JL)
      PSSFC(JL)=PSFL(JL)
      PAPRC(JL)=PAPRC(JL)+ZDIAGW*(PRFL(JL)+PSFL(JL))
      PAPRS(JL)=PAPRS(JL)+ZDIAGW*PSFL(JL)
      PSHEAT=PSHEAT+ZSHEAT(JL)
      PSRAIN=PSRAIN+PRAIN(JL)
      PSEVAP=PSEVAP-(PRFL(JL)+PSFL(JL))
      PSMELT=PSMELT+ZMELT(JL)
  310 CONTINUE
C
      PSEVAP=PSEVAP+PSRAIN
C
C
      RETURN
      END
