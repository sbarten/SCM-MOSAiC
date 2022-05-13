      SUBROUTINE CUFLX4
     *     (KIDIA, KFDIA,
     *     KLP2,     KLON,     KLEV,     KLEVP1,
     *     NSTEP,   NSTART,  TWODT,
     *     PQEN,     PQSEN,    PTENH,    PQENH,

C     LG- adding the tracers

     *     KTRAC,   PDTIME,
     *     PXTEN,    PXTENH,   PMFUXT,   PMFDXT,

C     LG-
C     *I US310195.48
     *     PXTU, PLU,ZXTD,

C     LG- end

     *     PAPHP1,     LDLAND,   PGEOH,
     *     KCBOT,    KCTOP,    KDTOP,
     *     KTYPE,    LDDRAF,   LDCUM,
     *     PMFU,     PMFD,     PMFUS,    PMFDS,
     *     PMFUQ,    PMFDQ,    PMFUL,    PLUDE,
     *     PDMFUP,   PDMFDP,   PRFL,     PRAIN,   PRFLCK,
     *     PTEN,     PSFL,     PDPMEL,   KTOPM2)
C     
C     M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C     
C     PURPOSE
C     -------
C     
C     THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
C     FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
C     
C     INTERFACE
C     ---------
C     THIS ROUTINE IS CALLED FROM *CUMASTR*.
C     
C     EXTERNALS
C     ---------
C     NONE
C     
c     jhc*CALL COMCON
c     jhc*CALL COMCUMF
c     jhc*CALL PARAM
c     jhc*CALL COMSDS
c     jhc*CALL COMCTL
c     jhc*CALL COMPH2
      INCLUDE 'comcon.h'
      INCLUDE 'comcumf.h'
      INCLUDE 'comph2.h'

C     LG- 

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

C     LG- end

C     
      REAL     PQEN(KLP2,KLEV),        PQSEN(KLP2,KLEV),
     *     PTENH(KLP2,KLEV),       PQENH(KLP2,KLEV),
     *     PAPHP1(KLP2,KLEVP1),      PGEOH(KLP2,KLEV)
C     
      REAL     PMFU(KLP2,KLEV),        PMFD(KLP2,KLEV),
     *     PMFUS(KLP2,KLEV),       PMFDS(KLP2,KLEV),
     *     PMFUQ(KLP2,KLEV),       PMFDQ(KLP2,KLEV),
     *     PDMFUP(KLP2,KLEV),      PDMFDP(KLP2,KLEV),
     *     PMFUL(KLP2,KLEV),       PLUDE(KLP2,KLEV),
     *     PRFL(KLP2),             PRAIN(KLP2),
     *     PRFLCK(KLP2,KLEV)
      REAL     PTEN(KLP2,KLEV),          PDPMEL(KLP2,KLEV),
     *     PSFL(KLP2)
      INTEGER  KCBOT(KLP2),            KCTOP(KLP2),
     *     KDTOP(KLP2),            KTYPE(KLP2)
      LOGICAL  LDDRAF(KLP2),           LDLAND(KLP2),
     *     LDCUM(KLP2)

C     LG- adding the tracers

      REAL   PXTEN(KLON,KLEV,KTRAC),   PXTENH(KLON,KLEV,KTRAC),
     *     PMFUXT(KLON,KLEV,KTRAC),  PMFDXT(KLON,KLEV,KTRAC)

C     LG-

C     LG- extra declarations for hydrocarbon species and some other
C     parameters

      REAL HS0(NTRACT),DHT(NTRACT),FRACL(KTRAC),
     *     XDEPCV(KLON,KLEV,KTRAC),PMX(NLON,KTRAC),
     *     PML(NLON,KTRAC),PMML(KLON,KTRAC)

      REAL DPM,DPML,TX

C     LG- end


*     I HF240591.187
      REAL PXTU(KLON,KLEV,KTRAC),PLU(KLP2,KLEV),RAINH(KLON),
     *     PXTD(KLON,KLEV,KTRAC)
      REAL PT2(KLON),PRHOA2(KLON),RDRAD(KLON),SO2X(KLON),O3X(KLON),
     *     H2O2X(KLON),HNO3X(KLON),SAERL(KLON),SIVL(KLON),H2O2L(KLON),
     *     ANO3ML(KLON),HPL(KLON),SOXPER(KLON),SOXO3(KLON)
      REAL RLWC(KLON),WATFR2(KLON),TOTWFLX(KLON)

C     LG- end

C     
      INCLUDE 'paramh.h'
      REAL     ZPSUBCL(JPHR)
c     jhc      POINTER (IZPSUBCL,ZPSUBCL)
C     
C     *             SPECIFY CONSTANTS
C     
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZCONS1=CPD/(ALF*G*ZTMST)
      ZCONS2=1./(G*ZTMST)
      ZCUCOV=0.05
      ZTMELP2=TMELT+2.
C     
C     LG-

C     LG- definition of constants used to calculate the Henry coefficients
C     for a selection of trace gasesof the hydrocarbon chemistry 
C     scheme. The required parameters are the Henry's law coefficient at
C     298 K in mole l-1 air mole l-1 water (see Geert-Jan Roelofs thesis,
C     page 43), HS0 and the parameter dH/R [K] with dH being the enthalpy 
C     change for the dissolution reaction and R is the universal gas constant.
C     For the applied relationship see page 43, equation II.4.7 of the thesis.

      CALL RESETR (HS0,NTRACT,0.)
      CALL RESETR (DHT,NTRACT,0.)

C     LG- the recalculation is: Henry coefficients are given in M atm-1, a gas
C     pressure of 1 atm at 298 K corresponds with a concentration of 
C     0.0409 (=R*T with R=0.082 atm M-1 K-1 and K is 298 K, see Seinfeld, page
C     198, the final definition of the Henry coeff. is in Mole l-1 air 
C     Mole l-1 water, thus taking the inverse value

      HS0(ich3o2h)=6.82E-3
      HS0(ich2o)=1.36E-5
      HS0(ipan)=8.18E-3
      IF (LCBM4_ECH) THEN
         HS0(ihcooh)=7.40E-6
         HS0(ich3co2h)=7.44E-6
      ENDIF

      DHT(ich3o2h)=-5316.
      DHT(ich2o)=-7217.
      DHT(ipan)=0.
      IF (LCBM4_ECH) THEN
         DHT(ihcooh)=-5629.
         DHT(ich3co2h)=-5894.
      ENDIF

*     I CUFLX.46

      CALL RESETR(DCONCCV,KLON*KLEV*KTRAC,0.)
      CALL RESETI(KCONBOT,KLON,0)
      CALL RESETR(CVWDFLX,KLON*(KTRAC+1),0.)
      CALL RESETR(RAINH,KLON,0.)
      CALL RESETR(XO3A,KLON*KLEV,0.)
      CALL RESETR(XO3B,KLON*KLEV,0.)
      CALL RESETR(XPERA,KLON*KLEV,0.)
      CALL RESETR(XPERB,KLON*KLEV,0.)

C     LG- resetting the Henry coeff.

      CALL RESETR(HENRY,KTRAC,0.)

C     -----------------------------------------------
C     Calculate heterogeneous chemistry in convective clouds
C     -----------------------------------------------

C     LG- end

C     
c     jhc      ILENX=KLP2
c     jhc      CALL ALLOCA(IZPSUBCL,ILENX,'CUFLX',99)
C     
C     *    1.0          DETERMINE FINAL CONVECTIVE FLUXES
C     ---------------------------------
C     
 100  CONTINUE
      ITOP=KLEV
      DO 110 JL=KIDIA,KFDIA
         ITOP=MIN(ITOP,KCTOP(JL))
         IF(.NOT.LDCUM(JL).OR.KDTOP(JL).LT.KCTOP(JL)) LDDRAF(JL)=.FALSE.
         IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
 110  CONTINUE
      KTOPM2=ITOP-2
C====================================================================
C     EVM     A.P. SIEBESMA modification on shallow convection
C     
C------------------------------------------------------------------
C     |
C     DETERMINATION OF KTOPM2 IN CASE OF HIGH RESOLUTION MODE	  |
C     |
C------------------------------------------------------------------
C     
      JL   = 1
      IF (LMFHRES.AND.(KTYPE(JL).LT.3)) THEN
C     
C     
C-----------------------------------------------
C     KTOPM2: LOWEST LEVEL WITH ZERO MASS FLUX  |
C     |
C     ORIGINAL TIEDTKE: KTOPM2 = ITOP-2		|
C------------------------------------------------
C     
C     
         ICTOP = KLEV
         DO 102 JK = 2,KCTOP(JL)
            IF ((PMFU(JL,JK).EQ.0.0).AND.(PMFU(JL,JK+1).GT.0.0)) THEN
               ICTOP = JK+1
            ENDIF
 102     CONTINUE
         KTOPM2 = ICTOP-1
C     
C     
      ENDIF
C     
C---------------END OF LMFHIRES  -------------
C     
C     
C====================================================================
      DO 120 JK=KTOPM2,KLEV
C     DIR$ IVDEP
         DO 115 JL=KIDIA,KFDIA
C     EVM  IF(LDCUM(JL).AND.JK.GE.KCTOP(JL)-1) THEN
            IF(LDCUM(JL).AND.JK.GE.(KTOPM2+1)) THEN
               PMFUS(JL,JK)=PMFUS(JL,JK)-PMFU(JL,JK)*
     1              (CPD*PTENH(JL,JK)+PGEOH(JL,JK))
               PMFUQ(JL,JK)=PMFUQ(JL,JK)-PMFU(JL,JK)*PQENH(JL,JK)
               IF(LDDRAF(JL).AND.JK.GE.KDTOP(JL)) THEN
                  PMFDS(JL,JK)=PMFDS(JL,JK)-PMFD(JL,JK)*
     1                 (CPD*PTENH(JL,JK)+PGEOH(JL,JK))
                  PMFDQ(JL,JK)=PMFDQ(JL,JK)-PMFD(JL,JK)*PQENH(JL,JK)
               ELSE
                  PMFD(JL,JK)=0.
                  PMFDS(JL,JK)=0.
                  PMFDQ(JL,JK)=0.
                  PDMFDP(JL,JK-1)=0.
               END IF
            ELSE
               PMFU(JL,JK)=0.
               PMFD(JL,JK)=0.
               PMFUS(JL,JK)=0.
               PMFDS(JL,JK)=0.
               PMFUQ(JL,JK)=0.
               PMFDQ(JL,JK)=0.
               PMFUL(JL,JK)=0.
               PDMFUP(JL,JK-1)=0.
               PDMFDP(JL,JK-1)=0.
               PLUDE(JL,JK-1)=0.
            END IF
 115     CONTINUE
C     

*     I CUFLX.93

         IF (LCHEM) THEN

            DO 1111 JL=1,NLON
               PT2(JL)=PTENH(JL,JK)
               PRHOA2(JL)=GRMASS(JL,JK)/GRVOL(JL,JK)
               IF (LDCUM(JL)) KCONBOT(JL)=KCBOT(JL)
 1111       CONTINUE
            DO 1112 JL=1,NLON
               IF (LDCUM(JL).AND.JK.GE.KCTOP(JL)-1.AND.JK.LE.KCBOT(JL)) THEN
c     determine fraction ice
                  TCELC=PTENH(JL,JK)-273.15
                  FICE=0.
                  IF (TCELC.LT.-30.) FICE=1.
                  IF (TCELC.GE.-30.0.AND.TCELC.LE.0.)
     *                 FICE=-1.*TCELC/30.
                  WATFR2(JL)=1.-FICE
                  TOTWFLX(JL)=PLU(JL,JK)*PMFU(JL,JK)+PLUDE(JL,JK)+PDMFUP(JL,JK)
                  RLWC(JL)=TOTWFLX(JL)*WATFR2(JL)/PMFU(JL,JK)*PRHOA2(JL)*1E3
               ELSE
                  RLWC(JL)=0.
               ENDIF
 1112       CONTINUE

C     LG- calculating the concentration changes due to convective transport
C     with cumulus clouds if LXTCONV=.TRUE., this statement is not introduced
C     in the beginning of the tracer calculations since for the subroutine
C     CLSCAV.f some of the previously calculated parameters are required

            IF (LXTCONV) THEN

               IAQCALC=0
               CALL AQCALC (RLWC,IAQCALC)
               IF (IAQCALC.EQ.0) GOTO 5555
               DO 1113 JL=1,NLON
                  IF (RLWC(JL).GT.1E-5) THEN
                     SCAVEFF=0.60
                     SAERL(JL)=WATFR2(JL)*SCAVEFF*PXTU(JL,JK,iso4)*PRHOA2(JL)
                     HPL(JL)=SAERL(JL)
                     SO2X(JL)=PXTU(JL,JK,iso2)*PRHOA2(JL)
                     O3X(JL)=PXTU(JL,JK,io3)*PRHOA2(JL)
                     H2O2X(JL)=PXTU(JL,JK,ih2o2)*PRHOA2(JL)
                     HNO3X(JL)=PXTU(JL,JK,ihno3)*PRHOA2(JL)
                     SIVL(JL)=0.
                     H2O2L(JL)=0.
                     SOXPER(JL)=0.
                     SOXO3(JL)=0.
                     ANO3ML(JL)=0.
                     RDRAD(JL)=1E-3
                     RAINH(JL)=RAINH(JL)+PDMFUP(JL,JK)

C     LG-   added some trace gases of the hydrocarbon chemistry scheme

                     DO JT=1,KTRAC

C     LG-    Calculation of the Henry coeffients for a selection of
C     trace gases of the hydrocarbon chemistry scheme, consistent 
C     with the definition used by Geert-Jan Roelofs (mol mol-1, see 
C     his thesis)

                        TX=1./PT2(JL)-1./298.

c     calculate rx-constants

                        HENRY(JT)=HS0(JT)*EXP(DHT(JT)*TX)

C     LG- end 

                        IF (HENRY(JT).GT.0.) THEN
                           PMX(JL,JT)=PXTU(JL,JK,JT)*PRHOA2(JL)
                           PMML(JL,JT)=0.
                        ENDIF

                     ENDDO

C     LG-   end


                  ENDIF
 1113          CONTINUE
c     
 5023          format (8e12.4)

               CALL SULFAQ (RLWC,RDRAD,SO2X,O3X,
     *              H2O2X,HNO3X,PMX,PT2,PRHOA2,ZTMST,
     *              SAERL,SIVL,H2O2L,ANO3ML,HPL,PML,PMML,SOXPER,
     *              SOXO3,0)
c     
               DO 1114 JL=1,NLON
                  IF (RLWC(JL).GT.1E-5) THEN
                     SCAVCO=WATFR2(JL)*PDMFUP(JL,JK)/TOTWFLX(JL)
                     ZMATS=0.3
                     FAC0=SCAVCO
                     FAC1=FAC0*ZMATS
                     ZMTOF=1000.*PMFU(JL,JK)
                     DHNO3L=FAC1*ANO3ML(JL)/PRHOA2(JL)
                     DHNO3L=MAX(0.,DHNO3L)
                     DHNO3L=MIN(PXTU(JL,JK,ihno3),DHNO3L)
                     DCONCCV(JL,JK,ihno3)=DHNO3L*ZMTOF
                     CVWDFLX(JL,ihno3)=CVWDFLX(JL,ihno3)+DCONCCV(JL,JK,ihno3)
                     PXTU(JL,JK,ihno3)=PXTU(JL,JK,ihno3)-DHNO3L
                     PMFUXT(JL,JK,ihno3)=PXTU(JL,JK,ihno3)*PMFU(JL,JK)
                     XPERA1=SOXPER(JL)/PRHOA2(JL)*FAC1
                     XPERB1=SOXPER(JL)/PRHOA2(JL)*(ZMATS-FAC1)
                     XPERA(JL,JK)=XPERA1*ZMTOF
                     XPERB(JL,JK)=XPERB1*ZMTOF
                     XO3A1=SOXO3(JL)/PRHOA2(JL)*FAC1
                     XO3B1=SOXO3(JL)/PRHOA2(JL)*(ZMATS-FAC1)
                     XO3A(JL,JK)=XO3A1*ZMTOF
                     XO3B(JL,JK)=XO3B1*ZMTOF
                     DCONC=XO3A1+XO3B1
                     DCONCCV(JL,JK,io3)=DCONC*ZMTOF
                     PXTU(JL,JK,io3)=PXTU(JL,JK,io3)-DCONC
                     PMFUXT(JL,JK,io3)=PXTU(JL,JK,io3)*PMFU(JL,JK)
                     DH2O2L=H2O2L(JL)*FAC1/PRHOA2(JL)
                     DH2O2L=MAX(0.,DH2O2L)
                     DH2O2L=MIN(PXTU(JL,JK,ih2o2),DH2O2L)
                     DCONC=DH2O2L+XPERA1+XPERB1
                     DCONCCV(JL,JK,ih2o2)=DCONC*ZMTOF
                     PXTU(JL,JK,ih2o2)=PXTU(JL,JK,ih2o2)-DCONC
                     PMFUXT(JL,JK,ih2o2)=PXTU(JL,JK,ih2o2)*PMFU(JL,JK)
                     CVWDFLX(JL,ih2o2)=CVWDFLX(JL,ih2o2)+DH2O2L*ZMTOF
                     DSIVL=SIVL(JL)*FAC1/PRHOA2(JL)
                     DSIVL=MAX(0.,DSIVL)
                     DSIVL=MIN(PXTU(JL,JK,iso2),DSIVL)
                     DCONC=DSIVL+XPERA1+XPERB1+XO3A1+XO3B1
                     DCONCCV(JL,JK,iso2)=DCONC*ZMTOF
                     PXTU(JL,JK,iso2)=PXTU(JL,JK,iso2)-DCONC
                     PMFUXT(JL,JK,iso2)=PXTU(JL,JK,iso2)*PMFU(JL,JK)
                     CVWDFLX(JL,iso2)=CVWDFLX(JL,iso2)+DSIVL*ZMTOF
                     DSAERL=FAC1*SAERL(JL)/PRHOA2(JL)
                     DSAERL=MAX(0.,DSAERL)
                     DSAERL=MIN(PXTU(JL,JK,iso4),DSAERL)
                     DCONC=DSAERL-XPERB1-XO3B1
                     DCONCCV(JL,JK,iso4)=DCONC*ZMTOF
                     PXTU(JL,JK,iso4)=PXTU(JL,JK,iso4)-DCONC
                     PMFUXT(JL,JK,iso4)=PXTU(JL,JK,iso4)*PMFU(JL,JK)
                     CVWDFLX(JL,iso4)=CVWDFLX(JL,iso4)+(DSAERL+XPERA1+XO3A1)*ZMTOF
                     DHPL=HPL(JL)*FAC1/PRHOA2(JL)
                     DHPL=MAX(0.,DHPL)
                     DCONCCV(JL,JK,ihplus)=DCONC*ZMTOF
                     CVWDFLX(JL,ihplus)=CVWDFLX(JL,ihplus)+DCONCCV(JL,JK,ihplus)
                     DHPL=FAC1*HPL(JL)/PRHOA2(JL)
                     DHPL=MAX(0.,DHPL)
                     DCONCCV(JL,JK,ihplus)=DHPL*ZMTOF
                     CVWDFLX(JL,ihplus)=CVWDFLX(JL,ihplus)+DCONCCV(JL,JK,ihplus)

C     LG-  added some trace gases of the hydrocarbon chemistry scheme

                     DO JT=1,KTRAC

                        IF (HENRY(JT).GT.0.) THEN

C     LG-    PMML is > 0 for dissociating trace gases such as HCOOH and 
C     then PMML is considered, otherwise the concentration of the 
C     dissolved gas is used
                           
                           IF (PMML(JL,JT).GT.0.) THEN         
                              DPML=FAC1*PMML(JL,JT)/PRHOA2(JL)
                           ELSE
                              DPML=FAC1*PML(JL,JT)/PRHOA2(JL)
                           ENDIF

                           DPML=MAX(0.,DPML)
                           DPML=MIN(PXTU(JL,JK,JT),DPML)
                           DCONCCV(JL,JK,JT)=DPML*ZMTOF
                           CVWDFLX(JL,JT)=CVWDFLX(JL,JT)+DCONCCV(JL,JK,JT)
                           PXTU(JL,JK,JT)=PXTU(JL,JK,JT)-DPML
                           PMFUXT(JL,JK,JT)=PXTU(JL,JK,JT)*PMFU(JL,JK)

                        ENDIF

                     ENDDO

C     LG-  end

                  ENDIF
 1114          CONTINUE
 5555          CONTINUE

C     LG- end IF (LXTCONV)

            ENDIF

C     LG- adding the tracers

            DO 1154 JT=1,KTRAC
               DO 1152 JL=KIDIA,KFDIA
                  IF(LDCUM(JL).AND.JK.GE.KCTOP(JL)-1) THEN
                     PMFUXT(JL,JK,JT)=PMFUXT(JL,JK,JT)
     *                    -PMFU(JL,JK)*PXTENH(JL,JK,JT)
                     IF(LDDRAF(JL).AND.JK.GE.KDTOP(JL)) THEN
                        PMFDXT(JL,JK,JT)=PMFDXT(JL,JK,JT)
     *                       -PMFD(JL,JK)*PXTENH(JL,JK,JT)
                     ELSE
                        PMFDXT(JL,JK,JT)=0.
                     ENDIF
                  ELSE
                     PMFUXT(JL,JK,JT)=0.
                     PMFDXT(JL,JK,JT)=0.
                  ENDIF
 1152          CONTINUE
 1154       CONTINUE

C     LG- endif (LCHEM)

         ENDIF

C     LG- end

C     
 120  CONTINUE
      DO 130 JK=KTOPM2,KLEV
C     DIR$ IVDEP
         DO 125 JL=KIDIA,KFDIA
            IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
               IKB=KCBOT(JL)
               ZZP=((PAPHP1(JL,KLEVP1)-PAPHP1(JL,JK))/
     1              (PAPHP1(JL,KLEVP1)-PAPHP1(JL,IKB)))
               ZZP=CVMGT(ZZP**2,ZZP,KTYPE(JL).EQ.3)
               PMFU(JL,JK)=PMFU(JL,IKB)*ZZP
               PMFUS(JL,JK)=PMFUS(JL,IKB)*ZZP
               PMFUQ(JL,JK)=PMFUQ(JL,IKB)*ZZP
               PMFUL(JL,JK)=PMFUL(JL,IKB)*ZZP
            END IF
 125     CONTINUE
C     

C     LG- adding the tracers

         DO 1254 JT=1,KTRAC
C     DIR$ IVDEP
            DO 1252 JL=KIDIA,KFDIA
               IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
                  IKB=KCBOT(JL)
                  ZZP=(PAPHP1(JL,KLEVP1)-PAPHP1(JL,JK))/
     *                 (PAPHP1(JL,KLEVP1)-PAPHP1(JL,IKB))
                  ZZP=CVMGT(ZZP**2,ZZP,KTYPE(JL).EQ.3)
                  PMFUXT(JL,JK,JT)=PMFUXT(JL,IKB,JT)*ZZP
               ENDIF
 1252       CONTINUE
 1254    CONTINUE

C     LG- end

C     
 130  CONTINUE
C     
C     
C     *    2.            CALCULATE RAIN/SNOW FALL RATES
C     *                  CALCULATE MELTING OF SNOW
C     *                  CALCULATE EVAPORATION OF PRECIP
C     -------------------------------
C     
 200  CONTINUE
      DO 210 JL=KIDIA,KFDIA
         PRFL(JL)=0.
         PSFL(JL)=0.
         PRAIN(JL)=0.
 210  CONTINUE
      DO 220 JK=KTOPM2,KLEV
         DO 215 JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               PRAIN(JL)=PRAIN(JL)+PDMFUP(JL,JK)
               IF(PTEN(JL,JK).GT.TMELT) THEN
                  PRFL(JL)=PRFL(JL)+PDMFUP(JL,JK)+PDMFDP(JL,JK)
                  IF(PSFL(JL).GT.0..AND.PTEN(JL,JK).GT.ZTMELP2) THEN
                     ZFAC=ZCONS1*(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))
                     ZSNMLT=MIN(PSFL(JL),ZFAC*(PTEN(JL,JK)-ZTMELP2))
                     PDPMEL(JL,JK)=ZSNMLT
                     PSFL(JL)=PSFL(JL)-ZSNMLT
                     PRFL(JL)=PRFL(JL)+ZSNMLT
                  END IF
               ELSE
                  PSFL(JL)=PSFL(JL)+PDMFUP(JL,JK)+PDMFDP(JL,JK)
               END IF
            END IF
C     EVM 971023
            PRFLCK(JL,JK)=MAX(0.,PRFL(JL))+MAX(0.,PSFL(JL))
 215     CONTINUE
 220  CONTINUE
      DO 230 JL=KIDIA,KFDIA
         PRFL(JL)=MAX(PRFL(JL),0.)
         PSFL(JL)=MAX(PSFL(JL),0.)
         ZPSUBCL(JL)=PRFL(JL)+PSFL(JL)
 230  CONTINUE
      DO 240 JK=KTOPM2,KLEV
         DO 235 JL=KIDIA,KFDIA
            IF(LDCUM(JL).AND.JK.GE.KCBOT(JL).AND.
     1           ZPSUBCL(JL).GT.1.E-20) THEN
               ZRFL=ZPSUBCL(JL)
               ZRNEW=(MAX(0.,SQRT(ZRFL/ZCUCOV)-
     1              CEVAPCU(JK)*(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))*
     1              MAX(0.,PQSEN(JL,JK)-PQEN(JL,JK))))**2*ZCUCOV
               ZRMIN=ZRFL-ZCUCOV*MAX(0.,0.8*PQSEN(JL,JK)-PQEN(JL,JK))
     1              *ZCONS2*(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))
               ZRNEW=MAX(ZRNEW,ZRMIN)
               ZRFLN=MAX(ZRNEW,0.)
               ZDRFL=MIN(0.,ZRFLN-ZRFL)
               PDMFUP(JL,JK)=PDMFUP(JL,JK)+ZDRFL
               ZPSUBCL(JL)=ZRFLN
C     EVM 971023
               PRFLCK(JL,JK)=ZRFLN
            END IF
 235     CONTINUE
 240  CONTINUE
      DO 250 JL=KIDIA,KFDIA
         ZDPEVAP=ZPSUBCL(JL)-(PRFL(JL)+PSFL(JL))
         PRFL(JL)=PRFL(JL)+ZDPEVAP*PRFL(JL)*
     1        (1./MAX(1.E-20,PRFL(JL)+PSFL(JL)))
         PSFL(JL)=PSFL(JL)+ZDPEVAP*PSFL(JL)*
     1        (1./MAX(1.E-20,PRFL(JL)+PSFL(JL)))
C     EVM 971023
         PRFLCK(JL,KLEV)=MAX(0.,PRFL(JL))+MAX(0.,PSFL(JL))
 250  CONTINUE
C     
c     jhc      CALL UNLOC('CUFLX',99)
C     
      RETURN
      END
