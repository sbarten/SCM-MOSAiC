      SUBROUTINE CBM4_ECH (NSTEP,NSTOP,NPRINT,NLEV,NLEVEL,
     *                     PM,PMLOC,PMZZ,KG3X,PTMST,
     *                     MAE,PPHOTCHEM,PZ,PT,PQ,PDP,PP,PRHOA)

C__________________________________________________________________
C     *CBM4_ECH* - routine for day- and nighttime chemistry
C
C     Interface: *CBM4_ECH* is called from *PRECHEM*
C__________________________________________________________________
C
      IMPLICIT NONE
 
      REAL G

      PARAMETER (G=9.80665)

      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'      

      INCLUDE 'parreact_cbm4_ech.h'

      INTEGER KG3X,NSTEP,IM,JKM,JL,JR,JR2T21,JR1T21,
     &        JL1T21,JL2T21,JR1,JR2,JK,NLEV,NLEVEL,DLEVEL,
     &        FJR,FJL,JL1,JL2,INDEX,JT,JG,ITER,IW,NSTOP,
     &        NPRINT,JJK,MAE,KTHEIGHT,NTRAC_OUT,NREACTION(NTRAC)

      REAL PM(NLON,NLEVT,NTRAC),PMZZ(NLON,NTRAC),
     *     PRHOA(NLON,NLEVT),PZ(NLON,NLEVT),PT(NLON,NLEVT),
     *     PQ(NLON,NLEVT),PDP(NLON,NLEVT),PP(NLON,NLEVT),
     *     PPHOTCHEM(NLON,NLEVT,MAE),PMLOC(NLON,NLEVT,KG3X)

      REAL ZPM(NLON,NLEVEL+1,NTRAC),ZPMLOC(NLON,NLEVEL+1,KG3X),
     &     ZPRHOA(NLON,NLEVEL+1),ZPT(NLON,NLEVEL+1),ZPQ(NLON,NLEVEL+1),
     &     ZPP(NLON,NLEVEL+1),ZPPHOTCHEM(NLON,NLEVEL+1,NTRACT)

      REAL RJNO2(NLON),RJHNO3(NLON),RJO3D(NLON),RJH2O2(NLON),
     *     RJMEPE(NLON),RJBCH2O(NLON),RJACH2O(NLON),RJN2O5(NLON),
     *     RJANO3(NLON),RJBNO3(NLON),RJHNO4(NLON),RJHONO(NLON),
     *     RNOO3(NLON),RHO2NO(NLON),RMO2NO(NLON),RNO2OH(NLON),
     *     ROHHNO3(NLON),RNO2O3(NLON),RNONO3(NLON),RNO2NO3(NLON),
     *     RN2O5(NLON),RHNO4OH(NLON),RHNO4M(NLON),RNO2HO2(NLON),
     *     RODM(NLON),RH2OOD(NLON),RMOO2(NLON),RO3HO2(NLON),
     *     RCOOH(NLON),RO3OH(NLON),RHPOH(NLON),RFRMOH(NLON),
     *     RCH4OH(NLON),ROHPCAT(NLON),ROHPFRM(NLON),RMO2HO2(NLON),
     *     RHO2OH(NLON),RHO2HO2(NLON),RN2O5AQ(NLON),ROHSO2(NLON),
     *     ROHDMS(NLON),RNO3DMS(NLON)
      REAL AIR(NLON),O2(NLON),H2O(NLON),
     *  CH40(NLON),CO0(NLON),HNO30(NLON),H2O20(NLON),CH3O2H0(NLON),
     *  ZNO0(NLON),ZNO20(NLON),ZNO30(NLON),ZN2O50(NLON),HNO40(NLON),
     *  OH0(NLON),HO20(NLON),O30(NLON),OD0(NLON),CH3O0(NLON),
     *  CH3O20(NLON),CH2O0(NLON),ODDN0(NLON),
     *  DMS0(NLON),SO20(NLON),SAER0(NLON),RADON0(NLON),
     *  CH4(NLON),CO(NLON),HNO3(NLON),H2O2(NLON),CH3O2H(NLON),
     *  ZNO(NLON),ZNO2(NLON),ZNO3(NLON),ZN2O5(NLON),HNO4(NLON),
     *  OH(NLON),HO2(NLON),O3(NLON),OD(NLON),CH3O(NLON),
     *  CH3O2(NLON),CH2O(NLON),
     *  DMS(NLON),SO2(NLON),SAER(NLON),RADON(NLON)

C LG- extra declarations

      REAL ZNOY0(NLON),ZNOY(NLON),TOTN0(NLON),TOTN(NLON),
     *     OLDNOY,NEWNOY,CHECK,RSP_LTIME
     
      REAL TOTMASSN_OLD,TOTMASSN_NEW,TOTMASS_OLD,TOTMASS_NEW

      REAL RJPAN(NLON),RJALD2(NLON),RJACET(NLON),RJMGLY(NLON),
     *     RJMEK(NLON),RJNITR(NLON),RH2OH(NLON),ROHOH(NLON),
     *     RFRMNO3(NLON),RMO2MO2(NLON),RALD2OH(NLON),RALD2NO3(NLON),
     *     RC23NO(NLON),RC23NO2(NLON),EQPAN,RPAN(NLON),RC23C23(NLON),
     *     RC23HO2(NLON),RC23MO2(NLON),RETHOH(NLON),RETHO3(NLON),
     *     RPAROH(NLON),RRORA(NLON),RRORB(NLON),RRORNO2(NLON),
     *     RRXPAR(NLON),ROLEOH(NLON),ROLEO3(NLON),RMGLYOH(NLON),
     *     RISOPOH(NLON),RISOPO3(NLON),RISOPNO2(NLON),RISPDOH(NLON),
     *     RISPDO3(NLON),RISPNO3(NLON),RMTHCOH(NLON),RMTHCO3(NLON),
     *     RMTHCNO3(NLON),RMVKOH(NLON),RMVKO3(NLON),RMC23NO(NLON),
     *     RMC23HO2(NLON),RMC23NO2(NLON),RMPAN(NLON),RMPANO3(NLON),
     *     RMPANOH(NLON),RMEKOH(NLON),RACETOH(NLON),
     *     RXO2NO(NLON),ROLENO3(NLON),RISOPNO3(NLON),RISPDNO3(NLON),
     *     RXO2XO2(NLON),RXO2HO2(NLON),RXO2NNO(NLON),RXO2N(NLON),
     *     RXO2NHO2(NLON),RXO2NXO2(NLON),RBXO2NNO(NLON),RBXO2N(NLON),
     *     RBXO2NHO2(NLON),RBXO2NXO2(NLON),RBXO2NXO2N(NLON),

C LG-  added trace gases

     *     RISONTROH(NLON),RHCOOHOH(NLON),RCH3CO2HOH(NLON),
     *     RC23MO2A(NLON),RC23MO2B(NLON),ROHNH3(NLON),RNONH2(NLON),
     *     RNO2NH2(NLON),RHO2NH2(NLON),RO2NH2(NLON),RO3NH2(NLON),
     *     RSO4NH3(NLON),RMATERPOH(NLON),RMBTERPOH(NLON),RSQTERPOH(NLON),
     *     RMATERPO3(NLON),RMBTERPO3(NLON),RSQTERPO3(NLON),
     *     RMATERPNO3(NLON),RMBTERPNO3(NLON),RSQTERPNO3(NLON),
     *     RNOOH(NLON),RNO2H2O(NLON), 
     *     RCH2OHO2H(NLON),RRCHOHO2H(NLON)
     
      REAL O3S0(NLON),ALD20(NLON),PAR0(NLON),OLE0(NLON),ETH0(NLON),
     *     PAN0(NLON),ACET0(NLON),ISOP0(NLON),MGLY0(NLON),ISOPRD0(NLON),
     *     METHAC0(NLON),MVK0(NLON),MEK0(NLON),MPAN0(NLON),NITR0(NLON),
     *     C2O30(NLON),XO20(NLON),ROR0(NLON),XO2N0(NLON),RXPAR0(NLON),
     *     BXO2N0(NLON),MC3O30(NLON),O3S(NLON),ALD2(NLON),
     *     PAR(NLON),OLE(NLON),ETH(NLON),PAN(NLON),ACET(NLON),
     *     ISOP(NLON),MGLY(NLON),ISOPRD(NLON),METHAC(NLON),MVK(NLON),
     *     MEK(NLON),MPAN(NLON),NITR(NLON),C2O3(NLON),XO2(NLON),
     *     ROR(NLON),XO2N(NLON),RXPAR(NLON),BXO2N(NLON),MC3O3(NLON),

C LG- added trace gases

     *     ISONTR0(NLON),ISONTR(NLON),HCOOH0(NLON),HCOOH(NLON),
     *     CH3CO2H0(NLON),CH3CO2H(NLON),NH20(NLON),NH2(NLON),
     *     NH30(NLON),NH3(NLON),NH40(NLON),NH4(NLON),ACID0(NLON),
     *     ACID(NLON),MATERP0(NLON),MBTERP0(NLON),SQTERP0(NLON),
     *     HONO0(NLON),CH2OHO2H0(NLON),RCHOHO2H0(NLON),
     *     MATERP(NLON),MBTERP(NLON),SQTERP(NLON),HONO(NLON),
     *     CH2OHO2H(NLON),RCHOHO2H(NLON)

      REAL PCH3O2,XLCH3O2,PC2O3,XLC2O3,PXO2,XLXO2,PROR,XLROR,PRXPAR,
     *     XLRXPAR,PXO2N,XLXO2N,PBXO2N,XLBXO2N,POLE,XLOLE,PMGLY,XLMGLY,
     *     PISOP,XLISOP,PISOPRD,XLISOPRD,PMETHAC,XLMETHAC,PMVK,XLMVK,
     *     PMEK,XLMEK,PMC3O3,XLMC3O3,PCH3O2H,XLCH3O2H,PPAN,XLPAN,
     *     PMPAN,XLMPAN,PETH,XLETH,PACET,XLACET,PNTR,XLNTR,
     *     PALD2,XLALD2,PPAR,XLPAR,PO3S,XLO3S,

C LG- added trace gases

     *     PISONTR,XLISONTR,PHCOOH,XLHCOOH,PCH3CO2H,XLCH3CO2H,
     *     PNH2,XLNH2,XLNH3,PMSA,PMATERP,XLMATERP,PMBTERP,XLMBTERP,
     *     PSQTERP,XLSQTERP,PHONO,XLHONO,PCH2OHO2H,XLCH2OHO2H, 
     *     PRCHOHO2H,XLRCHOHO2H
     
      REAL ZF3BOD2,ZFARR,RX1,RX2,ER,ZTREC,ZF3BOD,CT,C0,XP,XL,DLT,
     &     ZMH2O,PKL,DT,PTMST,DT2,DGHNO3,DGAIR,GAMICE,GAMWAT,ZRHOA,ZN2,
     &     ZT3REC,RX3,RX4,EQN2O5,EQHNO4,P1,R12,R21,XL1,P2,
     &     XL2,P3,XL3,X1,X2,X3,C1,C2,C3,Y2,XJT,R21T,R12T,R12TC,
     &     R21TC,XJTC,ACUB,BCUB,CCUB,CUBDET,DNO2,R57,R56,
     &     R65,R75,P5,XL5,R66,X5,P6,XL6,X6,C6,XL7,C7,Y1,R89,
     &     P8,XL8,X4,C5,R98,XL9,R1011,R1012,C10,R1211,R1112,
     &     P11,X11,C11,C12,C1112,XLOD,PHNO3,XLHNO3,PH2O2,XLH2O2,
     &     PCH2O,XLCH2O,PCO,XLSO2,PSAER,XLDMS,DODDN,ZFAC,
     
C LG- 12092004+ added     

     &     RTERPO3
      
      REAL DPM(NLON,NLEVT,NTRAC),DPMLOC(NLON,NLEVT,KG3X)
      CHARACTER*10 CHDUM
      REAL RN2O5L(NLON)
      REAL RXN2O5M(32,64,10),RN25T30(19,96,48)
      DIMENSION JL1T21(96),JR1T21(48),JL2T21(96),JR2T21(48)
      REAL FJLT21(96),FJRT21(48)
      SAVE RN25T30,RXN2O5M
      SAVE JL1T21,JR1T21,JL2T21,JR2T21,FJLT21,FJRT21

C LG- extra declarations 

      REAL O1D,OPP(NLEVT),POHHO2(NLON,NLEVT),PHO2OH(NLON,NLEVT),
     &     H2OLTR

      LOGICAL LWRITEREAC
      PARAMETER(LWRITEREAC=.TRUE.)

C LG- end

C - - -
C    Functions and molar weights
C

C LG- declaration of functions which are used in the EBI scheme

      ZFARR(RX1,ER,ZTREC)=RX1*EXP(ER*ZTREC)
      ZF3BOD(RX1,RX2)=RX1/(1+RX1/RX2)*0.6**
     *   (1./(1+ALOG10(RX1/RX2)**2))
      ZF3BOD2(RX1,RX2)=RX1/(1+RX1/RX2)*0.3**
     *   (1./(1+ALOG10(RX1/RX2)**2))
      CT(C0,XP,XL,DLT)=(C0-XP/XL)*EXP(-1.*XL*DLT)+XP/XL

C LG-

      ZMH2O=18.

      WRITE(NUNMDFL,*)'Start CBM4_ECH.f'

C LG- definition of number of extra levels for which the chemistry 
C     must be calculated

      DLEVEL=0
      IF (LXTMZZ) DLEVEL=DLEVEL+1	! chemistry at extra reference heigth zz

C LG-  

C - - - 
C    Read rx-rates for N2O5-reaction on aerosols
C    Taken from Dentener's MOGUNTIA
C    Values in 1000*s-1
C

      IF (NSTEP.EQ.NRESUM) THEN
c read field on MOGUNTIA grid
       DO 90 IM=1,IMON
        READ (119,'(A)') CHDUM
       DO 90 JKM=10,1,-1
       DO 90 JL=1,64
        READ (119,4999) (RXN2O5M(JR,JL,JKM),JR=1,32)
   90  CONTINUE
 4999  FORMAT (16F6.3)

c interpolate to T30
c get interpolation parameters latitude

C LG-  definining proper latitude

       DO 93 JR=ILAT,ILAT
        FJRT21(JR)=0.
        IF (MOD(JR,3).EQ.0) THEN
         FJRT21(JR)=0.
         JR2T21(JR)=JR*2/3
         JR1T21(JR)=JR2T21(JR)
        ELSE IF (MOD(JR,3).EQ.1) THEN          
         FJRT21(JR)=1.
         JR1T21(JR)=(JR-1)*2/3+1
         JR2T21(JR)=JR1T21(JR)
        ELSE IF (MOD(JR,3).EQ.2) THEN          
         FJRT21(JR)=.5
         JR1T21(JR)=(JR-2)*2/3+1
         JR2T21(JR)=JR1T21(JR)+1
        ENDIF
   93  CONTINUE

c get interpolation parameters longitude

C LG-  definining proper longitude

       DO 92 JL=ILON,ILON
        FJLT21(JL)=0
        IF (MOD(JL,3).EQ.0) THEN
         FJLT21(JL)=0.
         JL2T21(JL)=JL*2/3
         JL1T21(JL)=JL2T21(JL)
        ELSE IF (MOD(JL,3).EQ.1) THEN          
         FJLT21(JL)=1.
         JL1T21(JL)=(JL-1)*2/3+1
         JL2T21(JL)=JL1T21(JL)
        ELSE IF (MOD(JL,3).EQ.2) THEN          
         FJLT21(JL)=.5
         JL1T21(JL)=(JL-2)*2/3+1
         JL2T21(JL)=JL1T21(JL)+1
        ENDIF
   92  CONTINUE
       CALL RESETR (RN25T30,19*96*48,0.)
      ENDIF
      IF (NSTEP.EQ.NRESUM) THEN
        JR=ILAT
        FJR=FJRT21(JR)
        JR1=JR1T21(JR)
        JR2=JR2T21(JR)
       DO 91 JK=1,NLEV

C LG-  defining proper longitude

       DO 91 JL=ILON,ILON
        FJL=FJLT21(JL)
        JL1=JL1T21(JL)
        JL2=JL2T21(JL)

C LG-   the pressure is only 1-D in the longitude!!

        PKL=PP(1,JK)/100.
        INDEX=1
        IF (PKL.LE.50) THEN
         INDEX=1
        ELSE IF (PKL.GT.950) THEN
         INDEX=10
        ELSE
         INDEX=NINT(PKL/100.)
        ENDIF
        RN25T30(JK,JL,ILAT)=FJR*FJL*RXN2O5M(JR1,JL1,INDEX)+
     *    FJR*(1.-FJL)*RXN2O5M(JR1,JL2,INDEX)+
     *    (1.-FJR)*FJL*RXN2O5M(JR2,JL1,INDEX)+
     *    (1.-FJR)*(1.-FJL)*RXN2O5M(JR2,JL2,INDEX)
        RN25T30(JK,JL,ILAT)=RN25T30(JK,JL,ILAT)*1.E-3
   91  CONTINUE
      ENDIF

C - - -  
C  Prepare budgetcalculation

      DO 302 JT=1,NTRAC
      DO 302 JK=1,NLEVEL
      DO 302 JL=1,NLON
	IF(LBULKVEG.AND.JK.EQ.NLEV+1) THEN
         DPM(JL,JK,JT)=-1.*PM(JL,JK,JT)
	ELSEIF(LVEG_MLAY.AND.JK.GE.NLEV+1.OR.
     &         LSNOW_MLAY.AND.JK.GE.NLEV+1) THEN
         DPM(JL,JK,JT)=-1.*PM(JL,JK,JT)
        ELSE
         DPM(JL,JK,JT)=-1.*PM(JL,JK,JT)
	ENDIF	  
  302 CONTINUE

      DO 303 JG=1,KG3X
      DO 303 JK=1,NLEVEL
      DO 303 JL=1,NLON
        IF (LBULKVEG.AND.JK.EQ.NLEV+1) THEN
         DPMLOC(JL,JK,JG)=-1.*PMLOC(JL,JK,JG)	 
	ELSEIF(LVEG_MLAY.AND.JK.GE.NLEV+1.OR.
     &         LSNOW_MLAY.AND.JK.GE.NLEV+1) THEN
         DPMLOC(JL,JK,JG)=-1.*PMLOC(JL,JK,JG)
        ELSE
         DPMLOC(JL,JK,JG)=-1.*PMLOC(JL,JK,JG)
        ENDIF
  303 CONTINUE

C     -------
C     Some other constants; gamma values for N2O5 on clouds
C
      DT=PTMST
      DT2=DT**2
      DGHNO3=0.136
      DGAIR=0.133
      GAMICE=0.25
      GAMWAT=0.045

C LG- definition of N2O5 chemistry parameters which was originally
C     within the loop 1104

      DO 1004 JK=1,NLEV
      DO 1004 JL=1,NLON

C aerosol N2O5 -> 2 HNO3
C

C LG-   defining the proper longitude

        RN2O5AQ(JL)=RN25T30(JK,ILON,ILAT)

C cloud N2O5 -> 2 HNO3
C

C LG-   the proper longitude is defined

C LG-   09-12-97 the proper definition and assignment of ZLCOVER and
C       WATFRAC needs to be checked

        RN2O5L(JL)=ZCLCOVER(ILON,JK)/DT*
     *      (WATFRAC(ILON,JK)*GAMWAT+(1.-WATFRAC(ILON,JK))*GAMICE)        

 1004 CONTINUE

C LG- for additional chemistry calculations (at reference height zz and
C     within the bulk vegetation layer (both only possible for the 
C     "big leaf" approach and thus LBIOSPH=.FALSE.), the applied
C     parameter arrays are assigned the values of the first NLEVEL's 
C     for the first NLEVEL calculations and for the other DLEVEL 
C     calculations, the values representing level ZZ or the bulk 
C     vegetation layer are assigned  

      DO 99 JK=1,NLEVEL+DLEVEL
      DO 99 JL=1,NLON

C LG- assigning the tracer concentrations 

       DO JT=1,NTRAC
        IF (LBULKVEG.AND.JK.EQ.NLEV+1) THEN
         ZPM(JL,JK,JT)=PM(JL,JK,JT)
	ELSEIF (LVEG_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1.OR.
     &          LSNOW_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1) THEN
         ZPM(JL,JK,JT)=PM(JL,JK,JT)  
	ELSEIF (LXTMZZ.AND.JK.EQ.NLEVEL+1) THEN
         ZPM(JL,JK,JT)=PMZZ(JL,JT)
        ELSE
         ZPM(JL,JK,JT)=MAX(1.E-30,PM(JL,JK,JT)) ! avoiding numerical problems
	ENDIF
       ENDDO

C LG-  and the short-lived species

       DO JG=1,NG3X
        IF (LBULKVEG.AND.JK.EQ.NLEV+1) THEN
         ZPMLOC(JL,JK,JG)=PMLOC(JL,JK,JG)
	ELSEIF (LVEG_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1.OR.
     &          LSNOW_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1) THEN
         ZPMLOC(JL,JK,JG)=PMLOC(JL,JK,JG)
	ELSEIF (LXTMZZ.AND.JK.EQ.NLEVEL+1) THEN
         ZPMLOC(JL,JK,JG)=PMZZ(JL,JG)
        ELSE
         ZPMLOC(JL,JK,JG)=PMLOC(JL,JK,JG)
	ENDIF
       ENDDO

C LG-  assigning some other parameters, e.g. temperature etc.

       IF (LXTMZZ) THEN
        ZPRHOA(JL,JK)=PRHOA(JL,MIN(NLEV,JK))
        ZPQ(JL,JK)=PQ(JL,MIN(NLEV,JK))
        ZPT(JL,JK)=PT(JL,MIN(NLEV,JK))
        ZPP(JL,JK)=PP(JL,MIN(NLEV,JK))
       ELSE
        ZPRHOA(JL,JK)=PRHOA(JL,JK)
        ZPQ(JL,JK)=PQ(JL,JK)
        ZPT(JL,JK)=PT(JL,JK)
        ZPP(JL,JK)=PP(JL,JK)
       ENDIF

C LG-  assigning the photodissociation rates 

       DO JT=1,MAE
        IF (LXTMZZ.AND.JK.EQ.NLEVEL+1) THEN       
	 ZPPHOTCHEM(JL,JK,JT)=PPHOTCHEM(JL,NLEV,JT)
        ELSE 
         ZPPHOTCHEM(JL,JK,JT)=PPHOTCHEM(JL,JK,JT)
        ENDIF         
       ENDDO

   99 CONTINUE

C LG- end

C     --------
C     1. Start loop
C
  100 CONTINUE

C LG- the chemistry is default only calculated for the troposphere!

      IF (NSTEP.EQ.0)
     &  WRITE(NUNMDFL,'(1a)')
     &  ' The chemistry is only resolved within the troposphere'        

      KTHEIGHT=KTRHE(1)

      DO 101 JK=KTHEIGHT,NLEVEL+DLEVEL

C ** 1.1 ** prepare chemistry
C
      DO 1102 JL=1,NLON
        ZRHOA=ZPRHOA(JL,JK)
        AIR(JL)=ZRHOA*6.022045E23/ZMAIR
        O2(JL)=0.209476*AIR(JL)
        ZN2=0.78084*AIR(JL)
        H2O(JL)=ZPQ(JL,JK)*ZRHOA*AVO/ZMH2O  ! molecules cm-3
        H2O(JL)=MAX(H2O(JL),0.1)

C LG-   added the calculation of the available water in liters!

        H2OLTR=ZPQ(JL,JK)*ZRHOA*1.E-3 ! zpq in g H2O g-1 air and zrhoa 
	     ! in g air cm-3 gives g H2O cm-3. 1 g = 1e-3 kg = 1 liter

C LG-   01-2004: assigning the amount of water, H2O, to parameter ZH2O 
C       for diagnostics

        ZH2O(JL,JK)=1.E-17*H2O(JL)

 1102 CONTINUE

C - - -
C      Evaluation of reactionrates
C
C rjach2o: ch2o -> co+h2
C rjbch2o: ch2o -> co+2ho2
C rjano3:  no3 -> no2+o3
C rjbno3:  no3 -> no

      DO 1103 JL=1,NLON
        IF (ZPPHOTCHEM(JL,JK,7).GT.1E-20) THEN
          RJH2O2(JL)=ZPPHOTCHEM(JL,JK,1)
          RJHNO3(JL)=ZPPHOTCHEM(JL,JK,2)
          RJNO2(JL)=ZPPHOTCHEM(JL,JK,3)
          RJN2O5(JL)=ZPPHOTCHEM(JL,JK,4)
          RJACH2O(JL)=ZPPHOTCHEM(JL,JK,5)
          RJBCH2O(JL)=ZPPHOTCHEM(JL,JK,6)
          RJO3D(JL)=ZPPHOTCHEM(JL,JK,7)
          RJMEPE(JL)=ZPPHOTCHEM(JL,JK,8)
          RJHNO4(JL)=ZPPHOTCHEM(JL,JK,9)
          RJANO3(JL)=ZPPHOTCHEM(JL,JK,10)
          RJBNO3(JL)=ZPPHOTCHEM(JL,JK,11)
          RJPAN(JL)=ZPPHOTCHEM(JL,JK,12)
          RJALD2(JL)=ZPPHOTCHEM(JL,JK,13)
          RJACET(JL)=ZPPHOTCHEM(JL,JK,14)
          RJMGLY(JL)=ZPPHOTCHEM(JL,JK,15)
        ELSE
          RJNO2(JL)=0.
          RJHNO3(JL)=0.
          RJO3D(JL)=0.
          RJH2O2(JL)=0.
          RJMEPE(JL)=0.
          RJBCH2O(JL)=0.
          RJACH2O(JL)=0.
          RJN2O5(JL)=0.
          RJANO3(JL)=0.
          RJBNO3(JL)=0.
          RJHNO4(JL)=0.
          RJPAN(JL)=0.
          RJALD2(JL)=0.
          RJACET(JL)=0.
          RJMGLY(JL)=0.
        ENDIF
        RJMEK(JL)=3.6E-4*RJNO2(JL)
        RJNITR(JL)=4.8*RJHNO3(JL)

C LG-   12-2003, added the photodissociation of HONO (see personal
C       communication by Yvonne Trebs:
C       J(HONO)=0.189*J(NO2)+8.433*10^(-2)*J(NO2)^2  
C       Kraus and Hofzumahaus (1998)

        RJHONO(JL)=0.189*RJNO2(JL)+8.433E-2*RJNO2(JL)**2

C LG-   end

 1103 CONTINUE

      DO 1104 JL=1,NLON
        ZTREC=1./ZPT(JL,JK)
        ZT3REC=300./ZPT(JL,JK)
        RNOO3(JL)=ZFARR(2.E-12,-1400.,ZTREC)
        RHO2NO(JL)=ZFARR(3.5E-12,250.,ZTREC)
        RMO2NO(JL)=ZFARR(3.0E-12,280.,ZTREC)
           RX1=2.5E-30*ZT3REC**4.4*AIR(JL)
           RX2=1.6E-11*ZT3REC**1.7
        RNO2OH(JL)=ZF3BOD(RX1,RX2)
           RX1=ZFARR(7.2E-15,785.,ZTREC)
           RX3=ZFARR(1.9E-33,725.,ZTREC)
           RX4=ZFARR(4.1E-16,1440.,ZTREC)
        ROHHNO3(JL)=RX1+RX3*AIR(JL)/(1.+RX3*AIR(JL)/RX4)
        RNO2O3(JL)=ZFARR(1.2E-13,-2450.,ZTREC)
        RNONO3(JL)=ZFARR(1.5E-11,170.,ZTREC)
           RX1=2.2E-30*ZT3REC**3.9*AIR(JL)
           RX2=1.5E-12*ZT3REC**0.7
        RNO2NO3(JL)=ZF3BOD(RX1,RX2)
          EQN2O5=4.E-27*EXP(10930.*ZTREC)
        RN2O5(JL)=RNO2NO3(JL)/EQN2O5
        RHNO4OH(JL)=ZFARR(1.3E-12,380.,ZTREC)
           RX1=1.8E-31*ZT3REC**3.2*AIR(JL)
           RX2=4.7E-12*ZT3REC**1.4
        RNO2HO2(JL)=ZF3BOD(RX1,RX2)
          EQHNO4=2.1E-27*EXP(10900.*ZTREC)

C LG-   15092004- sensitivity analysis of role of thermal decomposition
C       of HNO4 for (nocturnal) HO2 production

        IF (RG.LT.1e-10) THEN
          RHNO4M(JL)=RNO2HO2(JL)/EQHNO4 !0.
	ELSE
	  RHNO4M(JL)=RNO2HO2(JL)/EQHNO4
	ENDIF

        RODM(JL)=0.2094*ZFARR(3.2E-11,70.,ZTREC)
     *        +0.7808*ZFARR(1.8E-11,110.,ZTREC)
        RH2OOD(JL)=2.2E-10
        RH2OH(JL)=0.5E-6*ZFARR(5.5E-12,-2000.,ZTREC)
        ROHOH(JL)=ZFARR(4.3E-12,-240.,ZTREC)
        RMOO2(JL)=ZFARR(3.9E-14,-900.,ZTREC)
        RO3HO2(JL)=ZFARR(1.1E-14,-500.,ZTREC)
        RCOOH(JL)=1.5E-13*(1.+0.6*ZPP(JL,JK)/101325.)
        RO3OH(JL)=ZFARR(1.6E-12,-940.,ZTREC)
        RHPOH(JL)=ZFARR(2.9E-12,-160.,ZTREC)
        RFRMOH(JL)=1.0E-11
        RFRMNO3(JL)=ZFARR(3.4E-13,-1900.,ZTREC)
        RCH4OH(JL)=ZFARR(2.45E-12,-1775.,ZTREC)
        ROHPCAT(JL)=0.7*ZFARR(3.8E-12,200.,ZTREC)
        ROHPFRM(JL)=0.3*ZFARR(3.8E-12,200.,ZTREC)
        RMO2HO2(JL)=ZFARR(3.8E-13,800.,ZTREC)
        RMO2MO2(JL)=ZFARR(9.1E-14,416.,ZTREC)
        RHO2OH(JL)=ZFARR(4.8E-11,250.,ZTREC)
        RHO2HO2(JL)=(ZFARR(2.3E-13,600.,ZTREC)+
     *               ZFARR(1.7E-33*AIR(JL),1000.,ZTREC))*
     *              (1.+H2O(JL)*ZFARR(1.4E-21,2200.,ZTREC))

C aerosol N2O5 -> 2 HNO3
C

C LG-   defining the proper longitude

        IF (JK.LE.NLEV) THEN
          RN2O5AQ(JL)=RN25T30(JK,ILON,ILAT)
        ELSE 
          RN2O5AQ(JL)=RN25T30(NLEV,ILON,ILAT)       
        ENDIF

C cloud N2O5 -> 2 HNO3
C

C LG-   the proper longitude is defined

        IF (JK.LE.NLEV) THEN
          RN2O5L(JL)=ZCLCOVER(ILON,JK)/DT*
     *      (WATFRAC(ILON,JK)*GAMWAT+(1.-WATFRAC(ILON,JK))*GAMICE) 
        ELSE 
          RN2O5L(JL)=ZCLCOVER(ILON,NLEV)/DT*
     *      (WATFRAC(ILON,NLEV)*GAMWAT+(1.-WATFRAC(ILON,NLEV))*GAMICE)       
        ENDIF

C O1D steady state
        RJO3D(JL)=RJO3D(JL)*H2O(JL)*RH2OOD(JL)/
     *     (H2O(JL)*RH2OOD(JL)+AIR(JL)*RODM(JL))

C organic reaction rates;
c taken from Duncan and Chameides (1998); Stockwell ea (1997);
c Houweling ea (1998)
       RALD2OH(JL)=ZFARR(7.E-12,250.,ZTREC)
       RALD2NO3(JL)=ZFARR(1.4E-12,-1900.,ZTREC)
       RC23NO(JL)=ZFARR(3.5E-11,-180.,ZTREC)
         RX1=9.7E-29*ZT3REC**5.6*AIR(JL)
         RX2=9.3E-12*ZT3REC**1.5
       RC23NO2(JL)=ZF3BOD2(RX1,RX2)
         EQPAN=8.62E-29*EXP(13954.*ZTREC)
       RPAN(JL)=RC23NO2(JL)/EQPAN
       RC23C23(JL)=ZFARR(2.5E-12,500.,ZTREC)
       RC23HO2(JL)=6.5E-12
       RC23MO2(JL)=6.5e-12
         RX1=1.E-28*ZT3REC**0.8*AIR(JL)
         RX2=8.8E-12
       RETHOH(JL)=ZF3BOD(RX1,RX2)
       RETHO3(JL)=ZFARR(9.14e-15,-2580.,ZTREC)
       RPAROH(JL)=8.1E-13
       RRORA(JL)=ZFARR(1.E15,-8000.,ZTREC)
       RRORB(JL)=1.6E3
       RRORNO2(JL)=1.5E-11
       RRXPAR(JL)=8E-11
       ROLEOH(JL)=ZFARR(5.2E-12,504.,ZTREC)
       ROLEO3(JL)=ZFARR(1.4E-14,-2100.,ZTREC)
       ROLENO3(JL)=7.7E-15
       RMGLYOH(JL)=1.7E-11
       RISOPOH(JL)=ZFARR(2.54E-11,410.,ZTREC)
       RISOPO3(JL)=ZFARR(12.3E-15,-2013.,ZTREC)
       RISOPNO3(JL)=ZFARR(4.E-12,-446.,ZTREC)
       RISOPNO2(JL)=1.5E-19
       RISPDOH(JL)=6.1E-11
       RISPDO3(JL)=4.2E-18
       RISPDNO3(JL)=1.E-13
       RMTHCOH(JL)=ZFARR(1.9E-11,176.,ZTREC)
       RMTHCO3(JL)=ZFARR(1.4E-15,-2114.,ZTREC)
       RMTHCNO3(JL)=ZFARR(1.5E-12,-1726.,ZTREC)
       RMVKOH(JL)=ZFARR(4.1E-12,453.,ZTREC)
       RMVKO3(JL)=ZFARR(7.5E-16,-1520.,ZTREC)
       RMC23NO(JL)=RC23NO(JL)
       RMC23HO2(JL)=RC23HO2(JL)
       RMC23NO2(JL)=RC23NO2(JL)
       RMPAN(JL)=RPAN(JL)
       RMPANO3(JL)=8.2E-18
       RMPANOH(JL)=3.6E-12
       RMEKOH(JL)=ZFARR(2.9E-13,413.,ZTREC)

! mz_lg_20060215+ check the mecca code for the MEK+OH reaction rate which 
!      appears to be a factor 10 larger compared to values of the CBM4 scheme
!       RMEKOH(JL)=ZFARR(1.3e-12,-25.,ZTREC) {&1207}

       RACETOH(JL)=ZFARR(1.7E-12,-600.,ZTREC)
       RXO2NO(JL)=8.1E-12
       RXO2XO2(JL)=ZFARR(1.7E-14,1300.,ZTREC)
       RXO2HO2(JL)=ZFARR(7.7E-14,1300.,ZTREC)
       RXO2NNO(JL)=8.1E-12
       RXO2N(JL)=ZFARR(1.7E-14,1300.,ZTREC)
       RXO2NHO2(JL)=ZFARR(7.7E-14,1300.,ZTREC)
       RXO2NXO2(JL)=ZFARR(3.5E-14,1300.,ZTREC)

C WP-  changed RBXO2NNO

       RBXO2NNO(JL)=ZFARR(0.5*4.2E-12,180.,ZTREC)

C WP-  end

       RBXO2N(JL)=ZFARR(1.7E-14,1300.,ZTREC)
       RBXO2NHO2(JL)=ZFARR(7.7E-14,1300.,ZTREC)
       RBXO2NXO2(JL)=0.25*8.1E-12
       RBXO2NXO2N(JL)=ZFARR(3.4E-14,1300.,ZTREC)

C LG-  sulfur chemistry reaction rates

           RX1=3.0E-31*ZT3REC**3.3*AIR(JL)
           RX2=1.5E-12
       ROHSO2(JL)=ZF3BOD(RX1,RX2)
       ROHDMS(JL)=ZFARR(1.2E-11,-260.,ZTREC)
       RNO3DMS(JL)=ZFARR(1.9E-13,500.,ZTREC)

C WP-  added isontr+oh

       RISONTROH(JL)=3.E-11

C WP-  end

C LG-  added hcooh+oh, and ch3coh2+oh

       RHCOOHOH(JL)=4.5E-13
       RCH3CO2HOH(JL)=ZFARR(4.E-13,200.,ZTREC)

C LG- updated reaction rates for c2o3 reactions

       RC23HO2(JL)=ZFARR(1.3E-13,1040.,ZTREC)
           RX1=9.8E-12
           RX2=ZFARR(2.2E6,-3870.,ZTREC)
       RC23MO2A(JL)=RX1/(1.+RX2)
       RC23MO2B(JL)=RX1*(1.+1./RX2)

C LG-  end       

C LG-  ammonia reaction rates

       ROHNH3(JL)=ZFARR(1.7E-12,-710.,ZTREC) !1.56e-13 at 298K
       RNONH2(JL)=ZFARR(3.8E-12,+450.,ZTREC) !1.72e-11
       RNO2NH2(JL)=ZFARR(2.1E-12,650.,ZTREC) !1.86e-11
       RHO2NH2(JL)=3.4E-11
       RO2NH2(JL)=6.0E-21
       RO3NH2(JL)=ZFARR(4.3E-12,-930.,ZTREC) !1.89e-13 at 298K

C LG-  see TM3 code for definition of the parameters

! knh3so4 is uptake coefficient on H2SO4. 1 uptake of NH3 consumes 1 acid molecule.  
! 
c  rr(jl,knh3so4)=het_nh3(jl)/1e-9/y(jl,iair) 

       RSO4NH3(JL)=0.

C LG-  added the monoterpenes destruction rates (Kostas Tsigaridis, personal 
C      communication, March, 2002)

       RMATERPOH(JL)=ZFARR(12.1E-12,444.,ZTREC)
       RMBTERPOH(JL)=ZFARR(23.8E-12,357.,ZTREC)
       RMATERPO3(JL)=ZFARR(1.01E-15,-732.,ZTREC)
       RMBTERPO3(JL)=1.5E-17
       RMATERPNO3(JL)=ZFARR(1.19E-12,490.,ZTREC)
       RMBTERPNO3(JL)=2.51E-12

C LG-  added the sesquiterpenes (personal communications with Boris Bonn, 
C      April, 2003)

       RSQTERPOH(JL)=2.9E-10
       RSQTERPO3(JL)=1170.E-17
       RSQTERPNO3(JL)=3.5E-11

C LG-  end
 
C LG-  12-2003, added some of the HONO chemistry
C      Rate constant: NO+OH+M --> HONO+M
C      9*10^(-12) cm^3 molecule^(-1) s^(-1) Seinfeld&Pandis (1998)

       RNOOH(JL)=9.E-12

C      Rate constant: 2NO2+H2O--> HONO+HNO3 (1st order reaction)
C      5.6*10^(-6) *(100/mixing height [m]) *s^(-1) Harrison et al. (1996)

C LG-  this seems a highly parameterized reaction rate to consider the
C      role of the heterogeneous production of HONO using the liquid water 
C      and PBL height. For a small PBL height the reaction rate is largest, 
C      likely expressing indirectly the presence of dew in the inversion layer
C      Since some of these parameters are actually available in the SCM, we 
C      could consider using a more explicit representation of this chemistry
C      component. The units of the reaction rate is in s-1 which indicates 
C      that in order to arrive at an HONO production rate in molecules cm-3 s-1
C      the reaction rate is only multiplied with the NO2 concentration. 
C      Obviously the parameterization is based on an assumption of the surface
C      layer moisture (H2O concentration, +/- 6e17 molecules cm-3), possibly
C      expressed also by the PBL height term

       IF (JK.GE.NLEV) THEN ! only considering the reaction in SL and canopy!
         RNO2H2O(JL)=5.6E-6*MIN(1.,100./PBLHGT(JL)) ! taking a max. value of 1 
       ELSE
         RNO2H2O(JL)=0.
       ENDIF

C LG-  see RC.49, pp174 of thesis of Tanja Winterrath. The rate is given in 
C      (mol l-1)-1 s-1, and which consequently must be recalculated to
C      molecules-1 s-1 by including avo and the amount of water available for 
C      heterogeneous chemistry

C ???? check further
C        RNO2H2O(JL)=ZFARR(8.4E7,-2900.,ZTREC)/(AVO*H2OLTR)
C 
C        IF (JK.EQ.NLEV) print *,'cbm4_ech',5.6E-6*MIN(1.,100./PBLHGT(JL)),
C      &    RNO2H2O(JL),H2OLTR,ZPQ(JL,JK),ZRHOA

C LG-  end

C LG-  added the decomposition rate of HMHP into HCOOH and H2O (see Peter Neebs
C      thesis, Tabel 3.15, pp 98) and that of HAHP into H2O2 and RCHO, see page 108
C      (Chapter 6, summary and comparison..., Valverde's thesis), k1, k2

       RCH2OHO2H(JL)=6.E-4

       IF (RG.LT.1e-10) THEN
         RRCHOHO2H(JL)=3.5E-3 ! 0. for sensitivity analysis, nocturnal H2O2 exchanges 
       ELSE
         RRCHOHO2H(JL)=3.5E-3 
       ENDIF

 1104 CONTINUE

C - - -
C      Set concentrations for time 0
C
      DO 1105 JL=1,NLON
        O30(JL)=ZPM(JL,JK,io3)
        CH40(JL)=ZPM(JL,JK,ich4)
        CO0(JL)=ZPM(JL,JK,ico)
        HNO30(JL)=ZPM(JL,JK,ihno3)
        H2O20(JL)=ZPM(JL,JK,ih2o2)
        CH3O2H0(JL)=ZPM(JL,JK,ich3o2h)
        CH2O0(JL)=ZPM(JL,JK,ich2o)
        O3S0(JL)=ZPM(JL,JK,io3s)
        ALD20(JL)=ZPM(JL,JK,iald2)
        PAR0(JL)=ZPM(JL,JK,ipar)
        OLE0(JL)=ZPM(JL,JK,iole)
        ETH0(JL)=ZPM(JL,JK,ieth)
        PAN0(JL)=ZPM(JL,JK,ipan)
        ACET0(JL)=ZPM(JL,JK,iacet)
        ISOP0(JL)=ZPM(JL,JK,iisop)
        MGLY0(JL)=ZPM(JL,JK,imgly)
        ISOPRD0(JL)=ZPM(JL,JK,iisoprd)
        METHAC0(JL)=ZPM(JL,JK,imethac)
        MVK0(JL)=ZPM(JL,JK,imvk)
        MEK0(JL)=ZPM(JL,JK,imek)
        MPAN0(JL)=ZPM(JL,JK,impan)
        NITR0(JL)=ZPM(JL,JK,intr)

C LG-   sulfur chemistry

        DMS0(JL)=ZPM(JL,JK,idms)
        SO20(JL)=ZPM(JL,JK,iso2)
        SAER0(JL)=ZPM(JL,JK,iso4)

C LG-   radon

	RADON0(JL)=ZPM(JL,JK,irad)

c WP-   added isontr

        ISONTR0(JL)=ZPM(JL,JK,iisontr)

c WP-   end

c LG-   added formic acid and acetic acid

        HCOOH0(JL)=ZPM(JL,JK,ihcooh)
        CH3CO2H0(JL)=ZPM(JL,JK,ich3co2h)

c LG-   added ammonia chemistry

        NH20(JL)=ZPM(JL,JK,inh2)
        NH30(JL)=ZPM(JL,JK,inh3)
        NH40(JL)=ZPM(JL,JK,inh4)

C LG-   initialisation of H+ concentration from sulfate and ammonium conc.

        ACID0(JL)=MAX(2.*SAER0(JL)-NH40(JL),0.)

C LG-   added the monoterpenes as a group

        MATERP0(JL)=ZPM(JL,JK,imaterp)
        MBTERP0(JL)=ZPM(JL,JK,imbterp)

C LG-   added the sesquiterpenes as a group

        SQTERP0(JL)=ZPM(JL,JK,isqterp)

C LG-   added HONO

        HONO0(JL)=ZPM(JL,JK,ihono)

C LG-   added hydroxy-methyl/alkylhydroxyperoxide

        CH2OHO2H0(JL)=ZPM(JL,JK,ich2oho2h)
        RCHOHO2H0(JL)=ZPM(JL,JK,irchoho2h)

C LG-   and short-lived species

        IF (NTRAC.EQ.NTRACT) THEN
         ZNO0(JL)=ZPM(JL,JK,ino)
         ZNO20(JL)=ZPM(JL,JK,ino2)
         ZNO30(JL)=MAX(ZPM(JL,JK,ino3),1.E-30)
         ZN2O50(JL)=MAX(ZPM(JL,JK,in2o5),1.E-30)
         HNO40(JL)=MAX(ZPM(JL,JK,ihno4),1.E-30)
         OH0(JL)=ZPM(JL,JK,ioh)
         HO20(JL)=ZPM(JL,JK,iho2)
         CH3O20(JL)=ZPM(JL,JK,ich3o2)
         C2O30(JL)=ZPM(JL,JK,ic2o3)
         XO20(JL)=ZPM(JL,JK,ixo2)
         ROR0(JL)=ZPM(JL,JK,iror)
         XO2N0(JL)=ZPM(JL,JK,ixo2n)
         RXPAR0(JL)=ZPM(JL,JK,irxpar)
         RXPAR0(JL)=MIN(RXPAR0(JL),PAR0(JL))
         BXO2N0(JL)=ZPM(JL,JK,ibxo2n)
         MC3O30(JL)=ZPM(JL,JK,imc3o3)
        ELSE

         IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
          ZNO0(JL)=ZPM(JL,JK,ino)
          ZNO20(JL)=ZPM(JL,JK,ino2)
          ZNO30(JL)=MAX(ZPM(JL,JK,ino3),1.E-30)
          ZN2O50(JL)=MAX(ZPM(JL,JK,in2o5),1.E-30)
          HNO40(JL)=MAX(ZPM(JL,JK,ihno4),1.E-30)
         ELSE
          ZNO0(JL)=ZPMLOC(JL,JK,ino)
          ZNO20(JL)=ZPMLOC(JL,JK,ino2)
          ZNO30(JL)=MAX(ZPMLOC(JL,JK,ino3),1.E-30)
          ZN2O50(JL)=MAX(ZPMLOC(JL,JK,in2o5),1.E-30)
          HNO40(JL)=MAX(ZPMLOC(JL,JK,ihno4),1.E-30)
         ENDIF

         OH0(JL)=ZPMLOC(JL,JK,ioh)
         HO20(JL)=ZPMLOC(JL,JK,iho2)
         CH3O20(JL)=ZPMLOC(JL,JK,ich3o2)
         C2O30(JL)=ZPMLOC(JL,JK,ic2o3)
         XO20(JL)=ZPMLOC(JL,JK,ixo2)
         ROR0(JL)=ZPMLOC(JL,JK,iror)
         XO2N0(JL)=ZPMLOC(JL,JK,ixo2n)
         RXPAR0(JL)=ZPMLOC(JL,JK,irxpar)
         RXPAR0(JL)=MIN(RXPAR0(JL),PAR0(JL))
         BXO2N0(JL)=ZPMLOC(JL,JK,ibxo2n)
         MC3O30(JL)=ZPMLOC(JL,JK,imc3o3)
        ENDIF
   
c LG-   end

C - - -
C      Set concentrations for time t
C
        O3(JL)=O30(JL)
        CH4(JL)=CH40(JL)
        CO(JL)=CO0(JL)
        HNO3(JL)=HNO30(JL)
        H2O2(JL)=H2O20(JL)
        CH3O2H(JL)=CH3O2H0(JL)
        CH2O(JL)=CH2O0(JL)
        O3S(JL)=O3S0(JL)
        ALD2(JL)=ALD20(JL)
        PAR(JL)=PAR0(JL)
        OLE(JL)=OLE0(JL)
        ETH(JL)=ETH0(JL)
        PAN(JL)=PAN0(JL)
        ACET(JL)=ACET0(JL)
        ISOP(JL)=ISOP0(JL)
        MGLY(JL)=MGLY0(JL)
        ISOPRD(JL)=ISOPRD0(JL)
        METHAC(JL)=METHAC0(JL)
        MVK(JL)=MVK0(JL)
        MEK(JL)=MEK0(JL)
        MPAN(JL)=MPAN0(JL)
        NITR(JL)=NITR0(JL)
        ZNO(JL)=ZNO0(JL)
        ZNO2(JL)=ZNO20(JL)
        ZNO3(JL)=ZNO30(JL)
        ZN2O5(JL)=ZN2O50(JL)
        HNO4(JL)=HNO40(JL)
        OH(JL)=OH0(JL)
        HO2(JL)=HO20(JL)
        CH3O2(JL)=CH3O20(JL)
        C2O3(JL)=C2O30(JL)
        XO2(JL)=XO20(JL)
        ROR(JL)=ROR0(JL)
        XO2N(JL)=XO2N0(JL)
        RXPAR(JL)=RXPAR0(JL)
        BXO2N(JL)=BXO2N0(JL)
        MC3O3(JL)=MC3O30(JL)

c
        ODDN0(JL)=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)+HNO4(JL)+
     *            PAN(JL)+MPAN(JL)+NITR(JL)

C LG-   sulfur chemistry

        DMS(JL)=DMS0(JL)
        SO2(JL)=SO20(JL)
        SAER(JL)=SAER0(JL)

C LG-   radon

	RADON(JL)=RADON0(JL)

C WP-   added isontr

        ISONTR(JL)=ISONTR0(JL)

C WP-   end

C WP-   added formic acid and acetic acid

        HCOOH(JL)=HCOOH0(JL)
	CH3CO2H(JL)=CH3CO2H0(JL)

C WP-   end

C LG-   ammonia chemistry
    
        NH2(JL)=NH20(JL)
        NH3(JL)=NH30(JL)
        NH4(JL)=NH40(JL)

C LG-   added the monoterpenes as a group

        MATERP(JL)=MATERP0(JL)
        MBTERP(JL)=MBTERP0(JL)

C LG-   added the sesquiterpenes as a group

        SQTERP(JL)=SQTERP0(JL)

C LG-   added HONO

        HONO(JL)=HONO0(JL)

C LG-   added hydroxy-methyl/alkylhydroxyperoxide

        CH2OHO2H(JL)=CH2OHO2H0(JL)
        RCHOHO2H(JL)=RCHOHO2H0(JL)

C LG-   end

C LG-   determining the ratio of NO and NO2 from JNO2, the O3 concentration
C       the reaction rate of the reaction between O3 and NO, this parameter
C       is being used within the subroutine bulkveg.f to estimate the subgrid
C       scale conversion of NO to NO2, this estimate can only be applied when
C       a photostationary equilibruim exists!
 
        IF (JK.GT.NLEV) THEN
         RATIONOX=RJNO2(JL)/(O3(JL)*RNOO3(JL))
        ENDIF

C LG-   end

 1105 CONTINUE

C
C  ** 1.2 ** CHEMISTRY STARTS HERE !!!!!!!
C

      DO 2104 JL=1,NLON

C LG- different iteration mechanism compared to the original ECHAM/CBM4
C     code, the number of iterations is determined by the relative
C     difference between the new and old concentration of NOy 

      ITER=0
      OLDNOY=ODDN0(JL)

 9999 CONTINUE

C --- First group: NO NO2 O3
          P1=RJBNO3(JL)*ZNO3(JL)

C LG- added the production of NO through the photolysis of HONO

     *     +RJHONO(JL)*HONO(JL)

          R12=0.2*RISOPNO2(JL)*ISOP(JL)
          R21=RHO2NO(JL)*HO2(JL)+RMO2NO(JL)*CH3O2(JL)+
     *      RC23NO(JL)*C2O3(JL)+RXO2NO(JL)*XO2(JL)+
     *      RMC23NO(JL)*MC3O3(JL)
          XL1=RNONO3(JL)*ZNO3(JL)+RXO2NNO(JL)*XO2N(JL)+
     *      RBXO2NNO(JL)*BXO2N(JL)

C LG- added the production of HONO from NO+OH

     *     +RNOOH(JL)*OH(JL)

          P2=RJHNO3(JL)*HNO3(JL)+RJN2O5(JL)*ZN2O5(JL)+
     *      RN2O5(JL)*ZN2O5(JL)+HNO4(JL)*
     *      (RJHNO4(JL)+RHNO4M(JL)+RHNO4OH(JL)*OH(JL))+
     *      2*RNONO3(JL)*ZNO3(JL)*ZNO(JL)+
     *      PAN(JL)*(RPAN(JL)+RJPAN(JL))+
     *      MPAN(JL)*(RMPAN(JL)+RJPAN(JL)+0.70*RMPANO3(JL)*O3(JL))+
     *      ZNO3(JL)*(RJANO3(JL)+ROLENO3(JL)*OLE(JL)+
     *       0.2*RISOPNO3(JL)*ISOP(JL))+
     *      RJNITR(JL)*NITR(JL)

C WP- added isontr:isontr+oh=no2

     *     +RISONTROH(JL)*ISONTR(JL)*OH(JL)

C WP- end

          XL2=RNO2OH(JL)*OH(JL)+RNO2NO3(JL)*ZNO3(JL)+
     *     RNO2HO2(JL)*HO2(JL)+RNO2O3(JL)*O3(JL)+
     *     RC23NO2(JL)*C2O3(JL)+RRORNO2(JL)*ROR(JL)+
     *     0.8*RISOPNO2(JL)*ISOP(JL)+RMC23NO2(JL)*MC3O3(JL)
          P3=RJANO3(JL)*ZNO3(JL)+ROHOH(JL)*OH(JL)*OH(JL)
	 
C LG- added acetic acid:o3 formation in reaction 
C     c2o3+ch3o2->0.3*(ch3co2h+o3) (the ratio is already accounted
C     for in the reaction rate)

     *    +RC23HO2(JL)*C2O3(JL)*HO2(JL)
	 
C LG- end

          XL3=RO3HO2(JL)*HO2(JL)+RO3OH(JL)*OH(JL)+
     *      RNO2O3(JL)*ZNO2(JL)+RJO3D(JL)+RETHO3(JL)*ETH(JL)+
     *      RISOPO3(JL)*ISOP(JL)+ROLEO3(JL)*OLE(JL)+
     *      RISPDO3(JL)*ISOPRD(JL)+RMTHCO3(JL)*METHAC(JL)+
     *      RMVKO3(JL)*MVK(JL)+RMPANO3(JL)*MPAN(JL)
          X1=ZNO0(JL)+P1*DT
          X2=ZNO20(JL)+P2*DT
          X3=O30(JL)+P3*DT
          C1=1.+XL1*DT
          C2=1.+XL2*DT
          C3=1.+XL3*DT
          Y2=RNOO3(JL)*DT/(C1*C3)
          XJT=RJNO2(JL)*DT
          R21T=R21*DT
          R12T=R12*DT
          R12TC=R12T/C2
          R21TC=R21T/C1
          XJTC=XJT/C2
C --- Oplossen onbekende x
          ACUB=-1.*Y2*(1.+R12TC+R21TC)
          BCUB=1.+R12TC+XJTC+R21TC+
     *           Y2*(R12TC*(X1-X2)+2.*R21TC*X1+X1+X3)
          CCUB=X2*(R12TC+XJTC)-X1*R21TC+Y2*X1*
     *           (X2*R12TC-X3-R21TC*X1)
          CUBDET=BCUB*BCUB-4.*ACUB*CCUB
          CUBDET=MAX(CUBDET,1.E-20)
          DNO2=(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB)
          ZNO2(JL)=(X2+DNO2)/C2
          ZNO(JL)=(X1-DNO2)/C1
          O3(JL)=(X3+XJTC*(X2+DNO2))/(C3+Y2*C3*(X1-DNO2))

C  --- Second group: HO2 OH HNO4
          R57=RJHNO4(JL)+RHNO4M(JL)
          R56=RCOOH(JL)*CO(JL)+RO3OH(JL)*O3(JL)+
     *       RHPOH(JL)*H2O2(JL)+RFRMOH(JL)*CH2O(JL)+
     *       RH2OH(JL)*AIR(JL)+RETHOH(JL)*ETH(JL)+
     *       ROLEOH(JL)*OLE(JL)+
     *       0.91*RISOPOH(JL)*ISOP(JL)+
     *       0.68*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL)+
     *       0.3*RMVKOH(JL)*MVK(JL)+0.11*RPAROH(JL)*PAR(JL)+
     *       0.4*RMPANOH(JL)*MPAN(JL)

C LG- 12092004+ the HO2 production term from the OH reactions

          PHO2OH(JL,JK)=R56

C LG- March 2002, adding the production of OH in the ozonolysis reactions of 
C     the terpenes, the yields are 0.76 and 0.33 for the ozonolysis of alpha 
C     and beta pinenes respectively (Kostas Tsigaridis, personal communication)
C     These reactions might be an important nocturnal source of OH

          R65=RHO2NO(JL)*ZNO(JL)+RO3HO2(JL)*O3(JL)

C LG- modified based on personal communications with Kostas, 23 July
C     2002, see also the change in the OH production term P6!

C LG- 12092004+ modified by adding an extra term that resembles the OH
C     production related to the terpene ozonolysis. This term was originally
C     included in the R65 term but the latter is also used to calculate the
C     HO2 production/destruction and the terpene ozonolysis reactions do not
C     result in an HO2 destruction, as it was initially suggested by including
C     the reactions in the R65 term!

          RTERPO3=
     *         0.76*RMATERPO3(JL)*O3(JL)
     *        +0.33*RMBTERPO3(JL)*O3(JL)

C LG- further modified based on personal communications with Boris Bonn, 
C     April, 2003. The yield of 0.22 is representative for the species alpha-
C     humulene, which reactions rates with OH, O3 and NO3 are also applied

     *        +0.22*RSQTERPO3(JL)*O3(JL)  ! check out this yield, a 50%
                       ! uncertainty would require a substantially lower
		       ! emission for a significant OH production 

C LG- end

C LG- assigning the OH production rate for interpretation of OH 
C     production rates  

C LG- 12092004+ the OH production term from the HO2 reactions

          POHHO2(JL,JK)=R65

          R75=RNO2HO2(JL)*ZNO2(JL)
          P5=2*RJBCH2O(JL)*CH2O(JL)+RMO2NO(JL)*CH3O2(JL)*ZNO(JL)+
     *      0.66*RMO2MO2(JL)*CH3O2(JL)*CH3O2(JL)+
     *      RJMEPE(JL)*CH3O2H(JL)+2*RJALD2(JL)*ALD2(JL)+
     *     C2O3(JL)*
     *       (RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL)+
     *       RC23MO2(JL)*CH3O2(JL))+
     *     O3(JL)*
     *       (0.44*ROLEO3(JL)*OLE(JL)+0.12*RETHO3(JL)*ETH(JL)+
     *       0.23*RISPDO3(JL)*ISOPRD(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+
     *       0.08*RMPANO3(JL)*MPAN(JL))+
     *     ZNO3(JL)*
     *       (0.8*RISOPNO3(JL)*ISOP(JL)+RISPDNO3(JL)*ISOPRD(JL)+
     *       0.5*RMTHCNO3(JL)*METHAC(JL)+RFRMNO3(JL)*CH2O(JL))+
     *     0.8*RISOPNO2(JL)*ISOP(JL)*ZNO2(JL)+
     *     RRORB(JL)*ROR(JL)+RJMEK(JL)*MEK(JL)+
     *     1.17*RRORA(JL)*ROR(JL)+
     *     RJNITR(JL)*NITR(JL)
     
C LG- added formic acid:hcooh+oh -> ho2

     *    +RHCOOHOH(JL)*HCOOH(JL)*OH(JL) 
     
C LG- end         

C LG- added acetic acid:production ho2 in reaction c2o3+ch3o2->ch2o+ho2

     *    +RC23MO2A(JL)*C2O3(JL)*CH3O2(JL) 
     
C LG- end       

          XL5=R65+R75+RMO2HO2(JL)*CH3O2(JL)+RHO2OH(JL)*OH(JL)+
     
C LG- added acetic acid: change of constant 0.21 in 1 (see paper Duncan),
C     There is production of 0.79 ho2 from 1 ho2, thus a net loss of
C     0.21 ho2, due to the removal of reaction 53 this yield is changed 
C     to 1

     *     1.*RC23HO2(JL)*C2O3(JL)+

C LG- end

     *      RXO2HO2(JL)*XO2(JL)+RXO2NHO2(JL)*XO2N(JL)+
     *      RBXO2NHO2(JL)*BXO2N(JL)+RMC23HO2(JL)*MC3O3(JL)
          R66=2.*RHO2HO2(JL)
          X5=HO20(JL)+P5*DT

C           P6=RJHNO3(JL)*HNO3(JL)+2*RJO3D(JL)*O3(JL)+
C      *     2.*RJH2O2(JL)*H2O2(JL)+RJMEPE(JL)*CH3O2H(JL)+
C      *     O3(JL)*
C      *       (0.4*RISPDO3(JL)*ISOPRD(JL)+0.1*ROLEO3(JL)*OLE(JL)+
C      *       0.2*RISOPO3(JL)*ISOP(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+
C      *       0.05*RMVKO3(JL)*MVK(JL)+0.04*RMPANO3(JL)*MPAN(JL))
C 
C C LG- removed acetic acid: removal of reaction 53 (see paper Duncan)
C 
C c     *     +0.79*RC23HO2(JL)*C2O3(JL)
C 
C C LG- end

C LG- modified based on personal communications with Kostas, July 2002

          P6=RJHNO3(JL)*HNO3(JL)+2.*RJO3D(JL)*O3(JL)+
     *     2.*RJH2O2(JL)*H2O2(JL)+RJMEPE(JL)*CH3O2H(JL)+
     *     O3(JL)*
     *       (0.4*RISPDO3(JL)*ISOPRD(JL)+0.1*ROLEO3(JL)*OLE(JL)+
     *       0.2*RISOPO3(JL)*ISOP(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+
     *       0.05*RMVKO3(JL)*MVK(JL)+0.04*RMPANO3(JL)*MPAN(JL))
     *      +0.76*RMATERPO3(JL)*MATERP(JL)*O3(JL)
     *      +0.33*RMBTERPO3(JL)*MBTERP(JL)*O3(JL)

C LG- futher modified based on personal communications with Boris Bonn, 
C     April 2003, to include the ozonolysis of the sesquiterpenes

     *      +0.22*RSQTERPO3(JL)*SQTERP(JL)*O3(JL)

C LG- and added the production of OH due to HONO photodissociation

     *      +RJHONO(JL)*HONO(JL)

C LG- end

          XL6=RHNO4OH(JL)*HNO4(JL)+RHO2OH(JL)*HO2(JL)+
     *       RNO2OH(JL)*ZNO2(JL)+ROHHNO3(JL)*HNO3(JL)+
     *       R56+RCH4OH(JL)*CH4(JL)+ROHPCAT(JL)*CH3O2H(JL)+
     *       2*ROHOH(JL)*OH(JL)+0.89*RPAROH(JL)*PAR(JL)+
     *       RMGLYOH(JL)*MGLY(JL)+0.09*RISOPOH(JL)*ISOP(JL)+
     *       0.32*RISPDOH(JL)*ISOPRD(JL)+
     *       0.5*RMTHCOH(JL)*METHAC(JL)+RMEKOH(JL)*MEK(JL)+
     *       0.7*RMVKOH(JL)*MVK(JL)+RALD2OH(JL)*ALD2(JL)+
     *       0.6*RMPANOH(JL)*MPAN(JL)+RACETOH(JL)*ACET(JL)

C LG- including oxidation of monoterpenes

     *      +RMATERPOH(JL)*MATERP(JL)
     *      +RMBTERPOH(JL)*MBTERP(JL)

C LG- futher modified based on personal communications with Boris Bonn, 
C     April 2003, to include the ozonolysis of the sesquiterpenes

     *      +RSQTERPOH(JL)*SQTERP(JL)

C LG- including HONO production through NO+OH

     *      +RNOOH(JL)*ZNO(JL)

C LG- end

          X6=OH0(JL)+P6*DT
          C6=1.+XL6*DT
          XL7=R57+RHNO4OH(JL)*OH(JL)
          C7=1.+XL7*DT
          Y1=R57/C7     ! R57=RJHNO4(JL)+RHNO4M(JL)
          Y2=R56/C6     ! R56=RCOOH(JL)*CO(JL)+RO3OH(JL)*O3(JL)+ etc.
          ACUB=R66*DT   ! R66=2.*RHO2HO2(JL)
          BCUB=1.+XL5*DT-DT2*(Y1*R75+Y2*R65) ! R75=RNO2HO2(JL)*ZNO2(JL)
	                                     ! R65=RHO2NO(JL)*ZNO(JL)+RO3HO2(JL)*O3(JL)
          CCUB=-1.*X5-DT*(Y1*HNO40(JL)+Y2*X6)
          CUBDET=BCUB*BCUB-4.*ACUB*CCUB
          CUBDET=MAX(CUBDET,1.E-20)
          HO2(JL)=(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB)
          OH(JL)=(X6+R65*HO2(JL)*DT)/C6
          HNO4(JL)=(HNO40(JL)+R75*DT*HO2(JL))/C7 ! R75=RNO2HO2(JL)*ZNO2(JL)

C LG- added the production of HONO including the reaction of NO and OH
C     and the reaction of NO2 and H2O. H2O is not included in the 
C     production term since this is already implictly considered (the
C     units of RNO2H2O is s-1 instead of cm3 molecules-1 s-1)

C LG- it still must be find out if the twice the NO2 concentration or
C     only once the NO2 concentration should be used (check the
C     parameterization). The HNO3 being formed is not gas-phase HNO3
C     but HNO3 present in the water/aerosol: 2NO2+H2O -> HONO+HNO3.
C     See also paper Harrisson et al, JGR, 1996. 

          PHONO=RNOOH(JL)*ZNO(JL)*OH(JL)+RNO2H2O(JL)*2.*ZNO2(JL)
	  XLHONO=RJHONO(JL)
          HONO(JL)=(HONO0(JL)+PHONO*DT)/(1.+XLHONO*DT)

C --- Third group: NO3 N2O5

          R89=RJN2O5(JL)+RN2O5(JL)
          P8=ROHHNO3(JL)*HNO3(JL)*OH(JL)+
     *        RNO2O3(JL)*ZNO2(JL)*O3(JL)+
     *        RMPANOH(JL)*MPAN(JL)*OH(JL)
          XL8=RJBNO3(JL)+RJANO3(JL)+RNONO3(JL)*ZNO(JL)+
     *        RNO2NO3(JL)*ZNO2(JL)+RALD2NO3(JL)*ALD2(JL)+
     *        RFRMNO3(JL)*CH2O(JL)+ROLENO3(JL)*OLE(JL)+
     *        RISOPNO3(JL)*ISOP(JL)+RISPDNO3(JL)*ISOPRD(JL)+
     *        RMTHCNO3(JL)*METHAC(JL)
          X4=ZNO30(JL)+P8*DT
          C5=1.+XL8*DT
          R98=RNO2NO3(JL)*ZNO2(JL)
          XL9=RJN2O5(JL)+RN2O5(JL)+RN2O5AQ(JL)+RN2O5L(JL)
          C6=1.+XL9*DT
          C7=(C5*C6-R89*R98*DT2)
          ZN2O5(JL)=(C5*ZN2O50(JL)+R98*DT*X4)/C7
          ZNO3(JL)=(C6*X4+R89*DT*ZN2O50(JL))/C7
          PCH3O2=RCH4OH(JL)*CH4(JL)*OH(JL)+
     *      OH(JL)*(ROHPCAT(JL)*CH3O2H(JL))+RJACET(JL)*ACET(JL)+
     *      0.5*RJMGLY(JL)*MGLY(JL)

C rc23mo2 not considered since P=L
          XLCH3O2=RMO2NO(JL)*ZNO(JL)+RMO2HO2(JL)*HO2(JL)+
     *      2*RMO2MO2(JL)*CH3O2(JL)
     
C LG- added acetic acid: destruction of ch3o2 in reaction with c2o3

     *     +RC23MO2A(JL)*C2O3(JL)+RC23MO2B(JL)*C2O3(JL) 
     
C LG- end        
     
          CH3O2(JL)=(CH3O20(JL)+PCH3O2*DT)/(1.+XLCH3O2*DT)

c --- Some others C2O3 XO2 etc
C C2O3

       PC2O3=ALD2(JL)*(RALD2OH(JL)*OH(JL)+RALD2NO3(JL)*ZNO3(JL))+
     *  OH(JL)*
     *   (RMGLYOH(JL)*MGLY(JL)+
     *    0.7*RMVKOH(JL)*MVK(JL)+0.5*RMEKOH(JL)*MEK(JL))+
     *  RJACET(JL)*ACET(JL)+PAN(JL)*(RPAN(JL)+RJPAN(JL))+
     *  RJMEK(JL)*MEK(JL)+0.5*RJMGLY(JL)*MGLY(JL)+
     *  O3(JL)*
     *    (0.17*RISPDO3(JL)*ISOPRD(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+
     *    0.70*RMPANO3(JL)*MPAN(JL))
       XLC2O3=RC23NO(JL)*ZNO(JL)+RC23NO2(JL)*ZNO2(JL)+
     * 2*RC23C23(JL)*C2O3(JL)+RC23HO2(JL)*HO2(JL)+
     
C LG- added acetic acid: using both rc23mo2a and b     
    
     *  RC23MO2A(JL)*CH3O2(JL)+RC23MO2B(JL)*CH3O2(JL)

C LG- end

       C2O3(JL)=(C2O30(JL)+PC2O3*DT)/(1.+XLC2O3*DT)
C XO2
       PXO2=RJALD2(JL)*ALD2(JL)+
     *  C2O3(JL)*
     *    (RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL))+

C LG- removed acetic acid: removal of reaction 53 (see paper Duncan)

c     *    0.79*RC23HO2(JL)*HO2(JL))+

C LG- end

     *  OH(JL)*
     *    (RACETOH(JL)*ACET(JL)+ROLEOH(JL)*OLE(JL)+
     *    RETHOH(JL)*ETH(JL)+0.87*RPAROH(JL)*PAR(JL)+
     *    RMGLYOH(JL)*MGLY(JL)+0.99*RISOPOH(JL)*ISOP(JL)+
     *    0.68*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL)+
     *    1.5*RMEKOH(JL)*MEK(JL)+RMVKOH(JL)*MVK(JL)+
     *    RMPANOH(JL)*MPAN(JL))+
     *  O3(JL)*
     *    (0.22*ROLEO3(JL)*OLE(JL)+0.2*RISPDO3(JL)*ISOPRD(JL)+
     *     0.2*RISOPO3(JL)*ISOP(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+
     *     0.05*RMVKO3(JL)*MVK(JL))+
     *  ZNO3(JL)*
     *    (0.91*ROLENO3(JL)*OLE(JL)+RISPDNO3(JL)*ISOPRD(JL)+
     *    RISOPNO3(JL)*ISOP(JL)+0.5*RMTHCNO3(JL)*METHAC(JL))+
     *  ISOP(JL)*RISOPNO2(JL)*ZNO2(JL)+0.76*RRORA(JL)*ROR(JL)+
     *  RJMEK(JL)*MEK(JL)
       XLXO2=RXO2NO(JL)*ZNO(JL)+2*RXO2XO2(JL)*XO2(JL)+
     *  RXO2HO2(JL)*HO2(JL)+RXO2NXO2(JL)*XO2N(JL)+
     *  RBXO2NXO2(JL)*BXO2N(JL)
       XO2(JL)=(XO20(JL)+PXO2*DT)/(1.+XLXO2*DT)
C shorties
       PROR=0.76*RPAROH(JL)*PAR(JL)*OH(JL)+0.03*RRORA(JL)*ROR(JL)
       XLROR=RRORA(JL)+RRORB(JL)+RRORNO2(JL)*ZNO2(JL)
       ROR(JL)=(ROR0(JL)+PROR*DT)/(1.+XLROR*DT)
       PRXPAR=3.39*RRORA(JL)*ROR(JL)+ROLEO3(JL)*OLE(JL)*O3(JL)+
     *  ROLENO3(JL)*OLE(JL)*ZNO3(JL)+
     *  OH(JL)*
     *    (0.1*RPAROH(JL)*PAR(JL)+ROLEOH(JL)*OLE(JL))
       XLRXPAR=RRXPAR(JL)*PAR(JL)
       RXPAR(JL)=(RXPAR0(JL)+PRXPAR*DT)/(1.+XLRXPAR*DT)
       PXO2N=0.13*RPAROH(JL)*PAR(JL)*OH(JL)+0.03*RRORA(JL)*ROR(JL)+
     *   0.09*ROLENO3(JL)*OLE(JL)*ZNO3(JL)
       XLXO2N=RXO2NNO(JL)*ZNO(JL)+2*RXO2N(JL)*XO2N(JL)+
     *   RXO2NHO2(JL)*HO2(JL)+RXO2NXO2(JL)*XO2(JL)+
     *   RBXO2NXO2N(JL)*BXO2N(JL)
       XO2N(JL)=(XO2N0(JL)+PXO2N*DT)/(1.+XLXO2N*DT)
       PBXO2N=0.09*RISOPOH(JL)*ISOP(JL)*OH(JL)
       XLBXO2N=RBXO2NNO(JL)*ZNO(JL)+2*RBXO2N(JL)*BXO2N(JL)+
     *  RBXO2NHO2(JL)*HO2(JL)+RBXO2NXO2(JL)*XO2(JL)+
     *  RBXO2NXO2N(JL)*XO2N(JL)
       BXO2N(JL)=(BXO2N0(JL)+PBXO2N*DT)/(1.+XLBXO2N*DT)
       POLE=0.
       XLOLE=ROLEOH(JL)*OH(JL)+ROLEO3(JL)*O3(JL)+ROLENO3(JL)*ZNO3(JL)
       OLE(JL)=(OLE0(JL)+POLE*DT)/(1.+XLOLE*DT)
       PMGLY=
     *   OH(JL)*
     *     (0.15*RISPDOH(JL)*ISOPRD(JL)+0.08*RMTHCOH(JL)*METHAC(JL)+
     *     0.3*RMVKOH(JL)*MVK(JL)+RACETOH(JL)*ACET(JL))+
     *   O3(JL)*
     *     (0.65*RISPDO3(JL)*ISOPRD(JL)+0.9*RMTHCO3(JL)*METHAC(JL)+
     *      0.95*RMVKO3(JL)*MVK(JL))
       XLMGLY=RMGLYOH(JL)*OH(JL)+RJMGLY(JL)
       MGLY(JL)=(MGLY0(JL)+PMGLY*DT)/(1.+XLMGLY*DT)
       PISOP=0.
       XLISOP=RISOPOH(JL)*OH(JL)+RISOPO3(JL)*O3(JL)+
     *   RISOPNO3(JL)*ZNO3(JL)+RISOPNO2(JL)*ZNO2(JL)
       ISOP(JL)=(ISOP0(JL)+PISOP*DT)/(1.+XLISOP*DT)
       PISOPRD=ISOP(JL)*
     *    (0.36*RISOPOH(JL)*OH(JL)+0.1*RISOPO3(JL)*O3(JL)+
     *     0.2*RISOPNO3(JL)*ZNO3(JL)+0.2*RISOPNO2(JL)*ZNO2(JL))
       XLISOPRD=RISPDOH(JL)*OH(JL)+RISPDO3(JL)*O3(JL)+
     *  RISPDNO3(JL)*ZNO3(JL)
       ISOPRD(JL)=(ISOPRD0(JL)+PISOPRD*DT)/(1.+XLISOPRD*DT)
       PMETHAC=ISOP(JL)*
     *    (0.23*RISOPOH(JL)*OH(JL)+0.39*RISOPO3(JL)*O3(JL))
       XLMETHAC=RMTHCOH(JL)*OH(JL)+RMTHCO3(JL)*O3(JL)+
     *   RMTHCNO3(JL)*ZNO3(JL)
       METHAC(JL)=(METHAC0(JL)+PMETHAC*DT)/(1.+XLMETHAC*DT)
       PMVK=ISOP(JL)*(0.32*RISOPOH(JL)*OH(JL)+0.16*RISOPO3(JL)*O3(JL))
       XLMVK=RMVKOH(JL)*OH(JL)+RMVKO3(JL)*O3(JL)
       MVK(JL)=(MVK0(JL)+PMVK*DT)/(1.+XLMVK*DT)
       PMEK=0.42*RMTHCOH(JL)*OH(JL)*METHAC(JL)+
     *  ISOPRD(JL)*(0.48*RISPDOH(JL)*OH(JL)+0.28*RISPDO3(JL)*O3(JL))
       XLMEK=RMEKOH(JL)*OH(JL)+RJMEK(JL)
       MEK(JL)=(MEK0(JL)+PMEK*DT)/(1.+XLMEK*DT)
       PMC3O3=ISOP(JL)*0.2*RISOPO3(JL)*O3(JL)+
     *   OH(JL)*
     *     (0.31*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL))+
     *   0.5*RMTHCNO3(JL)*METHAC(JL)*ZNO3(JL)+
     *   MPAN(JL)*(RMPAN(JL)+RJPAN(JL))
       XLMC3O3=RMC23NO(JL)*ZNO(JL)+RMC23NO2(JL)*ZNO2(JL)+
     *   RMC23HO2(JL)*HO2(JL)
       MC3O3(JL)=(MC3O30(JL)+PMC3O3*DT)/(1.+XLMC3O3*DT)

       PCH3O2H=RMO2HO2(JL)*CH3O2(JL)*HO2(JL)
       XLCH3O2H=(ROHPCAT(JL)+ROHPFRM(JL))*OH(JL)+RJMEPE(JL)
       CH3O2H(JL)=(CH3O2H0(JL)+PCH3O2H*DT)/(1.+XLCH3O2H*DT)

C LG-  082003 added the production and destruction terms of hydroxy-
C      methyl/alkylhydroxyperoxide

       PCH2OHO2H=0.105*RMBTERPO3(JL)*MBTERP(JL)*O3(JL)
       XLCH2OHO2H=RCH2OHO2H(JL)  ! LG- HMHP goes into HCOOH and H2O: 
                      ! The decomposition rate is taken from Peter Neebs thesis
       CH2OHO2H(JL)=(CH2OHO2H0(JL)+PCH2OHO2H*DT)/(1.+XLCH2OHO2H*DT)

C LG-  Table 2.7, Valverde's Thesis, alpha- and beta pinene destruction to mimic 
C      differences between endo and exocyclic alkene ozonolyis 

       PRCHOHO2H=0.05*RMATERPO3(JL)*MATERP(JL)*O3(JL)+
     *           0.23*RMBTERPO3(JL)*MBTERP(JL)*O3(JL)
       XLRCHOHO2H=RRCHOHO2H(JL) ! LG- all the HAHP goes to H2O2 (see below)
          ! LG- correct???? adding production of RCHO
       RCHOHO2H(JL)=(RCHOHO2H0(JL)+PRCHOHO2H*DT)/(1.+XLRCHOHO2H*DT)

C LG-  end

       PHNO3=RNO2OH(JL)*ZNO2(JL)*OH(JL)+
     *      2.*(RN2O5AQ(JL)+RN2O5L(JL))*ZN2O5(JL)+
     *      ZNO3(JL)*
     *        (RALD2NO3(JL)*ALD2(JL)+RFRMNO3(JL)*CH2O(JL)+
     *        RMTHCNO3(JL)*METHAC(JL)+0.8*RISOPNO3(JL)*ISOP(JL)+
     *        RISPDNO3(JL)*ISOPRD(JL))+
     *   0.8*RISOPNO2(JL)*ZNO2(JL)*ISOP(JL)

C WP- removed isontr:

c     *   +RBXO2NNO(JL)*BXO2N(JL)*ZNO(JL)

C WP- end
       
       XLHNO3=RJHNO3(JL)+OH(JL)*ROHHNO3(JL)
       HNO3(JL)=(HNO30(JL)+PHNO3*DT)/(1.+XLHNO3*DT)
       PH2O2=RHO2HO2(JL)*HO2(JL)*HO2(JL)

C LG-  added the production of H2O2 through the ozonolysis of the terpenes/
C      alkenes. The key production term is the further destruction of
C      RCHOHO2H

     *   +RRCHOHO2H(JL)*RCHOHO2H(JL) ! LG- assuming a 100% conversion of HAHP to H2O2

       XLH2O2=RJH2O2(JL)+RHPOH(JL)*OH(JL)
       H2O2(JL)=(H2O20(JL)+PH2O2*DT)/(1.+XLH2O2*DT)
       PCH2O=RMO2NO(JL)*CH3O2(JL)*ZNO(JL)+RJMEPE(JL)*CH3O2H(JL)+
     *   ROHPFRM(JL)*CH3O2H(JL)*OH(JL)+RJALD2(JL)*ALD2(JL)+
     *   1.33*RMO2MO2(JL)*CH3O2(JL)*CH3O2(JL)+
     *   C2O3(JL)*(RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL)+

C LG- removed acetic acid: removal of reaction 53 (see paper Duncan)

c     *     0.79*RC23HO2(JL)*HO2(JL)+

C LG- end

C LG- added acetic acid: updating rc23mo2a and b

     *   RC23MO2A(JL)*CH3O2(JL)+RC23MO2B(JL)*CH3O2(JL))+

C LG- end

     *   OH(JL)*(1.56*RETHOH(JL)*ETH(JL)+ROLEOH(JL)*OLE(JL)+
     *     0.63*RISOPOH(JL)*ISOP(JL)+0.14*RISPDOH(JL)*ISOPRD(JL)+
     *     0.08*RMTHCOH(JL)*METHAC(JL)+0.3*RMVKOH(JL)*MVK(JL)+
     *     0.4*RMPANOH(JL)*MPAN(JL)+0.5*RMEKOH(JL)*MEK(JL))+
     *   O3(JL)*(RETHO3(JL)*ETH(JL)+0.74*ROLEO3(JL)*OLE(JL)+
     *     0.6*RISOPO3(JL)*ISOP(JL)+0.39*RISPDO3(JL)*ISOPRD(JL)+
     *     0.7*RMPANO3(JL)*MPAN(JL)+0.2*RMTHCO3(JL)*METHAC(JL)+
     *     0.1*RMVKO3(JL)*MVK(JL))+
     *   ZNO3(JL)*(ROLENO3(JL)*OLE(JL)+0.33*RISPDNO3(JL)*ISOPRD(JL))+
     *   MC3O3(JL)*(RMC23NO(JL)*ZNO(JL)+2*RMC23HO2(JL)*HO2(JL))
       XLCH2O=RJACH2O(JL)+RJBCH2O(JL)+OH(JL)*RFRMOH(JL)+
     *    ZNO3(JL)*RFRMNO3(JL)
       CH2O(JL)=(CH2O0(JL)+PCH2O*DT)/(1.+XLCH2O*DT)
       PPAN=ZNO2(JL)*RC23NO2(JL)*C2O3(JL)+
     *    MPAN(JL)*(0.3*RMPANO3(JL)*O3(JL)+0.4*RMPANOH(JL)*OH(JL))
       XLPAN=RPAN(JL)+RJPAN(JL)
       PAN(JL)=(PAN0(JL)+PPAN*DT)/(1.+XLPAN*DT)
       PMPAN=RMC23NO2(JL)*MC3O3(JL)*ZNO2(JL)
       XLMPAN=RMPAN(JL)+RJPAN(JL)+RMPANO3(JL)*O3(JL)+
     *    RMPANOH(JL)*OH(JL)
       MPAN(JL)=(MPAN0(JL)+PMPAN*DT)/(1.+XLMPAN*DT)
       PETH=0.
       XLETH=RETHOH(JL)*OH(JL)+RETHO3(JL)*O3(JL)
       ETH(JL)=(ETH0(JL)+PETH*DT)/(1.+XLETH*DT)
       PACET=0.64*RRORA(JL)*ROR(JL)+0.6*RMPANOH(JL)*MPAN(JL)*OH(JL)+
     *     0.8*RJNITR(JL)*NITR(JL)
       XLACET=RACETOH(JL)*OH(JL)+RJACET(JL)
       ACET(JL)=(ACET0(JL)+PACET*DT)/(1.+XLACET*DT)
       PNTR=RRORNO2(JL)*ROR(JL)*ZNO2(JL)+
     *   RXO2NNO(JL)*XO2N(JL)*ZNO(JL)
       XLNTR=RJNITR(JL)
       NITR(JL)=(NITR0(JL)+PNTR*DT)/(1.+XLNTR*DT)

C WP- added isontr:

       PISONTR=RBXO2NNO(JL)*BXO2N(JL)*ZNO(JL)
       XLISONTR=RJNITR(JL)+RISONTROH(JL)*OH(JL)
       ISONTR(JL)=(ISONTR0(JL)+PISONTR*DT)/(1.+XLISONTR*DT)

C WP- end

C LG- March 2002, added the destruction of the monoterpenes

       PMATERP=0.
       XLMATERP=RMATERPOH(JL)*OH(JL)+RMATERPO3(JL)*O3(JL)+
     *   RMATERPNO3(JL)*ZNO3(JL)
       MATERP(JL)=(MATERP0(JL)+PMATERP*DT)/(1.+XLMATERP*DT)
       PMBTERP=0.
       XLMBTERP=RMBTERPOH(JL)*OH(JL)+RMBTERPO3(JL)*O3(JL)+
     *   RMBTERPNO3(JL)*ZNO3(JL)
       MBTERP(JL)=(MBTERP0(JL)+PMBTERP*DT)/(1.+XLMBTERP*DT)

C LG- April 2002, added the destruction of the sesquiterpenes

       PSQTERP=0.
       XLSQTERP=RSQTERPOH(JL)*OH(JL)+RSQTERPO3(JL)*O3(JL)+
     *   RSQTERPNO3(JL)*ZNO3(JL)
       SQTERP(JL)=(SQTERP0(JL)+PSQTERP*DT)/(1.+XLSQTERP*DT)

C LG- added formic/acetic acid:the formation and destruction of 
C     formic acid and acetic acid

       PHCOOH=0.56*0.5*RISOPO3(JL)*ISOP(JL)*O3(JL)
     *   +RCH2OHO2H(JL)*CH2OHO2H(JL)  ! LG- added the production through the
                                      ! decomposition of HMHP
       XLHCOOH=RHCOOHOH(JL)*OH(JL)
       HCOOH(JL)=(HCOOH0(JL)+PHCOOH*DT)/(1.+XLHCOOH*DT)
       PCH3CO2H=RC23HO2(JL)*C2O3(JL)*HO2(JL)+
     *           RC23MO2B(JL)*C2O3(JL)*CH3O2(JL)
       XLCH3CO2H=RCH3CO2HOH(JL)*OH(JL)
       CH3CO2H(JL)=(CH3CO2H0(JL)+PCH3CO2H*DT)/(1.+XLCH3CO2H*DT)

C LG- end

       CH4(JL)=CH40(JL)/(1.+RCH4OH(JL)*OH(JL)*DT)
       PCO=CH2O(JL)*XLCH2O+RJALD2(JL)*ALD2(JL)+
     *     O3(JL)*
     *       (0.44*RETHO3(JL)*ETH(JL)+0.33*ROLEO3(JL)*OLE(JL)+
     *       0.13*RMPANO3(JL)*MPAN(JL))+
     *     MGLY(JL)*(OH(JL)*RMGLYOH(JL)+1.5*RJMGLY(JL))
       CO(JL)=(CO0(JL)+PCO*DT)/(1.+RCOOH(JL)*OH(JL)*DT)
       XLALD2=RALD2OH(JL)*OH(JL)+RALD2NO3(JL)*ZNO3(JL)+RJALD2(JL)
       PALD2=
     *     OH(JL)*(0.22*RETHOH(JL)*ETH(JL)+0.11*RPAROH(JL)*PAR(JL)+
     *       ROLEOH(JL)*OLE(JL)+0.7*RMVKOH(JL)*MVK(JL)+
     *       RMEKOH(JL)*MEK(JL)+0.19*RISPDOH(JL)*ISOPRD(JL))+
     *     O3(JL)*(0.1*RMTHCO3(JL)*METHAC(JL)+0.5*ROLEO3(JL)*OLE(JL)+
     *       0.15*RISOPO3(JL)*ISOP(JL)+0.06*RISPDO3(JL)*ISOPRD(JL))+
     *     ZNO3(JL)*(ROLENO3(JL)*OLE(JL)+0.8*RISOPNO3(JL)*ISOP(JL)+
     *       0.33*RISPDNO3(JL)*ISOPRD(JL))+RJMEK(JL)*MEK(JL)+
     *     1.1*RRORA(JL)*ROR(JL)+0.8*RISOPNO2(JL)*ZNO2(JL)*ISOP(JL)+
     *     0.2*RJNITR(JL)*NITR(JL)
       ALD2(JL)=(ALD20(JL)+PALD2*DT)/(1.+XLALD2*DT)
       PPAR=0.
       XLPAR=RPAROH(JL)*OH(JL)+RRXPAR(JL)*RXPAR(JL)
       PAR(JL)=(PAR0(JL)+PPAR*DT)/(1.+XLPAR*DT)
C stratosferic ozone
       XLO3S=RO3HO2(JL)*HO2(JL)+RO3OH(JL)*OH(JL)+
     *       RJO3D(JL)
       PO3S=0.
       IF (JK.LT.KTRHE(JL)) THEN
          XLO3S=RO3HO2(JL)*HO2(JL)+RO3OH(JL)*OH(JL)+
     *         RNO2O3(JL)*ZNO2(JL)+RJO3D(JL)+RNOO3(JL)*ZNO(JL)
          PO3S=RJANO3(JL)*ZNO3(JL)+RJNO2(JL)*ZNO2(JL)
       ENDIF
       O3S(JL)=(O3S0(JL)+PO3S*DT)/(1.+XLO3S*DT)

C LG-  added sulfur chemistry:

       IF (LSULFCHEM) THEN
         XLSO2=ROHSO2(JL)*OH(JL)
         SO2(JL)=(SO20(JL)+(ROHDMS(JL)*OH(JL)+
     *          RNO3DMS(JL)*ZNO3(JL))*DMS(JL)*DT)/(1.+XLSO2*DT)
         PSAER=ROHSO2(JL)*OH(JL)*SO2(JL)
         SAER(JL)=SAER0(JL)+PSAER*DT
         XLDMS=ROHDMS(JL)*OH(JL)+RNO3DMS(JL)*ZNO3(JL)
         DMS(JL)=DMS0(JL)/(1.+XLDMS*DT)
       ENDIF

C LG-  and ammonia chemistry

       IF (LNH3CHEM) THEN
         ACID(JL)=(ACID0(JL)+(2.*XLSO2*SO2(JL))*DT)/
     $          (1.+RSO4NH3(JL)*NH3(JL)*DT)
         NH4(JL)=NH40(JL)+ACID(JL)*RSO4NH3(JL)*NH3(JL)*DT
         PNH2=OH(JL)*ROHNH3(JL)

C LG-    note that the destruction of NH3 by the reaction with OH is not
C        included here. In the paper by Dentener and Crutzen [1994], this
C        reaction is mentioned being an important source of N2O but also a
C        very uncertain source due to the uncertainty in the rate constants

         XLNH3=ACID(JL)*RSO4NH3(JL)+PNH2
         NH3(JL)=NH30(JL)/(1.+XLNH3*DT)
         XLNH2=RNONH2(JL)*ZNO(JL)+RNO2NH2(JL)*ZNO2(JL)+
     $      RHO2NH2(JL)*HO2(JL)+RO2NH2(JL)+RO3NH2(JL)*O3(JL)
         NH2(JL)=(NH20(JL)+NH3(JL)*PNH2*DT)/(1.+XLNH2*DT)

       ENDIF

C LG-  end

       NEWNOY=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)+HNO4(JL)+
     *        PAN(JL)+MPAN(JL)+NITR(JL)

C LG-  radon decay, Radioactive decay of 222^Ra; life time of 5.5 day

       RSP_LTIME=1./(5.5*86400.)
       RADON(JL)=RADON0(JL)*EXP(-PTMST*RSP_LTIME)

C LG-  end

C LG-  iteration counter

       ITER=ITER+1

       IF (ABS((NEWNOY-OLDNOY)/NEWNOY).LT.1.E-15.OR.ITER.GT.250) 
     *   GOTO 999
       OLDNOY=NEWNOY

       GOTO 9999

  999 CONTINUE  

 2104 CONTINUE 

C ** 1.3 ** --- Finish all
      DO 106 JL=1,NLON

        ZPM(JL,JK,inox)=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)+HNO4(JL)
        DODDN=ZPM(JL,JK,inox)+PAN(JL)+MPAN(JL)+NITR(JL)-ODDN0(JL)

        TOTN(JL)=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)
     &     +HNO4(JL)+PAN(JL)+MPAN(JL)+HNO3(JL)+NITR(JL)

        ZPM(JL,JK,io3)=O3(JL)
        ZPM(JL,JK,ich4)=CH4(JL)
        ZPM(JL,JK,ico)=CO(JL)
        ZPM(JL,JK,ihno3)=ZPM(JL,JK,ihno3)-DODDN

        ZPM(JL,JK,ih2o2)=H2O2(JL)

C         IF (JK.LT.NLEV) THEN
C          ZPM(JL,JK,ih2o2)=H2O2(JL)
C 	ELSE
C 	 print *,'cbm4_ech: H2O2 concentration NOT updated for level: ',jk
C         ENDIF

        ZPM(JL,JK,ich3o2h)=CH3O2H(JL)
        ZPM(JL,JK,ich2o)=CH2O(JL)
        ZPM(JL,JK,io3s)=O3S(JL)
        ZPM(JL,JK,io3s)=MAX(ZPM(JL,JK,io3s),0.)
        ZPM(JL,JK,iald2)=ALD2(JL)
        ZPM(JL,JK,ipar)=PAR(JL)
        ZPM(JL,JK,iole)=OLE(JL)
        ZPM(JL,JK,ieth)=ETH(JL)
        ZPM(JL,JK,ipan)=PAN(JL)
        ZPM(JL,JK,iacet)=ACET(JL)
        ZPM(JL,JK,iisop)=ISOP(JL)
        ZPM(JL,JK,imgly)=MGLY(JL)
        ZPM(JL,JK,iisoprd)=ISOPRD(JL)
        ZPM(JL,JK,imethac)=METHAC(JL)
        ZPM(JL,JK,imvk)=MVK(JL)
        ZPM(JL,JK,imek)=MEK(JL)
        ZPM(JL,JK,impan)=MPAN(JL)
        ZPM(JL,JK,intr)=NITR(JL)

C LG- sulfur chemistry

        ZPM(JL,JK,idms)=DMS(JL)
        ZPM(JL,JK,iso2)=SO2(JL)
        ZPM(JL,JK,iso4)=SAER(JL)

C LG- radon

        ZPM(JL,JK,irad)=RADON(JL)

C WP- added isontr

c        ZPM(JL,JK,iisontr)=ISONTR(JL)

C WP- end

C LG- added formic and acetic acid

        ZPM(JL,JK,ihcooh)=HCOOH(JL)
	ZPM(JL,JK,ich3co2h)=CH3CO2H(JL)

C LG- ammonia chemistry

        ZPM(JL,JK,inh2)=NH2(JL)
        ZPM(JL,JK,inh3)=NH3(JL)
        ZPM(JL,JK,inh4)=NH4(JL)

C LG- and the monoterpenes as a group

        ZPM(JL,JK,imaterp)=MATERP(JL)
        ZPM(JL,JK,imbterp)=MBTERP(JL)

C LG- and the sesquiterpenes as a group

        ZPM(JL,JK,isqterp)=SQTERP(JL)

C LG- and HONO

        ZPM(JL,JK,ihono)=HONO(JL)

C LG- and hydroxy-methyl/alkylhydroxyperoxide

        ZPM(JL,JK,ich2oho2h)=CH2OHO2H(JL)
        ZPM(JL,JK,irchoho2h)=RCHOHO2H(JL)

        IF (NTRAC.EQ.NTRACT) THEN

C LG-    For NTRAC = NTRACT, the defined spin-up time for chemistry 
C        is (not yet) being used. 

         IF (JK.EQ.NLEV.AND.NSTEP.EQ.1) THEN
          WRITE(*,'(1a)')
     &  ' The defined spin-up time (T_SPINUP) in namchem is NOT used',
     &  ' within CBM4_ECH.f !!'
          WRITE(*,'(1a)')' ENTER to continue'
          READ (*,*)
         ENDIF

         ZPM(JL,JK,ino)=ZNO(JL)
         ZPM(JL,JK,ino2)=ZNO2(JL)
         ZPM(JL,JK,ino3)=ZNO3(JL)
         ZPM(JL,JK,in2o5)=ZN2O5(JL)
         ZPM(JL,JK,ihno4)=HNO4(JL)
         ZPM(JL,JK,ioh)=OH(JL)
         ZPM(JL,JK,iho2)=HO2(JL)
         ZPM(JL,JK,ich3o2)=CH3O2(JL)
         ZPM(JL,JK,ic2o3)=C2O3(JL)
         ZPM(JL,JK,ixo2)=XO2(JL)
         ZPM(JL,JK,iror)=ROR(JL)
         ZPM(JL,JK,ixo2n)=XO2N(JL)
         ZPM(JL,JK,irxpar)=RXPAR(JL)
         ZPM(JL,JK,ibxo2n)=BXO2N(JL)          
         ZPM(JL,JK,imc3o3)=MC3O3(JL)
	ELSE

C LG-    For NTRAC << that of NTRACT, the defined spin-up time for chemistry 
C        is being used. 

         IF (JK.EQ.NLEV.AND.NSTEP.EQ.1) THEN
          WRITE(*,'(1a)')
     &  ' A defined spin-up time (T_SPINUP) in namchem is also used',
     &  ' within CBM4_ECH.f !!'
          WRITE(*,'(1a)')' ENTER to continue'
          READ (*,*)
	 ENDIF

         IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
          ZPM(JL,JK,ino)=ZNO(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZNO0(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPM(JL,JK,ino2)=ZNO2(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZNO20(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPM(JL,JK,ino3)=ZNO3(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZNO30(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPM(JL,JK,in2o5)=ZN2O5(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZN2O5(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPM(JL,JK,ihno4)=HNO4(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       HNO40(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPM(JL,JK,inox)=ZPM(JL,JK,ino)+ZPM(JL,JK,ino2)+
     &      ZPM(JL,JK,ino3)+2*ZPM(JL,JK,in2o5)+ZPM(JL,JK,ihno4)
         ELSE
          ZPMLOC(JL,JK,ino)=ZNO(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZNO0(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPMLOC(JL,JK,ino2)=ZNO2(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZNO20(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPMLOC(JL,JK,ino3)=ZNO3(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZNO30(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPMLOC(JL,JK,in2o5)=ZN2O5(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       ZN2O5(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPMLOC(JL,JK,ihno4)=HNO4(JL)*MIN(1.,NSTEP/(T_SPINUP/PTMST))+
     &       HNO40(JL)*(1.-MIN(1.,NSTEP/(T_SPINUP/PTMST)))
          ZPM(JL,JK,inox)=ZPMLOC(JL,JK,ino)+ZPMLOC(JL,JK,ino2)+
     &      ZPMLOC(JL,JK,ino3)+2*ZPMLOC(JL,JK,in2o5)+ZPMLOC(JL,JK,ihno4)
         ENDIF

         ZPMLOC(JL,JK,ioh)=OH(JL)
         ZPMLOC(JL,JK,iho2)=HO2(JL)
         ZPMLOC(JL,JK,ich3o2)=CH3O2(JL)
         ZPMLOC(JL,JK,ic2o3)=C2O3(JL)
         ZPMLOC(JL,JK,ixo2)=XO2(JL)
         ZPMLOC(JL,JK,iror)=ROR(JL)
         ZPMLOC(JL,JK,ixo2n)=XO2N(JL)
         ZPMLOC(JL,JK,irxpar)=RXPAR(JL)
         ZPMLOC(JL,JK,ibxo2n)=BXO2N(JL)
         ZPMLOC(JL,JK,imc3o3)=MC3O3(JL)
	ENDIF

C LG- end

  106 CONTINUE

C LG- not considering the concentrations at zz for budgets and
C     and other I/O

      IF (LXTMZZ.AND.JK.EQ.NLEVEL+1) GOTO 101

C LG- writing the budgets of the chemical reactions

      DO 107 JL=1,NLON
       IW=IWHERE(JL,JK)
       ZFAC=PTMST*GRVOL(JL,JK)
        IW=IWHERE(JL,JK)
        ZFAC=PTMST*GRVOL(JL,JK)
        BRX(IW,1)=BRX(IW,1)+RJNO2(JL)*ZNO2(JL)*ZFAC
        BRX(IW,2)=BRX(IW,2)+RNOO3(JL)*ZNO(JL)*O3(JL)*ZFAC
        BRX(IW,3)=BRX(IW,3)+RHO2NO(JL)*HO2(JL)*ZNO(JL)*ZFAC
        BRX(IW,4)=BRX(IW,4)+RMO2NO(JL)*CH3O2(JL)*ZNO(JL)*ZFAC
        BRX(IW,5)=BRX(IW,5)+RNO2OH(JL)*ZNO2(JL)*OH(JL)*ZFAC
        BRX(IW,6)=BRX(IW,6)+RJHNO3(JL)*HNO3(JL)*ZFAC
        BRX(IW,7)=BRX(IW,7)+ROHHNO3(JL)*OH(JL)*HNO3(JL)*ZFAC
        BRX(IW,8)=BRX(IW,8)+RJN2O5(JL)*ZN2O5(JL)*ZFAC
        BRX(IW,9)=BRX(IW,9)+RJBNO3(JL)*ZNO3(JL)*ZFAC
        BRX(IW,10)=BRX(IW,10)+RJANO3(JL)*ZNO3(JL)*ZFAC
        BRX(IW,11)=BRX(IW,11)+RNO2O3(JL)*ZNO2(JL)*O3(JL)*ZFAC
        BRX(IW,12)=BRX(IW,12)+RNONO3(JL)*ZNO(JL)*ZNO3(JL)*ZFAC
        BRX(IW,13)=BRX(IW,13)+RNO2NO3(JL)*ZNO2(JL)*ZNO3(JL)*ZFAC
        BRX(IW,14)=BRX(IW,14)+RN2O5(JL)*ZN2O5(JL)*ZFAC
        BRX(IW,15)=BRX(IW,15)+(RN2O5AQ(JL)+RN2O5L(JL))*ZN2O5(JL)*ZFAC
        BRX(IW,16)=BRX(IW,16)+ROHOH(JL)*OH(JL)*OH(JL)*ZFAC
        BRX(IW,17)=BRX(IW,17)+RFRMNO3(JL)*CH2O(JL)*ZNO3(JL)*ZFAC
        BRX(IW,20)=BRX(IW,20)+RJO3D(JL)*O3(JL)*ZFAC
        BRX(IW,25)=BRX(IW,25)+RJBCH2O(JL)*CH2O(JL)*ZFAC
        BRX(IW,27)=BRX(IW,27)+RO3HO2(JL)*O3(JL)*HO2(JL)*ZFAC
        BRX(IW,28)=BRX(IW,28)+RCOOH(JL)*CO(JL)*OH(JL)*ZFAC
        BRX(IW,29)=BRX(IW,29)+RO3OH(JL)*O3(JL)*OH(JL)*ZFAC
        BRX(IW,31)=BRX(IW,31)+RFRMOH(JL)*CH2O(JL)*OH(JL)*ZFAC
        BRX(IW,32)=BRX(IW,32)+RCH4OH(JL)*CH4(JL)*OH(JL)*ZFAC
        BRX(IW,38)=BRX(IW,38)+RJACH2O(JL)*CH2O(JL)*ZFAC
        BRX(IW,41)=BRX(IW,41)+RC23NO(JL)*ZNO(JL)*C2O3(JL)*ZFAC
        BRX(IW,42)=BRX(IW,42)+RMC23NO(JL)*MC3O3(JL)*ZNO(JL)*ZFAC
        BRX(IW,43)=BRX(IW,43)+RXO2NO(JL)*XO2(JL)*ZNO(JL)*ZFAC
        BRX(IW,44)=BRX(IW,44)+RBXO2NNO(JL)*BXO2N(JL)*ZNO(JL)*ZFAC
        BRX(IW,45)=BRX(IW,45)+RALD2NO3(JL)*ALD2(JL)*ZNO3(JL)*ZFAC
        BRX(IW,46)=BRX(IW,46)+ROLENO3(JL)*OLE(JL)*ZNO3(JL)*ZFAC
        BRX(IW,47)=BRX(IW,47)+RISOPNO3(JL)*ISOP(JL)*ZNO3(JL)*ZFAC
        BRX(IW,48)=BRX(IW,48)+RISPDNO3(JL)*ISOPRD(JL)*ZNO3(JL)*ZFAC
        BRX(IW,49)=BRX(IW,49)+RMTHCNO3(JL)*METHAC(JL)*ZNO3(JL)*ZFAC
        BRX(IW,50)=BRX(IW,50)+RETHO3(JL)*ETH(JL)*O3(JL)*ZFAC
        BRX(IW,51)=BRX(IW,51)+ROLEO3(JL)*OLE(JL)*O3(JL)*ZFAC
        BRX(IW,52)=BRX(IW,52)+RISOPO3(JL)*ISOP(JL)*O3(JL)*ZFAC
        BRX(IW,53)=BRX(IW,53)+RISPDO3(JL)*ISOPRD(JL)*O3(JL)*ZFAC
        BRX(IW,54)=BRX(IW,54)+RMTHCO3(JL)*METHAC(JL)*O3(JL)*ZFAC
        BRX(IW,55)=BRX(IW,55)+RMVKO3(JL)*MVK(JL)*O3(JL)*ZFAC
        BRX(IW,56)=BRX(IW,56)+RISOPNO2(JL)*ISOP(JL)*ZNO2(JL)*ZFAC
        BRX(IW,57)=BRX(IW,57)+RRORNO2(JL)*ROR(JL)*ZNO2(JL)*ZFAC
        BRX(IW,58)=BRX(IW,58)+ZFAC*(
     *  ISOP(JL)*(0.37*RISOPOH(JL)*OH(JL)+0.6*RISOPO3(JL)*O3(JL)+
     *    2.4*RISOPNO2(JL)*ZNO2(JL)+2.4*RISOPNO3(JL)*ZNO3(JL))+
     *  ISOPRD(JL)*(0.87*RISPDOH(JL)*OH(JL)+
     *    1.08*RISPDO3(JL)*O3(JL)+4.01*RISPDNO3(JL)*ZNO3(JL))+
     *  METHAC(JL)*(0.7*RMTHCO3(JL)*O3(JL)+2*RMTHCNO3(JL)*ZNO3(JL))+
     *  MC3O3(JL)*(3*RMC23NO(JL)*ZNO(JL)+2*RMC23HO2(JL)*HO2(JL))+
     *  MPAN(JL)*(RMPANOH(JL)*OH(JL)+1.17*RMPANO3(JL)*O3(JL))+
     *  1.05*RMVKO3(JL)*MVK(JL)*O3(JL)+0.5*RMEKOH(JL)*MEK(JL)*OH(JL)+
     *  C2O3(JL)*(RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL)+
     *    1.21*RC23HO2(JL)*HO2(JL)+RC23MO2(JL)*CH3O2(JL))+
     *  ROR(JL)*(1.158*RRORB(JL)+1.158*RRORNO2(JL)*ZNO2(JL)-
     *    .22516*RRORA(JL))
     *  -2.8*RJNITR(JL)*NITR(JL)+0.58*RETHO3(JL)*O3(JL)
     *  +0.93*ROLEO3(JL)*O3(JL)
     *  )

C LG-   sulfur chemistry

        IF (LSULFCHEM) THEN

         BRX(IW,40)=BRX(IW,40)+DMS(JL)*OH(JL)*ROHDMS(JL)*ZFAC
         BRX(IW,41)=BRX(IW,41)+DMS(JL)*ZNO3(JL)*RNO3DMS(JL)*ZFAC
         BRX(IW,42)=BRX(IW,42)+SO2(JL)*OH(JL)*ROHSO2(JL)*ZFAC

        ENDIF

C LG-   end

        IF (LWRITEREAC) THEN

C WP-   calculation of the Ozone Producing Potential for all levels.
C       The ozone producing potential is a linear estimate of the ozone
C       production/loss due to chemistry. It represents all the reactions
C       that disturb the fast photochemical equilibrium between NO2/NO and O3.
C       It is a multiplication of species Y, the rate constant K and the [O3].
C       This gives an estimate of P(O3) in molec cm-3 s-1.
C       The fast chemical equilibrium with NO2 is not included, and neither is
C       the formation of reservoir species, since these reactions are reversible 
C       in time. 

        IF (NSTEP.EQ.0.AND.JK.EQ.KTHEIGHT) THEN

C LG-     opening files for chemical production/destruction of O3 and
C         some other trace gases

          OPEN(UNIT=NUNREAC,FILE='/data/ganzevl/racmo/output/reaction.out',STATUS='UNKNOWN')
          OPEN(UNIT=NUNREAC2,FILE='/data/ganzevl/racmo/output/reactionOH.out',STATUS='UNKNOWN')
          OPEN(UNIT=NUNREAC3,FILE='/data/ganzevl/racmo/output/reactionHO2.out',STATUS='UNKNOWN')

          WRITE(NUNREAC,'(2a)')
     &      'Chemical production and destruction terms for a selection ',
     &      'of trace gases'
 	  WRITE(NUNREAC,*)NSTOP,PTMST,NPRINT,NLEVEL

          WRITE(NUNREAC2,'(2a)')
     &      'Chemical production and destruction terms for ',
     &      'OH'
 	  WRITE(NUNREAC2,*)NSTOP,PTMST,NPRINT,NLEVEL

          WRITE(NUNREAC3,'(2a)')
     &      'Chemical production and destruction terms for ',
     &      'HO2'
 	  WRITE(NUNREAC3,*)NSTOP,PTMST,NPRINT,NLEVEL

          IF (LBIOSPH.AND.LAGRIAN) THEN
	   WRITE (NUNREAC,*)  (PZ(JL,JJK),JJK=NLEVEL,1,-1)
	   WRITE (NUNREAC2,*)  (PZ(JL,JJK),JJK=NLEVEL,1,-1)
	   WRITE (NUNREAC3,*)  (PZ(JL,JJK),JJK=NLEVEL,1,-1)
	  ELSE 
	   WRITE (NUNREAC,*)  (HGHT(JJK),JJK=NLEVEL,1,-1)
	   WRITE (NUNREAC2,*)  (HGHT(JJK),JJK=NLEVEL,1,-1)
	   WRITE (NUNREAC3,*)  (HGHT(JJK),JJK=NLEVEL,1,-1)
          ENDIF

C LG-     writing of the different reaction/photolysis rates numbers/names as 
C         they are defined in the parameter list of the hydrocarbon chemistry code
C         so that the different reactions can be coupled to the proper reaction
C         names, this list must resemble that of the listed reactions in the 
C         write statement of the actual prodcution/destruction terms,
C         the photolysis rates numbers all start with an j** whereas the reaction
C         rate numbers all start with with nr**** (since r*** is already used
C         for the reaction rates as a function of the longitude)

          NTRAC_OUT=7
	  NREACTION(io3)=17  ! mz_lg_20060316+ modified
	  NREACTION(ich2o)=32
	  NREACTION(iisop)=4
          NREACTION(ino_pr)=12
	  NREACTION(ino2_pr)=32
	  NREACTION(ih2o2)=4
	  NREACTION(ich3o2h)=4

	  NREACTION(ioh_pr)=49

	  NREACTION(iho2_pr)=48

	  WRITE(NUNREAC,'(3i4)') ICHEMTYPE,NTRACT,NTRAC_OUT
	  WRITE(NUNREAC,'(50i4)')io3,ich2o,iisop,ino_pr,ino2_pr,
     &          ih2o2,ich3o2h 
          WRITE(NUNREAC,'(50i4)')NREACTION(io3),NREACTION(ich2o),
     &          NREACTION(iisop),NREACTION(ino_pr),NREACTION(ino2_pr),
     &          NREACTION(ih2o2),NREACTION(ich3o2h)

C LG-     OH, to other outputfile!

	  WRITE(NUNREAC2,'(3i4)') ICHEMTYPE,NTRACT,1
	  WRITE(NUNREAC2,'(50i4)')ioh_pr  
          WRITE(NUNREAC2,'(50i4)')NREACTION(ioh_pr)

C LG-     HO2, to other outputfile!

	  WRITE(NUNREAC3,'(3i4)') ICHEMTYPE,NTRACT,1
	  WRITE(NUNREAC3,'(50i4)')iho2_pr  

          WRITE(NUNREAC3,'(50i4)')NREACTION(iho2_pr)

C LG-     writing the reaction/photodissociation numbers which are defined
C         in the parameter file parreact.h and are used to define the legends in
C         IDL scripts by coupling these numbers (positive for photodiss., and
C         negative for reaction rates) to three characterstrings in the IDL
C         scripts.

C LG-     ozone

	  WRITE(NUNREAC,'(100i4)')NRHO2NO,NRMO2NO,
     *     NRC23NO,NRXO2NO,NROHOH,NRC23HO2,NRO3OH,NRO3HO2,NROLEO3,
     *     NRISOPO3,NRISPDO3,NRMVKO3,NRMPANO3,NRH2OOD,NRMATERPO3,
     *     NRMBTERPO3,NRSQTERPO3

C LG-     formaldehyde

	  WRITE(NUNREAC,'(100i4)')NRMO2NO,JMEPE,
     *     NROHPFRM,JALD2,NRMO2MO2,NRC23NO,NRC23C23,
     *     NRC23HO2,NRC23MO2,NRETHOH,NROLEOH,NRISOPOH,
     *     NRISPDOH,NRMTHCOH,NRMVKOH,NRMPANOH,NRMEKOH,
     *     NRETHO3,NROLEO3,NRISOPO3,NRISPDO3,NRMPANO3,
     *     NRMTHCO3,NRMVKO3,NROLENO3,NRISPDNO3,NRMC23NO,
     *     NRMC23HO2,JACH2O,JBCH2O,NRFRMOH,NRFRMNO3

C LG-     isoprene

	  WRITE(NUNREAC,'(100i4)')NRISOPOH,NRISOPO3,NRISOPNO3,
     *          NRISOPNO2

C LG-     NO

	  WRITE(NUNREAC,'(100i4)')JNO2,JBNO3,NRISOPNO2,NRHO2NO,
     *          NRMO2NO,NRC23NO,NRXO2NO,NRMC23NO,NRNONO3,NRXO2NNO,
     *          NRBXO2NNO,NRNOO3

C LG-     NO2

	  WRITE(NUNREAC,'(100i4)')NRHO2NO,NRMO2NO,NRC23NO,NRXO2NO,
     *          NRMC23NO,NRNOO3,JHNO3,JN2O5,NRN2O5,JHNO4,
     *          NRHNO4M,NRHNO4OH,NRNONO3,NRPAN,JPAN,NRMPAN,JMPAN,
     *          NRMPANO3,JANO3,NROLENO3,NRISOPNO3,JNTR,NRISONTROH,
     *          JNO2,NRNO2OH,NRNO2NO3,NRNO2HO2,NRNO2O3,NRC23NO2,
     *          NRRORNO2,NRISOPNO2,NRMC23NO2

C LG-     H2O2

	  WRITE(NUNREAC,'(100i4)')NRHO2HO2,NRRCHOHO2H,NRHPOH,JH2O2

C LG-     CH3O2H

	  WRITE(NUNREAC,'(100i4)')NRMO2HO2,NROHPCAT,NROHPFRM,JMEPE

C LG-     OH to other outputfile!

	  WRITE(NUNREAC2,'(100i4)')JHNO3,JO3D,JH2O2,JMEPE,
     *      NRISPDO3,NROLEO3,NRISOPO3,NRMTHCO3,NRMVKO3,NRMPANO3,
     *      NRMATERPO3,NRMBTERPO3,NRSQTERPO3,JHONO,NR65,NRHNO4OH,NRHO2OH,
     *      NRNO2OH,NROHHNO3,NRCOOH,NRO3OH,NRHPOH,NRFRMOH,NRH2OH,NRETHOH,
     *      NROLEOH,NRISOPOH,NRISPDOH,NRMTHCOH,NRMVKOH,NRPAROH,NRMPANOH,
     *      NRCH4OH,NROHPCAT,NROHOH,NRPAROH,NRMGLYOH,NRISOPOH,NRISPDOH,
     *      NRMTHCOH,NRMEKOH,NRMVKOH,NRALD2OH,NRMPANOH,NRACETOH,
     *      NRMATERPOH,NRMBTERPOH,NRSQTERPOH,NRNOOH

C LG-     HO2 to other outputfile!

C LG-     15092004, I5 format since we have more then 100 reactions!

	  WRITE(NUNREAC3,'(100i5)')JBCH2O,JMEPE,JALD2,JNTR,
     *      NRMO2NO,NRMO2MO2,NRC23NO,NRC23C23,NRC23MO2,NROLEO3,NRETHO3,
     *      NRISPDO3,NRMTHCO3,NRMPANO3,NRISOPNO3,NRISPDNO3,NRMTHCNO3,
     *      NRFRMNO3,NRISOPNO2,NRRORB,NRRORA,NRHCOOHOH,NRC23MO2A,
     *      NRHNO4M,NRCOOH,NRO3OH,NRHPOH,NRFRMOH,NRH2OH,NRETHOH,NROLEOH,
     *      NRISOPOH,NRISPDOH,NRMTHCOH,NRMVKOH,NRPAROH,NRMPANOH,
     *      NRHO2NO,NRO3HO2,NRNO2HO2,NRMO2HO2,NRHO2OH,NRC23HO2,NRXO2HO2,
     *      NRXO2NHO2,NRBXO2NHO2,NRMC23HO2,NRHO2HO2     

        ENDIF 
    	
C LG-   writing some integration specific data with frequency nprint

        IF(MOD(NSTEP,NPRINT).EQ.0) THEN

         IF (JK.EQ.KTHEIGHT) THEN
           WRITE(NUNREAC,'(a14)') LDATLTIME
           WRITE(NUNREAC,*) NSTEP,JDAY,GMT,LTIME,HC,TRPHGT,PBLHGT,
     *                      KTHEIGHT

           WRITE(NUNREAC2,'(a14)') LDATLTIME
           WRITE(NUNREAC2,*) NSTEP,JDAY,GMT,LTIME,HC,TRPHGT,PBLHGT,
     *                       KTHEIGHT

           WRITE(NUNREAC3,'(a14)') LDATLTIME
           WRITE(NUNREAC3,*) NSTEP,JDAY,GMT,LTIME,HC,TRPHGT,PBLHGT,
     *                       KTHEIGHT

       	 ENDIF

         IF (NTRAC.EQ.NTRACT) THEN

C WP-    O1d is needed for an estimate of the loss reaction h2o+o1d
	
         O1D=ZPM(JL,JK,io3)*ZPPHOTCHEM(JL,JK,7)/(RODM(JL)*AIR(JL))	
       
C WP-    Ozone Producing Potential calculation, this is actual the sum
C        of all the reactions (being calculated within the IDL script
C        being the total chemical tendency!), The reaction between NO and
C        O3 and the photodissociation of NO2 are not considered assuming
C        a photostationary equilibrium 

         OPP(JK)=

C LG-           production terms 

     * 		ZPM(JL,JK,ino)*RHO2NO(JL)*ZPM(JL,JK,iho2)	!,'NO+HO2'
     * 		+ZPM(JL,JK,ino)*RMO2NO(JL)*ZPM(JL,JK,ich3o2)	!,'NO+CH3O2'
     *		+ZPM(JL,JK,ino)*RC23NO(JL)*ZPM(JL,JK,ic2o3)	!,'NO+C2O3'
     * 		+ZPM(JL,JK,ino)*RXO2NO(JL)*ZPM(JL,JK,ixo2)	!,'NO+XO2'
     * 		+ZPM(JL,JK,ioh)*ROHOH(JL)*ZPM(JL,JK,ioh)	!,'OH+OH'
     * 		+ZPM(JL,JK,ic2o3)*RC23HO2(JL)*ZPM(JL,JK,iho2)	!,'C23+HO2'

C LG-           and destruction terms

     * 		-ZPM(JL,JK,io3)*RO3OH(JL)*ZPM(JL,JK,ioh)	!,'O3+OH'
     * 		-ZPM(JL,JK,io3)*RO3HO2(JL)*ZPM(JL,JK,iho2)	!,'O3+HO2'
     * 		-ZPM(JL,JK,io3)*ROLEO3(JL)*ZPM(JL,JK,iole)	!,'O3+OLE'
     * 		-ZPM(JL,JK,io3)*RISOPO3(JL)*ZPM(JL,JK,iisop)	!,'O3+ISOP'
     * 		-ZPM(JL,JK,io3)*RISPDO3(JL)*ZPM(JL,JK,iisoprd)	!,'O3+ISOPRD'
     * 		-ZPM(JL,JK,io3)*RMVKO3(JL)*ZPM(JL,JK,imvk)	!,'O3+MVK'
     * 		-ZPM(JL,JK,io3)*RMPANO3(JL)*ZPM(JL,JK,impan)	!,'O3+MPAN'
     * 		-RH2OOD(JL)*H2O(JL)*O1D				!,'H2O+O1D'
     *          -ZPM(JL,JK,io3)*RMATERPO3(JL)*ZPM(JL,JK,imaterp)!,'O3+aterp'
     *          -ZPM(JL,JK,io3)*RMBTERPO3(JL)*ZPM(JL,JK,imbterp)!,'O3+bterp'
     *          -ZPM(JL,JK,io3)*RSQTERPO3(JL)*ZPM(JL,JK,isqterp)!,'O3+sqterp'

C LG- end

         WRITE(NUNREAC,*)          
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     * 		ZPM(JL,JK,ino)*RHO2NO(JL)*ZPM(JL,JK,iho2),	!,'NO+HO2'
     * 		ZPM(JL,JK,ino)*RMO2NO(JL)*ZPM(JL,JK,ich3o2),	!,'NO+CH3O2'
     *		ZPM(JL,JK,ino)*RC23NO(JL)*ZPM(JL,JK,ic2o3),	!,'NO+C2O3'
     * 		ZPM(JL,JK,ino)*RXO2NO(JL)*ZPM(JL,JK,ixo2),	!,'NO+XO2'
     * 		ZPM(JL,JK,ioh)*ROHOH(JL)*ZPM(JL,JK,ioh),	!,'OH+OH'
     * 		ZPM(JL,JK,ic2o3)*RC23HO2(JL)*ZPM(JL,JK,iho2),	!,'C23+HO2'
     * 		-ZPM(JL,JK,io3)*RO3OH(JL)*ZPM(JL,JK,ioh),	!,'O3+OH'
     * 		-ZPM(JL,JK,io3)*RO3HO2(JL)*ZPM(JL,JK,iho2),	!,'O3+HO2'
     * 		-ZPM(JL,JK,io3)*ROLEO3(JL)*ZPM(JL,JK,iole),	!,'O3+OLE'
     * 		-ZPM(JL,JK,io3)*RISOPO3(JL)*ZPM(JL,JK,iisop),	!,'O3+ISOP'
     * 		-ZPM(JL,JK,io3)*RISPDO3(JL)*ZPM(JL,JK,iisoprd),	!,'O3+ISOPRD'
     * 		-ZPM(JL,JK,io3)*RMVKO3(JL)*ZPM(JL,JK,imvk),	!,'O3+MVK'
     * 		-ZPM(JL,JK,io3)*RMPANO3(JL)*ZPM(JL,JK,impan),	!,'O3+MPAN'
     * 		-RH2OOD(JL)*H2O(JL)*O1D,			!,'H2O+O1D'
     *          -ZPM(JL,JK,io3)*RMATERPO3(JL)*ZPM(JL,JK,imaterp), !,'O3+aterp'
     *          -ZPM(JL,JK,io3)*RMBTERPO3(JL)*ZPM(JL,JK,imbterp), !,'O3+bterp'
     *          -ZPM(JL,JK,io3)*RSQTERPO3(JL)*ZPM(JL,JK,isqterp)  !,'O3+sqterp'

C LG-    production terms of formaldehyde

         WRITE(NUNREAC,*)
     *      RMO2NO(JL)*ZPM(JL,JK,ich3o2)*ZPM(JL,JK,ino),
     *      RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *      ROHPFRM(JL)*ZPM(JL,JK,ich3o2h)*ZPM(JL,JK,ioh),
     *      RJALD2(JL)*ZPM(JL,JK,iald2),
     *      1.33*RMO2MO2(JL)*ZPM(JL,JK,ich3o2)*ZPM(JL,JK,ich3o2),
     *      RC23NO(JL)*ZPM(JL,JK,ic2o3)*ZPM(JL,JK,ino),
     *      2*RC23C23(JL)*ZPM(JL,JK,ic2o3)*ZPM(JL,JK,ic2o3),
     *      0.79*RC23HO2(JL)*ZPM(JL,JK,iho2),
     *      RC23MO2(JL)*ZPM(JL,JK,ich3o2),
     *      1.56*RETHOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,ieth),
     *      ROLEOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,iole),
     *      0.63*RISOPOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,iisop),
     *      0.14*RISPDOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,iisoprd),
     *      0.08*RMTHCOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,imethac),
     *      0.3*RMVKOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,imvk),
     *      0.4*RMPANOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,impan),
     *      0.5*RMEKOH(JL)*ZPM(JL,JK,ioh)*ZPM(JL,JK,imek),
     *      RETHO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,ieth),
     *      0.74*ROLEO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iole),
     *      0.6*RISOPO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iisop),
     *      0.39*RISPDO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iisoprd),
     *      0.7*RMPANO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,impan),
     *      0.2*RMTHCO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,imethac),
     *      0.1*RMVKO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,imvk),
     *      ROLENO3(JL)*ZPM(JL,JK,ino3)*ZPM(JL,JK,iole),
     *      0.33*RISPDNO3(JL)*ZPM(JL,JK,ino3)*ZPM(JL,JK,iisoprd),
     *      RMC23NO(JL)*ZPM(JL,JK,imc3o3)*ZPM(JL,JK,ino),
     *      2*RMC23HO2(JL)*ZPM(JL,JK,imc3o3)*ZPM(JL,JK,iho2),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ich2o)*RJACH2O(JL),
     *          -ZPM(JL,JK,ich2o)*RJBCH2O(JL),
     *          -ZPM(JL,JK,ich2o)*RFRMOH(JL)*ZPM(JL,JK,ioh),        
     *          -ZPM(JL,JK,ich2o)*RFRMNO3(JL)*ZPM(JL,JK,ino3)

C LG-    production terms of Isoprene

         WRITE(NUNREAC,*),

C LG-    and destruction terms

     *          -ZPM(JL,JK,iisop)*RISOPOH(JL)*ZPM(JL,JK,ioh),
     *          -ZPM(JL,JK,iisop)*RISOPO3(JL)*ZPM(JL,JK,io3),
     *          -ZPM(JL,JK,iisop)*RISOPNO3(JL)*ZPM(JL,JK,ino3),
     *          -ZPM(JL,JK,iisop)*RISOPNO2(JL)*ZPM(JL,JK,ino2)

C LG-    production terms of NO

         WRITE(NUNREAC,*),
     *          RJNO2(JL)*ZPM(JL,JK,ino2),
     *          RJBNO3(JL)*ZPM(JL,JK,ino3),
     *          ZPM(JL,JK,ino2)*0.2*RISOPNO2(JL)*ZPM(JL,JK,iisop),

C LG-    and destruction terms

     *          -ZPM(JL,JK,iho2)*RHO2NO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ich3o2)*RMO2NO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ic2o3)*RC23NO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ixo2)*RXO2NO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,imc3o3)*RMC23NO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ino)*RNONO3(JL)*ZPM(JL,JK,ino3),
     *          -ZPM(JL,JK,ixo2n)*RXO2NNO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ibxo2n)*RBXO2NNO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ino)*RNOO3(JL)*ZPM(JL,JK,io3)

C LG-    production terms of NO2

         WRITE(NUNREAC,*),
     *          ZPM(JL,JK,iho2)*RHO2NO(JL)*ZPM(JL,JK,ino),
     *          ZPM(JL,JK,ich3o2)*RMO2NO(JL)*ZPM(JL,JK,ino),
     *          ZPM(JL,JK,ic2o3)*RC23NO(JL)*ZPM(JL,JK,ino),
     *          ZPM(JL,JK,ixo2)*RXO2NO(JL)*ZPM(JL,JK,ino),
     *          ZPM(JL,JK,imc3o3)*RMC23NO(JL)*ZPM(JL,JK,ino),
     *          ZPM(JL,JK,ino)*RNOO3(JL)*ZPM(JL,JK,io3),
     *          RJHNO3(JL)*ZPM(JL,JK,ihno3),
     *          RJN2O5(JL)*ZPM(JL,JK,in2o5),
     *          RN2O5(JL)*ZPM(JL,JK,in2o5),
     *          RJHNO4(JL)*ZPM(JL,JK,ihno4),
     *          RHNO4M(JL)*ZPM(JL,JK,ihno4),
     *          ZPM(JL,JK,ihno4)*RHNO4OH(JL)*ZPM(JL,JK,ioh),
     *          ZPM(JL,JK,ino)*2.*RNONO3(JL)*ZPM(JL,JK,ino3),
     *          RPAN(JL)*ZPM(JL,JK,ipan),
     *          RJPAN(JL)*ZPM(JL,JK,ipan),
     *          RMPAN(JL)*ZPM(JL,JK,impan),
     *          RJPAN(JL)*ZPM(JL,JK,impan),
     *          ZPM(JL,JK,impan)*0.70*RMPANO3(JL)*ZPM(JL,JK,io3),
     *          RJANO3(JL)*ZPM(JL,JK,ino3),
     *          ZPM(JL,JK,iole)*ROLENO3(JL)*ZPM(JL,JK,ino3),
     *          ZPM(JL,JK,iisop)*0.2*RISOPNO3(JL)*ZPM(JL,JK,ino3),
     *          RJNITR(JL)*ZPM(JL,JK,intr),
     *          ZPM(JL,JK,iisontr)*RISONTROH(JL)*ZPM(JL,JK,ioh),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ino2)*RJNO2(JL),
     *          -ZPM(JL,JK,ino2)*RNO2OH(JL)*ZPM(JL,JK,ioh),
     *          -ZPM(JL,JK,ino2)*RNO2NO3(JL)*ZPM(JL,JK,ino3),
     *          -ZPM(JL,JK,ino2)*RNO2HO2(JL)*ZPM(JL,JK,iho2),
     *          -ZPM(JL,JK,ino2)*RNO2O3(JL)*ZPM(JL,JK,io3),
     *          -ZPM(JL,JK,ic2o3)*RC23NO2(JL)*ZPM(JL,JK,ino2),
     *          -ZPM(JL,JK,iror)*RRORNO2(JL)*ZPM(JL,JK,ino2),
     *          -ZPM(JL,JK,iisop)*0.8*RISOPNO2(JL)*ZPM(JL,JK,ino2),
     *          -ZPM(JL,JK,imc3o3)*RMC23NO2(JL)*ZPM(JL,JK,ino2)

C LG-    production terms of H2O2

         WRITE(NUNREAC,*),RHO2HO2(JL)*ZPM(JL,JK,iho2)*ZPM(JL,JK,iho2),
     *                    RRCHOHO2H(JL)*ZPM(JL,JK,irchoho2h),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ih2o2)*RHPOH(JL)*ZPM(JL,JK,ioh),
     *          -ZPM(JL,JK,ih2o2)*RJH2O2(JL)

C LG-    production terms of CH3O2H

         WRITE(NUNREAC,*),RMO2HO2(JL)*ZPM(JL,JK,ich3o2)*ZPM(JL,JK,iho2),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ich3o2h)*ROHPCAT(JL)*ZPM(JL,JK,ioh),
     *          -ZPM(JL,JK,ich3o2h)*ROHPFRM(JL)*ZPM(JL,JK,ioh),
     *          -ZPM(JL,JK,ich3o2h)*RJMEPE(JL)

C LG-    production terms of OH

         WRITE(NUNREAC2,*),
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     *          RJHNO3(JL)*ZPM(JL,JK,ihno3),
     *          2.*RJO3D(JL)*ZPM(JL,JK,io3),
     *          2.*RJH2O2(JL)*ZPM(JL,JK,ih2o2),
     *          RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *          ZPM(JL,JK,iisoprd)*0.4*RISPDO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,iole)*0.1*ROLEO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,iisop)*0.2*RISOPO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imethac)*0.1*RMTHCO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imvk)*0.05*RMVKO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,impan)*0.04*RMPANO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imaterp)*0.76*RMATERPO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imbterp)*0.33*RMBTERPO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,isqterp)*0.22*RSQTERPO3(JL)*ZPM(JL,JK,io3),
     *          RJHONO(JL)*ZPM(JL,JK,ihono),
     *          POHHO2(JL,JK)*ZPM(JL,JK,iho2),

C LG-    and destruction terms
     *         -ZPM(JL,JK,ioh)*RHNO4OH(JL)*ZPM(JL,JK,ihno4),
     *         -ZPM(JL,JK,ioh)*RHO2OH(JL)*ZPM(JL,JK,iho2),
     *         -ZPM(JL,JK,ioh)*RNO2OH(JL)*ZPM(JL,JK,ino2),
     *         -ZPM(JL,JK,ioh)*ROHHNO3(JL)*ZPM(JL,JK,ihno3),
     *         -ZPM(JL,JK,ioh)*RCOOH(JL)*ZPM(JL,JK,ico),
     *         -ZPM(JL,JK,ioh)*RO3OH(JL)*ZPM(JL,JK,io3),
     *         -ZPM(JL,JK,ioh)*RHPOH(JL)*ZPM(JL,JK,ih2o2),
     *         -ZPM(JL,JK,ioh)*RFRMOH(JL)*ZPM(JL,JK,ich2o),
     *         -ZPM(JL,JK,ioh)*RH2OH(JL)*ZPRHOA(JL,JK)*6.022045E23/ZMAIR,
     *         -ZPM(JL,JK,ioh)*RETHOH(JL)*ZPM(JL,JK,ieth),
     *         -ZPM(JL,JK,ioh)*ROLEOH(JL)*ZPM(JL,JK,iole),
     *         -ZPM(JL,JK,ioh)*0.91*RISOPOH(JL)*ZPM(JL,JK,iisop),
     *         -ZPM(JL,JK,ioh)*0.68*RISPDOH(JL)*ZPM(JL,JK,iisoprd),
     *         -ZPM(JL,JK,ioh)*0.5*RMTHCOH(JL)*ZPM(JL,JK,imethac),
     *         -ZPM(JL,JK,ioh)*0.3*RMVKOH(JL)*ZPM(JL,JK,imvk),
     *         -ZPM(JL,JK,ioh)*0.11*RPAROH(JL)*ZPM(JL,JK,ipar),
     *         -ZPM(JL,JK,ioh)*0.4*RMPANOH(JL)*ZPM(JL,JK,impan),
     *         -ZPM(JL,JK,ioh)*RCH4OH(JL)*ZPM(JL,JK,ich4),
     *         -ZPM(JL,JK,ioh)*ROHPCAT(JL)*ZPM(JL,JK,ich3o2h),
     *         -ZPM(JL,JK,ioh)*ROHOH(JL)*ZPM(JL,JK,ioh),
     *         -ZPM(JL,JK,ioh)*0.89*RPAROH(JL)*ZPM(JL,JK,ipar),
     *         -ZPM(JL,JK,ioh)*RMGLYOH(JL)*ZPM(JL,JK,imgly),
     *         -ZPM(JL,JK,ioh)*0.09*RISOPOH(JL)*ZPM(JL,JK,iisop),
     *         -ZPM(JL,JK,ioh)*0.32*RISPDOH(JL)*ZPM(JL,JK,iisoprd),
     *         -ZPM(JL,JK,ioh)*0.5*RMTHCOH(JL)*ZPM(JL,JK,imethac),
     *         -ZPM(JL,JK,ioh)*RMEKOH(JL)*ZPM(JL,JK,imek),
     *         -ZPM(JL,JK,ioh)*0.7*RMVKOH(JL)*ZPM(JL,JK,imvk),
     *         -ZPM(JL,JK,ioh)*RALD2OH(JL)*ZPM(JL,JK,iald2),
     *         -ZPM(JL,JK,ioh)*0.6*RMPANOH(JL)*ZPM(JL,JK,impan),
     *         -ZPM(JL,JK,ioh)*RACETOH(JL)*ZPM(JL,JK,iacet),
     *         -ZPM(JL,JK,ioh)*RMATERPOH(JL)*ZPM(JL,JK,imaterp),
     *         -ZPM(JL,JK,ioh)*RMBTERPOH(JL)*ZPM(JL,JK,imbterp),
     *         -ZPM(JL,JK,ioh)*RSQTERPOH(JL)*ZPM(JL,jk,isqterp),
     *         -ZPM(JL,JK,ioh)*RNOOH(JL)*ZPM(JL,JK,ino)

C LG-    production terms of HO2

         WRITE(NUNREAC3,*),
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     *          2.*RJBCH2O(JL)*ZPM(JL,JK,ich2o),
     *          RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *          2.*RJALD2(JL)*ZPM(JL,JK,iald2),
     *          RJNITR(JL)*ZPM(JL,JK,intr),
     *          RMO2NO(JL)*ZPM(JL,JK,ich3o2)*ZPM(JL,JK,ino),
     *          0.66*RMO2MO2(JL)*ZPM(JL,JK,ich3o2)*ZPM(JL,JK,ich3o2),
     *          ZPM(JL,JK,ic2o3)*RC23NO(JL)*ZPM(JL,JK,ino),     
     *          ZPM(JL,JK,ic2o3)*2.*RC23C23(JL)*ZPM(JL,JK,ic2o3),
     *          ZPM(JL,JK,ic2o3)*RC23MO2(JL)*ZPM(JL,JK,ich3o2),
     *          ZPM(JL,JK,io3)*0.44*ROLEO3(JL)*ZPM(JL,JK,iole),
     *          ZPM(JL,JK,io3)*0.12*RETHO3(JL)*ZPM(JL,JK,ieth),
     *          ZPM(JL,JK,io3)*0.23*RISPDO3(JL)*ZPM(JL,JK,iisoprd),
     *          ZPM(JL,JK,io3)*0.1*RMTHCO3(JL)*ZPM(JL,JK,imethac),
     *          ZPM(JL,JK,io3)*0.08*RMPANO3(JL)*ZPM(JL,JK,impan),
     *          ZPM(JL,JK,ino3)*0.8*RISOPNO3(JL)*ZPM(JL,JK,iisop),
     *          ZPM(JL,JK,ino3)*RISPDNO3(JL)*ZPM(JL,JK,iisoprd),
     *          ZPM(JL,JK,ino3)*0.5*RMTHCNO3(JL)*ZPM(JL,JK,imethac),
     *          ZPM(JL,JK,ino3)*RFRMNO3(JL)*ZPM(JL,JK,ich2o),
     *          RISOPNO2(JL)*ZPM(JL,JK,iisop)*ZPM(JL,JK,ino2),
     *          RRORB(JL)*ZPM(JL,JK,iror),
     *          1.17*RRORA(JL)*ZPM(JL,JK,iror),
     *          RHCOOHOH(JL)*ZPM(JL,JK,ihcooh)*ZPM(JL,JK,ioh),
     *          RC23MO2A(JL)*ZPM(JL,JK,ic2o3)*ZPM(JL,JK,ich3o2), 
     *          RHNO4M(JL)*ZPM(JL,JK,ihno4),
     *          ZPM(JL,JK,ioh)*RCOOH(JL)*ZPM(JL,JK,ico),
     *          ZPM(JL,JK,ioh)*RO3OH(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,ioh)*RHPOH(JL)*ZPM(JL,JK,ih2o2),
     *          ZPM(JL,JK,ioh)*RFRMOH(JL)*ZPM(JL,JK,ich2o),
     *          ZPM(JL,JK,ioh)*RH2OH(JL)*ZPRHOA(JL,JK)*6.022045E23/ZMAIR,
     *          ZPM(JL,JK,ioh)*RETHOH(JL)*ZPM(JL,JK,ieth),
     *          ZPM(JL,JK,ioh)*ROLEOH(JL)*ZPM(JL,JK,iole),
     *          ZPM(JL,JK,ioh)*0.91*RISOPOH(JL)*ZPM(JL,JK,iisop),
     *          ZPM(JL,JK,ioh)*0.68*RISPDOH(JL)*ZPM(JL,JK,iisoprd),
     *          ZPM(JL,JK,ioh)*0.5*RMTHCOH(JL)*ZPM(JL,JK,imethac),
     *          ZPM(JL,JK,ioh)*0.3*RMVKOH(JL)*ZPM(JL,JK,imvk),
     *          ZPM(JL,JK,ioh)*0.11*RPAROH(JL)*ZPM(JL,JK,ipar),
     *          ZPM(JL,JK,ioh)*0.4*RMPANOH(JL)*ZPM(JL,JK,impan),

C LG-    and destruction terms

     *         -ZPM(JL,JK,iho2)*RHO2NO(JL)*ZPM(JL,JK,ino),
     *         -ZPM(JL,JK,iho2)*RO3HO2(JL)*ZPM(JL,JK,io3),
     *         -ZPM(JL,JK,iho2)*RNO2HO2(JL)*ZPM(JL,JK,ino2),
     *         -ZPM(JL,JK,iho2)*RMO2HO2(JL)*ZPM(JL,JK,ich3o2),
     *         -ZPM(JL,JK,iho2)*RHO2OH(JL)*ZPM(JL,JK,ioh),
     *         -ZPM(JL,JK,iho2)*RC23HO2(JL)*ZPM(JL,JK,ic2o3),
     *         -ZPM(JL,JK,iho2)*RXO2HO2(JL)*ZPM(JL,JK,ixo2),
     *         -ZPM(JL,JK,iho2)*RXO2NHO2(JL)*ZPM(JL,JK,ixo2n),
     *         -ZPM(JL,JK,iho2)*RBXO2NHO2(JL)*ZPM(JL,JK,ibxo2n),
     *         -ZPM(JL,JK,iho2)*RMC23HO2(JL)*ZPM(JL,JK,imc3o3),
     *         -ZPM(JL,JK,iho2)*2.*RHO2HO2(JL)*ZPM(JL,JK,iho2)

         ELSEIF (NTRAC.LT.NTRACT.AND.ino.EQ.1) THEN 

C LG-    for NTRAC < NTRACT and ino=1
 
C WP-    O1d is needed for an estimate of the loss reaction h2o+o1d
	
         O1D=ZPM(JL,JK,io3)*ZPPHOTCHEM(JL,JK,7)/(RODM(JL)*AIR(JL))	
       
C WP-    Ozone Producing Potential calculation, this is actual the sum
C        of all the reactions (being calculated within the IDL script
C        being the total chemical tendency!), The reaction between NO and
C        O3 and the photodissociation of NO2 are not considered assuming
C        a photostationary equilibrium 

         OPP(JK)=

C LG-           production terms 

     * 	ZPMLOC(JL,JK,ino)*RHO2NO(JL)*ZPMLOC(JL,JK,iho2)		!,'NO+HO2'
     * 	+ZPMLOC(JL,JK,ino)*RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2)	!,'NO+CH3O2'
     *	+ZPMLOC(JL,JK,ino)*RC23NO(JL)*ZPMLOC(JL,JK,ic2o3)	!,'NO+C2O3'
     * 	+ZPMLOC(JL,JK,ino)*RXO2NO(JL)*ZPMLOC(JL,JK,ixo2)	!,'NO+XO2'
     * 	+ZPMLOC(JL,JK,ioh)*ROHOH(JL)*ZPMLOC(JL,JK,ioh)		!,'OH+OH'
     * 	+ZPMLOC(JL,JK,ic2o3)*RC23HO2(JL)*ZPMLOC(JL,JK,iho2)	!,'C23+HO2'

C LG-           and destruction terms

     * 	-ZPM(JL,JK,io3)*RO3OH(JL)*ZPMLOC(JL,JK,ioh)		!,'O3+OH'
     * 	-ZPM(JL,JK,io3)*RO3HO2(JL)*ZPMLOC(JL,JK,iho2)		!,'O3+HO2'
     * 	-ZPM(JL,JK,io3)*ROLEO3(JL)*ZPM(JL,JK,iole)		!,'O3+OLE'
     * 	-ZPM(JL,JK,io3)*RISOPO3(JL)*ZPM(JL,JK,iisop)		!,'O3+ISOP'
     * 	-ZPM(JL,JK,io3)*RISPDO3(JL)*ZPM(JL,JK,iisoprd)		!,'O3+ISOPRD'
     * 	-ZPM(JL,JK,io3)*RMVKO3(JL)*ZPM(JL,JK,imvk)		!,'O3+MVK'
     * 	-ZPM(JL,JK,io3)*RMPANO3(JL)*ZPM(JL,JK,impan)		!,'O3+MPAN'
     * 	-RH2OOD(JL)*H2O(JL)*O1D					!,'H2O+O1D'

C LG- end

         WRITE(NUNREAC,*)          
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     * 	ZPMLOC(JL,JK,ino)*RHO2NO(JL)*ZPMLOC(JL,JK,iho2),	!,'NO+HO2'
     * 	ZPMLOC(JL,JK,ino)*RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2),	!,'NO+CH3O2'
     *	ZPMLOC(JL,JK,ino)*RC23NO(JL)*ZPMLOC(JL,JK,ic2o3),	!,'NO+C2O3'
     * 	ZPMLOC(JL,JK,ino)*RXO2NO(JL)*ZPMLOC(JL,JK,ixo2),	!,'NO+XO2'
     * 	ZPMLOC(JL,JK,ioh)*ROHOH(JL)*ZPMLOC(JL,JK,ioh),		!,'OH+OH'
     * 	ZPMLOC(JL,JK,ic2o3)*RC23HO2(JL)*ZPMLOC(JL,JK,iho2),	!,'C23+HO2'
     * 	-ZPM(JL,JK,io3)*RO3OH(JL)*ZPMLOC(JL,JK,ioh),		!,'O3+OH'
     * 	-ZPM(JL,JK,io3)*RO3HO2(JL)*ZPMLOC(JL,JK,iho2),		!,'O3+HO2'
     * 	-ZPM(JL,JK,io3)*ROLEO3(JL)*ZPM(JL,JK,iole),		!,'O3+OLE'
     * 	-ZPM(JL,JK,io3)*RISOPO3(JL)*ZPM(JL,JK,iisop),		!,'O3+ISOP'
     * 	-ZPM(JL,JK,io3)*RISPDO3(JL)*ZPM(JL,JK,iisoprd),		!,'O3+ISOPRD'
     * 	-ZPM(JL,JK,io3)*RMVKO3(JL)*ZPM(JL,JK,imvk),		!,'O3+MVK'
     * 	-ZPM(JL,JK,io3)*RMPANO3(JL)*ZPM(JL,JK,impan),		!,'O3+MPAN'
     * 	-RH2OOD(JL)*H2O(JL)*O1D					!,'H2O+O1D'

C LG-    production terms of formaldehyde

         WRITE(NUNREAC,*)
     *      RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,ino),
     *      RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *      ROHPFRM(JL)*ZPM(JL,JK,ich3o2h)*ZPMLOC(JL,JK,ioh),
     *      RJALD2(JL)*ZPM(JL,JK,iald2),
     *      1.33*RMO2MO2(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,ich3o2),
     *      RC23NO(JL)*ZPMLOC(JL,JK,ic2o3)*ZPMLOC(JL,JK,ino),
     *      2*RC23C23(JL)*ZPMLOC(JL,JK,ic2o3)*ZPMLOC(JL,JK,ic2o3),
     *      0.79*RC23HO2(JL)*ZPMLOC(JL,JK,iho2),
     *      RC23MO2(JL)*ZPMLOC(JL,JK,ich3o2),
     *      1.56*RETHOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,ieth),
     *      ROLEOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,iole),
     *      0.63*RISOPOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,iisop),
     *      0.14*RISPDOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,iisoprd),
     *      0.08*RMTHCOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,imethac),
     *      0.3*RMVKOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,imvk),
     *      0.4*RMPANOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,impan),
     *      0.5*RMEKOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,imek),
     *      RETHO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,ieth),
     *      0.74*ROLEO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iole),
     *      0.6*RISOPO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iisop),
     *      0.39*RISPDO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iisoprd),
     *      0.7*RMPANO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,impan),
     *      0.2*RMTHCO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,imethac),
     *      0.1*RMVKO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,imvk),
     *      ROLENO3(JL)*ZPMLOC(JL,JK,ino3)*ZPM(JL,JK,iole),
     *      0.33*RISPDNO3(JL)*ZPMLOC(JL,JK,ino3)*ZPM(JL,JK,iisoprd),
     *      RMC23NO(JL)*ZPMLOC(JL,JK,imc3o3)*ZPMLOC(JL,JK,ino),
     *      2*RMC23HO2(JL)*ZPMLOC(JL,JK,imc3o3)*ZPMLOC(JL,JK,iho2),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ich2o)*RJACH2O(JL),
     *          -ZPM(JL,JK,ich2o)*RJBCH2O(JL),
     *          -ZPM(JL,JK,ich2o)*RFRMOH(JL)*ZPMLOC(JL,JK,ioh),        
     *          -ZPM(JL,JK,ich2o)*RFRMNO3(JL)*ZPMLOC(JL,JK,ino3)

C LG-    production terms of Isoprene

         WRITE(NUNREAC,*),

C LG-    and destruction terms

     *          -ZPM(JL,JK,iisop)*RISOPOH(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,iisop)*RISOPO3(JL)*ZPM(JL,JK,io3),
     *          -ZPM(JL,JK,iisop)*RISOPNO3(JL)*ZPMLOC(JL,JK,ino3),
     *          -ZPM(JL,JK,iisop)*RISOPNO2(JL)*ZPMLOC(JL,JK,ino2)

C LG-    production terms of NO

         WRITE(NUNREAC,*),
     *          RJNO2(JL)*ZPMLOC(JL,JK,ino2),
     *          RJBNO3(JL)*ZPMLOC(JL,JK,ino3),
     *          ZPMLOC(JL,JK,ino2)*0.2*RISOPNO2(JL)*ZPM(JL,JK,iisop),

C LG-    and destruction terms

     *          -ZPMLOC(JL,JK,iho2)*RHO2NO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ich3o2)*RMO2NO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ic2o3)*RC23NO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ixo2)*RXO2NO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,imc3o3)*RMC23NO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ino)*RNONO3(JL)*ZPMLOC(JL,JK,ino3),
     *          -ZPMLOC(JL,JK,ixo2n)*RXO2NNO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ibxo2n)*RBXO2NNO(JL)*ZPMLOC(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ino)*RNOO3(JL)*ZPM(JL,JK,io3)

C LG-    production terms of NO2

         WRITE(NUNREAC,*),
     *          ZPMLOC(JL,JK,iho2)*RHO2NO(JL)*ZPMLOC(JL,JK,ino),
     *          ZPMLOC(JL,JK,ich3o2)*RMO2NO(JL)*ZPMLOC(JL,JK,ino),
     *          ZPMLOC(JL,JK,ic2o3)*RC23NO(JL)*ZPMLOC(JL,JK,ino),
     *          ZPMLOC(JL,JK,ixo2)*RXO2NO(JL)*ZPMLOC(JL,JK,ino),
     *          ZPMLOC(JL,JK,imc3o3)*RMC23NO(JL)*ZPMLOC(JL,JK,ino),
     *          ZPMLOC(JL,JK,ino)*RNOO3(JL)*ZPM(JL,JK,io3),
     *          RJHNO3(JL)*ZPM(JL,JK,ihno3),
     *          RJN2O5(JL)*ZPMLOC(JL,JK,in2o5),
     *          RN2O5(JL)*ZPMLOC(JL,JK,in2o5),
     *          RJHNO4(JL)*ZPMLOC(JL,JK,ihno4),
     *          RHNO4M(JL)*ZPMLOC(JL,JK,ihno4),
     *          ZPMLOC(JL,JK,ihno4)*RHNO4OH(JL)*ZPMLOC(JL,JK,ioh),
     *          ZPMLOC(JL,JK,ino)*2.*RNONO3(JL)*ZPMLOC(JL,JK,ino3),
     *          RPAN(JL)*ZPM(JL,JK,ipan),
     *          RJPAN(JL)*ZPM(JL,JK,ipan),
     *          RMPAN(JL)*ZPM(JL,JK,impan),
     *          RJPAN(JL)*ZPM(JL,JK,impan),
     *          ZPM(JL,JK,impan)*0.70*RMPANO3(JL)*ZPM(JL,JK,io3),
     *          RJANO3(JL)*ZPMLOC(JL,JK,ino3),
     *          ZPM(JL,JK,iole)*ROLENO3(JL)*ZPMLOC(JL,JK,ino3),
     *          ZPM(JL,JK,iisop)*0.2*RISOPNO3(JL)*ZPMLOC(JL,JK,ino3),
     *          RJNITR(JL)*ZPM(JL,JK,intr),
     *          ZPM(JL,JK,iisontr)*RISONTROH(JL)*ZPMLOC(JL,JK,ioh),

C LG-    and destruction terms

     *          -ZPMLOC(JL,JK,ino2)*RJNO2(JL),
     *          -ZPMLOC(JL,JK,ino2)*RNO2OH(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPMLOC(JL,JK,ino2)*RNO2NO3(JL)*ZPMLOC(JL,JK,ino3),
     *          -ZPMLOC(JL,JK,ino2)*RNO2HO2(JL)*ZPMLOC(JL,JK,iho2),
     *          -ZPMLOC(JL,JK,ino2)*RNO2O3(JL)*ZPM(JL,JK,io3),
     *          -ZPMLOC(JL,JK,ic2o3)*RC23NO2(JL)*ZPMLOC(JL,JK,ino2),
     *          -ZPMLOC(JL,JK,iror)*RRORNO2(JL)*ZPMLOC(JL,JK,ino2),
     *          -ZPM(JL,JK,iisop)*0.8*RISOPNO2(JL)*ZPMLOC(JL,JK,ino2),
     *          -ZPMLOC(JL,JK,imc3o3)*RMC23NO2(JL)*ZPMLOC(JL,JK,ino2)

C LG-    production terms of H2O2

         WRITE(NUNREAC,*),RHO2HO2(JL)*ZPMLOC(JL,JK,iho2)*ZPMLOC(JL,JK,iho2),
     *                    RRCHOHO2H(JL)*ZPM(JL,JK,irchoho2h),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ih2o2)*RHPOH(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ih2o2)*RJH2O2(JL)

C LG-    production terms of CH3O2H

         WRITE(NUNREAC,*),RMO2HO2(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,iho2),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ich3o2h)*ROHPCAT(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ich3o2h)*ROHPFRM(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ich3o2h)*RJMEPE(JL)

C LG-    production terms of OH

         WRITE(NUNREAC2,*),
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     *          RJHNO3(JL)*ZPM(JL,JK,ihno3),
     *          2.*RJO3D(JL)*ZPM(JL,JK,io3),
     *          2.*RJH2O2(JL)*ZPM(JL,JK,ih2o2),
     *          RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *          ZPM(JL,JK,iisoprd)*0.4*RISPDO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,iole)*0.1*ROLEO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,iisop)*0.2*RISOPO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imethac)*0.1*RMTHCO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imvk)*0.05*RMVKO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,impan)*0.04*RMPANO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imaterp)*0.76*RMATERPO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,imbterp)*0.33*RMBTERPO3(JL)*ZPM(JL,JK,io3),
     *          ZPM(JL,JK,isqterp)*0.22*RSQTERPO3(JL)*ZPM(JL,JK,io3),
     *          RJHONO(JL)*ZPM(JL,JK,ihono),
     *          POHHO2(JL,JK)*ZPMLOC(JL,JK,iho2),

C LG-    and destruction terms
     *         -ZPMLOC(JL,JK,ioh)*RHNO4OH(JL)*ZPMLOC(JL,JK,ihno4),
     *         -ZPMLOC(JL,JK,ioh)*RHO2OH(JL)*ZPMLOC(JL,JK,iho2),
     *         -ZPMLOC(JL,JK,ioh)*RNO2OH(JL)*ZPMLOC(JL,JK,ino2),
     *         -ZPMLOC(JL,JK,ioh)*ROHHNO3(JL)*ZPM(JL,JK,ihno3),
     *         -ZPMLOC(JL,JK,ioh)*RCOOH(JL)*ZPM(JL,JK,ico),
     *         -ZPMLOC(JL,JK,ioh)*RO3OH(JL)*ZPM(JL,JK,io3),
     *         -ZPMLOC(JL,JK,ioh)*RHPOH(JL)*ZPM(JL,JK,ih2o2),
     *         -ZPMLOC(JL,JK,ioh)*RFRMOH(JL)*ZPM(JL,JK,ich2o),
     *         -ZPMLOC(JL,JK,ioh)*RH2OH(JL)*ZPRHOA(JL,JK)*6.022045E23/ZMAIR,
     *         -ZPMLOC(JL,JK,ioh)*RETHOH(JL)*ZPM(JL,JK,ieth),
     *         -ZPMLOC(JL,JK,ioh)*ROLEOH(JL)*ZPM(JL,JK,iole),
     *         -ZPMLOC(JL,JK,ioh)*0.91*RISOPOH(JL)*ZPM(JL,JK,iisop),
     *         -ZPMLOC(JL,JK,ioh)*0.68*RISPDOH(JL)*ZPM(JL,JK,iisoprd),
     *         -ZPMLOC(JL,JK,ioh)*0.5*RMTHCOH(JL)*ZPM(JL,JK,imethac),
     *         -ZPMLOC(JL,JK,ioh)*0.3*RMVKOH(JL)*ZPM(JL,JK,imvk),
     *         -ZPMLOC(JL,JK,ioh)*0.11*RPAROH(JL)*ZPM(JL,JK,ipar),
     *         -ZPMLOC(JL,JK,ioh)*0.4*RMPANOH(JL)*ZPM(JL,JK,impan),
     *         -ZPMLOC(JL,JK,ioh)*RCH4OH(JL)*ZPM(JL,JK,ich4),
     *         -ZPMLOC(JL,JK,ioh)*ROHPCAT(JL)*ZPM(JL,JK,ich3o2h),
     *         -ZPMLOC(JL,JK,ioh)*ROHOH(JL)*ZPMLOC(JL,JK,ioh),
     *         -ZPMLOC(JL,JK,ioh)*0.89*RPAROH(JL)*ZPM(JL,JK,ipar),
     *         -ZPMLOC(JL,JK,ioh)*RMGLYOH(JL)*ZPM(JL,JK,imgly),
     *         -ZPMLOC(JL,JK,ioh)*0.09*RISOPOH(JL)*ZPM(JL,JK,iisop),
     *         -ZPMLOC(JL,JK,ioh)*0.32*RISPDOH(JL)*ZPM(JL,JK,iisoprd),
     *         -ZPMLOC(JL,JK,ioh)*0.5*RMTHCOH(JL)*ZPM(JL,JK,imethac),
     *         -ZPMLOC(JL,JK,ioh)*RMEKOH(JL)*ZPM(JL,JK,imek),
     *         -ZPMLOC(JL,JK,ioh)*0.7*RMVKOH(JL)*ZPM(JL,JK,imvk),
     *         -ZPMLOC(JL,JK,ioh)*RALD2OH(JL)*ZPM(JL,JK,iald2),
     *         -ZPMLOC(JL,JK,ioh)*0.6*RMPANOH(JL)*ZPM(JL,JK,impan),
     *         -ZPMLOC(JL,JK,ioh)*RACETOH(JL)*ZPM(JL,JK,iacet),
     *         -ZPMLOC(JL,JK,ioh)*RMATERPOH(JL)*ZPM(JL,JK,imaterp),
     *         -ZPMLOC(JL,JK,ioh)*RMBTERPOH(JL)*ZPM(JL,JK,imbterp),
     *         -ZPMLOC(JL,JK,ioh)*RSQTERPOH(JL)*ZPM(JL,jk,isqterp),
     *         -ZPMLOC(JL,JK,ioh)*RNOOH(JL)*ZPMLOC(JL,JK,ino)

C LG-    production terms of HO2

         WRITE(NUNREAC3,*),
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     *          2.*RJBCH2O(JL)*ZPM(JL,JK,ich2o),
     *          RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *          2.*RJALD2(JL)*ZPM(JL,JK,iald2),
     *          RJNITR(JL)*ZPM(JL,JK,intr),
     *          RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,ino),
     *          0.66*RMO2MO2(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,ich3o2),
     *          ZPMLOC(JL,JK,ic2o3)*RC23NO(JL)*ZPMLOC(JL,JK,ino),     
     *          ZPMLOC(JL,JK,ic2o3)*2.*RC23C23(JL)*ZPMLOC(JL,JK,ic2o3),
     *          ZPMLOC(JL,JK,ic2o3)*RC23MO2(JL)*ZPMLOC(JL,JK,ich3o2),
     *          ZPM(JL,JK,io3)*0.44*ROLEO3(JL)*ZPM(JL,JK,iole),
     *          ZPM(JL,JK,io3)*0.12*RETHO3(JL)*ZPM(JL,JK,ieth),
     *          ZPM(JL,JK,io3)*0.23*RISPDO3(JL)*ZPM(JL,JK,iisoprd),
     *          ZPM(JL,JK,io3)*0.1*RMTHCO3(JL)*ZPM(JL,JK,imethac),
     *          ZPM(JL,JK,io3)*0.08*RMPANO3(JL)*ZPM(JL,JK,impan),
     *          ZPMLOC(JL,JK,ino3)*0.8*RISOPNO3(JL)*ZPM(JL,JK,iisop),
     *          ZPMLOC(JL,JK,ino3)*RISPDNO3(JL)*ZPM(JL,JK,iisoprd),
     *          ZPMLOC(JL,JK,ino3)*0.5*RMTHCNO3(JL)*ZPM(JL,JK,imethac),
     *          ZPMLOC(JL,JK,ino3)*RFRMNO3(JL)*ZPM(JL,JK,ich2o),
     *          RISOPNO2(JL)*ZPM(JL,JK,iisop)*ZPMLOC(JL,JK,ino2),
     *          RRORB(JL)*ZPMLOC(JL,JK,iror),
     *          1.17*RRORA(JL)*ZPMLOC(JL,JK,iror),
     *          RHCOOHOH(JL)*ZPM(JL,JK,ihcooh)*ZPMLOC(JL,JK,ioh),
     *          RC23MO2A(JL)*ZPMLOC(JL,JK,ic2o3)*ZPMLOC(JL,JK,ich3o2), 
     *          RHNO4M(JL)*ZPMLOC(JL,JK,ihno4),
     *          ZPMLOC(JL,JK,ioh)*RCOOH(JL)*ZPM(JL,JK,ico),
     *          ZPMLOC(JL,JK,ioh)*RO3OH(JL)*ZPM(JL,JK,io3),
     *          ZPMLOC(JL,JK,ioh)*RHPOH(JL)*ZPM(JL,JK,ih2o2),
     *          ZPMLOC(JL,JK,ioh)*RFRMOH(JL)*ZPM(JL,JK,ich2o),
     *          ZPMLOC(JL,JK,ioh)*RH2OH(JL)*ZPRHOA(JL,JK)*6.022045E23/ZMAIR,
     *          ZPMLOC(JL,JK,ioh)*RETHOH(JL)*ZPM(JL,JK,ieth),
     *          ZPMLOC(JL,JK,ioh)*ROLEOH(JL)*ZPM(JL,JK,iole),
     *          ZPMLOC(JL,JK,ioh)*0.91*RISOPOH(JL)*ZPM(JL,JK,iisop),
     *          ZPMLOC(JL,JK,ioh)*0.68*RISPDOH(JL)*ZPM(JL,JK,iisoprd),
     *          ZPMLOC(JL,JK,ioh)*0.5*RMTHCOH(JL)*ZPM(JL,JK,imethac),
     *          ZPMLOC(JL,JK,ioh)*0.3*RMVKOH(JL)*ZPM(JL,JK,imvk),
     *          ZPMLOC(JL,JK,ioh)*0.11*RPAROH(JL)*ZPM(JL,JK,ipar),
     *          ZPMLOC(JL,JK,ioh)*0.4*RMPANOH(JL)*ZPM(JL,JK,impan),

C LG-    and destruction terms

     *         -ZPMLOC(JL,JK,iho2)*RHO2NO(JL)*ZPMLOC(JL,JK,ino),
     *         -ZPMLOC(JL,JK,iho2)*RO3HO2(JL)*ZPM(JL,JK,io3),
     *         -ZPMLOC(JL,JK,iho2)*RNO2HO2(JL)*ZPMLOC(JL,JK,ino2),
     *         -ZPMLOC(JL,JK,iho2)*RMO2HO2(JL)*ZPMLOC(JL,JK,ich3o2),
     *         -ZPMLOC(JL,JK,iho2)*RHO2OH(JL)*ZPMLOC(JL,JK,ioh),
     *         -ZPMLOC(JL,JK,iho2)*RC23HO2(JL)*ZPMLOC(JL,JK,ic2o3),
     *         -ZPMLOC(JL,JK,iho2)*RXO2HO2(JL)*ZPMLOC(JL,JK,ixo2),
     *         -ZPMLOC(JL,JK,iho2)*RXO2NHO2(JL)*ZPMLOC(JL,JK,ixo2n),
     *         -ZPMLOC(JL,JK,iho2)*RBXO2NHO2(JL)*ZPMLOC(JL,JK,ibxo2n),
     *         -ZPMLOC(JL,JK,iho2)*RMC23HO2(JL)*ZPMLOC(JL,JK,imc3o3),
     *         -ZPMLOC(JL,JK,iho2)*2.*RHO2HO2(JL)*ZPMLOC(JL,JK,iho2)

         ELSEIF (NTRAC.LT.NTRACT.AND.ino.GT.1) THEN 
	
C LG-    for NTRAC < NTRACT and ino > 1
 
C WP-    O1d is needed for an estimate of the loss reaction h2o+o1d
	
         O1D=ZPM(JL,JK,io3)*ZPPHOTCHEM(JL,JK,7)/(RODM(JL)*AIR(JL))	
       
C WP-    Ozone Producing Potential calculation, this is actual the sum
C        of all the reactions (being calculated within the IDL script
C        being the total chemical tendency!), The reaction between NO and
C        O3 and the photodissociation of NO2 are not considered assuming
C        a photostationary equilibrium 

         OPP(JK)=

C LG-           production terms 

     * 	ZPM(JL,JK,ino)*RHO2NO(JL)*ZPMLOC(JL,JK,iho2)		!,'NO+HO2'
     * 	+ZPM(JL,JK,ino)*RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2)	!,'NO+CH3O2'
     *	+ZPM(JL,JK,ino)*RC23NO(JL)*ZPMLOC(JL,JK,ic2o3)	!,'NO+C2O3'
     * 	+ZPM(JL,JK,ino)*RXO2NO(JL)*ZPMLOC(JL,JK,ixo2)	!,'NO+XO2'
     * 	+ZPMLOC(JL,JK,ioh)*ROHOH(JL)*ZPMLOC(JL,JK,ioh)		!,'OH+OH'
     * 	+ZPMLOC(JL,JK,ic2o3)*RC23HO2(JL)*ZPMLOC(JL,JK,iho2)	!,'C23+HO2'

C LG-           and destruction terms

     * 	-ZPM(JL,JK,io3)*RO3OH(JL)*ZPMLOC(JL,JK,ioh)		!,'O3+OH'
     * 	-ZPM(JL,JK,io3)*RO3HO2(JL)*ZPMLOC(JL,JK,iho2)		!,'O3+HO2'
     * 	-ZPM(JL,JK,io3)*ROLEO3(JL)*ZPM(JL,JK,iole)		!,'O3+OLE'
     * 	-ZPM(JL,JK,io3)*RISOPO3(JL)*ZPM(JL,JK,iisop)		!,'O3+ISOP'
     * 	-ZPM(JL,JK,io3)*RISPDO3(JL)*ZPM(JL,JK,iisoprd)		!,'O3+ISOPRD'
     * 	-ZPM(JL,JK,io3)*RMVKO3(JL)*ZPM(JL,JK,imvk)		!,'O3+MVK'
     * 	-ZPM(JL,JK,io3)*RMPANO3(JL)*ZPM(JL,JK,impan)		!,'O3+MPAN'
     * 	-RH2OOD(JL)*H2O(JL)*O1D					!,'H2O+O1D'

C LG- end

         WRITE(NUNREAC,*)          
     *          NSTEP,JK,ZPRHOA(JL,JK), 
     * 	ZPM(JL,JK,ino)*RHO2NO(JL)*ZPMLOC(JL,JK,iho2),	!,'NO+HO2'
     * 	ZPM(JL,JK,ino)*RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2),	!,'NO+CH3O2'
     *	ZPM(JL,JK,ino)*RC23NO(JL)*ZPMLOC(JL,JK,ic2o3),	!,'NO+C2O3'
     * 	ZPM(JL,JK,ino)*RXO2NO(JL)*ZPMLOC(JL,JK,ixo2),	!,'NO+XO2'
     * 	ZPMLOC(JL,JK,ioh)*ROHOH(JL)*ZPMLOC(JL,JK,ioh),		!,'OH+OH'
     * 	ZPMLOC(JL,JK,ic2o3)*RC23HO2(JL)*ZPMLOC(JL,JK,iho2),	!,'C23+HO2'
     * 	-ZPM(JL,JK,io3)*RO3OH(JL)*ZPMLOC(JL,JK,ioh),		!,'O3+OH'
     * 	-ZPM(JL,JK,io3)*RO3HO2(JL)*ZPMLOC(JL,JK,iho2),		!,'O3+HO2'
     * 	-ZPM(JL,JK,io3)*ROLEO3(JL)*ZPM(JL,JK,iole),		!,'O3+OLE'
     * 	-ZPM(JL,JK,io3)*RISOPO3(JL)*ZPM(JL,JK,iisop),		!,'O3+ISOP'
     * 	-ZPM(JL,JK,io3)*RISPDO3(JL)*ZPM(JL,JK,iisoprd),		!,'O3+ISOPRD'
     * 	-ZPM(JL,JK,io3)*RMVKO3(JL)*ZPM(JL,JK,imvk),		!,'O3+MVK'
     * 	-ZPM(JL,JK,io3)*RMPANO3(JL)*ZPM(JL,JK,impan),		!,'O3+MPAN'
     * 	-RH2OOD(JL)*H2O(JL)*O1D					!,'H2O+O1D'

C LG-    production terms of formaldehyde

         WRITE(NUNREAC,*)
     *      RMO2NO(JL)*ZPMLOC(JL,JK,ich3o2)*ZPM(JL,JK,ino),
     *      RJMEPE(JL)*ZPM(JL,JK,ich3o2h),
     *      ROHPFRM(JL)*ZPM(JL,JK,ich3o2h)*ZPMLOC(JL,JK,ioh),
     *      RJALD2(JL)*ZPM(JL,JK,iald2),
     *      1.33*RMO2MO2(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,ich3o2),
     *      RC23NO(JL)*ZPMLOC(JL,JK,ic2o3)*ZPM(JL,JK,ino),
     *      2*RC23C23(JL)*ZPMLOC(JL,JK,ic2o3)*ZPMLOC(JL,JK,ic2o3),
     *      0.79*RC23HO2(JL)*ZPMLOC(JL,JK,iho2),
     *      RC23MO2(JL)*ZPMLOC(JL,JK,ich3o2),
     *      1.56*RETHOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,ieth),
     *      ROLEOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,iole),
     *      0.63*RISOPOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,iisop),
     *      0.14*RISPDOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,iisoprd),
     *      0.08*RMTHCOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,imethac),
     *      0.3*RMVKOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,imvk),
     *      0.4*RMPANOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,impan),
     *      0.5*RMEKOH(JL)*ZPMLOC(JL,JK,ioh)*ZPM(JL,JK,imek),
     *      RETHO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,ieth),
     *      0.74*ROLEO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iole),
     *      0.6*RISOPO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iisop),
     *      0.39*RISPDO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,iisoprd),
     *      0.7*RMPANO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,impan),
     *      0.2*RMTHCO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,imethac),
     *      0.1*RMVKO3(JL)*ZPM(JL,JK,io3)*ZPM(JL,JK,imvk),
     *      ROLENO3(JL)*ZPM(JL,JK,ino3)*ZPM(JL,JK,iole),
     *      0.33*RISPDNO3(JL)*ZPM(JL,JK,ino3)*ZPM(JL,JK,iisoprd),
     *      RMC23NO(JL)*ZPMLOC(JL,JK,imc3o3)*ZPM(JL,JK,ino),
     *      2*RMC23HO2(JL)*ZPMLOC(JL,JK,imc3o3)*ZPMLOC(JL,JK,iho2),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ich2o)*RJACH2O(JL),
     *          -ZPM(JL,JK,ich2o)*RJBCH2O(JL),
     *          -ZPM(JL,JK,ich2o)*RFRMOH(JL)*ZPMLOC(JL,JK,ioh),        
     *          -ZPM(JL,JK,ich2o)*RFRMNO3(JL)*ZPM(JL,JK,ino3)

C LG-    production terms of Isoprene

         WRITE(NUNREAC,*),

C LG-    and destruction terms

     *          -ZPM(JL,JK,iisop)*RISOPOH(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,iisop)*RISOPO3(JL)*ZPM(JL,JK,io3),
     *          -ZPM(JL,JK,iisop)*RISOPNO3(JL)*ZPM(JL,JK,ino3),
     *          -ZPM(JL,JK,iisop)*RISOPNO2(JL)*ZPM(JL,JK,ino2)

C LG-    production terms of NO

         WRITE(NUNREAC,*),
     *          RJNO2(JL)*ZPM(JL,JK,ino2),
     *          RJBNO3(JL)*ZPM(JL,JK,ino3),
     *          ZPM(JL,JK,ino2)*0.2*RISOPNO2(JL)*ZPM(JL,JK,iisop),

C LG-    and destruction terms

     *          -ZPMLOC(JL,JK,iho2)*RHO2NO(JL)*ZPM(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ich3o2)*RMO2NO(JL)*ZPM(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ic2o3)*RC23NO(JL)*ZPM(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ixo2)*RXO2NO(JL)*ZPM(JL,JK,ino),
     *          -ZPMLOC(JL,JK,imc3o3)*RMC23NO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ino)*RNONO3(JL)*ZPM(JL,JK,ino3),
     *          -ZPMLOC(JL,JK,ixo2n)*RXO2NNO(JL)*ZPM(JL,JK,ino),
     *          -ZPMLOC(JL,JK,ibxo2n)*RBXO2NNO(JL)*ZPM(JL,JK,ino),
     *          -ZPM(JL,JK,ino)*RNOO3(JL)*ZPM(JL,JK,io3)

C LG-    production terms of NO2

         WRITE(NUNREAC,*),
     *          ZPMLOC(JL,JK,iho2)*RHO2NO(JL)*ZPM(JL,JK,ino),
     *          ZPMLOC(JL,JK,ich3o2)*RMO2NO(JL)*ZPM(JL,JK,ino),
     *          ZPMLOC(JL,JK,ic2o3)*RC23NO(JL)*ZPM(JL,JK,ino),
     *          ZPMLOC(JL,JK,ixo2)*RXO2NO(JL)*ZPM(JL,JK,ino),
     *          ZPMLOC(JL,JK,imc3o3)*RMC23NO(JL)*ZPM(JL,JK,ino),
     *          ZPM(JL,JK,ino)*RNOO3(JL)*ZPM(JL,JK,io3),
     *          RJHNO3(JL)*ZPM(JL,JK,ihno3),
     *          RJN2O5(JL)*ZPM(JL,JK,in2o5),
     *          RN2O5(JL)*ZPM(JL,JK,in2o5),
     *          RJHNO4(JL)*ZPM(JL,JK,ihno4),
     *          RHNO4M(JL)*ZPM(JL,JK,ihno4),
     *          ZPM(JL,JK,ihno4)*RHNO4OH(JL)*ZPMLOC(JL,JK,ioh),
     *          ZPM(JL,JK,ino)*2.*RNONO3(JL)*ZPM(JL,JK,ino3),
     *          RPAN(JL)*ZPM(JL,JK,ipan),
     *          RJPAN(JL)*ZPM(JL,JK,ipan),
     *          RMPAN(JL)*ZPM(JL,JK,impan),
     *          RJPAN(JL)*ZPM(JL,JK,impan),
     *          ZPM(JL,JK,impan)*0.70*RMPANO3(JL)*ZPM(JL,JK,io3),
     *          RJANO3(JL)*ZPM(JL,JK,ino3),
     *          ZPM(JL,JK,iole)*ROLENO3(JL)*ZPM(JL,JK,ino3),
     *          ZPM(JL,JK,iisop)*0.2*RISOPNO3(JL)*ZPM(JL,JK,ino3),
     *          RJNITR(JL)*ZPM(JL,JK,intr),
     *          ZPM(JL,JK,iisontr)*RISONTROH(JL)*ZPMLOC(JL,JK,ioh),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ino2)*RJNO2(JL),
     *          -ZPM(JL,JK,ino2)*RNO2OH(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ino2)*RNO2NO3(JL)*ZPM(JL,JK,ino3),
     *          -ZPM(JL,JK,ino2)*RNO2HO2(JL)*ZPMLOC(JL,JK,iho2),
     *          -ZPM(JL,JK,ino2)*RNO2O3(JL)*ZPM(JL,JK,io3),
     *          -ZPMLOC(JL,JK,ic2o3)*RC23NO2(JL)*ZPM(JL,JK,ino2),
     *          -ZPMLOC(JL,JK,iror)*RRORNO2(JL)*ZPM(JL,JK,ino2),
     *          -ZPM(JL,JK,iisop)*0.8*RISOPNO2(JL)*ZPM(JL,JK,ino2),
     *          -ZPMLOC(JL,JK,imc3o3)*RMC23NO2(JL)*ZPM(JL,JK,ino2)
 
C LG-    production terms of H2O2

         WRITE(NUNREAC,*),RHO2HO2(JL)*ZPMLOC(JL,JK,iho2)*ZPMLOC(JL,JK,iho2),
     *                    RRCHOHO2H(JL)*ZPM(JL,JK,irchoho2h),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ih2o2)*RHPOH(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ih2o2)*RJH2O2(JL)

C LG-    production terms of CH3O2H

         WRITE(NUNREAC,*),RMO2HO2(JL)*ZPMLOC(JL,JK,ich3o2)*ZPMLOC(JL,JK,iho2),

C LG-    and destruction terms

     *          -ZPM(JL,JK,ich3o2h)*ROHPCAT(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ich3o2h)*ROHPFRM(JL)*ZPMLOC(JL,JK,ioh),
     *          -ZPM(JL,JK,ich3o2h)*RJMEPE(JL)

         ENDIF

       ENDIF

       ENDIF ! ENDIF (LWRITEREAC) THEN

C LG-  end

  107 CONTINUE
C
  101 CONTINUE

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNREAC)

C LG- assigning the calculated values for the different levels to the 
C     proper array names

      DO 200 JK=1,NLEVEL+DLEVEL
      DO 200 JL=1,NLON

C LG- assigning the tracer concentrations

       DO JT=1,NTRAC
	 IF(LBULKVEG.AND.JK.EQ.NLEV+1) THEN
     	   PM(JL,JK,JT)=ZPM(JL,JK,JT)
	 ELSEIF(LVEG_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1.OR.
     &          LSNOW_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1) THEN
     	   PM(JL,JK,JT)=ZPM(JL,JK,JT)
         ELSEIF (LXTMZZ.AND.JK.EQ.NLEVEL+1) THEN
           PMZZ(JL,JT)=ZPM(JL,JK,JT)
	 ELSE
	   PM(JL,JK,JT)=ZPM(JL,JK,JT)
	 ENDIF
       ENDDO

C LG-  and the short lived species

       DO JG=1,KG3X
	 IF(LBULKVEG.AND.JK.EQ.NLEV+1) THEN
     	   PMLOC(JL,JK,JG)=ZPMLOC(JL,JK,JG)
	 ELSEIF(LVEG_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1.OR.
     &          LSNOW_MLAY.AND.JK.GE.NLEV+1.AND.JK.LT.NLEVEL+1) THEN
     	   PMLOC(JL,JK,JG)=ZPMLOC(JL,JK,JG)
         ELSEIF (LXTMZZ.AND.JK.EQ.NLEVEL+1) THEN
           PMZZ(JL,JG+NTRAC)=ZPMLOC(JL,JK,JG)
	 ELSE
	   PMLOC(JL,JK,JG)=ZPMLOC(JL,JK,JG)
	 ENDIF
       ENDDO

 200  CONTINUE

C LG- end

C
C      END budgetcalculation
C

      DO 401 JT=1,NTRAC
      DO 401 JK=1,NLEVEL
      DO 401 JL=1,NLON
	IF(LBULKVEG.AND.JK.EQ.NLEV+1) THEN
         DPM(JL,JK,JT)=(PM(JL,JK,JT)+DPM(JL,JK,JT))*GRVOL(JL,JK)
	ELSEIF(LVEG_MLAY.AND.JK.GE.NLEV+1.OR.
     &         LSNOW_MLAY.AND.JK.GE.NLEV+1) THEN
         DPM(JL,JK,JT)=(PM(JL,JK,JT)+DPM(JL,JK,JT))*GRVOL(JL,JK)
        ELSE
         DPM(JL,JK,JT)=(PM(JL,JK,JT)+DPM(JL,JK,JT))*GRVOL(JL,JK)
 	ENDIF
  401 CONTINUE

      DO 402 JG=1,KG3X
      DO 402 JK=1,NLEVEL
      DO 402 JL=1,NLON
	IF(LBULKVEG.AND.JK.EQ.NLEV+1) THEN
         DPMLOC(JL,JK,JG)=(PMLOC(JL,JK,JG)+DPMLOC(JL,JK,JG))*
     &                     GRVOL(JL,JK)
	ELSEIF(LVEG_MLAY.AND.JK.GE.NLEV+1.OR.
     &         LSNOW_MLAY.AND.JK.GE.NLEV+1) THEN
         DPMLOC(JL,JK,JG)=(PMLOC(JL,JK,JG)+DPMLOC(JL,JK,JG))*
     &                     GRVOL(JL,JK)
        ELSE
         DPMLOC(JL,JK,JG)=(PMLOC(JL,JK,JG)+DPMLOC(JL,JK,JG))*
     &                     GRVOL(JL,JK)
	ENDIF
  402 CONTINUE

      CALL BUDCALC(NLEV,NLEVEL,DPM,DPMLOC,KG3X,IBCHEM)

C LG- 

      WRITE(NUNMDFL,*)'End CBM4_ECH.f'

      RETURN
      END

C LG- function from echam5 (mo_mz_mecca_kpp_mem.f90) to calculate
C     the termolecular reaction rate for NO+OH+M -> HONO

      REAL FUNCTION atk_3 (ztemp,zairfac,pa1,pa2,pb1,pb2,pfc)
! calculate third body reactions according to Atkinson '92

! ztemp: temperature in [K]
! airfac : concentration of air molecules in [molecules air/cm^3]

      REAL  :: za0,zb0,zx2
      REAL  :: pa1,pa2,pb1,pb2,pfc,ztemp,zairfac
 
      INTRINSIC  LOG10

      za0       = pa1*zairfac*(ztemp/300.)**pa2
      zb0       = pb1*(ztemp/300.)**pb2
      zx2       = pfc

      atk_3     = (za0/(1.0+za0/zb0))*(zx2**(1./(1.0+log10(za0/zb0)* 
     &   log10(za0/zb0))))

      END FUNCTION atk_3
