      SUBROUTINE XTTROPO(                                                                                                                                                     
     *      KSTEP, NSTART, KLON, KLEV,  ! mz_lg_20050829+ bugfixer, NSTART

C LG- adding the total number of levels

     *      KLEVEL,

C LG- end

     *      KROW,  KTRAC,
     *      PTM1, PQM1, PAPM1, PAPHM1, DAPM1, DAPHM1, 
     *      PGEOM1, PXTM1, KTROPO)                                                                                

C                                                                                               
C   *XTTROPO* CALCULATES THE TROPOPAUSE HEIGHT                                                   
C                                                                                                 
C    JOHANN FEICHTER         UNI-HH                                                              
C    REWRITTEN  U. SCHLESE    DKRZ - HAMBURG    SEP 93                                          
C    REWRITTEN FOR MULTI-TASKING   U. SCHLESE  JUNE 94                                           
C   PURPOSE                                                                                      
C   *XTTROPO* DETERMINES THE HIGHEST SIGMA-LEVEL IN THE TROPOSPHERE                             
C   "KTROPO"                                                                                  
C   DEPENDING ON THE VERTICAL GRADIENT OF THE POTENTIAL TEMPERATURE                            
C   AND CALCULATES ZONAL MASS-BUDGETS OF THE DIFFERENT TRACERS                                  
C                                                                                               
C   INTERFACE                                                                                 
C                                                                                             
C   *XTTROPO* IS CALLED FROM *PHYSC*                                                             
C                                                                                               
C   NO EXTERNALS                                                                                
C                                                                                               
C                                                                                                
C*    *PARAMETERS* CONTROLLING ARRAY SIZES.                                                      
C                                                                                                 
C       ----------------------------------------------------------------                         
C                                                                                                 

C LG- 

      IMPLICIT NONE

C LG- including common blocks

      INCLUDE 'comcon.h'


*/********************* budget calculations ********************* 
*/
*C XTTROPO
*I XTTROPO.33
C LG- The file PARCOMS contains the COMMON block BUD

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

C LG- declarations

      INTEGER KLON,KROW,KTRAC,KSTEP,NSTART,KLEV,
     *        KLEVELS,KLEVELSVEG,ITOP,ITOPM1,KTP1,ISTRAT,ITROP,IPBL,
     *        JK,JT,JJ

      REAL ZQG,ZFAC,ZLIMDTHDZ,ZKAPPA,ZTHETAN,ZDZ,ZTHETAP,ZDTHETA

      INTEGER KTROPO(KLON)                                                                       

      INTEGER JL,JGL

      REAL FIELD(24),BUDW(NLATT30)
      
      REAL PTM1(KLON,KLEV), PQM1(KLON,KLEV), PAPM1(KLON,KLEV),                                   
     *     PGEOM1(KLON,KLEV), PAPHM1(KLON,KLEV), PXTM1(KLON,NLEVT,KTRAC)                                                                
      REAL ZMA(KLON,NLEVT), ZMT(KLON,NLEVT,KTRAC),                                                
     *     ZTFAC(KLON,KLEV),  ZSFAC(KLON,KLEV)                                                   

C LG- more declarations
	
      REAL ZT(KLON,NLEVT),ZAPHM1(KLON,NLEVT+1),
     *     ZAPM1(KLON,NLEVT),ZHP(KLON,NLEVT+1),DZ(KLON,NLEVT),
     *     ZDP(KLON,NLEVT),ZQ(KLON,NLEVT),ZRHOA(KLON,NLEVT),Z(KLON,NLEVT)

C LG- pressure levels within the canopy, kept constant throughout the
C     model intergration at the value calculated for the first timestep 
C     from the defined reference heights as a function of the canopy 
C     height and the number of vegetation layers.  

      REAL DAPM1(NNLON,NLEVT),DAPHM1(NNLON,NLEVT+1)
      
      REAL TVIRT

      INTEGER KLEVEL,KLEVELP1

      SAVE BUDW

C LG- 

      WRITE(NUNMDFL,*)'Start XTTROPO.f'

C LG- determiming the number of levels for which tracer fields must
C     be calculated

      IF (LBIOSPH.AND.LAGRIAN) THEN
       KLEVELS=NLEVT
       KLEVELSVEG=NLEVV
      ELSEIF (LBULKVEG) THEN
       KLEVELS=KLEVEL+1
       KLEVELSVEG=1
      ELSEIF (LVEG_MLAY) THEN
       KLEVELS=KLEVEL+NLEVV_ML
       KLEVELSVEG=NLEVV_ML
      ELSEIF (LSNOW_MLAY) THEN
       KLEVELS=KLEVEL+NLEVS_ML
       KLEVELSVEG=NLEVS_ML
      ELSE
       KLEVELS=KLEVEL
       KLEVELSVEG=NLEVVEG
      ENDIF

C LG- reseting of PXTM1 at the first timestep
 
      IF (KSTEP.EQ.NSTART) THEN
       CALL RESETR(PXTM1,KLON*KLEVELS*NTRAC,0.)                                                                                                                                 
      ENDIF

      ZQG=1./G                                                                                  

C LG- the term BUDW(IROW) is normally calculated in ECHAM and put into
C     a common block, which is not available in RACMO and the calculated 
C     field of BUDW in ECHAM T30 is therefore written to an output file
C     which is read in this programm (for BUDW, see ECHAM source code, 
C     subroutine INIGAU, the values are defined for 24 longitudes and then 
C     copied to an array of 48 longitudes

      IF (KSTEP.EQ.NSTART) THEN
       DO JGL=1,24
        READ(118,*)JL,FIELD(JL)
        BUDW(2*JGL-1)=FIELD(JL)
        BUDW(2*JGL)=FIELD(JL)
       ENDDO
      ENDIF

C LG- ZFAC is a resolution dependent parameter which is used to determine
C     the total amount of mass of each specific grid square.

      ZFAC=BUDW(IROW)*510.0644719E+12*ZQG                                                        

      ITOP=KLEV-5                                                                                
      ITOPM1=ITOP-1                                                                              
      ZLIMDTHDZ=0.28E-02                                                                        
      ZKAPPA=RD/CPD                                                                              
C                                                                                                
      KTP1=KTRAC+1                                                                               
C                                                                                               
C                                                                                                

      DO 120 JK=ITOPM1,2,-1                                                                      
       DO 110 JL=1,KLON                                                                          
        ZTHETAP=PTM1(JL,JK+1)*(1000./PAPM1(JL,JK+1))**ZKAPPA                                     
        ZTHETAN=PTM1(JL,JK-1)*(1000./PAPM1(JL,JK-1))**ZKAPPA                                     
        ZDZ=(PGEOM1(JL,JK-1)-PGEOM1(JL,JK+1))*ZQG                                               
        ZDTHETA=(ZTHETAN-ZTHETAP)/ZDZ                                                            
        IF(ZDTHETA.LT.ZLIMDTHDZ) THEN                                                            
          KTROPO(JL)=JK 

C LG-
*I XTTROPO.59
C tropopauze determined after d(th)/dz
          KTRHE(JL)=KTROPO(JL)

C LG-     determining the tropopause height, being the height of the
C         layer to which the tropopause height has been assigned

          TRPHGT=PGEOM1(JL,KTRHE(JL))/9.8

C LG-     end

        END IF                                                                                   
  110  CONTINUE                                                                                 
  120 CONTINUE                                                                                   
C                                                                                                
                                                                                                

C LG-
*I XTTROPO.63

C LG- in the next statements the parameter IWHERE is assigned. This
C     is normally used for the distinction of 8 different regions in 
C     the horizontal as well as the vertical in ECHAM. This parameter is 
C     used in routines such as BUDGET.f to calculate the budgets and 
C     contributions of the different processes in this budgets for 4 equal 
C     areas of the globe (90N-30N,30N-0,0-30S,30S-90S), see the 
C     explanation below. For the actual version of the 1-D model (14-10-97), 
C     the stratosphere, the free troposhere, the PBL and the biosphere 
C     are distinguised, with IWHERE being 1 for the biosphere, 2 for the
C     PBL, 3 for the tropopshere, 4 for the stratosphere.

C LG- originally distinguished compartiments

C - - -
C  Determine localization budget bit-wise:
C  Bit 0 0:NH           1:SH
C  Bit 1 0:troposphere  1:stratosphere
C  Bit 2 0:>30o(midlat) 1:<30o (tropics)
C

C LG- only the stratosphere is distinguished here whereas
C     the other three compartments are assigned in PRECHEM.f 

C - - -
C  Determine localization budget bit-wise:
C  Bit 0 0:biosphere    1:Planetary Boundary Layer
C  Bit 1 0:troposphere  1:stratosphere

c      ISHEM=MOD(KROW+1,2)
c      ITROPIC=1
c      IF (KROW.LE.32) ITROPIC=0
c      DO 141 JK=1,KLEV
c      DO 141 JL=1,KLON
c        ISTRAT=0
c        IF (JK.LT.KTRHE(JL)) ISTRAT=1
c        IWHERE(JL,JK)=1+1*ISHEM+2*ISTRAT+4*ITROPIC
c  141 CONTINUE

C LG-

      DO 142 JK=1,KLEVELS
      DO 142 JL=1,KLON
        IPBL=0
        ITROP=0
        ISTRAT=0
        IF (JK.LE.KLEV) IPBL=1
        IF (JK.LT.KPBLHE(JL)) ITROP=1
        IF (JK.LT.KTRHE(JL)) ISTRAT=1
        IWHERE(JL,JK)=1+1*IPBL+1*ITROP+1*ISTRAT
  142 CONTINUE

C - - -
C  Determine localization budget bit-wise:
C  Bit 0 0:Biosphere    1:Planetary Boundary Layer
C  Bit 1 0:Troposphere  1:Stratosphere

C LG- defining the reference height of the half pressure level of
C     the surface layer which is used to determine the pressure levels
C     within the canopy

      DO 127 JK=1,KLEV
      DO 127 JL=1,KLON
       ZDP(JL,JK)=PAPHM1(JL,JK+1)-PAPHM1(JL,JK)
       ZT(JL,JK)=PTM1(JL,JK)
       ZQ(JL,JK)=PQM1(JL,JK)
       TVIRT=ZT(JL,JK)*(1.+0.607717*ZQ(JL,JK))
       ZRHOA(JL,JK)=PAPM1(JL,JK)*ZMAIR*1E-6/(TVIRT*ZGASC)

C LG-  Calculation of depth of layers and altitude,
C      the canopy height is accounted for by increasing the
C      model's altitude with the displacement height, which is about
C      2/3 the canopy height. 

       IF (LAGRIAN.AND.ILSTYPE.GT.1) HCCORR=HC ! mz_lg_20051212+

       DZ(JL,JK)=ZDP(JL,JK)/(ZRHOA(JL,JK)*1.E3*G)
       Z(JL,JK)=PGEOM1(JL,JK)/G+HCCORR
       ZHP(JL,JK)=Z(JL,JK)+0.5*DZ(JL,JK)

  127 CONTINUE

C LG- defining different physical parameters for the vegetation/snow layers,
C     e.g. pressure in order to determine the grid mass of the grids within the
C     canopy.

      IF (LBIOSPH.OR.LBULKVEG.OR.LVEG_MLAY.OR.LSNOW_MLAY) THEN

      DO 128 JK=KLEV+1,KLEVELS
      DO 128 JL=1,KLON  

        JJ=KLEVELSVEG+1-(JK-KLEV)

        ZT(JL,JK)=ZT(JL,KLEV)
        ZQ(JL,JK)=ZQ(JL,KLEV)

C LG-   reference height full pressure and half pressure level, these
C       are only calculated for the first timestep, then from these
C       defined reference heights, the pressure levels within the canopy
C       are calculated and then these are saved throughout the model
C       intergration to determine the reference heights 

        IF (KSTEP.EQ.0.OR.(LAGRIAN.AND.ILSTYPE.GT.1)) THEN ! mz_lg_20051212+

         IF (LDTHICK) THEN

C LG-     determining coefficients to define the height coordinates within
C         the canopy with an increasing thickness with increasing height

          CTHICK_CANLAY(JJ)=
     *     (1.-(FLOAT(KLEVELSVEG-JJ)/FLOAT(KLEVELSVEG-1)))*ALOG(HCCORR)+
     *     (FLOAT(KLEVELSVEG-JJ)/FLOAT(KLEVELSVEG-1))*ALOG(MINTHICK)

          CTHICK_CANLAY(JJ-1)=
     *     (1.-(FLOAT(KLEVELSVEG-(JJ-1))/FLOAT(KLEVELSVEG-1)))*ALOG(HCCORR)+
     *     (FLOAT(KLEVELSVEG-(JJ-1))/FLOAT(KLEVELSVEG-1))*ALOG(MINTHICK)

C LG-     vertical coordinate system within the canopy with increasing
C         thickness of the layers with increasing height, the lowest layer
C         has a thickness of MINTHICK being defined in the parameter file
C         parchem.h 

          ZHP(JL,JK)=EXP(CTHICK_CANLAY(JJ))
	  Z(JL,JK)=(EXP(CTHICK_CANLAY(JJ))+
     &       MAX(0.,EXP(CTHICK_CANLAY(JJ-1))))/2.

          IF (JK.EQ.KLEVELS) 
     &      Z(JL,JK)=ZHP(JL,JK)/2.

         ELSE

C LG-     equidistant levels within the canopy 

          ZHP(JL,JK)=(HCCORR/KLEVELSVEG)*JJ
	  Z(JL,JK)=(HCCORR/KLEVELSVEG)*JJ-(1./2.)*(HCCORR/KLEVELSVEG)

         ENDIF

C LG-    pressure difference for full pressure level

         DAPM1(JL,JK)=PAPM1(JL,KLEV)*G/(RD*ZT(JL,JK))*
     *                (Z(JL,KLEV)-Z(JL,JK))

C LG-    pressure difference for half pressure level

         DAPHM1(JL,JK)=PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK))*
     *                (ZHP(JL,KLEV)-ZHP(JL,JK))

C LG-    pressure at canopy floor

         IF (JK.EQ.KLEVELS) THEN
          ZHP(JL,KLEVELS+1)=0.
          DAPHM1(JL,KLEVELS+1)=PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK))*
     *                (ZHP(JL,KLEV)-ZHP(JL,KLEVELS+1))
         ENDIF

        ENDIF

C LG-   pressure and reference height at full pressure level

        ZAPM1(JL,JK)=PAPM1(JL,KLEV)+DAPM1(JL,JK)
        Z(JL,JK)=Z(JL,KLEV)-DAPM1(JL,JK)/
     *      (PAPM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
 
C LG-   pressure and reference height at half pressure level

        ZAPHM1(JL,JK)=PAPHM1(JL,KLEV)+DAPHM1(JL,JK)
        ZHP(JL,JK)=ZHP(JL,KLEV)-DAPHM1(JL,JK)/
     *         (PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))

C LG-   pressure at canopy floor

        IF (JK.EQ.KLEVELS) THEN
         ZAPHM1(JL,JK+1)=PAPHM1(JL,KLEV)+DAPHM1(JL,JK+1)
         ZHP(JL,JK+1)=ZHP(JL,KLEV)-DAPHM1(JL,JK+1)/
     *            (PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
        ENDIF

  128 CONTINUE
 
      DO 129 JK=KLEV+1,KLEVELS
      DO 129 JL=1,KLON
        TVIRT=ZT(JL,JK)*(1.+0.607717*ZQ(JL,JK))
        ZRHOA(JL,JK)=ZAPM1(JL,JK)*28.9644E-6/(TVIRT*8.3144)
        ZMA(JL,JK)=(ZAPHM1(JL,JK+1)-ZAPHM1(JL,JK))*ZFAC
  129 CONTINUE

C LG- end IF (LBIOSPH.OR.LBULKVEG.OR.LVEG_MLAY)

      ENDIF

      DO 130 JK=1,KLEV
      DO 130 JL=1,KLON
       TVIRT=PTM1(JL,JK)*(1.+0.607717*PQM1(JL,JK))
       ZRHOA(JL,JK)=PAPM1(JL,JK)*28.9644E-6/(PTM1(JL,JK)*8.3144)
       ZMA(JL,JK)=(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))*ZFAC
  130 CONTINUE

      DO 140 JK=1,KLEVELS                                                                        
      DO 140 JL=1,KLON                                                                         

C LG-
*I XTTROPO.67
C --- air mass contained in one gridbox (g)
       GRMASS(JL,JK)=ZMA(JL,JK)*1000.
       GRVOL(JL,JK)=GRMASS(JL,JK)/ZRHOA(JL,JK) ! [g / g cm-3] -> cm-3
  140 CONTINUE                                                                                   

C                                                                                               
      DO 170 JT=1,KTRAC                                                                         
      DO 160 JK=1,KLEVEL                                                                           
      DO 150 JL=1,KLON                                                                                 
      ZMT(JL,JK,JT)=ZMA(JL,JK)*PXTM1(JL,JK,JT)                                                   
  150 CONTINUE                                                                                  
  160 CONTINUE                                                                                   
  170 CONTINUE                                                                                   

C                                                                                                
C                                                                                               
      DO 190 JK=1,KLEV                                                                           
      DO 180 JL=1,KLON                                                                           
      IF(JK.GE.KTROPO(JL)) THEN                                                                 
       ZTFAC(JL,JK)=1.                                                                           
       ZSFAC(JL,JK)=0.                                                                           
      ELSE                                                                                       
       ZTFAC(JL,JK)=0.                                                                           
       ZSFAC(JL,JK)=1.                                                                           
      ENDIF                                                                                     
  180 CONTINUE                                                                                   
  190 CONTINUE                                                                                
C                                                                                               

C LG-
c     DO 200 JT=1,KTRAC                                                                          
c         TROPM(KROW,JT)=TROPM(KROW,JT)                                                          
c    *              +SDOT(KLON*KLEV,ZMT(1,1,JT),1,ZTFAC,1)                                      
c         STRATM(KROW,JT)=STRATM(KROW,JT)                                                         
c    *              +SDOT(KLON*KLEV,ZMT(1,1,JT),1,ZSFAC,1)                                       
c  200 CONTINUE                                                                                   
                                                                                                
c         TROPM(KROW,KTP1)=TROPM(KROW,KTP1)                                                      
c    *              +SDOT(KLON*KLEV,ZMA,1,ZTFAC,1)                                              
c         STRATM(KROW,KTP1)=STRATM(KROW,KTP1)                                                    
c    *              +SDOT(KLON*KLEV,ZMA,1,ZSFAC,1)                                              

C                                                                                               
C  -------------------------------------------------------                                      
C                                                                                               
      WRITE(NUNMDFL,'(1a,i3)')
     * ' End calculation trop. height, level from the top: ',
     *   KTRHE(1)

      RETURN                                                                                    
      END                                                                                     
