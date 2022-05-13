      SUBROUTINE CUPR
     *    (KHOR,     KHOR2,    KLON,
     *     KSTART,   KSTOP,    KLEN,
     *     KLEV,     KLEVP1,   KLEVM1,   KSTEP,     KPRINT,
     *     KCBOT,    KCTOP,    KTYPE,    KTOPM2,
     *     PTEN,     PQEN,     PQSEN,    PUEN,      PVEN,
     *     PVERV,    PGEO,     PAP,      PAPH,      PGEOH,
     *     PTENH,    PQENH,    PQSENH,   KLWMIN,
     *     PTU,      PQU,      PTD,      PQD,
     *     PUU,      PVU,      PUD,      PVD,
     *     PMFU,     PMFD,     PMFUS,    PMFDS,
     *     PMFUQ,    PMFDQ,    PDMFUP,   PDMFDP,   PRFLCK,
     *     PMFUL,    PLU,      PLUDE,    KLAB,
CEVM---    CONDENSATION RATES
     *     TEMFCD ,   QEMFCD , XEMFCD
     *   )
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE
C          -------
C
C          THIS ROUTINE PRINTS CHARACTERISTICS OF THE MASSFLUX SCHEME
C	   
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CUMASTR*.
C
C
      INCLUDE 'comcon.h'
      INCLUDE 'comio.h'
C
      REAL     PTEN(KHOR,KLEV),        PQEN(KHOR,KLEV),
     *         PUEN(KHOR2,KLEV),       PVEN(KHOR2,KLEV),
     *         PQSEN(KHOR,KLEV),       PVERV(KHOR,KLEV),
     *         PGEO(KHOR,KLEV),        PGEOH(KHOR,KLEV),
     *         PAPH(KHOR,KLEVP1),      PAP(KHOR,KLEV),
     *         PTENH(KHOR,KLEV),       PQENH(KHOR,KLEV),      
     *         PQSENH(KHOR,KLEV)     

C
      REAL     PTU(KHOR,KLEV),         PQU(KHOR,KLEV),
     *         PTD(KHOR,KLEV),         PQD(KHOR,KLEV),
     *         PUU(KHOR,KLEV),         PUD(KHOR,KLEV),
     *         PVU(KHOR,KLEV),         PVD(KHOR,KLEV),
     *         PMFU(KHOR,KLEV),        PMFD(KHOR,KLEV),
     *         PMFUS(KHOR,KLEV),       PMFDS(KHOR,KLEV),
     *         PMFUQ(KHOR,KLEV),       PMFDQ(KHOR,KLEV),
     *         PDMFUP(KHOR,KLEV),      PDMFDP(KHOR,KLEV),
     *         PMFUL (KHOR,KLEV),
     *         PLU(KHOR,KLEV),         PLUDE(KHOR,KLEV)
      INTEGER  KLAB(KHOR,KLEV),        KLWMIN(KHOR),
     * 	       KCBOT(KHOR),           
     *         KCTOP(KHOR)
C
      REAL PRFLCK(KHOR,KLEV)
C
      REAL TEMFCD(KLEV),QEMFCD(KLEV),XEMFCD(KLEV)
C
C     -------------------------------------------------------
C     DECLARATION OF WORKING ARRAYS
C
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
C
      LOGICAL
     *        LOFLAG (JPHR)
      REAL
     *        ZWMAX  (JPHR)
      COMMON/COMPR4/ISUMTYPE(0:3),ZBOTCL(3),ZTOPCL(3)
C--------------------------------------------------------
C     DECLARE SOME LOCAL FIELDS
C
      REAL  TVIRT    (MLEV),  TVIRTH    (MLEV),
     *	    TUVIRT   (MLEV),  TUVIRTH   (MLEV),
     *      THETAU   (MLEV),  THETAUH   (MLEV),
     *	    THETA    (MLEV),  THETAH    (MLEV),
     *      TURBQTEND(MLEV),  TURBTTEND (MLEV),
     *      CONDTTEND(MLEV),  CONDQTEND (MLEV),	
     *	    TURBTTOT (MLEV),  TURBQTOT  (MLEV),
     *      EXNER    (MLEV)
C
C
      DO JK=1,KLEV
        TVIRT    (JK) = 0. 
	TVIRTH   (JK) = 0.
        TUVIRT   (JK) = 0.
	TUVIRTH  (JK) = 0.
        THETAU   (JK) = 0.  
	THETAUH  (JK) = 0.
        THETA    (JK) = 0.
	THETAH   (JK) = 0.
        TURBQTEND(JK) = 0.
	TURBTTEND(JK) = 0.
        CONDTTEND(JK) = 0.	
	CONDQTEND(JK) = 0.	
        TURBTTOT (JK) = 0.	
	TURBQTOT (JK) = 0.
	EXNER    (JK) = 0.	
      ENDDO      
      DAYL    = 86400
      CONVTMP = DAYL
      CONVHUM = DAYL*1000.

C------------------------------------------------------------------
C  1.0   DETERMINE SOME FIELDS AT FULL AND HALF LEVELS
C------------------------------------------------------------------

      P0 = PAPH(1,KLEVP1)
      DO JK=2,KLEV
        EXNER  (JK) = (PAPH(1,JK)/P0)**(RD/CPD)
	THETAH (JK) = (PTENH(1,JK) + PGEOH(1,JK)/CPD)
	THETAUH(JK) = (PTU(1,JK)   + PGEOH(1,JK)/CPD)
      ENDDO
      DO JK=2,KLEV
         EXNER( JK) = (PAP(1,JK)/P0)**(RD/CPD)
         THETAU(JK) = (THETAUH(JK) + THETAUH(JK+1))/2.  
         THETA (JK) = (THETAH (JK) + THETAH (JK+1))/2.
      ENDDO
C
      IF (KSTEP.EQ.0) THEN
        WRITE (NUNCUM1,'(A/A)')
     $  '#----  H A L F   L E V E L   F I E L D S ---- ','#'
         WRITE (NUNCUM1,'(A)')
     $   '# HALF-LEVEL FIELDS GENERATED IN MASS FLUX SCHEME'
      END IF
C
      IF (KTYPE .NE. 0 ) THEN
C
      WRITE (NUNCUM1,'(A/A,I8,A)')
     $    '#'
     $   ,'#',KSTEP,' --- TIMESTEP ---'
      WRITE (NUNCUM1,'(2A/2A/2A/2A/2A)')
     1 '#-----------------------------------------------------------'
     1,'------------------'
     2,'# LV   HGHT  PRES       S/CP            TEMP           Q_VAP   '
     2,'Q_LIQ  Q_SAT'
     3,'#                    AV      CUM     AV      CUM     AV    CUM'
     3,'   CUM    AV'                
     4,'#     (M)   (MB)   (---------- KELVIN ----------)(----------' 
     4,' G/KG ----------)'
     1,'#-----------------------------------------------------------'
     1,'------------------'
C
      DO JLEV=2,KLEV
       IF (PMFU(KHOR,JLEV).EQ.0.0) THEN
	 THETAUH  (JLEV) = 0.0
         PTU (KHOR,JLEV) = 0.0
         PQU (KHOR,JLEV) = 0.0 
         PLU (KHOR,JLEV) = 0.0
       ENDIF
          WRITE (NUNCUM1,'(F4.1,F7.0,F6.1,4F8.2,4F6.2)')
     $             FLOAT  (JLEV)-0.5
     $            ,PGEOH  (KHOR,JLEV)/G
     $            ,PAPH   (KHOR,JLEV)/100.
     $            ,THETAH (JLEV)
     $            ,THETAUH(JLEV)
     $            ,PTENH  (KHOR,JLEV)
     $            ,PTU    (KHOR,JLEV)
     $            ,PQENH  (KHOR,JLEV)*1000.0
     $            ,PQU    (KHOR,JLEV)*1000.0
     $            ,PLU    (KHOR,JLEV)*1000.0
     $            ,PQSENH (KHOR,JLEV)*1000.0
       ENDDO
C      
       WRITE (NUNCUM1,'(A/A/A)')
     $   '#--------------------------------------------------------'
     $     ,'#CL_BASE: LEV   HGHT   PRES  CL_TOP: LEV   HGHT   PRES'
     $     ,'#               (M)    (MB)                 (M)   (MB)'
       WRITE (NUNCUM1,'(A,F4.1,F8.0,F7.1,A,F4.1,F8.0,F7.1)')
     $      '        '
     $            ,FLOAT(KCBOT(KHOR))-0.5 
     $            ,PGEOH(KHOR,KCBOT(KHOR))/G
     $            ,PAPH (KHOR,KCBOT(KHOR))/100.
     $     ,'        '
     $            ,FLOAT(KTOPM2)-0.5 
     $            ,PGEOH(KHOR,KTOPM2)/G
     $            ,PAPH (KHOR,KTOPM2)/100.
       WRITE (NUNCUM1,'(A)')
     $   '#--------------------------------------------------------'
C
C     ENDIF (KTYPE)
      ENDIF
C
C--------------------------------------------------------------------
C    2.0  PRINT TURBULENT FLUXES AND CONDENSATION RATES
C--------------------------------------------------------------------

      IF (KSTEP.EQ.0) THEN
        WRITE (NUNCUM2,'(A/A)')
     $ '#-- F L U X E S  A N D   C O N D E N S A T I O N   R A T E S --'
     $,'#'
      END IF
C
      IF (KTYPE .NE. 0) THEN
C
        WRITE (NUNCUM2,'(A,I8,A)')   
     $   '#',KSTEP,' --- TIMESTEP --- '
        IF (.NOT.LDDRAF) THEN
          WRITE (NUNCUM2,'(A/A)') 
     $    '#         ----   N O   D O W N D R A F T  ---- '   
     $   ,'#         TOTAL MASS FLUX =   UPWARD  MASS FLUX' 
          WRITE (NUNCUM2,'(A20,I4)')
     &    '#TYPE OF CONVECTION =',KTYPE
          WRITE (NUNCUM2,'(A/A/A/2A/2A)')
     $     '#-------------------------------------------------------'
     $    ,'#  ---     U P D R A F T   M A S S - F L U X     ---    '
     $    ,'#                                                       '
     $    ,'#LEV   HGHT  ILAB    M_UP    M*S_UP    M*Q_UP   M*L_UP  '
     $    ,'HGHT  G_p-e_p  '
     $    ,'#      (M)       (KG/M^2*S)  (------ W/M^2*S ------)    '
     $    ,'(M)  (W/M^2*S)'
C
          WRITE (NUNCUM2,'(F4.1,F7.0,I5,F10.4,3F9.2,F7.0,F9.2)')
     $             (FLOAT (JLEV)-0.5 
     $             ,PGEOH (KHOR,JLEV)/G
     $             ,KLAB  (KHOR,JLEV)
     $             ,PMFU  (KHOR,JLEV)
     $             ,PMFUS (KHOR,JLEV)
     $             ,PMFUQ (KHOR,JLEV)*ALV
     $             ,PMFUL (KHOR,JLEV)*ALV
     $             ,PGEO  (KHOR,JLEV)/G
     $		   ,PDMFUP(KHOR,JLEV)*ALV
     $             ,JLEV=1,KLEV)
C
          WRITE (NUNCUM2,'(A)')
     $     '#---------------------------------------------------------'
        ELSE
          WRITE (NUNCUM2,'(A/A/A/A/A)')
     $     '#-------------------------------------------------------'
     $    ,'#--  D O W N   D R A F T   M A S S - F L U X     ---    '
     $    ,'#                                                       '
     $    ,'#LEV  HGHT  ILAB   M_DOWN  M*S_DOWN  M*Q_DOWN M*L_DOWN'
     $    ,'#      (M)       (KG/M^2*S)    (------ W/M^2*S ------)'
C
          WRITE (NUNCUM2,'(F4.1,F8.0,I5,F10.3,2F9.2)')
     $               (FLOAT (JLEV)-0.5 
     $               ,PGEOH (KHOR,JLEV)/G
     $               ,KLAB  (KHOR,JLEV)
     $               ,PMFD  (KHOR,JLEV)
     $               ,PMFDS (KHOR,JLEV)
     $               ,PMFDQ (KHOR,JLEV)*ALV
     $               ,JLEV=KCTOP(KHOR)-2,KLEV)
          WRITE (NUNCUM2,'(A/A/A/A/A/A)')
     $     '#-------------------------------------------------------'
     $    ,'#      ---  T OT A L    M A S S - F L U X     ---    '
     $    ,'#                                                       '
     $    ,'#LEV  HGHT  ILAB    M_FLX   M*S_FLX   M*Q_FLX  M*L_FLX '
     $    ,'#      (M)       (KG/M^2*S)    (------ W/M^2*S ------)'
C
          WRITE (NUNCUM2,'(F4.1,F8.0,F10.3,3F9.2)')
     $               (FLOAT (JLEV)-0.5 
     $               ,PGEOH (KHOR,JLEV)/G
     $               ,KLAB  (KHOR,JLEV)
     $               ,PMFD  (KHOR,JLEV) + PMFU  (KHOR,JLEV)
     $               ,PMFDS (KHOR,JLEV) + PMFUS (KHOR,JLEV)
     $               ,(PMFDQ (KHOR,JLEV) + PMFUQ (KHOR,JLEV))*ALV
     $               ,PMFUL (KHOR,JLEV)*ALV
     $               ,JLEV=KCTOP(KHOR)-2,KLEV)
C       ENDIF (.NOT.LDDRAF)
        ENDIF
C
C     ENDIF (KTYPE)
      ENDIF
C
      DO 250 JK=2,KLEV
C
      IF(JK.LT.KLEV) THEN
         DO 220 JL=KSTART,KSTOP
            LLO1=(PTEN(JL,JK)-TMELT).GT.0.
            ZALV=CVMGT(ALV,ALS,LLO1)
            TURBTTEND(JK) = (G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     1                     (PMFUS(JL,JK+1)-PMFUS(JL,JK)+
     1                      PMFDS(JL,JK+1)-PMFDS(JL,JK))
            CONDTTEND(JK) =  (G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     1                      (-ZALV*(PMFUL(JL,JK+1)-PMFUL(JL,JK)-
Cacp  Volgende regel ontbrak nog:
     +                       PLUDE(JL,JK)-
     1                      (PDMFUP(JL,JK)+PDMFDP(JL,JK))))
            TURBQTEND(JK) = (G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1                       (PMFUQ(JL,JK+1)-PMFUQ(JL,JK)+
     1                        PMFDQ(JL,JK+1)-PMFDQ(JL,JK))
            CONDQTEND(JK) = (G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1                       (PMFUL(JL,JK+1)-PMFUL(JL,JK)-
Cacp  Volgende regel ontbrak nog:
     +                        PLUDE(JL,JK)-
     1                       (PDMFUP(JL,JK)+PDMFDP(JL,JK)))
  220 CONTINUE
      ELSE
      IF (JK.EQ.KLEV) THEN
         DO 230 JL=KSTART,KSTOP
            LLO1=(PTEN(JL,JK)-TMELT).GT.0.
            ZALV=CVMGT(ALV,ALS,LLO1)
            TURBTTEND(JK) = -(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     1                     (PMFUS(JL,JK)+PMFDS(JL,JK))
            CONDTTEND(JK) =  (G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     1                      ZALV*(PMFUL(JL,JK)
Cacp  Volgende regel ontbrak nog:
     +                       +PLUDE(JL,JK)
     1                      +PDMFUP(JL,JK)+PDMFDP(JL,JK))
            TURBQTEND(JK) = -(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1                       (PMFUQ(JL,JK)+PMFDQ(JL,JK))
            CONDQTEND(JK) = -(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     1                        (PMFUL(JL,JK)
Cacp  Volgende regel ontbrak nog:
     +                       +PLUDE(JL,JK)
     1			      +PDMFUP(JL,JK)+PDMFDP(JL,JK))
  230    CONTINUE
      END IF
      END IF

  250 CONTINUE
C
      IF (KTYPE .NE. 0) THEN
C
        WRITE (NUNCUM3,'(A/A/A/A/A/A/A)')
     $   '#-------------------------------------------------'
     $  ,'#--- C U M U L U S    T E N D E N C I E S  -------'
     $  ,'#                                                -'
     $  ,'#               |  TEMPERATURE  |   SPEC. HUM.   -'
     $  ,'#LEV   HEIGHT   | TURB     COND | TURB      COND -'
     $  ,'#       (M)     |-----[K/DAY]--------[G/KG/DAY]---'
     $  ,'#-------------------------------------------------'
        WRITE (NUNCUM3,'(A/A,I8,A/A)')
     $  '#','#',KSTEP,' --- TIMESTEP ---','#'
        WRITE (NUNCUM3,'(I4,F8.0,4F9.3)')
     $             (JLEV 
     $             ,PGEO (KHOR,JLEV)/G
     $             ,TURBTTEND(JLEV)*CONVTMP
     $             ,CONDTTEND(JLEV)*CONVTMP
     $             ,TURBQTEND(JLEV)*CONVHUM
     $             ,CONDQTEND(JLEV)*CONVHUM
     $             ,JLEV=1,KLEV)
C
C     ENDIF (KTYPE)
      ENDIF
C
CEVM
C     MOVE TENDENCIES DUE TO CONDENSATION INTO xxMFCD FIELDS
C     SUBTRACT xxMFCD TENDENCIES FROM xxMFTB TENDENCIES (IN *PHECHAM*)
C     IN THIS WAY TURBULENCE AND CONDENSATION COMPONENTS WITHIN
C     THE MASS-FLUX SCHEME CAN BE DISCRIMINATED IN *WRPRFL*
C
      DO 310 JK=1,KLEV
	TEMFCD(JK) =  CONDTTEND(JK)
	QEMFCD(JK) =  CONDQTEND(JK)
	XEMFCD(JK) = -CONDQTEND(JK)
 310  CONTINUE
C
C       -----------------------------------------
C       WRITE STATISTICS AT MOD(KSTEP,KPRINT) = 0
C
C       RESET ZBOTCL,ZTOPCL & ISUMTYPE IF MOD(KSTEP-1,KPRINT) = 0
C
	IF (KSTEP.EQ.0 .OR. MOD(KSTEP-1,KPRINT).EQ.0) THEN
	  DO JT=0,3
	    ZBOTCL(JT) = 0.
	    ZTOPCL(JT) = 0.
	  ENDDO
	  DO JT=0,3
	    ISUMTYPE(JT) = 0
	  ENDDO
	ENDIF
C
C       COUNT NUMBER OF OCCURRENCES OF KTYPE-CONVECTION
C
        IF (KTYPE .GE. 0 .AND. KTYPE .LE. 3) THEN
	  ISUMTYPE(KTYPE) = ISUMTYPE(KTYPE) + 1
	ENDIF
C
C       DO CONDITIONAL SAMPLING OF BOTTOM AND TOP OF CONVECTION
C
        IF (KTYPE .GT. 0 .AND. KTYPE .LE. 3) THEN
	  DO JK=KLEV,2,-1
	    IF (KLAB(KHOR,JK+1).EQ.1. AND. KLAB(KHOR,JK).EQ.2) THEN
	      ZBOTCL(KTYPE) = ZBOTCL(KTYPE) + FLOAT(JK)
            ENDIF
	    IF (KLAB(KHOR,JK+1).EQ.2. AND. KLAB(KHOR,JK).EQ.0) THEN
	      ZTOPCL(KTYPE) = ZTOPCL(KTYPE) + FLOAT(JK+1)
            ENDIF
	  ENDDO
	ENDIF
C
C       PRINT OUT
C
        IF (KSTEP .EQ. 0) THEN
          WRITE (NUNCUM4,'(A/A/2A/2A/2A/2A/2A)')
     $     '#------------------------------------------------'
     $    ,'#--- C O N V E C T I O N    S T A T I S T I C S -'
     $    ,'#-----------------------------'
     $    ,'----------------------------------'
     $    ,'#TIME|     NUMBER OF EVENTS  |'
     $    ,'    DEEP  |   SHALLOW| MID-LEVEL -'
     $    ,'#STEP|   NO  DEEP SHLLW MIDLV|'
     $    ,'  BOT  TOP|  BOT  TOP|  BOT  TOP -'
     $    ,'#                             '
     $    ,'   UNITS ARE CENTI-MODEL_LEVELS  -'
     $    ,'#-----------------------------'
     $    ,'----------------------------------'
	ENDIF
        IF (KSTEP .NE. 0 . AND. MOD(KSTEP,KPRINT) .EQ. 0) THEN
	  DO JT=1,3
	    IF ( ISUMTYPE(JT) .GT. 0) THEN
	      ZBOTCL(JT) = ZBOTCL(JT) / ISUMTYPE(JT)
	      ZTOPCL(JT) = ZTOPCL(JT) / ISUMTYPE(JT)
            ELSE
	      ZBOTCL(JT) = -0.01
	      ZTOPCL(JT) = -0.01
	    ENDIF
	  ENDDO
	  WRITE (NUNCUM4,'(I5,4I6,1X,2I5,1X,2I5,1X,2I5)')
     $         KSTEP
     $        ,(ISUMTYPE(JT),JT=0,3)
     $        ,(NINT(ZBOTCL(JT)*100.)
     $         ,NINT(ZTOPCL(JT)*100.),JT=1,3)
	ENDIF
C
       RETURN
      END
