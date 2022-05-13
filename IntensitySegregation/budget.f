      SUBROUTINE BUDGET (NSTOP,DTIME,NSTEP,NLEV,NLEVEL,PM,PMLOC,
     *                   ZRHOA,ZP,Z,ZDZ,KG3X,IS,IBS,ISHOR,ISBULK)

      IMPLICIT NONE

C LG- 
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

C LG- more declarations, IS is the frequency of the output of
C     profiles and time series whereas IBS is the frequency of the
C     the output of the budgets of the tracers

      INTEGER KG3X,JT,JK,NLEV,JL,IW,NSTEP,IS,IBS,NRESUM,NSTOP,
     *        NLEV1,NLEV2,NLEVPR,NLEVEL,IDAY,JJ,I

      REAL PM(NLON,NLEVT,NTRAC),PMLOC(NLON,NLEVT,KG3X),
     *     PMLOCBULK(KG3X),
     *   PMTB(NDIM,NTRAC),PMLOCTB(NDIM,NG3X),ATVOL(NDIM),ATWEI(NDIM),
     *   ZRHOA(NLON,NLEVT),ZP(NLON,NLEVT),Z(NLON,NLEVT),ZDZ(NLON,NLEVT)

C ACP-
      REAL ISHOR(NBINREAC,NLEVT),ISBULK(NBINREAC)
C ACP-end
     
C LG- interpolation parameters
   
      REAL ZPM(NLON,NLEVT,NTRAC),DH,DZ
      REAL ZPMBULK(NTRAC)
     
      SAVE PMTB,PMLOCTB,ATVOL,ATWEI

      REAL DTIME

      CHARACTER*70 FNAME
      
C LG-

      WRITE(NUNMDFL,*)'Start of BUDGET.f',NSTEP

C LG- the budgets of the tracers and the contribution of the different
C     processes in the total budgets are written away to an output file
C     which is afterwards being manipulated (see budget.f) to yield
C     the budgets in readable format. The budgets are written away with
C     intervals defined by the parameter IBS. The budgets are for a column
C     with horizontal dimensions resembling the horizontal resolution for
C     the specific latitude/longitude in the ECHAM T30 model (~ 300 * 300 km).  
C     In order to specify the budgets for a column with a surface area of 1 m2
C     a correction should be applied.  
 
      IF (NSTEP.EQ.0)
     *  OPEN(UNIT=NUNBUD,FILE='budget.out',STATUS='UNKNOWN')

      IF (MOD(NSTEP,IBS).EQ.0.OR.NSTEP.EQ.NSTOP) THEN

       CALL RESETR (PMTB,NDIM*NTRAC,0.)
       CALL RESETR (PMLOCTB,NDIM*KG3X,0.)
       CALL RESETR (ATVOL,NDIM,0.)
       CALL RESETR (ATWEI,NDIM,0.)

       DO 100 JT=1,NTRAC
       DO 100 JK=1,NLEVEL
       DO 100 JL=1,NLON
         IW=IWHERE(JL,JK)
         PMTB(IW,JT)=PMTB(IW,JT)
     *      +PM(JL,JK,JT)*GRVOL(JL,JK)
  100  CONTINUE
       DO 101 JT=1,KG3X
       DO 101 JK=1,NLEVEL
       DO 101 JL=1,NLON
         IW=IWHERE(JL,JK)
         PMLOCTB(IW,JT)=PMLOCTB(IW,JT)
     *      +PMLOC(JL,JK,JT)*GRVOL(JL,JK)
  101  CONTINUE
       DO 102 JK=1,NLEVEL
       DO 102 JL=1,NLON
         IW=IWHERE(JL,JK)
         ATWEI(IW)=ATWEI(IW)+GRMASS(JL,JK)
         ATVOL(IW)=ATVOL(IW)+GRVOL(JL,JK)
  102  CONTINUE

C LG-  writing the month, the total number of timesteps and the timestep

       IF (NSTEP.EQ.NRESUM) THEN
        WRITE (NUNBUD,*) IMON,JDAY,LTIME,DTIME,NSTOP

C LG-  writing of location of the selected grid square/column

        WRITE (NUNBUD,*) DLAT,DLON 
       ENDIF
       
       WRITE (NUNBUD,*) NSTEP,IBS
       WRITE (NUNBUD,*) ATWEI,ATVOL
       WRITE (NUNBUD,*) PMTB,BXT,PMLOCTB,BG3,PMTBM,PMLOCTBM,BRX

       IF (NSTEP+1.NE.NRESUM) THEN
	CALL RESETR (BXT,NDIM*NFIELD*NTRAC,0.)
	CALL RESETR (BG3,NDIM*NFIELD*NG3X,0.)
	CALL RESETR (BRX,NDIM*NREAC,0.)
	CALL RESETR (PMTBM,NDIM*NTRAC,0.)
	CALL RESETR (PMLOCTBM,NDIM*NG3X,0.)
       ENDIF

      ENDIF

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNBUD)

      IF (MOD(NSTEP,IS).EQ.0) THEN

C LG-  writing of data for plotting time series for different levels
C      One should be careful with the MAX statements used in the write
C      statements since this removes negative data which might be 
C      introduced due to errors in the model. However, since the plotting
C      program IDL can not deal with very small values (about 1e-50) this
C      MAX statement is required here!!!!!!!!

       NLEV1=1
       NLEV2=NLEVEL
       IF (LAGRIAN.AND.ILSTYPE.GT.1) NLEV2=NLEVT
       NLEVPR=NLEV2-NLEV1+1

       IF (NSTEP.EQ.0) THEN
        FNAME='output/1d_chem.out'
        OPEN(UNIT=NUN1DCHEM,FILE=FNAME,STATUS='UNKNOWN')

        WRITE(NUN1DCHEM,*) NSTOP,NLEVPR,IS,DTIME
        WRITE(NUN1DCHEM,*) NLEV1,NLEV2 

       ENDIF

       DO 300 JL=1,NLON

       DO 899 JK=1,NLEVEL
        JJ=NLEVEL+1-JK

C LG-   calculation of actual height

        IF (JK.LE.NLEV) THEN
	 HGHT(JK)=Z(JL,JK)
	ELSE
         IF (NLEVVEG.GT.0) THEN
          HGHT(JK)=(HCCORR/NLEVVEG)*JJ-0.5*(HCCORR/NLEVVEG)
         ELSE
          HGHT(JK)=0.
         ENDIF
	ENDIF

  899  CONTINUE

       IF (NSTEP.EQ.0) THEN
         IF (LAGRIAN.AND.LBIOSPH) THEN
	  WRITE (NUN1DCHEM,*) (Z(JL,JK),JK=NLEV2,NLEV1,-1)
	 ELSE 
	  WRITE (NUN1DCHEM,*) (HGHT(JK),JK=NLEV2,NLEV1,-1)
         ENDIF
         WRITE (NUN1DCHEM,*) (ZP(JL,JK),JK=NLEV2,NLEV1,-1)
       ENDIF

       WRITE(NUN1DCHEM,*)NSTEP,JDAY,GMT,LTIME,HC,TRPHGT,PBLHGT

C LG-  for the an experiment with an actual canopy height less than
C      the default maximum canopy height for forest (30 meters), the 
C      concentration fields are copied in such a way that the concentration 
C      fields, written to the files used for graphical interpretation, are  
C      representative for a reference height relatively to a zero height
C      of the canopy floor, bare soil, ocean surface. So for a actual canopy
C      height of 15 m, resembling three layers for a 6-layer canopy
C      representation, the concentration at the reference height of 64 m is 
C      corrected with 15 meters (HCMAX-HC)

       DO 999 JT=1,NTRAC
       DO 999 JK=1,NLEVT
         ZPM(JL,JK,JT)=PM(JL,JK,JT)
         IF (LBIOSPH) THEN
           DH=(NLEVV-NLEVVEG)*(HCMAX/NLEVV)
           IF (JK.LT.NLEVT) THEN
             DZ=Z(JL,JK)-Z(JL,JK+1)
	   ELSE
	     DZ=Z(JL,NLEVT)
	   ENDIF
	   ZPM(JL,JK,JT)=PM(JL,JK,JT)*MAX(0.,((DZ-DH)/DZ))+
     *      PM(JL,MIN(NLEVT,JK+1),JT)*MIN(1.,(DH/DZ))
           IF (DH/DZ.GT.1) 
     *     ZPM(JL,JK,JT)=PM(JL,MIN(NLEV,JK-NLEVVEG),JT)

C LG-      not for methane

           IF (JT.EQ.3) ZPM(JL,JK,JT)=PM(JL,JK,JT)

         ENDIF

  999  CONTINUE

C LG-  writing the tracer concentrations

       DO 301 JT=1,NTRAC
        IF (JT.NE.30.AND.JT.NE.31) THEN
         WRITE (NUN1DCHEM,5005) NSTEP,JT,
     *     (MAX(1.E-10,ZPM(JL,JK,JT)/ZRHOA(JL,JK)),
     *       JK=NLEV2,NLEV1,-1)
	 ZPMBULK(JT) = 0.
	 DO JK=KPBLHE(1),NLEV
	   ZPMBULK(JT) = ZPMBULK(JT) + ZDZ(JL,JK)*ZPM(JL,JK,JT)/ZRHOA(JL,JK)
         ENDDO
	 ZPMBULK(JT) = ZPMBULK(JT)/(Z(JL,KPBLHE(1))+0.5*ZDZ(JL,KPBLHE(1)))
        ELSE
          WRITE (NUN1DCHEM,5005) NSTEP,JT,(MAX(1.E-10,ZPM(JL,JK,JT)),
     *       JK=NLEV2,NLEV1,-1) 
	 ZPMBULK(JT) = 0.
	 DO JK=KPBLHE(1),NLEV
	   ZPMBULK(JT) = ZPMBULK(JT) + ZDZ(JL,JK)*ZPM(JL,JK,JT)
         ENDDO
	 ZPMBULK(JT) = ZPMBULK(JT)/(Z(JL,KPBLHE(1))+0.5*ZDZ(JL,KPBLHE(1)))
         ENDIF 
  301  CONTINUE

       DO 302 JT=1,KG3X
         IF (JT.NE.1.AND.JT.NE.2) THEN 
          WRITE (NUN1DCHEM,5005) NSTEP,JT,
     *     (MAX(1.E-10,PMLOC(JL,JK,JT)/ZRHOA(JL,JK)),
     *       JK=NLEV2,NLEV1,-1)
	  PMLOCBULK(JT) = 0.
	  DO JK=KPBLHE(1),NLEV
	    PMLOCBULK(JT) = PMLOCBULK(JT) + ZDZ(JL,JK)*PMLOC(JL,JK,JT)/
     *	    ZRHOA(JL,JK)
          ENDDO
	  PMLOCBULK(JT) = PMLOCBULK(JT)/(Z(JL,KPBLHE(1))+0.5*ZDZ(JL,KPBLHE(1)))
         ELSE
          WRITE (NUN1DCHEM,5005) NSTEP,JT,(MAX(1.E-10,PMLOC(JL,JK,JT)),
     *       JK=NLEV2,NLEV1,-1) 
	  PMLOCBULK(JT) = 0.
	  DO JK=KPBLHE(1),NLEV
	    PMLOCBULK(JT) = PMLOCBULK(JT) + ZDZ(JL,JK)*PMLOC(JL,JK,JT)
          ENDDO
	  PMLOCBULK(JT) = PMLOCBULK(JT)/(Z(JL,KPBLHE(1))+0.5*ZDZ(JL,KPBLHE(1)))
         ENDIF        
  302  CONTINUE

  300  CONTINUE 

      ENDIF

      IF (NSTEP.EQ.NSTOP) CLOSE(NUN1DCHEM)

 5005 FORMAT (2I5,24(1pE10.2))

C ACP-output of bulk quantities and profiles of intensities of segregation
        IF (NSTEP.EQ.0) THEN
          FNAME='output/bulk.out'
          OPEN(UNIT=NUNBULK,FILE=FNAME,STATUS='UNKNOWN')

          WRITE(NUNBULK,*) NSTOP,IS,DTIME

        ENDIF
        WRITE(NUNBULK,*)NSTEP,JDAY,GMT,LTIME,HC,TRPHGT,PBLHGT

        DO 305 JT=1,NTRAC
        IF (JT.NE.30.AND.JT.NE.31) THEN
         WRITE (NUNBULK,5005) NSTEP,JT,
     *     MAX(1.E-10,ZPMBULK(JT))
         ELSE
          WRITE (NUNBULK,5005) NSTEP,JT,MAX(1.E-10,ZPMBULK(JT))
         ENDIF 
  305  CONTINUE

       DO 306 JT=1,KG3X
         IF (JT.NE.1.AND.JT.NE.2) THEN 
          WRITE (NUNBULK,5005) NSTEP,JT+NTRAC,
     *     MAX(1.E-10,PMLOCBULK(JT))
         ELSE
          WRITE (NUNBULK,5005) NSTEP,JT+NTRAC,
     *	    MAX(1.E-10,PMLOCBULK(JT))
         ENDIF        
  306  CONTINUE

          DO 304 I=1,NBINREAC
            WRITE (NUNBULK,5005) NSTEP,I,
     *         SIGN(1.,ISBULK(I))*MAX(1.E-10,ABS(ISBULK(I)))
  304   CONTINUE

        IF (NSTEP.EQ.NSTOP) CLOSE(NUNBULK)

        IF (NSTEP.EQ.0) THEN
          FNAME='output/1d_ishor.out'
          OPEN(UNIT=NUN1DISHOR,FILE=FNAME,STATUS='UNKNOWN')

          WRITE(NUN1DISHOR,*) NSTOP,NLEVPR,IS,DTIME
          WRITE(NUN1DISHOR,*) NLEV1,NLEV2 

        ENDIF

        IF (NSTEP.EQ.0) THEN
          IF (LAGRIAN.AND.LBIOSPH) THEN
	    WRITE (NUN1DISHOR,*) (Z(JL,JK),JK=NLEV2,NLEV1,-1)
  	  ELSE 
	    WRITE (NUN1DISHOR,*) (HGHT(JK),JK=NLEV2,NLEV1,-1)
          ENDIF
          WRITE (NUN1DISHOR,*) (ZP(JL,JK),JK=NLEV2,NLEV1,-1)
        ENDIF

        WRITE(NUN1DISHOR,*)NSTEP,JDAY,GMT,LTIME,HC,TRPHGT,PBLHGT
          DO 303 I=1,NBINREAC
            WRITE (NUN1DISHOR,5005) NSTEP,I,
     *         (SIGN(1.,ISHOR(I,JK))*MAX(1.E-10,ABS(ISHOR(I,JK))),
     *           JK=NLEV2,NLEV1,-1)
  303   CONTINUE

        IF (NSTEP.EQ.NSTOP) CLOSE(NUN1DISHOR)
C ACP-end
      WRITE(NUNMDFL,*)'End BUDGET.f'

      RETURN
      END
