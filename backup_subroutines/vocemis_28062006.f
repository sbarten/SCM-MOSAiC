      SUBROUTINE VOCEMIS(NSTEP,NSTOP,NPRINT,DTIME,EMISFACT,ISOPEM_INT,
     &                   MTERPEM_INT,OVOCEM_INT)

c ----------------------------------------------------------------------
c     This program calculates the emission of volatile organic compounds 
c     from vegetation as a function of: biome (Leaf Area Index) and
c     temperature and Photosynthetically Active Radiation (PAR). The model 
c     considers the extinction of PAR within the canopy as a function
c     of the Leaf Area Index which is derived from the Olson ecosystems
c     database (1992), which discerns 72 ecosystems and 
c     their characteristics (see CDROM and paper by Guenther et al., 1995).
c     This Olson database is also applied to distinguish between different
c     biomes which show distinct different standard emission factors
c     (the emission rate taken at a standard temperature and for a standard
c     amount of PAR, e.g. 30 degrees C and 1000 umol m-2 s-1). The
c     units of the GEIA database which contains monthly average
c     emission fluxes is mg C m-2 month-1 and in order to compare
c     this model with these data the same units are applied where 
c     possible. This model version applies the model of Weiss and Norman 
c     (1985) to calculate the extinction of PAR as a function of the the 
c     Leaf Area Index, the distribution of the LAI (Leaf Area Density), 
c     the fraction of leaves and the orientation of these leaves. This in
c     contrast with the original model used for the GEIA emission inventory
c     which applies the formulas by Norman, 1982. 
c
c     Laurens Ganzeveld 1997
c ----------------------------------------------------------------------
c 

      IMPLICIT NONE

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

      CHARACTER HEAD*80
      CHARACTER*70 FNAME,FNAME1
      CHARACTER*50 FNAME2,FNAME3
      CHARACTER*40 dir
      CHARACTER*2 MON(12)    
      CHARACTER*15 VEG(7)
      CHARACTER*15 VEGTYPE(3)
      CHARACTER*70 CHDUM

      LOGICAL LSUNRISE

      INTEGER NSTEP, NPRINT, OLSCLASS, MAXCLASS, VEGETATION

      PARAMETER (OLSCLASS=59,MAXCLASS=73)

      INTEGER IHEAD,IIMONTH,IMONTH,NSTOP,I,J,IC,
     &        IS,IM,K,II,JJ,JK,JT,IP,N,ISPEC,IIP,
     &        NLEVEL,NLEV1,NLEV2,NLEVPR,NSTEP_SUNRISE

      INTEGER OLSNUMB(0:OLSCLASS)

      REAL FLUXSHADE(NLEVV,NSPEC_EMISVOC),FLUXSUN(NLEVV,NSPEC_EMISVOC),
     &     FLUXSUM,CT,CT1,CT2,TS,TM,RECALC,R,CLSUN,CLSHADE,ALPHA,BETA,
     &     CL1,DTIME,EL,FLUXTOT,FOLDENSLAY(NLEVV),VOCOLSINP(0:OLSCLASS,12),
     &     VOCPARAM(0:MAXCLASS,12),RES,PARDIR,PARDIF(NLEVV),
     &     EMISFACT(NSPEC_EMISVOC),ISOPEM_INT,MTERPEM_INT,
     &     OVOCEM_INT,DM_FSL,HEIGHT(NLEVT),PARSL,RBVD_OLD,
     &     DTIME_SPINUP_EMIS,FEMIS_SUNRISE

      PARAMETER (DTIME_SPINUP_EMIS=3600.)

      DATA MON /'01','02','03','04','05','06','07','08','09',
     &  '10','11','12'/

      DATA VEG /'Decid. forest','Conifer. forest','Mixed forest',
     &          'Wetland forest','Scrub woods','Mixed woods/crop',
     &          'other'/
c
c  -- assigning of values of the used constants, (see Guenther et al., 1993, JGR)
c     The term RECALC is the recalculation factor for getting the net short
c     wave radiation/PAR im umol m-2 s-1 instead of W m-2. This term is taken
c     as the average of the recalculation factor for clear sky (4.24) and
c     diffuse conditions (4.57). See Ecological Physics by J. Hage, D96-6, IMAU
c     and the official reference is (likely): Grace, J., Plant-Atmosphere
c     relationships, Chapman & Hall
c
      DATA ALPHA/0.0027/,CL1/1.066/,CT1/95000./,
     &     CT2/230000./,TS/303./,TM/314./,R/8.314/, 
     &     EL/0.42/,BETA/0.09/,RECALC/4.405/

c  -- resolution of data

      RES=0.5

  100 FORMAT (A80)
  101 FORMAT (10I8)
  102 FORMAT (4(E12.4,1X))
  104 FORMAT (4(I8,1X))
      
      IF (NSTEP.EQ.0) THEN 

C WP     adding the directory for reading vegetation input

	 DIR='../../echam/'
	 IIP=index(dir,' ')
C
	 OPEN (2,FILE=
     &     dir(1:iip-1)//'input/veg/VOC_Olson_92.dat'
     &     ,STATUS='OLD')

	 DO J=1,8
           READ (2,'(70A)') CHDUM
	 ENDDO

	 DO I=1,OLSCLASS
           READ (2,*) OLSNUMB(I),(VOCOLSINP(I,J),J=1,12)
	 ENDDO
	 CLOSE(UNIT=2)

c  --    The input parameters are copied to an array of 73 classes 
c        in order to have assigned specific values to all the Olson
c        classification numbers

	 I=1
	 DO II=0,MAXCLASS-1
	  DO J=1,12
            VOCPARAM(II,J)=0.0
            IF (II.EQ.OLSNUMB(I)) THEN
              VOCPARAM(II,J)=VOCOLSINP(I,J)
              IF (J.EQ.12) I=I+1
            ENDIF
	  ENDDO
	 ENDDO

      ENDIF ! endif (nstep.eq.0)

C LG- including some first-order estimate of a spinning-up of the
C     emissions after sunset. The assumption is being made that
C     it takes some time for the vegetation to produce the VOC's
C     to be emitted.

      IF (RBVD.GT.1.E-2.AND.RBVD_OLD.LE.1.E-2.AND..NOT.LSUNRISE) 
     &  LSUNRISE=.TRUE.

      LSUNRISE=.FALSE. ! outcomment to activate

      IF (LSUNRISE) THEN
         IF (NSTEP_SUNRISE*DTIME.LE.DTIME_SPINUP_EMIS) THEN
           NSTEP_SUNRISE=NSTEP_SUNRISE+1
	 ELSE
	   NSTEP_SUNRISE=0.
           LSUNRISE=.FALSE.
         ENDIF
      ELSE
         NSTEP_SUNRISE=0
      ENDIF
      IF (NSTEP_SUNRISE.GT.0) THEN
         FEMIS_SUNRISE=MIN(1.,(NSTEP_SUNRISE*DTIME)/DTIME_SPINUP_EMIS)
         print *,' VOC emissions: spinup for C5H8 emissions considered!'
         print *,' FEMIS_SUNRISE is: ',FEMIS_SUNRISE
      ELSE 
         FEMIS_SUNRISE=1.
      ENDIF

C LG- initialisation if integrated isoprene emission rate 
C     which is used in the hydrocarbon chemistry mode without
C     considering the interactions within the biosphere

      ISOPEM_INT=0.
      MTERPEM_INT=0.
      OVOCEM_INT=0.
      DM_FSL=0.
           
C LG- when LTVEG is FALSE, then the canopy temperature resembles the
C     surface temperature, otherwise the explicitly calculated canopy
C     temperature is being used

      CT=EXP((CT1*(TVEG(1)-TS))/
     &       (R*TS*TVEG(1)))/
     &     (1.+EXP((CT2*(TVEG(1)-TM))/
     &       (R*TS*TVEG(1))))

c  -- A number of layers within the canopy are distinguished
c     and for each layer the extinction of PAR is determined.
c     The amount of total biomass is distributed over these
c     layers, expressed by the Leaf Area Density (LAD) and 
c     combined with the LAI to calculate the emission from each 
c     layer and the total emission from the biome. NLEVV is the
c     top layer!!! The direct PAR is calculated from the direct visible
c     radiation and the zenith angle and combined with the fraction
c     of sunlit leaves and the total biomass yielding the emission flux 
c     of the fraction directly effected by the sun. The diffuse PAR is
c     a function of the location within the canopy and the fraction of
c     shaded leaves (1.-FSL)  

      IF (.NOT.LBIOSPH) THEN
        NLEVEL=NLEVV
      ELSE
        NLEVEL=NLEVVEG
      ENDIF

C LG- for LRAD_1LAY=.TRUE., the emissions are calculated considering
C     the vertical radiation profiles for only one vegetation layer

      IF (LRAD_1LAY) NLEVEL=1

      PARDIR=RBVD*RECALC
      PARSL=PARDIR+RVD(NLEVEL)*RECALC
           
      DO II=1,NLEVEL
        JJ=NLEVEL+1-II

        PARDIF(II)=RVD(II)*RECALC
        IF (RBVD.EQ.0.01) THEN
          PARDIR=0.
          PARDIF(II)=0.
        ENDIF

        FOLDENSLAY(II)=DM*LAD(II)

C LG-   calculating the intergrated amount of biomass exposed to
C       sunlight in order to make the comparison for different LAD
C       profiles 

        DM_FSL=DM_FSL+FOLDENSLAY(II)*FSL(II)

        CLSUN=(ALPHA*CL1*PARDIR)/
     &      (SQRT(1.+ALPHA**2*PARDIR**2))
        CLSHADE=(ALPHA*CL1*PARDIF(II))/
     &      (SQRT(1.+ALPHA**2*PARDIF(II)**2))

c
c  --   assigning a biome-related emission-rate
c       (mg/ug g-1 h-1) with g-1 stands for the 
c       amount of biomass
c
c  --   The emission factors by Guenther et al., 1995 are in
c       ug C g-1 h-1 and the foliar density DM is g m-2
c
        EMISFACT(1)=FEMIS_SUNRISE*EMISFACT_VOC(1)

C LG-   introduce a dependence of the emission factor as a function 
C       of the height in the canopy to mimic the difference in the
C       emission potential between sunlit and shaded leaves

C         IF (NSTEP.EQ.0.AND.JJ.EQ.1) THEN
C           print *,' VOC emissions: assumed vertical profile in emission factor'
C           print *,' see lines +/- 240 for modifications'
C           print *,' ENTER TO CONTINUE'
C           READ (*,*)
C 	ENDIF
C         IF (JJ.EQ.1) EMISFACT(1)=1.5*EMISFACT_VOC(1)
C         IF (JJ.EQ.2) EMISFACT(1)=1.25*EMISFACT_VOC(1)
C         IF (JJ.EQ.3) EMISFACT(1)=0.75*EMISFACT_VOC(1)
C         IF (JJ.EQ.4) EMISFACT(1)=0.5*EMISFACT_VOC(1)

c
c  --   make flux EMISFACT dependent on the radiation (PAR)
c       and temperature (see Guenther et al., 1993, JGR)
c       The PAR is calculated from the net short wave radiation
c       at the surface (See subroutine CALCPAR) 

        FLUXSUN(II,1)=EMISFACT(1)*CLSUN*CT*
     &          FOLDENSLAY(II)*FSL(II)
        FLUXSHADE(II,1)=EMISFACT(1)*CLSHADE*CT*
     &          FOLDENSLAY(II)*(1.-FSL(II))

c    -- The emission is in ug C m-2 hr-1 and must
c       be recalculated to the emission in [kg C m-2 s-1]

        ISOPEM(II)=(FLUXSUN(II,1)+FLUXSHADE(II,1))*1.E-9/(3600.)

C LG-   calculation of integrated isoprene emission rate for
C       bulk approach

        ISOPEM_INT=ISOPEM_INT+ISOPEM(II)

C LG------------------------------------------------------------------------
C       The calculated isoprene emission flux is a grid average emission
C       flux since the used input parameters controlling the emission, e.g.
C       the dry matter and and the LAI through controlling the vertical 
C       profiles of radiation, are also grid average values. Therefore, 
C       applying this emission flux in the bulkveg routine, in which the
C       fractions surface cover are also considered, the emission flux needs
C       to be increased by dividing it through the sum of the vegetation and
C       wet skin fraction, so that the total intergrated emission flux of
C       all the four surface cover fractions resembles the flux calculated
C       within this routine 
C LG-------------------------------------------------------------------------

c  --   emission of monoterpenes 

        EMISFACT(2)=EMISFACT_VOC(2)

c  --   The emission is in ug C m-2 hr-1 and must
c       be recalculated to the emission in [kg C m-2 s-1]

        MTERPEM(II)=EMISFACT(2)*EXP(BETA*(TVEG(1)-TS))*
     &             FOLDENSLAY(II)*1.E-9/(3600.)

!         IF (NSTEP.EQ.0.AND.II.EQ.NLEVEL) THEN
! 	  WRITE(*,'(1a)')' VOCEMIS.f: modified repres. of terpene emiss.!'
!           WRITE(*,'(1a)')' Also including the light dependence'
!           READ (*,*)
!         ENDIF
! 
!         FLUXSUN(II,2)=EMISFACT(2)*CLSUN*CT*
!      &          FOLDENSLAY(II)*FSL(II)
!         FLUXSHADE(II,2)=EMISFACT(2)*CLSHADE*CT*
!      &          FOLDENSLAY(II)*(1.-FSL(II))
! 
! c  --   The emission is in ug C m-2 hr-1 and must
! c       be recalculated to the emission in [kg C m-2 s-1]
! 
!         MTERPEM(II)=(FLUXSUN(II,2)+FLUXSHADE(II,2))*1.E-9/(3600.)

C LG-   calculation of integrated emission rate for
C       bulk approach

        MTERPEM_INT=MTERPEM_INT+MTERPEM(II)

c  --   and OVOC's

        EMISFACT(3)=EMISFACT_VOC(3)

        OVOCEM(II)=EMISFACT(3)*EXP(BETA*(TSURFACE(1)-TS))*
     &             FOLDENSLAY(II)*1.E-9/(3600.)

!         IF (NSTEP.EQ.0.AND.II.EQ.NLEVEL) THEN
! 	  WRITE(*,'(1a)')' VOCEMIS.f: modified repres. of OVOCs emiss.!'
!           WRITE(*,'(1a)')' Also including the light dependence'
!           READ (*,*)
!         ENDIF
! 
!         FLUXSUN(II,3)=EMISFACT(3)*CLSUN*CT*
!      &          FOLDENSLAY(II)*FSL(II)
!         FLUXSHADE(II,3)=EMISFACT(3)*CLSHADE*CT*
!      &          FOLDENSLAY(II)*(1.-FSL(II))
! 
! c  --   The emission is in ug C m-2 hr-1 and must
! c       be recalculated to the emission in [kg C m-2 s-1]
! 
!         OVOCEM(II)=(FLUXSUN(II,3)+FLUXSHADE(II,3))*1.E-9/(3600.)

C LG-   calculation of integrated emission rate for
C       bulk approach

        OVOCEM_INT=OVOCEM_INT+OVOCEM(II)

C LG-   opening of file and writing of radiation data for checking
C       the influence of the canopy structure on the hydrocarbon emissions
C       through the radiation profiles

        IF (NSTEP.EQ.0.AND.II.EQ.1) THEN
          FNAME='/data/ganzevl/racmo/output/radiation_veg.out'
          OPEN(UNIT=NUNRADV,FILE=FNAME,STATUS='UNKNOWN')

          NLEV1=NLEVATM+1
          NLEV2=NLEVATM+NLEVEL
          NLEVPR=NLEV2-NLEV1+1

          WRITE (NUNRADV,*) NSTOP,NLEVPR,NPRINT,DTIME
          WRITE (NUNRADV,*) NLEV1,NLEV2

C LG-     calculation of actual height for big leaf approach, bulk or
C         two-layer vegetation model, otherwise the HGHT has already been 
C         assigned an explicit value

          IF (.NOT.LBIOSPH) THEN
            DO JK=NLEV1,NLEV2
	      JJ=NLEV2+1-JK
              HEIGHT(JK)=(HC/NLEVPR)*JJ-0.5*(HC/NLEVPR)
            ENDDO
	  ENDIF

	  WRITE (NUNRADV,*) (HEIGHT(JK),JK=NLEV2,NLEV1,-1)
        ENDIF

        IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN

         IF (II.EQ.1) THEN
          WRITE(NUNRADV,'(a14)') LDATLTIME
	  WRITE(NUNRADV,*)NSTEP,JDAY,GMT,LTIME,HC,
     *                                   RG,RBVD,CT,CLSUN
         ENDIF
	 WRITE(NUNRADV,'(I4.4,1X,8F10.3,5E10.3)')
     *          II,PARDIR,PARDIF(II),PARDIR+PARDIF(II),
     *          (PARDIR/RECALC)*FSL(II)/MAX(1.E-20,RG),
     *          (PARDIR*FSL(II)+PARDIF(II)*(1.-FSL(II)))/
     *           MAX(1.E-20,PARSL),CLSHADE,FOLDENSLAY(II),
     *           FSL(II),ISOPEM(II)/(1.E-9/3600.),
     *           FLUXSUN(II,1),FLUXSHADE(II,1),
     *           MTERPEM(II)/(1.E-9/3600.),
     *           OVOCEM(II)/(1.E-9/3600.)

         IF (II.EQ.NLEVEL) THEN
           IF (DM.GT.0.) THEN 
             WRITE(NUNRADV,*)DM_FSL/DM
	   ELSE
	     WRITE(NUNRADV,*)DM_FSL
           ENDIF

C LG-         writing the total intergrated emission flux for comparison

           WRITE(NUNRADV,*)ISOPEM_INT/(1.E-9/3600.)
           WRITE(NUNRADV,*)MTERPEM_INT/(1.E-9/3600.)
           WRITE(NUNRADV,*)OVOCEM_INT/(1.E-9/3600.)

	 ENDIF

        ENDIF

C LG-   end loop DO II=1,NLEVEL

      ENDDO

C LG- setting the old global radiation

      RBVD_OLD=RBVD

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNRADV)

      RETURN
      END




