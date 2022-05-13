      SUBROUTINE VEGETATION(NSTEP,NSTOP,NPRINT,NLEV,NLEVEL,IRESET_VEG,
     *     SINLAT,SINLON,COSLAT,COSLON,FORESTM,ALBEDOM)

c     ------------------------------------------------------
c     this routine reads/assigns the vegetation database, 
c     LAI data and assigns the LAD profiles to each 
c     gridsquare and scales this to the proper resolution
c     ------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

      INTEGER ZILON,ZILAT,IILON,IILAT,IIILON,IIILAT,
     &     ILONINIT,ILATINIT,OLSCLASS,MAXCLASS,NNVEG,
     &     NPARAM,NRESET,IRESET_VEG,IRESETINP,IRESET_OLS

      PARAMETER (IILON=720,IILAT=360,OLSCLASS=59,MAXCLASS=73,
     &     NNVEG=73,NPARAM=20)
      
      INTEGER JK,I,II,J,IRESET,IP,IIP,JJ,NSTEP,NSTOP,NPRINT,
     &     NLEV,NLEVEL,NLEVELS,JL

      INTEGER IVEG (2160,1080),GVI(IILON,IILAT),GVIMAX(IILON,IILAT),
     &     OLSVEG(IILON,IILAT),OLSNUMB(0:OLSCLASS)
      
      LOGICAL LREAD,LRESET,L_OLSON_RESET

      REAL OLSINP(0:OLSCLASS,NPARAM),PARAM(0:MAXCLASS,NPARAM)

      REAL
     *     SINLAT (NLON), COSLAT(NLON),
     *     SINLON (NLON), COSLON(NLON)

      REAL API,FORESTM(NLON),ALBEDOM(NLON),ZLAT,ZLON,SLW,
     *     HCOLS,Z0MOLS,ZHP(NLEVT),Z(NLEVT),MIN_HC,Z0MCALC,
     *     DISPCALC,LAI_AVG,HC_AVG,Z0M_AVG,DISP_AVG,
     *     EMISFACT_VOC_AVG(NSPEC_EMISVOC),
     *     DM_AVG,SLW_AVG
      
      CHARACTER*255 CHDUM
      CHARACTER*255 DUMMY
      CHARACTER*70 FNAME
      CHARACTER*70 DIR
      CHARACTER*2 MON(12)
      CHARACTER*2 NUMBER(0:MAXCLASS),NUMBOLS
      CHARACTER*250 OLSNAME(0:MAXCLASS),OLSON

      SAVE IVEG,OLSVEG,GVI,GVIMAX,OLSNUMB,OLSINP,PARAM,
     &     ILONINIT,ILATINIT    ! mz_lg_20051211+

      API=3.1415927

C     LG- definition of month

      DATA MON /'01','02','03','04','05','06','07','08','09',
     &     '10','11','12'/

C     WP- adding the directory for reading vegetation input

      DIR='../../echam/'
      IIP=INDEX(DIR,' ')

C     WP- end
      
C     LG- added the switch set to true when the Olson ecosystem class is modified      
      
      L_OLSON_RESET=.FALSE. 

      WRITE(NUNMDFL,'(2a)')
     *     ' Assigning land cover type and its structure'
      IF (.NOT.LBIOSPH.AND.LVOCEMIS)
     *     WRITE(NUNMDFL,'(2a)')' (required for the calculation of ',
     *     'the hydrocarbon emissions)' 

      IF (.NOT.LAGRIAN) THEN
         WRITE(*,'(2a)')
     *        ' Assigning land cover type and its structure'
         IF (.NOT.LBIOSPH.AND..NOT.LBULKVEG.AND..NOT.LVEG_MLAY.
     *        AND.LVOCEMIS)
     *        WRITE(*,'(2a)')' (required for the calculation of ',
     *        'the hydrocarbon emissions)'
         WRITE(*,*)
      ENDIF

      IF (NSTEP.EQ.0.OR.(LCHMONTH.AND.LRESETINP)) THEN
         IRESET=0
         NRESET=0
         IRESETINP=0
         LRESETINP=.TRUE.       ! setting LRESETINP default to TRUE. (see also 889)
                                ! so that for a change in the month the input is
                                ! updated. Only when the input parameters are reset
                                ! and the user indicates that the initial parameters
                                ! have to be used throughout the simulation,
                                ! lresetinp is set to .FALSE.
         LREAD=.TRUE.
      ENDIF

C     LG- calculation of the column and row numbers of the 0.5 * 0.5 grid
C     coordinates resembling the position of the RACMO column for the
C     calculation of the proper vegetation characteristics
      
      ZLAT=ASIN(SINLAT(1))*180./API
      ZILAT=(90.-ZLAT)/(2*90./IILAT)+1
      ZLON=ASIN(SINLON(1))*180./API
      IF (COSLON(1).LT.0.AND.SINLON(1).GT.0) THEN
         ZLON=180.-ZLON
      ELSEIF (COSLON(1).LT.0.AND.SINLON(1).LT.0) THEN
         ZLON=-180.-ZLON
      ENDIF
      ZILON=(ZLON+180.)/(2*180./IILON)+1

C     LG- assigning filename of file with vegetation data. It is 
C     important to use the criteria NSTEP.LE.1 instead of NSTEP
C     .EQ.0 here otherwise the initial longitude and latitude are
C     not properly defined.
      
!     mz_lg_20051211+ modified by putting iloninit and ilatinit in SAVE statement

      IF (NSTEP.EQ.0.OR.(LCHMONTH.AND.LRESETINP)) THEN
         ILONINIT=ZILON
         ILATINIT=ZILAT
      ENDIF

!     iloninit=261
!     ilatinit=172

      print *,'vegetation.f, modify; iloninit',nstep,lchmonth,lresetinp,zilon,iloninit

      FNAME=dir(1:iip-1)//'case/input/veg_'
      IP=INDEX(FNAME,' ')
      FNAME=FNAME(1:IP-1)//MON(IMON)
      IP=INDEX(FNAME,' ')
      WRITE(DUMMY,'(i3.3)')
     +     ILONINIT
      FNAME=FNAME(1:IP-1)//DUMMY
      IP=INDEX(FNAME,' ')
      WRITE(DUMMY,'(i3.3)')
     +     ILATINIT
      FNAME=FNAME(1:IP-1)//DUMMY
      
      WRITE(*,'(2a)')' vegetation.f: Reading file with the name: ',
     *     FNAME

      OPEN(UNIT=2,FILE=FNAME,FORM='FORMATTED',
     +     STATUS='UNKNOWN')

C     LG- reading of file with vegetation data derived from the 
C     Olson ecosystem database and the satellite data, if this
C     file is not available for the location and the specific month
C     this file is created (ERR=999)

      REWIND(2) 

      READ(2,'(256A)',END=999,ERR=999) DUMMY
      READ(2,'(256A)',END=999,ERR=999) DUMMY

 200  READ(2,'(2i5,1x,i9,1x,f9.2,1x,i9,1x,8(f9.2,1x),3(i9,1x),a150)',
     +     END=201,ERR=999) IIILON,IIILAT,OLSVEG(IIILON,IIILAT),
     +     LAI,IPROF,HCOLS,Z0MOLS,DISP,EMISFACT_VOC(1),EMISFACT_VOC(2),
     +     EMISFACT_VOC(3),DM,SLW,INOCLASS(1),INOCLASS(2),PC3C4TYPE,
     +     OLSNAME(OLSVEG(IIILON,IIILAT))

      IF (IIILON.EQ.ZILON.AND.IIILAT.EQ.ZILAT) GOTO 201

      GOTO 200

 201  CONTINUE

      IF (NSTEP.EQ.NSTOP) CLOSE(2)

C     LG- if the read longitude/latitude does not resemble the longitude/
C     latitude of the available "preprocessed" surface cover data,
C     this is then done for the right coordinates 

      IF (IIILON.NE.ZILON.OR.IIILAT.NE.ZILAT) GOTO 999

      GOTO 9999

 999  CONTINUE

      WRITE(*,'(2a,2i4)')' Determining the land cover and its ',
     +     'characteristics for the Olson grid square: ',ZILON,ZILAT

      CLOSE(2)

      IF (LREAD) THEN

 1111    WRITE(*,'(2a)')' Reading files with vegetation data for ',
     +        'determining surface cover'

C     LG-  reading file with specific characteristics of the 73 ecosystems
C     discerned in the '92 database

         OPEN (10,FILE=
     &        dir(1:iip-1)//'input/veg/VOC_NO_Olson_92.dat'
     &        ,STATUS='OLD')

         DO J=1,14
            READ (10,'(255A)') CHDUM
         ENDDO

         DO I=1,OLSCLASS
            READ (10,*) OLSNUMB(I),(OLSINP(I,J),J=1,NPARAM)
         ENDDO

         CLOSE(10)

c     --  The input parameters are copied to an array of 73 classes 
c     in order to have assigned specific values to all the Olson
c     classification numbers

         I=1
         DO II=0,MAXCLASS-1
            DO J=1,NPARAM
               PARAM(II,J)=0.0
               IF (II.EQ.OLSNUMB(I)) THEN
                  PARAM(II,J)=OLSINP(I,J)
                  IF (J.EQ.NPARAM) I=I+1
               ENDIF
            ENDDO
         ENDDO

C     LG-  reading file with the ecosystem category number and the 
C     ecosystem name and a short description

         OPEN (11,FILE=
     &        dir(1:iip-1)//'input/veg/Olson_73.dat',STATUS='OLD')

         DO I=0,MAXCLASS
            READ (11,'(A,1X,A)') NUMBER(I),OLSNAME(I)
         ENDDO

         CLOSE(11)

c     --  call subroutine to read the vegetation data and to scale this
c     to the proper resolution, IVEG are the high resolution input data
c     (0.167 deg, 2160/1080) whereas OLSVEG contains the 0.5 deg resolution
c     vegetation data (see VOC_geia_85.f in ganzevl/Olson)
c     
         IF (.NOT.L_OLSON_RESET) ! LG- modified
     &        CALL READVEG(IVEG,OLSVEG)
c     
c     --   call subroutine to read GVI (Global Vegetation Index data 
c     which are used to determine the foliar density and the LAI of 
c     each grid square

         CALL READGVI(IMON,GVI,GVIMAX)

         LREAD=.FALSE.

      ENDIF
c     
c     --  call subroutine in which the foliar density and LAI are calculated
c     from the GVI data according to Guenther et al.
c     
      CALL CALCVEG(IMON,ZILON,ZILAT,PARAM,OLSVEG,GVI,GVIMAX,LAIMAX,
     &     DM,LAI)

C     LG- assigning the proper LAD profile, SLW value, canopy height and 
C     surface roughness

      IPROF=PARAM(OLSVEG(ZILON,ZILAT),12)
      SLW=PARAM(OLSVEG(ZILON,ZILAT),10)      
      HCOLS=PARAM(OLSVEG(ZILON,ZILAT),13)
      Z0MOLS=PARAM(OLSVEG(ZILON,ZILAT),17)

C     LG- assigning the isoprene, monoterpene and other VOC's emission factor
      
      EMISFACT_VOC(1)=PARAM(OLSVEG(ZILON,ZILAT),1)
      EMISFACT_VOC(2)=PARAM(OLSVEG(ZILON,ZILAT),2)

C     LG- for the other VOC's only the global total emission flux per ecosystem
C     is available and to estimate the ecosystem specific emission factor in
C     mug C g-1 hr-1, the ratio's of the global total emission fluxes of the
C     OVOC's and monoterpenes are used for the scaling 

      IF (PARAM(OLSVEG(ZILON,ZILAT),4).GT.0)
     &     EMISFACT_VOC(3)=
     &     (PARAM(OLSVEG(ZILON,ZILAT),5)/PARAM(OLSVEG(ZILON,ZILAT),4))*
     &     EMISFACT_VOC(2)

C     LG- assigning the NO emission classes, which are used in the routine
C     NOXEMIS

      INOCLASS(1)=0
      INOCLASS(2)=0.

      INOCLASS(1)=PARAM(OLSVEG(ZILON,ZILAT),15)
      INOCLASS(2)=PARAM(OLSVEG(ZILON,ZILAT),16)

C     LG- assigning the C3/C4 type, required for the calculation of the
C     stomatal exchange by the physiological model by Ronda et al., 2001

      PC3C4TYPE=PARAM(OLSVEG(ZILON,ZILAT),20)

!     mz_lg_20051212+ defining minimum LAI in case there is vegetation, with a
!     canopy height but an inferred LAI of zero. This happens due to different
!     resolutions in the NDVI satellite data and Olson ecosystem database

      IF (LAGRIAN.AND.OLSVEG(ZILON,ZILAT).GT.0.AND.
     &     HCOLS.GT.1e-10.AND.LAI.LT.1.E-10) THEN
         WRITE(2,'(2a)')'LAI set to value of 0.1 since canopy height > 0 ',
     &        ' whereas the inferred LAI from satellite data and Olson =0'
         LAI=0.1
         DM=LAI*MAX(1.E-20,SLW)
      ENDIF

      IF ((.NOT.LAGRIAN.AND.IRESET_VEG.EQ.1).OR.
     &     (L_OLSON_RESET)) GOTO 1001

C     C LG- writing outputfield to check the global vegetation field
C     
C     OPEN(UNIT=3,FILE='/data/ganzevl/racmo/output/noclass1.out',
C     &      FORM='FORMATTED',STATUS='UNKNOWN',RECORDTYPE='STREAM')
C     WRITE(3,'(72i4)') ((PARAM(OLSVEG(I,J),14),I=1,720),J=360,1,-1)
C     CLOSE(3)
C     
C     C LG- end

C     LG- writing to output file which is used for repetitive runs for
C     same location and month, these data are only written away for the
C     first simulation for the specific month and location

      OPEN(UNIT=2,FILE=FNAME,FORM='FORMATTED',
     +     STATUS='OLD')

C     LG- in order to add data and not to overwrite data!

 20   READ(2,*,END=21)
      GOTO 20
 21   BACKSPACE 2
      
      IF (NSTEP.EQ.0.OR.LCHMONTH) THEN
         WRITE(2,'(2a)')'Surface cover and some of its characteristics ',
     +        'derived from Olson database and satellite data'
         WRITE(2,'(2(a5),14(a10),1a)')
     +        'ILON','ILAT','Olson no.','LAI','LAD prof.',
     +        'HC','z0','DISP','C5H8-emis','C10-emis','OVOC-emis','DM','SLW',
     +        'NO-class1','NO-class2','C3/C4',' Olson ecosystem type descr.'
      ENDIF
      
      OLSON=OLSNAME(OLSVEG(ZILON,ZILAT))
      IP=INDEX(OLSON,'  ')      
      OLSON=OLSON(1:IP-1)

      WRITE(2,'(2i5,1x,i9,1x,f9.2,1x,i9,1x,8(f9.2,1x),3(i9,1x),1a)')
     +     ZILON,ZILAT,OLSVEG(ZILON,ZILAT),LAI,
     +     IPROF,HCOLS,Z0MOLS,DISP,EMISFACT_VOC(1),EMISFACT_VOC(2),
     +     EMISFACT_VOC(3),DM,SLW,INOCLASS(1),INOCLASS(2),PC3C4TYPE,
     +     OLSON(1:IP-1)

      IF (NSTEP.EQ.NSTOP) CLOSE(2)

 9999 CONTINUE

C     LG- modification of surface cover parameters for specific domain,
C     the modified parameters *_RESET, are defined in oned.f

      IF (IRESET_VEG.EQ.1) THEN ! mz_lg_20040826+ modified IF (IRESET_VEG)
         OLSVEG(ZILON,ZILAT)=OLSON_RESET
         IF (NRESET.GT.1) THEN
            IPROF=PARAM(OLSVEG(ZILON,ZILAT),12)
            SLW=PARAM(OLSVEG(ZILON,ZILAT),10)      
            HCOLS=PARAM(OLSVEG(ZILON,ZILAT),13)
            Z0MOLS=PARAM(OLSVEG(ZILON,ZILAT),17)

C     LG-   assigning the isoprene, monoterpene, OVOC's and NOx emission factors
            
            EMISFACT_VOC(1)=PARAM(OLSVEG(ZILON,ZILAT),1)
            EMISFACT_VOC(2)=PARAM(OLSVEG(ZILON,ZILAT),2)
            EMISFACT_VOC(3)=PARAM(OLSVEG(ZILON,ZILAT),3)
            INOCLASS(1)=PARAM(OLSVEG(ZILON,ZILAT),15)
            INOCLASS(2)=PARAM(OLSVEG(ZILON,ZILAT),16)
            PC3C4TYPE=PARAM(OLSVEG(ZILON,ZILAT),20)
         ENDIF
         LAI=LAI_RESET

C     LG-  recalculation of foliar density, scaled with the LAI,
C     08-2002, in case of a lagragian analyses the dry matter is corrected
C     for the vegetation fraction in a grid square. It still must checked
C     if this should also be done for non-lagrian simulations. The reason
C     why this has not be done yet is that the inferred DM should already
C     reflect an area average measurement, reflecting the vegetation 
C     fraction implicitly. We could have also corrected the predefined LAI
C     but this parameter often reflects local value and not grid-average
C     values. The dry matter is used to determine the isoprene emission 
C     flux and for a vegetation fraction close to zero, the isoprene 
C     emission flux should be close to zero.

         DM=LAI*MAX(1.E-20,SLW)
         IF (LAGRIAN) DM=VEGFRAC(1)*LAI*MAX(1.E-20,SLW)

         HCOLS=HC_RESET
         NRESET=NRESET+1
         IF (NRESET.EQ.1) GOTO 1111
      ENDIF

 1001 CONTINUE

      OLSONVEG=OLSVEG(ZILON,ZILAT)
      LAI=MAX(1.E-10,LAI)

C     LG- determining canopy height and surface roughness from these 
C     parameters defined for each ecosystem discerned in the Olson 
C     ecosystem database. The displacement height is taken as a 
C     constant fraction of the canopy height (0.66). The minimal 
C     canopy height is set to 1e-10 unless HCOLS > 0, for the 
C     simulations in the lagragian mode, the minimal canopy height
C     is dependent on the number of layers, the maximum canopy height 
C     and the verticl coordinate system. For a equidistant coordinate
C     system and four canopy layers, the minimum canopy height is 7.5 m
C     for a maximum canopy height of 30 m whereas for the logaritmic
C     coordinate system the minimal canopy height resembles the value
C     of MINTHICK (minimal layer thickness)

      MIN_HC=1.E-10

      IF (LBIOSPH.AND.LAGRIAN.AND.ILSTYPE.GT.1) THEN
         IF (LDTHICK) THEN
            IF (NSTEP.EQ.0) THEN
               WRITE(*,'(1a,f4.1)')
     &              ' The minimal resolved canopy height is: ',MINTHICK
               WRITE(*,'(2a)')' For a higher vert. resolution, increase the ',
     &              'number of layers (NLEVV parchem.h)'
               WRITE(*,'(1a)')' Enter to continue'
               READ (*,*)
            ENDIF
            IF (HCOLS.GT.0.) MIN_HC=MINTHICK
         ELSE
            IF (NSTEP.EQ.0) THEN
               WRITE(*,'(1a,f4.1)')
     &              ' The minimal resolved canopy height is: ',MAX_CANHEIGHT/NLEVV
               WRITE(*,'(2a)')' For a higher vert. resolution, increase the ',
     &              'number of layers (NLEVV parchem.h)'
               WRITE(*,'(2a)')' or set the switch LDTHICK in namchem ',
     &              'to .TRUE. for logaritmic vert. coord. system.'
               WRITE(*,'(1a)')' Enter to continue'
               READ (*,*)
            ENDIF
            IF (HCOLS.GT.0.) MIN_HC=MAX_CANHEIGHT/NLEVV 
         ENDIF 
      ENDIF

      HC=MAX(MIN_HC,HCOLS)
      HCCORR=HC

C     LG- The maximum canopy height is set to a maximum being defined in 
C     parchem.h. HCMAX is set to zero for the "big leaf" approach
      
      IF (.NOT.LBIOSPH.AND..NOT.LBULKVEG.AND..NOT.LVEG_MLAY) THEN
         HCMAX=0.
      ELSEIF (LAGRIAN.AND.ILSTYPE.GT.1.) THEN
         HCMAX=MAX_CANHEIGHT
         HCCORR=HCMAX
      ENDIF

C     LG- end

      NLEVVEG=0
      IF (LBULKVEG.AND.HC.GT.1.E-10) THEN
         NLEVVEG=1
      ELSEIF (LVEG_MLAY.AND.HC.GT.1.E-10) THEN
         NLEVVEG=NLEVV_ML
      ENDIF

C     LG- for the multilayer, bulk vegetation or two layer vegetation mode 
C     calculations, the canopy height is related to the number of layers in 
C     order to ensure minimal thickness's of the layers

      IF (LBIOSPH.OR.(LVEG_MLAY.AND..NOT.LAGRIAN)) THEN ! mz_lg_20060109+

         IF (.NOT.LDTHICK) THEN

            IF (HCCORR/NLEVV.GE.MINTHICK) THEN
               HCMAX=HC
               NLEVVEG=MIN(NLEVV,INT(HCCORR/MINTHICK))
            ELSEIF (LBIOSPH) THEN
               WRITE(*,'(1a)')
     &              ' Too thin vegetation layers (numerical problems). '
               WRITE(*,'(1a,i3)')
     &              ' The actual number of vegetation layers is set to: ',
     &              INT(HCCORR/MINTHICK)
               NLEVVEG=MIN(NLEVV,INT(HCCORR/MINTHICK))
	       NLEVEL=NLEV+NLEVVEG
            ENDIF

         ELSEIF (LDTHICK) THEN

C     LG-     determining the reference height full pressure and half pressure level
C     for the logarithmic vertical coordinate system. This is done in order
C     to calculate the number of vegetation levels ensuring a minimal 
C     thickness of all the layers resembling the value of the parameter
C     MINTHICK, which is being set in the file parchem.h

            NLEVVEG=MIN(NLEVV,INT(HCCORR/MINTHICK))
            NLEVEL=NLEV+NLEVVEG

C     LG-     for lagrian simulations over land and for the biosphere mode
C     the reference height at full pressure level for the maximum 
C     number of levels must be determined

            IF (LAGRIAN.AND.ILSTYPE.GT.1.AND.LBIOSPH) THEN
               NLEVVEG=NLEVV
               NLEVEL=NLEV+NLEVVEG
            ENDIF

 98         CONTINUE

            DO 127 JK=NLEV+1,NLEVEL

               JJ=NLEVVEG+1-(JK-NLEV)

C     LG-       determining coefficients to define the height coordinates within
C     the canopy with an increasing thickness with increasing height

               CTHICK_CANLAY(JJ)=
     *              (1.-(FLOAT(NLEVVEG-JJ)/FLOAT(NLEVVEG-1)))*ALOG(HCCORR)+
     *              (FLOAT(NLEVVEG-JJ)/FLOAT(NLEVVEG-1))*ALOG(MINTHICK)

               CTHICK_CANLAY(JJ-1)=
     *              (1.-(FLOAT(NLEVVEG-(JJ-1))/FLOAT(NLEVVEG-1)))*ALOG(HCCORR)+
     *              (FLOAT(NLEVVEG-(JJ-1))/FLOAT(NLEVVEG-1))*ALOG(MINTHICK)

C     LG-       vertical coordinate system within the canopy with increasing
C     thickness of the layers with increasing height, the lowest layer
C     has a thickness of MINTHICK being defined in the parameter file
C     parchem.h 

               ZHP(JK)=EXP(CTHICK_CANLAY(JJ))
               Z(JK)=(EXP(CTHICK_CANLAY(JJ))+
     &              MAX(0.,EXP(CTHICK_CANLAY(JJ-1))))/2.

               IF (JK.EQ.NLEVEL) THEN
                  ZHP(JK+1)=0.  
                  Z(JK)=ZHP(JK)/2.
               ENDIF

               IF ((ZHP(JK)-Z(JK))*2.LT.MINTHICK) THEN
                  WRITE(*,'(1a)')
     &                 ' Too thin vegetation layer within log. vert. coord. system' 
                  WRITE(*,'(1a,i3)')
     &                 ' The actual number of vegetation layers is reset to: ',
     &                 NLEVVEG-1
                  NLEVVEG=NLEVVEG-1
                  NLEVEL=NLEVEL-1
                  GOTO 98
               ENDIF

 127        CONTINUE

C     LG-     end correction of number of levels for logarithmic vertical coordinate
C     system

         ENDIF

C     LG-   The canopy height assigned to the different ecosystem types is corrected
C     such that it resembles the full pressure height of the top canopy
C     layer

         DO 128 JK=NLEV+1,NLEVEL

            JJ=NLEVVEG+1-(JK-NLEV)

C     LG-     defining the height for the equidistant coord. system

            IF (.NOT.LDTHICK) THEN
               ZHP(JK)=(HCCORR/NLEVVEG)*JJ
               Z(JK)=(HCCORR/NLEVVEG)*JJ-(1./2.)*(HCCORR/NLEVVEG)
            ENDIF

 128     CONTINUE
         IF (LBIOSPH.AND.HC.GT.Z(NLEV+1)) HC=ZHP(NLEV+1) ! mz_lg_20060110+ lbiosph

C     LG- end if (LBIOSPH)

      ENDIF

C     LG- call subroutine for calculation of the surface roughness and 
C     displacement heigth from the LAI and canopy height according
C     to Raupach, "Simplified expressions for vegetation roughness.."
C     in Boundary Layer Meteorology, 71, 211-216, 1994
C     
      CALL CALCZ0D(LAI,HC,Z0MCALC,DISPCALC)
      DISP=DISPCALC

      IF (HC.GT.1.E-10.AND.Z0MCALC.GE.0.15*HC)
     &     WRITE(*,'(1a,f6.2,1a)')
     &     ' The surface roughness is set to a maximum of 0.15 * HC: ',
     &     0.15*HC,' [m]'

      Z0M=MIN(0.15*HC,MAX(1.E-10,Z0MCALC))

C     LG- printing of initialised vegetation characteristics as a function
C     of the location. 

 110  CONTINUE 

      IF (NSTEP.EQ.0) THEN

         WRITE(*,'(2a)')' The model will be initialised with the next ',
     *        'land cover characteristics :'

         OLSON=OLSNAME(OLSVEG(ZILON,ZILAT))
         IP=INDEX(OLSON,'  ')      
         OLSON=OLSON(1:IP-1)

         WRITE(*,'(1a,i3,1x,a)')' Land cover type acc. to Olson [1992] :',
     *        OLSVEG(ZILON,ZILAT),OLSON(1:IP-1)

         IF (OLSVEG(ZILON,ZILAT).GT.0) THEN
            WRITE(*,'(1a,f4.1)')' Leaf Area Index: ',LAI 
            IF (NLEVVEG.GT.0)
     *           WRITE(*,'(1a,i3)')' The number of vegetation layers: ',
     *           NLEVVEG
            WRITE(*,'(1a,i2)')' The Leaf Area Density profile: ',IPROF 
            WRITE(*,'(1a,f4.1)')' Canopy height [m]: ',HC
            WRITE(*,'(1a,f4.1)')' Displacement height [m]: ',DISP 
            WRITE(*,'(1a,f6.1)')
     *           ' Surface roughness [cm] (calc. from LAI and canopy height): ',
     *           Z0M*100.
            WRITE(*,'(1a,f6.1)')
     *           ' Foliar density (hydrocarbon emiss.) [g m-2]: ',DM
            IF (LVOCEMIS) THEN
               WRITE(*,'(1a,f6.1)')
     *              ' Isoprene emission factor [ug C m-2 hr-1]: ',EMISFACT_VOC(1)    
               WRITE(*,'(1a,f6.1)')
     *              ' Monoterpene emission factor [ug C m-2 hr-1]: ',EMISFACT_VOC(2)  
               WRITE(*,'(1a,f6.1)')
     *              ' OVOCs emission factor [ug C m-2 hr-1]: ',EMISFACT_VOC(3)  
            ENDIF
            IF (LAGS) 
     *           WRITE(*,'(1a,I6.1)')
     *           ' C3/C4 Plant type (1=C3, 2=C4, 3=coniferous trees): ',
     *           PC3C4TYPE    
            IF (.NOT.LAGRIAN) THEN
               WRITE(*,'(2a)')' Do you want to (re)set some of the ',
     *              'vegetation characteristics manually (1)' 
               READ (*,*)IRESET

               IF (IRESET.EQ.1) THEN

 123              WRITE(*,'(1a)')' TYPE 1 for also RESETTING OLSON CLASS'
                  READ *,IRESET_OLS
                  IF (IRESET_OLS.LT.0.OR.IRESET_OLS.GT.1) GOTO 123
                  IF (IRESET_OLS.EQ.1) THEN 

C     LG-         reading file with the ecosystem category number and the 
C     ecosystem name and a short description

                     OPEN (11,FILE=
     &                    dir(1:iip-1)//'input/veg/Olson_73.dat',STATUS='OLD')
                     DO I=0,MAXCLASS
                        READ (11,'(A,1X,A)') NUMBER(I),OLSNAME(I)
                     ENDDO
                     CLOSE(11)

 124                 WRITE(*,'(1a)')' What Olson ecosystem class do you want to select? '
                     WRITE(*,'(1a)')' Give the NUMBER of the next list'
                     DO I=0,MAXCLASS
                        WRITE(*,'(A2,1X,A70)') NUMBER(I),OLSNAME(I)
                     ENDDO
                     READ *,OLSON_RESET
                     WRITE(NUMBOLS,'(i2.2)')
     +                    OLSON_RESET
                     DO I=0,MAXCLASS
                        IF (NUMBOLS.EQ.NUMBER(I)) L_OLSON_RESET=.TRUE.
                     ENDDO            
                     IF (.NOT.L_OLSON_RESET) GOTO 124
                     OLSVEG(ZILON,ZILAT)=OLSON_RESET
                     GOTO 999
                  ENDIF

 222              WRITE(*,'(1a)')' What Leaf Area Index [0-12.5] ? '
                  WRITE(*,'(2a)')' (LAI=amount of biomass)? ',
     *                 'typical value: max. of 10 for trop. forest, 1 for grass) '
                  READ *,LAI
                  IF (LAI.LT.0.OR.LAI.GT.12.5) GOTO 222

C     LG-       recalculation of foliar density, scaled with the LAI

                  DM=LAI*MAX(1.E-20,SLW)
                  WRITE(*,'(1a,f6.1)')
     *                 ' The rescaled foliar density in [g m-2] is : ',DM
                  
 223              WRITE(*,'(1a)')
     *                 ' Which Leaf Area Density (LAD) profile ? (=vert. distr. LAI) '

                  WRITE(*,'(1a)')' 1=agricultural vegetation.'
                  WRITE(*,'(1a)')' 2=decideous forest'
                  WRITE(*,'(1a)')' 3=top profile (nearly all biomass in top)'
                  WRITE(*,'(1a)')' 4=average'
                  WRITE(*,'(1a)')' 5=trop. rainforest Cuieiras, Brazil'
                  WRITE(*,'(1a)')' 6=trop. rainforest Jaru, Brazil'
                  WRITE(*,'(1a)')' 7=avg of Cuieiras and Jaru'
                  WRITE(*,'(1a)')' 8=fruit orchard, UK'

                  WRITE(*,'(2a)')' (See comments in vegetation.f for more ',
     *                 'extensive information about profiles)'
                  READ *,IPROF
                  IF (IPROF.LT.1.OR.IPROF.GT.8) GOTO 223

                  NRESET=0
                  LRESET=.FALSE.

C     LG-       modified 2002-02-07

                  IF (LVOCEMIS.AND..NOT.LBULKVEG.AND..NOT.LVEG_MLAY.AND..NOT.LBIOSPH) 
     &                 NLEVVEG=NLEVV
                  
 101              CONTINUE

 224              WRITE(*,'(1a,f4.1,1a)')
     &                 ' What must be the canopy height [m] [0-',MAX_CANHEIGHT,' m]'
                  READ *,HC
                  IF (HC.LT.0.OR.HC.GT.MAX_CANHEIGHT) GOTO 224
                  HC=MAX(1.E-10,HC)
                  HCCORR=HC
                  
C     LG-       for the multilayer calculations the canopy height is 
C     related to the number of layers in order to ensure minimal 
C     thickness's of the layers

                  IF (LBIOSPH.OR.LVOCEMIS) THEN

                     IF (.NOT.LDTHICK) THEN

                        IF (HC/NLEVVEG.GE.MINTHICK) THEN
                           HCMAX=HC
                           HCCORR=HC
                           NLEVELS=NLEV+NLEVVEG

                           DO 1288 JK=NLEV+1,NLEVELS

                              JJ=NLEVVEG+1-(JK-NLEV)

C     LG-               defining the height for the equidistant coord. system

                              ZHP(JK)=(HCCORR/NLEVVEG)*JJ
                              Z(JK)=(HCCORR/NLEVVEG)*JJ-(1./2.)*(HCCORR/NLEVVEG)
 1288                      CONTINUE
                           ZHP(NLEVELS+1)=0.
                        ELSEIF (LBIOSPH) THEN
                           WRITE(*,'(1a)')
     &                          ' Too thin vegetation layers (numerical problems). '
                           WRITE(*,'(1a,i3)')
     &                          ' The actual number of vegetation layers is reset to: ',
     &                          INT(HC/MINTHICK)
                           NLEVVEG=MIN(NLEVV,INT(HC/MINTHICK))
                           NRESET=NRESET+1
                           IF (NRESET.EQ.2) GOTO 112
                           WRITE(*,'(1a,f4.1,1a)')
     &                          ' or increase the canopy height to at least: ',
     &                          NLEVV*MINTHICK,' [m]'
                           GOTO 101
                        ENDIF

                     ELSEIF (LDTHICK) THEN

                        NLEVVEG=MIN(NLEVV,INT(HC/MINTHICK))
                        NLEVEL=NLEV+NLEVVEG
                        MIN_HC=HC

 102                    CONTINUE

C     LG-           for NRESET=0 the minimum canopy is being calculated for which
C     for all the defined vegetation layers NLEVV, for the logaritmic
C     vertical coordinate system, the layers have at least a thickness
C     with the value MINTHICK

                        IF (NRESET.EQ.0) THEN
                           DO 129 JK=NLEV+1,NLEV+NLEVV
                              
                              JJ=NLEVV+1-(JK-NLEV)

C     LG-               determining coefficients to define the height coordinates within
C     the canopy with an increasing thickness with increasing height
                              
                              CTHICK_CANLAY(JJ)=
     *                             (1.-(FLOAT(NLEVV-JJ)/FLOAT(NLEVV-1)))*ALOG(MIN_HC)+
     *                             (FLOAT(NLEVV-JJ)/FLOAT(NLEVV-1))*ALOG(MINTHICK)

                              CTHICK_CANLAY(JJ-1)=
     *                             (1.-(FLOAT(NLEVV-(JJ-1))/FLOAT(NLEVV-1)))*ALOG(MIN_HC)+
     *                             (FLOAT(NLEVV-(JJ-1))/FLOAT(NLEVV-1))*ALOG(MINTHICK)

C     LG-               vertical coordinate system within the canopy with increasing
C     thickness of the layers with increasing height, the lowest layer
C     has a thickness of MINTHICK being defined in the parameter file
C     parchem.h 

                              ZHP(JK)=EXP(CTHICK_CANLAY(JJ))
                              Z(JK)=(EXP(CTHICK_CANLAY(JJ))+
     &                             MAX(0.,EXP(CTHICK_CANLAY(JJ-1))))/2.
                              IF (JK.EQ.NLEV+NLEVV) THEN 
                                 ZHP(JK+1)=0.
                                 Z(JK)=ZHP(JK)/2.
                              ENDIF

                              IF ((ZHP(JK)-Z(JK))*2.LT.MINTHICK) THEN
                                 MIN_HC=MIN_HC+0.1
                                 GOTO 102
                              ENDIF

 129                       CONTINUE

                        ENDIF

C     LG-           end determining minimum canopy height

 111                    CONTINUE

                        DO 130 JK=NLEV+1,NLEVEL

                           JJ=NLEVVEG+1-(JK-NLEV)

C     LG-             determining coefficients to define the height coordinates within
C     the canopy with an increasing thickness with increasing height

                           CTHICK_CANLAY(JJ)=
     *                          (1.-(FLOAT(NLEVVEG-JJ)/FLOAT(NLEVVEG-1)))*ALOG(HC)+
     *                          (FLOAT(NLEVVEG-JJ)/FLOAT(NLEVVEG-1))*ALOG(MINTHICK)

                           CTHICK_CANLAY(JJ-1)=
     *                          (1.-(FLOAT(NLEVVEG-(JJ-1))/FLOAT(NLEVVEG-1)))*ALOG(HC)+
     *                          (FLOAT(NLEVVEG-(JJ-1))/FLOAT(NLEVVEG-1))*ALOG(MINTHICK)

C     LG-             vertical coordinate system within the canopy with increasing
C     thickness of the layers with increasing height, the lowest layer
C     has a thickness of MINTHICK being defined in the parameter file
C     parchem.h 

                           ZHP(JK)=EXP(CTHICK_CANLAY(JJ))
                           Z(JK)=(EXP(CTHICK_CANLAY(JJ))+
     &                          MAX(0.,EXP(CTHICK_CANLAY(JJ-1))))/2.

                           IF (JK.EQ.NLEVEL) THEN 
                              ZHP(JK+1)=0.
                              Z(JK)=ZHP(JK)/2.
                           ENDIF

                           IF ((ZHP(JK)-Z(JK))*2.LT.MINTHICK) THEN
                              LRESET=.TRUE.
                              NLEVVEG=NLEVVEG-1
                              NLEVEL=NLEVEL-1
                              GOTO 111
                           ENDIF

 130                    CONTINUE

                        IF (LRESET) THEN
                           WRITE(*,'(1a)')
     &                          ' Too thin vegetation layer within log. vert. coord. system' 
                           WRITE(*,'(1a,i3)')
     &                          ' The actual number of vegetation layers is reset to: ',
     &                          NLEVVEG
                           NRESET=NRESET+1
                           WRITE(*,'(1a,f4.1,1a)')
     &                          ' or increase the canopy height to at least: ',MIN_HC,' [m]'
                           IF (NRESET.EQ.2) GOTO 112
                           GOTO 101
                        ENDIF

                     ENDIF

 112                 CONTINUE

C     LG-       end if (LBIOSPH.OR.LVOCEMIS)

                  ENDIF

C     LG-       recalculation of displacement height and surface roughness

C     LG-       call subroutine for calculation of the surface roughness and 
C     displacement heigth from the LAI and canopy height according
C     to Raupach, "Simplified expressions for vegetation roughness.."
C     in Boundary Layer Meteorology, 71, 211-216, 1994
C     
                  CALL CALCZ0D(LAI,HC,Z0MCALC,DISPCALC)
                  DISP=DISPCALC

                  WRITE(*,'(1a,f6.1,1a)')
     *                 ' The rescaled displacement height acc. to Raupach [1994] : ',
     *                 DISP,' [m]'

                  WRITE(*,'(1a,f6.1,1a)')
     *                 ' And the rescaled surface roughness is : ',Z0MCALC,' [m]'

                  WRITE(*,'(1a)')
     *                 ' What must be the surface roughness (< 0.15*HC, in cm !!) ?'
                  READ *,Z0M
                  IF (HC.GT.1.E-10) Z0M=MIN(0.15*HC,Z0M/100.)

                  IF (LVOCEMIS) THEN
                     WRITE(*,'(1a)')
     *                    ' What must be the isoprene emission factor [ug C g-1 hr-1]?'
                     READ *,EMISFACT_VOC(1)
                     WRITE(*,'(1a)')
     *                    ' What must be the monoterpene emission factor [ug C g-1 hr-1]?'
                     READ *,EMISFACT_VOC(2)
                     WRITE(*,'(1a)')
     *                    ' What must be the OVOC emission factor [ug C g-1 hr-1]?'
                     READ *,EMISFACT_VOC(3)
                  ENDIF

 888              CONTINUE
                  IF (LAGS) THEN
                     WRITE(*,'(1a)')
     *                    ' What must be the C3/C4 type [1=C3, 2=C4, 3 for coniferous trees]?'
                     READ *,PC3C4TYPE
                  ENDIF
                  IF (PC3C4TYPE.GT.3.OR.PC3C4TYPE.LT.1) GOTO 888

                  GOTO 110

               ENDIF

            ENDIF

         ENDIF

      ENDIF

C     LG- 20031218+ new feature that in case of resetting the initial surface
C     cover parameters, which is for example done for a site evaluation,
C     that the parameter LRESETINP is set to true. In case that this 
C     parameter is set to true, the user is asked to indicate if he/she
C     wants for a change in the month that the updated surface cover
C     properties (and others such as emissions) are assigned to the 
C     actually applied parameters to constrain the SCM

      IF (NSTEP.EQ.0.AND..NOT.LAGRIAN) THEN
 889     CONTINUE
         WRITE(*,'(1a)')' Do you want to apply the initial surface cover param.',
     *        ' throughout the simulation, also for a changing month (YES=1,NO=0): ?'
         WRITE(*,'(1a)')' (not only surface cover param. but also emissions, etc.)'
         READ *,IRESETINP
         IF (IRESETINP.GT.1.OR.IRESETINP.LT.0) GOTO 889
         IF (IRESETINP.EQ.1) THEN
            LRESETINP=.FALSE.
         ENDIF
      ENDIF

C     LG- The maximum canopy height is set to a maximum being defined in 
C     parchem.h. HCMAX is set to zero for the "big leaf" approach
      
      IF (.NOT.LBIOSPH.AND..NOT.LBULKVEG.AND..NOT.LVEG_MLAY) THEN
         HCMAX=0.
      ELSEIF (LAGRIAN.AND.ILSTYPE.GT.1) THEN
         HCMAX=MAX_CANHEIGHT
         HCCORR=HCMAX
      ENDIF

C     LG- end

C     LG- defining the number of layers for different approaches. For the
C     biosphere mode the maximu number of levels is restricted by the
C     canopy height and the minimal canopy thickness. The radiation profiles
C     are for this specific case then only resolved for these layers since
C     the multilayer dry deposition calculations require these radiation
C     parameters. For the bulkveg and bigleaf model calculations, the 
C     radiation profiles are resolved for the maximum number of canopy 
C     layers in order to use the optimal resolution of the vertical radiation
C     profiles used for the calculations of the VOC emissions! 

      IF (.NOT.LBIOSPH.AND..NOT.LBULKVEG.AND..NOT.LVEG_MLAY) THEN
         NLEVELS=NLEVV
         NLEVVEG=0
      ELSEIF (LBULKVEG) THEN
         NLEVELS=NLEVV
         NLEVVEG=1
      ELSEIF (LVEG_MLAY) THEN
         NLEVELS=NLEVV
         NLEVVEG=NLEVV_ML
      ELSEIF (LBIOSPH.AND.LAGRIAN.AND.ILSTYPE.GT.1) THEN
         NLEVELS=NLEVV
      ELSE
         NLEVELS=NLEVVEG
      ENDIF

C     -- calling LADPROF in which the LAD profiles are assigned
C     to the different biomes/ecosystems, the LAD profile is determined
C     using the optimal vertical resolution within the canopy due the fact
C     that this LAD profile strongly controls the VOC emissions through
C     the radiation regime within the canopy

      print *,'hc,ladprof',hc
      read (*,*)
      IF (HC.GT.1.E-10) THEN
         CALL LADPROF(IPROF,NLEV,NLEVELS,NLEVT,LDTHICK,
     &        LAGRIAN,LBIOSPH,LVEG_MLAY,LVOCEMIS,ZHP,Z,HC,HCCORR,LAD)

         OPEN(UNIT=2,FILE='/data/ganzevl/racmo/output/LAD.out',STATUS='UNKNOWN')
         
         WRITE(2,'(1a)')'Distribution of LAI within canopy (LAI*LAD)'
         IF (LDTHICK) THEN
            WRITE(2,'(1a)')'Logarithmic vertical coord. system'
         ELSE
            WRITE(2,'(1a)')'Equidistant vertical coord. system'
         ENDIF

         IF (IPROF.EQ.1) WRITE(2,'(2a)')
     *        'LAD profile representative. for agric. vegetation'
         IF (IPROF.EQ.2) WRITE(2,'(2a)')
     *        'LAD profile representative. for decideous forest (e.g., Oak)'
         IF (IPROF.EQ.3) WRITE(2,'(2a)')
     *        'LAD profile with most of biomass in the top of the canopy'
         IF (IPROF.EQ.4) WRITE(2,'(2a)')
     *        'Average LAD profile'
         IF (IPROF.EQ.5) WRITE(2,'(2a)')
     *        'LAD profile for Tropical Rainforest, Cuieiras (paper Kruijt 1999)'
         IF (IPROF.EQ.6) WRITE(2,'(2a)')
     *        'LAD profile for Tropical Rainforest, Jaru (paper Kruijt 1999)'
         IF (IPROF.EQ.7) WRITE(2,'(2a)')
     *        'LAD profile, average of Cuieiras and Jaru (paper Kruijt 1999)'
         IF (IPROF.EQ.8) WRITE(2,'(2a)')
     *        'LAD profile for fruit orchard (paper Walton et al. 1997)'

         WRITE(2,'(a5,5a12)')
     &        'Level','RACMO-lev','z[m]','z/h[-]','LAD[-]','LAD*LAI'      

         DO JK=1,NLEVELS
            JJ=NLEV+NLEVELS+1-JK
            WRITE(2,'(i5.5,9x,i3,4f12.3)')JK,JJ,Z(JJ),Z(JJ)/HCCORR,
     &           LAD(JK),LAI*LAD(JK)
         ENDDO
         
         CLOSE(2)
         
      ENDIF

c     -- writing of some surface cover characteristics

      IF (NSTEP.EQ.0) THEN
         OPEN (UNIT=NUNSURF,FILE='/data/ganzevl/racmo/output/surfcov.out',
     *        STATUS='UNKNOWN')
         WRITE(NUNSURF,'(2a)')'Surface cover and some of its ',
     +        'characteristics derived from Olson database and satellite data'
         WRITE(NUNSURF,'(a5,2a8,13(a10),1a)')
     +        'istep','Long','Lat','LAI','LAD-prof.','HC','z0','DISP',
     +        'C5H8-emis','C10-emis','OVOC-emis','DM','SLW',
     +        'NO-class1','NO-class2','C3/C4',' Olson ecosystem type descr.'

C     LG-  opening file to write average surface cover properties

         OPEN (UNIT=NUNSURF2,FILE='/data/ganzevl/racmo/output/surfcov_avg.out',
     *        STATUS='UNKNOWN')
         WRITE(NUNSURF2,'(2a)')'Average surface cover and some of its ',
     +        'characteristics derived from Olson database and satellite data'
         WRITE(NUNSURF2,'(a5,2a8,9(a10))')
     +        'istep','Long','Lat','LAI','HC','z0','DISP',
     +        'C5H8-emis','C10-emis','OVOC-emis','DM','SLW'

      ENDIF

      OLSON=OLSNAME(OLSVEG(ZILON,ZILAT))
      IP=INDEX(OLSON,'  ')      
      OLSON=OLSON(1:IP-1)

      WRITE(NUNSURF,'(I5.5,2f8.2,1x,f9.1,1x,i9,1x,8(f9.2,1x),
     +     3(i9,1x),1a)') 
     +     NSTEP,ZLON,ZLAT,LAI,IPROF,HC,Z0M,DISP,EMISFACT_VOC(1),EMISFACT_VOC(2),
     +     EMISFACT_VOC(3),DM,SLW,INOCLASS(1),INOCLASS(2),PC3C4TYPE,OLSON(1:IP-1)

C     LG- writing the average surface cover parameters, ofcourse only the
C     "contineous" data can be averaged and not the discrete data such as
C     as the emission classes and vegetation classes

      LAI_AVG=LAI_AVG+LAI
      HC_AVG=HC_AVG+HC
      Z0M_AVG=Z0M_AVG+Z0M
      DISP_AVG=DISP_AVG+DISP
      EMISFACT_VOC_AVG(1)=EMISFACT_VOC_AVG(1)+EMISFACT_VOC(1)
      EMISFACT_VOC_AVG(2)=EMISFACT_VOC_AVG(2)+EMISFACT_VOC(2)
      EMISFACT_VOC_AVG(3)=EMISFACT_VOC_AVG(3)+EMISFACT_VOC(3)
      DM_AVG=DM_AVG+DM
      SLW_AVG=SLW_AVG+SLW

      IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN
         IF (NSTEP.GT.0.) THEN
            LAI_AVG=LAI_AVG/NPRINT
            HC_AVG=HC_AVG/NPRINT
            Z0M_AVG=Z0M_AVG/NPRINT
            DISP_AVG=DISP_AVG/NPRINT
            EMISFACT_VOC_AVG(1)=EMISFACT_VOC_AVG(1)/NPRINT
            EMISFACT_VOC_AVG(2)=EMISFACT_VOC_AVG(2)/NPRINT
            EMISFACT_VOC_AVG(3)=EMISFACT_VOC_AVG(3)/NPRINT
            DM_AVG=DM_AVG/NPRINT
            SLW_AVG=SLW_AVG/NPRINT
         ENDIF

         WRITE(NUNSURF2,'(I5.5,2f8.2,1x,f9.1,1x,8(f9.2,1x))') 
     +        NSTEP,ZLON,ZLAT,LAI_AVG,HC_AVG,Z0M_AVG,DISP_AVG,
     +        EMISFACT_VOC_AVG(1),EMISFACT_VOC_AVG(2),EMISFACT_VOC_AVG(2),
     +        DM_AVG,SLW_AVG

         LAI_AVG=0.
         HC_AVG=0.
         Z0M_AVG=0.
         DISP_AVG=0.
         EMISFACT_VOC_AVG(1)=0.
         EMISFACT_VOC_AVG(2)=0.
         EMISFACT_VOC_AVG(3)=0.
         DM_AVG=0.
         SLW_AVG=0.
      ENDIF

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNSURF)
      IF (NSTEP.EQ.NSTOP) CLOSE(NUNSURF2)

C     LG- in order to put the emissions and dry deposition in the layers
C     bordering the surface layer, the LAD profiles are shifted upward
C     if the LAD in one of the vegetation layers equals zero. This
C     is done to deal with the changing canopy height in the Lagragian
C     simulations

      NLEVEL=NLEV

      IF (LBIOSPH) THEN
         NLEVEL=NLEV+NLEVVEG
      ENDIF

      JJ=1

 997  CONTINUE

      IF (LBIOSPH.AND.LAGRIAN.AND.LAD(JJ)*LAI.LE.1.E-10) THEN
         NLEVEL=NLEVEL-1

C     LG-   actual number of vegetation layers

         NLEVVEG=NLEVEL-NLEV
         LAD(JJ)=0.
         IF (NLEVEL.EQ.NLEV) GOTO 998
         JJ=JJ+1
         GOTO 997       
      ENDIF

 998  CONTINUE 
      
      WRITE(NUNMDFL,'(1a,i3)')
     &     ' The actual number of vegetation layers is: ',NLEVVEG

      IF (NSTEP.EQ.0.AND.NLEVVEG.AND.LBIOSPH.AND..NOT.LDTHICK) THEN  
         WRITE(*,'(1a,f6.2,1a)')
     &        ' The thickness of the vegetation layers is',HC/NLEVVEG,' m'
      ENDIF

      WRITE(NUNMDFL,'(1a)')' End reading and scaling land cover data'

      RETURN
      
      END

      SUBROUTINE READVEG(IVEG,OLSVEG)

c ------------------------------------------------------
c     this routine reads the vegetation database with
c     a high resolution and scales this to the to proper 
c     resolution
c ------------------------------------------------------

      IMPLICIT NONE
   
      INTEGER ZILON, ZILAT, IILON, IILAT, I1, I2, J1, J2,
     &        I, J, II, JJ, iip

      PARAMETER (ZILON=2160, ZILAT=1080, IILON=720, IILAT=360)

      INTEGER*1 IVEG(ZILON,ZILAT)

      INTEGER OLSVEG(IILON,IILAT)
      
      REAL ECHRES,XRES,YRES,RESIN,WEIGHT,RESOUT,a

      CHARACTER*40 dir
      
C WP- adding the directory for reading vegetation input

      DIR='../../echam/'
      IIP=INDEX(DIR,' ')

C WP- end

c  -- 0.5 x 0.5 output resolution

      ECHRES=0.5

c  -- The resolution of the input file

      RESOUT=0.16666667

      XRES=ECHRES
      YRES=ECHRES

      RESIN=1./RESOUT

c
c  -- Opening and file which contains the vegetation classification     
c 
      OPEN(UNIT=10,FILE=dir(1:iip-1)//'input/veg/OWE14D.IMG',
     &   FORM='unformatted',STATUS='OLD',RECORDTYPE='STREAM')

      READ(10) IVEG

      CLOSE(10)

      DO J=1,IILAT
        DO I=1,IILON

c    -- definition of area of high resolution data covering
c       the ECHAM input grid

        I1=INT(RESIN*(I-1)*XRES)+1.
        I2=INT(RESIN*(((I-1)*XRES+XRES)+RESOUT))
        J1=INT(RESIN*(J-1)*YRES)+1.
        J2=INT(RESIN*(((J-1)*YRES+YRES)+RESOUT))

        IF (I1.EQ.(RESIN*(I-1)*XRES)) I1=I1+1
        IF (J1.EQ.(RESIN*(J-1)*YRES)) J1=J1+1
        IF(I1.EQ.0) I1=1
        IF(I2.GT.360*RESIN) I2=360*RESIN
        IF(J1.EQ.0) J1=1
        IF(J2.GT.180*RESIN) J2=180*RESIN

        DO II=I1,I2
          DO JJ=J1,J2
             WEIGHT=1.

c         -- weighting calculation of edges/corners

             IF(II-1.LT.RESIN*(I-1)*XRES) 
     &          WEIGHT=(II/RESIN-(I-1)*XRES)*RESIN
             IF(II.GT.RESIN*I*XRES) 
     &          WEIGHT=(I*XRES-(II-1)/RESIN)*RESIN
             IF(JJ-1.LT.RESIN*(J-1)*YRES) 
     &          WEIGHT=(JJ/RESIN-(J-1)*YRES)*RESIN
             IF(JJ.GT.RESIN*J*YRES) 
     &          WEIGHT=(J*YRES-(JJ-1)/RESIN)*RESIN
             IF(II-1.LT.RESIN*(I-1)*XRES.AND.JJ-1.LT.
     &          RESIN*(J-1)*YRES)
     &          WEIGHT=(II/RESIN-(I-1)*XRES)*RESIN*
     &          (JJ/RESIN-(J-1)*YRES)*RESIN
             IF(II.GT.RESIN*I*XRES.AND.JJ-1.LT.RESIN*(J-1)*YRES)
     &          WEIGHT=(I*XRES-(II-1)/RESIN)*RESIN*
     &          (JJ/RESIN-(J-1)*YRES)*RESIN
             IF(II-1.LT.RESIN*(I-1)*XRES.AND.JJ.GT.RESIN*J*YRES)
     &          WEIGHT=(II/RESIN-(I-1)*XRES)*RESIN*
     &          (J*YRES-(JJ-1)/RESIN)*RESIN
             IF(II.GT.RESIN*I*XRES.AND.JJ.GT.RESIN*J*YRES)
     &         WEIGHT=(I*XRES-(II-1)/RESIN)*RESIN*
     &         (J*YRES-(JJ-1)/RESIN)*RESIN

             IF (WEIGHT.EQ.1) OLSVEG(I,J)=IVEG(II,JJ)

           ENDDO
          ENDDO
        ENDDO
      ENDDO

      WRITE(*,'(1a)')' End reading and scaling vegetation type data'

      END

      SUBROUTINE READGVI(IMONTH,GVI05X05,GVIMAX05X05)

c ----------------------------------------------------
c     reading of GVI input data which are derived from
c     Global Ecosystem database CDROM. The fields are already
c     preprocessed and scaled to the 0.5 X 0.5 horizontal 
c     grid resolution in order to reduce the I/O 
c ----------------------------------------------------

      IMPLICIT NONE

      INTEGER IILON, IILAT,IMONTH, I, II, J, K, IP, iip

      PARAMETER (IILON=720, IILAT=360)

      CHARACTER*75 FNAME,FNAME1
      CHARACTER*2 MON(12)    

      INTEGER GVI05X05(IILON,IILAT),GVIMAX05X05(IILON,IILAT)
   
      DATA MON /'01','02','03','04','05','06','07','08','09',
     &  '10','11','12'/

      CHARACTER*40 dir

C WP  adding the directory for reading vegetation input

      DIR='../../echam/'
      IIP=index(dir,' ')
C

c  -- reading input data

      DO II=1,2
       IF (II.EQ.1) THEN

C LG-   reading the preprocessed GVI and maximum GVI fields 

        FNAME=DIR(1:IIP-1)//'input/veg/GVI_AVHRR_05X05_'
        IP=INDEX(FNAME,' ')
        FNAME=FNAME(1:IP-1)//MON(IMONTH)

        OPEN(UNIT=2,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
        READ(2,'(36i3)') ((GVI05X05(I,J),I=1,IILON),J=IILAT,1,-1)
        CLOSE(2)
       ELSE 
        FNAME=
     &  DIR(1:IIP-1)//'input/veg/GVIMAX_AVHRR_05X05'

        OPEN(UNIT=2,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
        READ(2,'(36i3)') ((GVIMAX05X05(I,J),I=1,IILON),J=IILAT,1,-1)
        CLOSE(2)
       ENDIF

      ENDDO
       
      WRITE(*,'(1a)')' End reading GVI data'

      END


      SUBROUTINE CALCVEG(IMONTH,LON,LAT,PARAM,OLSVEG,GVI,GVIMAX,
     &                   LAIMAX,DM,LAI)

c ----------------------------------------------------
c     Calculation of foliar density and LAI from the GVI 
c     data according to Guenther et al. 1995
c ----------------------------------------------------

      IMPLICIT NONE

      INTEGER ZILON, ZILAT, MAXCLASS
      REAL DMMAX

      PARAMETER (ZILON=720, ZILAT=360, MAXCLASS=73,DMMAX=1375.)

      CHARACTER*70 FNAME,FNAME1
      CHARACTER*2 MON(12)    

      INTEGER IMONTH,K,IP,II,LON,LAT,J

      INTEGER GVI(ZILON,ZILAT),GVIMAX(ZILON,ZILAT),OLSVEG(ZILON,ZILAT)

      REAL NPP,G2,DR,SLW,
     &    PARAM(0:MAXCLASS,11),DP,LAIMAX,DM,LAI 
   
      DATA MON /'01','02','03','04','05','06','07','08','09',
     &  '10','11','12'/

c  -- calculation foliar density and LAI

      LAI=0.0
      NPP=PARAM(OLSVEG(LON,LAT),7)
      G2=PARAM(OLSVEG(LON,LAT),8)
      DR=PARAM(OLSVEG(LON,LAT),9)
      SLW=PARAM(OLSVEG(LON,LAT),10)
      DP=DR*NPP

      IF (GVI(LON,LAT).GT.GVIMAX(LON,LAT)) print *,' GVI > GVImax !' 

      IF (GVI(LON,LAT).LT.G2) THEN
        DM=0.0
      ELSE
        DM=MAX(0.,DP*(EXP(ALOG(2.)*((GVI(LON,LAT)-G2)/
     &        MAX(1.,(GVIMAX(LON,LAT)-G2))))-1.))
      ENDIF

      IF (OLSVEG(LON,LAT).GT.0) LAI=DM/MAX(1.E-20,SLW)

C ============================================================================
C LG- modification of the calculation of LAI from the dry matter, in order
C     to reduce the maximum values which seem to be much too large. In order
C     to somehow include the saturation effect, an asymptotic curve has been 
C     applied to reach the maximum value, given by LAIMAX (07-2002). This 
C     cut-off introduces an inconsitency between the foliar density, specific
C     leaf weight and LAI. It still must be adressed if the problem is related
C     to DM (what are maximum estimates??), the specific leaf weight or any
C     other parameter???
C ===========================================================================

C LG- see spreadsheet with different functions to introduce cut-off in 
C     inferred LAI. Now the cut-off is starting at the half the maximum LAI
C     value, taking the minimum of the linear and the exponential value.

C LG- to check shape of function, data are written to file.

C       OPEN(UNIT=3,FILE='/data/ganzevl/racmo/output/LAI_cut-off.out',STATUS='UNKNOWN')
C       DO J=1,28
C       DM=J*50.
C       SLW=125.
C       LAI=DM/SLW
C       IF (LAI.GT.0.5*LAIMAX)
C      &     LAI=MIN(LAI,0.5*LAIMAX*(1.+((DM/SLW-0.5*LAIMAX)/
C      &        (DMMAX/SLW-0.5*LAIMAX))**0.5))  
C       WRITE(3,'(2F8.2)')DM/SLW,LAI
C       ENDDO
C       CLOSE(3)
C       READ (*,*)

      IF (OLSVEG(LON,LAT).GT.0.AND.LAI.GT.0.5*LAIMAX) THEN 
        LAI=MIN(LAI,0.5*LAIMAX*(1.+((DM/SLW-0.5*LAIMAX)/
     &        (DMMAX/SLW-0.5*LAIMAX))**0.5))  
      ENDIF

      END

      SUBROUTINE LADPROF(IPROF,NLEV,NLEVV,NLEVT,LDTHICK,
     &                   LAGRIAN,LBIOSPH,LVEG_MLAY,LVOCEMIS,
     &                   ZHP,Z,HC,HCCORR,LAD)

c ----------------------------------------------------
c     reading of LAD input data taken from the DDIM model.
c     and interpolation from the 21-layer values to the no. of
c     layers used in this study.
c ----------------------------------------------------

      IMPLICIT NONE

      INTEGER NLEV,NLEVV,NLEVT,NLAYER,I,JK,II,J,JJ,JJJ,JJJJ,IIP
      
      INTEGER IPROF

      REAL LAD(NLEVV),ZLAD(NLEVV),ZZLAD(NLEVV),LADINP(21),
     &     RLAYER,SUMLAD,SUMLAD2,HC,HCCORR,ZHP(NLEVT+1),Z(NLEVT),
     &     ZHP_EQUI(NLEVT),Z_INP(21),ZHP_INP(22),ZDZ_INP,ZDZ_EQUI,
     &     ZDZ,ZH(21)

      LOGICAL LAGRIAN,LBIOSPH,LVEG_MLAY,LVOCEMIS,LDTHICK ! mz_lg_20051224+
  
      CHARACTER*70 FNAME
      CHARACTER*40 dir
      CHARACTER*40 DUMMY
      
C WP  adding the directory for reading vegetation input

      DIR='../../echam/'
      IIP=index(dir,' ')
C

c  -- reading input data

c  -- PADPROF4.20 contains an average profile which is 
c     constructed from the other three profiles, whereas 
c     PADPROF1.20 should be used for agricultural crops,
c     PADPROF2.20 for decideous forest and PADPROF3.20 for
c     coniferous forest. The default profile is the average profile

      FNAME=dir(1:iip-1)//'input/veg/PADPROF4.20'
      IF(IPROF.EQ.1) 
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF1.20'
      IF(IPROF.EQ.2) 
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF2.20'
      IF(IPROF.EQ.3) 
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF3.20'

C LG- Observed lad profiles for the Cuieiras and Jaru site, see paper
C     "Turbulence above and within two Amazon rain forest canopies" by
C     Kruijt et al., 1999. These profiles are specifically defined
C     for the canopy heights of the measurement sites. Profile number 
C     7 is the average profile for the two sites. The LAD profiles are
C     defined in m2 m-3 which is this already the absolute amount of
C     biomass for each point. Taking the average of the profiles and
C     multiplying the average with the observed canopy height gives the
C     total amount of biomass, for Guieiras the average LAD is +/- 0.16
C     which yields together with the canopy height of 30-35 m an LAI
C     of 5.5-6. In order to use the profiles, they are normalized by the
C     total intergrated value so that the total amount of LAI does not
C     change 

      IF(IPROF.EQ.5)
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF5.20'
      IF(IPROF.EQ.6)
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF6.20'
      IF(IPROF.EQ.7)
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF7.20'
      IF(IPROF.EQ.8)
     & FNAME=dir(1:iip-1)//'input/veg/PADPROF8.20'

      OPEN(UNIT=2,FILE=FNAME,STATUS='UNKNOWN') 
      IF (IPROF.LT.5) THEN
       J=1
 100   CONTINUE
       READ(2,*,ERR=200,END=300) LADINP(J)
       J=J+1
       GOTO 100
 200   WRITE(*,'(1a)')' ERROR reading LAD input file' 
 300   CONTINUE 
       NLAYER=J-1   

C LG-  determining the normalized height ZH 

       DO J=1,NLAYER
        ZH(J)=(HC/NLAYER*J)/MAX(1.E-10,HC)
       ENDDO

      ELSEIF (IPROF.GE.5) THEN
       READ(2,*)DUMMY
       READ(2,*)DUMMY
       READ(2,*)DUMMY
       J=1
       SUMLAD=0.
 400   CONTINUE
       READ(2,*,ERR=500,END=600) LADINP(J),ZH(J)
       SUMLAD=SUMLAD+LADINP(J)
       J=J+1
       GOTO 400
 500   WRITE(*,'(1a)')' ERROR reading LAD input file' 
 600   CONTINUE 
       NLAYER=J-1    

C LG-  normalising with the summed LAD

       DO J=1,NLAYER
        LADINP(J)=LADINP(J)/SUMLAD
        SUMLAD2=SUMLAD2+LADINP(J)
       ENDDO
       
      ENDIF
      CLOSE(2)

C LG- for a selected number of vegetation layers being larger than the
C     defined number of layers of the Leaf Area Density profiles, the
C     program is stopped. A smaller selection of vegetation layers is
C     is then required or a higher resolution of the Leaf Area Density
C     is required.

      IF (NLAYER.LT.NLEVV) THEN
       WRITE(*,'(3a)')
     &    ' STOP in vegetation.f, the resolution of the Leaf',
     &    ' Area Density profile is smaller then the number of',
     &    ' vegetation layers'
       WRITE(*,'(1a,i4)')
     &    ' Set NLEVV to a maximum of :',MIN(NLAYER,MIN(20,NLEVV))
       STOP
      ENDIF

C LG- in case of simulations with the big leaf approach in the lagragian 
C     mode, the height levels are not defined and thus are calculated here!

      IF (LVOCEMIS.AND..NOT.LBIOSPH.AND..NOT.LVEG_MLAY) THEN
       DO 128 JK=NLEV+1,NLEV+NLEVV

        JJ=NLEVV+1-(JK-NLEV)

C LG-   defining the height for the equidistant coord. system

        ZHP(JK)=(HCCORR/NLEVV)*JJ
	Z(JK)=(HCCORR/NLEVV)*JJ-(1./2.)*(HCCORR/NLEVV)
 128   CONTINUE
      ENDIF

! mz_lg_20051224+ in case of lagrian simulations with lveg_mlay
      IF (LAGRIAN.AND.LVEG_MLAY) THEN
       DO 129 JK=NLEV+1,NLEV+NLEVV

        JJ=NLEVV+1-(JK-NLEV)

C LG-   defining the height for the equidistant coord. system

        ZHP(JK)=(HC/NLEVV)*JJ
	Z(JK)=(HC/NLEVV)*JJ-(1./2.)*(HC/NLEVV)
 129   CONTINUE
      ENDIF
! mz_lg_20051224-

      SUMLAD=0

C LG- determining the LAD profile for the equidistant vertical coord.
C     system from the different resolution for the LAD input data and the
C     equidistant vertical coordinate system

      DO 21 J=1,NLAYER
       JJ=NLAYER+1-J
       IF (J.GT.1) THEN
        ZHP_INP(J)=((ZH(JJ+1)+ZH(JJ))/2.)*HC
       ELSE
        ZHP_INP(J)=HC
       ENDIF
       Z_INP(J)=ZH(JJ)*HC
 21   CONTINUE

      ZHP_INP(NLAYER+1)=0.

      DO 22 J=NLEVV,1,-1
       LAD(J)=0.
 22   CONTINUE

      SUMLAD=0.
      I=1
      II=NLAYER+1-I

      JJJ=0

      DO 23 J=NLEVV,1,-1

        JJ=NLEV+NLEVV+1-J
        ZLAD(J)=0.

C LG-   in case of the canopy height being smaller than the reference
C       height of the half pressure level the calculation of LAI in the 
C       layer is not done and the loop continues to the next layer until 
C       the layer which does designate the real canopy top

        IF (HC.LT.ZHP(JJ))  THEN
         JJJ=JJJ+1
         GOTO 23
	ENDIF

C LG-   in case of the height of the full pressure level of the vertical
C       coord. system being smaller than that of the input data

 222    CONTINUE
        ZDZ_INP=ZHP_INP(I)-ZHP_INP(I+1)	

        IF (ZHP(JJ+1).LT.ZHP_INP(I).AND.ZHP(JJ+1).GE.ZHP_INP(I+1)) THEN
         ZLAD(J)=ZLAD(J)+LADINP(II)*((ZHP_INP(I)-ZHP(JJ+1))/ZDZ_INP)                   
	ELSEIF (ZHP(JJ+1).LT.ZHP_INP(I).AND.ZHP(JJ).LT.ZHP_INP(I)) THEN
         ZLAD(J)=ZLAD(J)+LADINP(II)*((ZHP(JJ)-ZHP_INP(I+1))/ZDZ_INP)     
         I=I+1
	 II=II-1
	 IF (II.GT.1) GOTO 222
        ELSE
         ZLAD(J)=ZLAD(J)+LADINP(II)*
     &     MIN(1.,((ZHP_INP(I)-ZHP(JJ+1))/ZDZ_INP))
         I=I+1
	 II=II-1
	 IF (II.GT.1) GOTO 222
        ENDIF

        IF (J.EQ.1.AND.Z_INP(I).LT.ZHP(JJ)) 
     &   ZLAD(J)=ZLAD(J)+LADINP(II)  

 223    CONTINUE

        LAD(J+JJJ)=ZLAD(J)

c    -- The total sum of LAD of all the layers must be one!

        SUMLAD=SUMLAD+LAD(J+JJJ)

 23   CONTINUE

C LG- determining the LAD profile for the logarithmic vertical coord.
C     system from the different thicknesses for this system and the
C     equidistant vertical coordinate system

      IF (LDTHICK) THEN
       DO 24 J=NLEVV,1,-1
        LAD(J)=0.
        JJ=J+NLEV
        JJJJ=NLEVV+1-J
	ZHP_EQUI(JJ)=(HC/NLEVV)*JJJJ
 24    CONTINUE

       ZHP_EQUI(NLEV+NLEVV+1)=0.

       SUMLAD=0.
       ZDZ_EQUI=(HC/NLEVV)
       I=NLEV+1
       II=NLEVV+1-(I-NLEV)

       JJJ=0

       DO 25 J=NLEVV,1,-1

        JJ=NLEV+NLEVV+1-J
        ZZLAD(J)=0.

C LG-   in case of the canopy height being smaller than the reference
C       height of the full pressure level the calculation of LAI in the 
C       layer is not done and the loop continues to the next layer until 
C       the layer which does designate the real canopy top

        IF (HC.LE.Z(NLEV+NLEVV+1-J))  THEN
         JJJ=JJJ+1
         GOTO 25
	ENDIF

C LG-   in case of the height of the half pressure level of the logaritmic
C       coord. system being smaller than that of the equidistant height

        IF (ZHP(JJ+1).GT.ZHP_EQUI(I+1)) THEN
         ZZLAD(J)=ZZLAD(J)+ZLAD(II)*((ZHP(JJ)-ZHP(JJ+1))/ZDZ_EQUI)
        ELSE

C LG-    in case of the height of the half pressure level of the logaritmic
C        coord. system being larger than that of the equidistant height

         IF (ZHP(JJ+1).LT.ZHP_EQUI(I+1).AND.ZHP(JJ+1).GE.ZHP_EQUI(I)) THEN
          ZZLAD(J)=ZZLAD(J)+ZLAD(II)*
     &       ((ZHP_EQUI(I)-ZHP(JJ+1))/ZDZ_EQUI)
	  II=II-1
	 ENDIF	 
 225     CONTINUE
         ZZLAD(J)=ZZLAD(J)+ZLAD(II)*
     &    ((MIN(ZHP(JJ),ZHP_EQUI(I))-MAX(ZHP(JJ+1),ZHP_EQUI(I+1)))/
     &      ZDZ_EQUI)

         IF (II.EQ.1) GOTO 226 

         IF (ZHP(JJ+1).LT.ZHP_EQUI(I+1)) THEN
	  I=I+1
          II=II-1
          GOTO 225
	 ENDIF
        ENDIF

 226    CONTINUE

        LAD(J+JJJ)=ZZLAD(J)

c    -- The total sum of LAD of all the layers must be one!

        SUMLAD=SUMLAD+LAD(J+JJJ)

 25    CONTINUE

      ENDIF

C LG- setting the LAD values to zero for the lowest layers within the
C     canopy for canopy height smaller than the maximum canopy height

      DO J=1,JJJ
	LAD(J)=1.E-20
      ENDDO

C LG- for a summed LAD larger than zero, implying that there is a defined
C     canopy height and thus vegetation, the summed LAD over the whole 
C     canopy must be 1 !

      IF (SUMLAD.GT.0.AND.SUMLAD.LT.0.99.OR.SUMLAD.GT.1.01) THEN
       WRITE(*,'(1a,F4.2)')
     &    ' Total int. LAD (must be 1) = ',SUMLAD
       WRITE(*,'(1a)')
     &     ' The LAD profile is rescaled with the intergrated LAD'
       IF (.NOT.LAGRIAN) THEN
        WRITE(*,'(1a)')' Enter to continue'
        READ (*,*)
       ENDIF	
       SUMLAD2=0.
       DO J=1,NLEVV
        LAD(J)=LAD(J)/SUMLAD
       ENDDO
      ENDIF
          
      END

