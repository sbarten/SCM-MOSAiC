      SUBROUTINE EMISS (NSTEP,NSTOP,NPRINT,NLEV,
     *     NLEVEL,PM,PMLOC,PRHOA,PTMST,PDP,PP)

C__________________________________________________________________
C     
C     ** EMISS ** for emission calculation of species
C     
C     Interface: *EMISS* is called from *PRECHEM*
C__________________________________________________________________
C     

      IMPLICIT NONE

      REAL ZAVO,G,PRCNT

      PARAMETER (ZAVO=6.022045E23,G=9.80665)      
      
C     LG- 

      REAL XMN,XMC,XMNO,XMCO,XMISOP,XMMTERP,XMSQTERP,XMS,XMRAD,
     &     XMCO2,XMNH3

C     LG- moler mass of some species. For the monoterpenes this is based on
C     the molecular structure C10H16 whereas for sesquiterpenes this is 
C     reflecting C15H24 (diterpenes=4*C5H8=C20H32)

      PARAMETER (XMN=14.,XMC=12.,XMNO=32.,XMCO=28.,XMISOP=68.,
     &     XMMTERP=136.,XMSQTERP=204.,XMS=32.,XMRAD=222.,
     &     XMCO2=44.,XMNH3=17.)       

      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'
      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

! mz_lg_20050118+ added

      INCLUDE 'daycent_SCM.h'

! mz_lg_20050118-

      REAL PRHOA(NLON,NLEVT),PDP(NLON,NLEVT),PP(NLON,NLEVT)
      REAL PM(NLON,NLEVT,NTRAC),PMLOC(NLON,NLEVT,NG3X)
      REAL ZCOLMSS(NLON)

C     LG- declaration of some more parameters

      REAL XNOEMIS(NLON,NLEVT),EMIS(NLON,NTRACT),
     &     EMISFACT(NSPEC_EMISVOC),ISOPEM_INT,MTERPEM_INT,OVOCEM_INT,
     &     NO,NOMHR

      REAL XNOXEM(NLONT30,NLATT30,3),COEM(NLONT30,NLATT30)
      SAVE XNOXEM,COEM

      REAL DMSEM(NLONT30,NLATT30,1),SO2EM(NLONT30,NLATT30,2)
      REAL*8 SVOLC(NLONT30,NLATT30),IVTOP(NLONT30,NLATT30)
      SAVE DMSEM,SO2EM,SVOLC,IVTOP

C     LG- more declarations

      INTEGER IW,NSTEP,JL,JR,JK,JT,JJ,NSTOP,NPRINT,
     &     NLEV,NLEVEL,DLEV,IP,IIP,INO_CLASS
      
      CHARACTER*70 FNAME
      CHARACTER*40 DIR
      CHARACTER*2 MON(12)
      
      REAL FFAC,FPG3,PTMST,EMIFAC,EMRECALC,COEMX,XNOXEMT,EMTOT,
     &     MODROW,OLDC,ZEMNOFL,TOTNOFL,COMOLEC,EMNO,EMISOP,EMMATERP,
     &     EMMBTERP,EMMTTERP,EMSQTERP,SFAC1,SFAC2,BBNS,SEMTOT,EMIFAC2,
     &     ZEMI1,ZEMI2,ZEMI3,EMTOT10,EMTOT11,EMTOT15,EMTOT16,RADEMX,
     &     CULT_INDEX,FERTFLX,SO2MOLEC,CO2EMX,NH3EMX,CH3CLEMX,CHCL3EMX,
     &     FT,HONOEMX,NOXEMX,JHNO3S_HONO,JHNO3S_NOX,HNO3S,EMHNO3_ACCUM,
     &     DDHNO3_ACCUM,CH3OHEMX,ACETEMX,ACETALDEMX
      
      SAVE HNO3S ! saving the accumulated mass of HNO3 deposited at the
                 ! the surface

      REAL CRF_AVG,CULT_AVG,NOMOLEC_AVG,NOEM_AVG,FERTFLX_AVG,
     &     XNOXEMT_AVG,CRFFLX_AVG
      REAL EMISFACT_AVG(NSPEC_EMISVOC),ISOPEM_INT1_AVG,ISOPEM_INT2_AVG,
     &     ISOPMOLEC1_AVG,ISOPMOLEC2_AVG,MTERPMOLEC_AVG,OVOCMOLEC_AVG

! mz_lg_20050602+ extra declarations for testing CH4 emission calculation

!      REAL zxtm1,pxtm1_ch4,pxtte_ch4,ch4_emflux,zdensair,zdz

C LG- definition of month

      DATA MON /'01','02','03','04','05','06','07','08','09',
     &  '10','11','12'/

! mz _lg_20060312+ for reading VOC emission output file

      INTEGER NUNVOCEMINP, NNSTEP
      REAL VOCINP
      CHARACTER*255 DUMMY,LLDATLTIME
      LOGICAL LREAD_VOCEMIS

      PARAMETER (LREAD_VOCEMIS=.FALSE.,NUNVOCEMINP=222)

C     LG- Default: if the biosphere is considered than the biogenic 
C     NO emissions and VOC emissions are used

      IF (LBIOSPH.AND..NOT.LNOEMIS.OR.LBIOSPH.AND..NOT.LVOCEMIS) THEN
         IF (NSTEP.EQ.0)
     &        WRITE(NUNMDFL,'(1a)')
     &        'LBIOSPH=.TRUE., -> LNOEMIS and LVOCEMIS set to TRUE.'
         LNOEMIS=.TRUE. ! LDAYCENT=.TRUE. (alternative option)
         LVOCEMIS=.TRUE.
      ENDIF
      
      CALL RESETR (ZCOLMSS,NLON,0.)

C     LG- reseting the parameters EMIS, these parameters are used in the
C     subroutines calcxtmzz.f and bulkveg.f. This parameter should only 
C     being assigned a value for biogenic emissions within the canopy or
C     through the interface of the surface and atmospheric surface layer,
C     thus not for example for the anthropogenically produced SO2 and 
C     emitted over the whole surface layer or even in the layer aloft
      
      CALL RESETR (EMIS_TOT,NLON*NTRACT,0.)
      CALL RESETR (EMIS_VEG,NLON*NTRACT,0.)
      CALL RESETR (EMIS_WS,NLON*NTRACT,0.)
      CALL RESETR (EMIS_BS,NLON*NTRACT,0.)
      CALL RESETR (EMIS_SN,NLON*NTRACT,0.)
      CALL RESETR (EMIS_WAT,NLON*NTRACT,0.)

C     LG- writing the emission scaling factors for the different species to
C     the screen to avoid possible misinterpretations

      IF (NSTEP.EQ.0) THEN

         WRITE(*,'(1a)')
     &        ' The emission scaling factors have the values: '
         WRITE(*,'(1a,f5.2)')' CO	: ',FEMCO
         WRITE(*,'(1a,f5.2)')' NOx	: ',FEMNOX
         WRITE(*,'(1a,f5.2)')' NO	: ',FEMNO
         WRITE(*,'(1a,f5.2)')' ISOP	: ',FEMISOP
         WRITE(*,'(1a,f5.2)')' MATERP	: ',FEMMATERP
         WRITE(*,'(1a,f5.2)')' MBTERP	: ',FEMMBTERP
         WRITE(*,'(1a,f5.2)')' MTTERP	: ',FEMMTTERP
         WRITE(*,'(1a,f5.2)')' OVOC	: ',FEMOVOC
         WRITE(*,'(1a,f5.2)')' SO2	: ',FEMSO2
         WRITE(*,'(1a,f5.2)')' DMS	: ',FEMDMS
         WRITE(*,'(1a,f5.2)')' RADON	: ',FEMRAD
         WRITE(*,'(1a,f5.2)')' CO2	: ',FEMCO2
         WRITE(*,'(1a,f5.2)')' NH3	: ',FEMNH3
         WRITE(*,'(1a,f5.2)')' CH3CL	: ',FEMCH3CL
         WRITE(*,'(1a,f5.2)')' CHCL3	: ',FEMCHCL3
         WRITE(*,'(1a,f5.2)')' CH3OH	: ',FEMCH3OH
         WRITE(*,'(1a,f5.2)')' Acetone	: ',FEMACET
         WRITE(*,'(1a,f5.2)')' Acetald.	: ',FEMACETALD
         WRITE(*,'(1a)')' Enter to continue'
         READ (*,*)
      ENDIF

C     LG- end

      WRITE(NUNMDFL,*)'Start EMISS.f'
      
C     LG- setting the switch for considering the emission within the
C     vertical diffusion routines, default the parameter IDD_TURB is
C     set to one for all the tracers (April 2000, only option for LVEG_MLAY
C     is TRUE) 

      IF ((LVEG_MLAY.OR.LSNOW_MLAY).AND.LXTVDIFF.AND.LEM_TURB) THEN ! mz_lg_20051218+
         DO JT=1,NTRAC
            IEM_TURB(JT)=1
         ENDDO
      ENDIF

C     WP  adding the directory for reading vegetation input

      DIR='../../echam/'
      IIP=INDEX(DIR,' ')
C     

C     --- 
C     1=anthr  2=biogene     3=biomass burning
C     

C     NOx in molecules m-2 s-1 (xnoxem)
C     CO  in molec CO m-2 s-1 (coem)
C     SO2 in molec S m-2 s-1 (so2em)
C     DMS in molec S m-2 s-1 (dmsem)
C     SO2 volcanic in kg S m-2 s-1

C     LG- radon being added

      DO JL=1,NLON

C       RADON in atom m-2 s-1, for the emission rate see the paper by 
C       Trumbore et al., 1990, the emission rate is representrative for 
C       tropical rain forest soils (22-12-98), solely emissions from the
C       soils.

        RADEM=FEMRAD*(1.-WTFRAC(1))*0.30*1.E4

C LG-   20031218+ added an estimate of the HONO and NO2 emission flux related
C       to the photo-dissociation of HNO3 accumulated at the surface through
C       dry deposition (Paper Zhou et al., GRL, 2003/2004). The parameter
C       FSLSUM is the summated fraction of sunlit leaves calculated from 
C       the FSL and the LAD in each layer and then finally scaled here with
C       LAI to get the total surface area exposed directly to sunlight.

        JHNO3S_HONO=2.5E-5 ! Table 2: 50% humidity
        JHNO3S_NOX=2.2E-5  ! Table 2: 50% humidity

        IF (NSTEP.EQ.0) THEN
          IF (JHNO3S_HONO.GT.0.) THEN
	    PRINT *,'EMISS line 192: biogenic emissions of HONO/NOx !'
          ELSE
            PRINT *,'EMISS line 192: deactivated emissions of HONO/NOx !'
          ENDIF
	  PRINT *,'ENTER to continue'
	  READ (*,*)
	ENDIF

C LG-   accumulated amount of HNO3 at surface

        IF (RAINTM(JL).GT.1.E-10) THEN

C LG-     It is assumed that the all the accumulated HNO3 at that part of
C         surface that interceps the rain is removed. The parameter CTHRU,
C         determined in ec4_surf as a function of the relative contribution
C         of convective rainfall to the total rainfall, is applied to
C         estimate this HNO3 removal. The effect of rainfall interception is
C         already partly accounted for by only calculating the emission flux
C         for the dry vegetation fraction but it can be assumed that the
C         rainfall is falling randomly within the grid square at an area that
C         was dry previously. To avoid the contineous accumulation of HNO3
C         at the dry vegetation fraction, whenever it rains it is assumed 
C         that a for a fraction resembling CTHRU the HNO3 is washed off  

	  HNO3S=(1.-CTHRU)*HNO3S  
	ELSE
	  HNO3S=HNO3S+DDHNO3*PTMST
        ENDIF
        HONO_JHNO3SEM=(1.-WSFRAC(JL))*  ! no emission from wet surface
     &         VEGFRAC(JL)*LAI*FSLSUM*JHNO3S_HONO*HNO3S
        NOX_JHNO3SEM=(1.-WSFRAC(JL))*
     &         VEGFRAC(JL)*LAI*FSLSUM*JHNO3S_NOX*HNO3S

! mz_lg_20051215+ source of NOx through HNO3 photolysis< is this a real
!       feature also in the snow-pack???
        IF (LSNOW_MLAY.AND.SNFRAC(JL).GT.0.) THEN
          IF (NSTEP.EQ.0.AND.JHNO3S_HONO.GT.0.) THEN
            PRINT *,'EMISS line 241: emissions of HONO/NOx from snow !'
	    PRINT *,'ENTER to continue'
	    READ (*,*)
          ENDIF
          HONO_JHNO3SEM=SNFRAC(JL)*JHNO3S_HONO*HNO3S
          NOX_JHNO3SEM=SNFRAC(JL)*JHNO3S_NOX*HNO3S
        ENDIF
! mz_lg_20051215-

        EMHNO3_ACCUM=EMHNO3_ACCUM+HONO_JHNO3SEM+NOX_JHNO3SEM
        DDHNO3_ACCUM=DDHNO3_ACCUM+MAX(1.E-20,DDHNO3*PTMST)
      ENDDO

C     LG- end

C     LG- reading of CO and NOx emission files for the T30 resolution

      IF (NSTEP.EQ.0.OR.(LCHMONTH.AND.LRESETINP)) THEN

C LG-    reading in total CO emissions

         FNAME=dir(1:iip-1)//'input/emis/CO_total_T30.'
         IP=INDEX(FNAME,' ')
         FNAME=FNAME(1:IP-1)//MON(IMON)

         WRITE(*,'(2a)')' emiss.f: Reading emission files, e.g., CO: ',
     *     FNAME

         OPEN(UNIT=NUNCOEM,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNCOEM,*) ((COEM(JL,JR),JL=1,NLONT30),JR=1,47,2),
     *        ((COEM(JL,JR),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNCOEM)

C LG-    reading in anthropogenic NOx emissions

         FNAME=dir(1:iip-1)//'input/emis/NO_anthr_T30.00'
         OPEN(UNIT=NUNNOXANT,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNNOXANT,*) ((XNOXEM(JL,JR,1),JL=1,NLONT30),JR=1,47,2),
     *        ((XNOXEM(JL,JR,1),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNNOXANT)

C LG-    reading in biomass burning NOx emissions

         FNAME=dir(1:iip-1)//'input/emis/NO_biobur_T30.'
         IP=INDEX(FNAME,' ')
         FNAME=FNAME(1:IP-1)//MON(IMON)
         OPEN(UNIT=NUNNOXBB,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNNOXBB,*) ((XNOXEM(JL,JR,2),JL=1,NLONT30),JR=1,47,2),
     *        ((XNOXEM(JL,JR,2),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNNOXBB)

C LG-    reading in biogenic NOx emissions

         FNAME=dir(1:iip-1)//'input/emis/NO_biog_T30.'
         IP=INDEX(FNAME,' ')
         FNAME=FNAME(1:IP-1)//MON(IMON)
         OPEN(UNIT=NUNNOXBIO,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNNOXBIO,*) ((XNOXEM(JL,JR,3),JL=1,NLONT30),JR=1,47,2),
     *        ((XNOXEM(JL,JR,3),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNNOXBIO)

C     LG- reading of sulfur emission files

C     SO2 anthr.: index 2: surface layer  index 1: 100m layer
C     note: in the filenames it's just the other way around
C     surface level contains also DMS

C     LG- A good reference to the the biogenic emission of sulfur (tropical forest) 
C     is the paper by Andreae et al., JGR, 95, 16813-, 1990 (DMS, H2S, MeSH), 
C     DMS is emitted by both the soil and the vegetation and has a relatively 
C     long chemical lifetime. There is a continueous soil source
C     whereas Andreae et al. scale the vegetation DMS emission flux using
C     the isoprene emission flux as a function of the radiation and 
C     temperature. 

C LG-    reading in anthropogenic SO2 emissions

         FNAME=dir(1:iip-1)//'input/emis/SO2_level1_T30.'
         IP=INDEX(FNAME,' ')
         FNAME=FNAME(1:IP-1)//MON(IMON)
         OPEN(UNIT=NUNSO2L1,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNSO2L1,*) ((SO2EM(JL,JR,2),JL=1,NLONT30),JR=1,47,2),
     *        ((SO2EM(JL,JR,2),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNSO2L1)

C LG-    reading in anthropogenic SO2 emissions

         FNAME=dir(1:iip-1)//'/input/emis/SO2_level2_T30.'
         IP=INDEX(FNAME,' ')
         FNAME=FNAME(1:IP-1)//MON(IMON)
         OPEN(UNIT=NUNSO2L2,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNSO2L2,*) ((SO2EM(JL,JR,1),JL=1,NLONT30),JR=1,47,2),
     *        ((SO2EM(JL,JR,1),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNSO2L2)

C LG-    reading in DMS emissions

         FNAME=dir(1:iip-1)//'input/emis/DMS_level1_T30.'
         IP=INDEX(FNAME,' ')
         FNAME=FNAME(1:IP-1)//MON(IMON)
         OPEN(UNIT=NUNDMSEM,FILE=FNAME,FORM='formatted',STATUS='unknown')
         READ (NUNDMSEM,*) ((DMSEM(JL,JR,1),JL=1,NLONT30),JR=1,47,2),
     *        ((DMSEM(JL,JR,1),JL=1,NLONT30),JR=48,2,-2)
         CLOSE (NUNDMSEM)

C     LG-   opening of unformatted CRAY input file

         OPEN(UNIT=2,FILE=dir(1:iip-1)//'input/emis/Volcano_spiro',
     *        FORM='unformatted',STATUS='unknown',RECORDTYPE='stream')

         READ (2) SVOLC,IVTOP

         CLOSE(2)

C     LG-   end

c     total volcanic 7.8 Tg S yr-1
         DO JR=1,NLATT30
            DO JL=1,NLONT30
               SVOLC(JL,JR)=SVOLC(JL,JR)*7.8*1.013
            ENDDO
	 ENDDO

C C LG-  writing of emission fields
C     
C       OPEN (UNIT=NUNEMIS,
C      &    FILE='/data/ganzevl/racmo/output/COemis_T30.out',FORM='FORMATTED',
C      &          STATUS='UNKNOWN')
C       WRITE(NUNEMIS,'(1a,i4)')
C      &  'emission rates for the month: ',IMON
C 
C c  --  writing of the emission fields to check the global fields
C     
C       WRITE(NUNEMIS,'(96e10.3)') 
C      &   ((COEM(JL,JR),JL=1,NLONT30),JR=1,NLATT30,2)
C       WRITE(NUNEMIS,'(96e10.3)') 
C      &   ((COEM(JL,JR),JL=1,NLONT30),JR=NLATT30,2,-2)
C     
C       CLOSE(NUNEMIS)

      ENDIF

C     -- Calling of routine in which the VOC emission is calculated. 
C     The routine is only called if the hydrocarbon chemistry is 
C     used. The subroutine calculates the isoprene and monoterpenes
C     emission as a function of the radiation intensity, temperature 
C     and the canopy structure, which determines the radiation 
C     extinction profile. For the bulk approach, the total
C     summed emission fluxes are calculated and emitted in the
C     surface layer whereas for the calculations, explictly considering
C     the biosphere, the emission fluxes are in each of the distinguished 
C     layers within the canopy. The emission rate (ISOP) is in [kg C m-2 s-1]
      
      DO JK=1,NLEVV
         ISOPEM(JK)=0.
         MTERPEM(JK)=0.
         OVOCEM(JK)=0.
      ENDDO
      ISOPEM_INT=0. 
      MTERPEM_INT=0. 
      OVOCEM_INT=0. 

      IF(LVOCEMIS) THEN
         IF (.NOT.LVOCEMIS_MEGAN) THEN  ! mz_lg_20051229+ modified
           CALL VOCEMIS(NSTEP,NSTOP,NPRINT,PTMST,EMISFACT,ISOPEM_INT,
     &                  MTERPEM_INT,OVOCEM_INT)
         ELSE
           CALL VOCEMIS_MEGAN(NSTEP,NSTOP,NPRINT,PTMST,EMISFACT,ISOPEM_INT,
     &                        MTERPEM_INT,OVOCEM_INT)
         ENDIF

c     --   recalculation of emission flux to differents units for
c     comparison and writing of data to screen/file, 1.E3 it to 
c     recalculate from kg to g, 1/XMC to recalculate to mol C and the term 
c     1/5 is to correct for the 5 C molecules. In order to recalc from mol 
c     isoprene cm-2 s-1 to molecules its multiplied with the avogadro number

         ISOPMOLEC=FEMISOP*ISOPEM_INT*
     *        1.E3*(1./XMC)*(1./5.)*ZAVO

         WRITE(NUNMDFL,'(1a,f8.3)')
     *        ' The bulk isoprene flux in 10**11 molec. cm-2 s-1 is:',
     *        1.E-11*1.E-4*ISOPMOLEC

         MTERPMOLEC=(FEMMATERP+FEMMBTERP+FEMMTTERP)*MTERPEM_INT*
     *        1.E3*(1./XMC)*(1./10.)*ZAVO

         WRITE(NUNMDFL,'(1a,f8.3)')
     *        ' The bulk monoterpene flux in 10**11 molec. cm-2 s-1 is:',
     *        1.E-11*1.E-4*MTERPMOLEC

         OVOCMOLEC=FEMOVOC*OVOCEM_INT*
     *        1.E3*(1./XMC)*(1./15.)*ZAVO ! C15 has been assumed to represent
                                ! OVOC (sesquiterpenes?)

         WRITE(NUNMDFL,'(1a,f8.3)')
     *        ' The bulk OVOC flux in 10**11 molec. cm-2 s-1 is:',
     *        1.E-11*1.E-4*OVOCMOLEC

      ENDIF

C     -- Calling of routine in which the NO emission is read/
C     calculated. The NO is emitted in the model's lowest layer

      IF (LNOEMIS) THEN
         CALL NOXEMIS_YL95(NSTEP,NSTOP,INO_CLASS,CULT_INDEX,FERTFLX,FT) 

c     --   recalculation of emission flux to differents units for
c     comparison and writing of data to screen/file, see the isoprene
c     emission recalculation for the different terms

         NOMOLEC=CRF_YL95*FEMNO*NOEM*
     *        1.E3*(1./XMN)*(1./1.)*ZAVO !molec. m-2 s-1
         NOMHR=CRF_YL95*FEMNO*NOEM*1.E9*3600.

         WRITE(NUNMDFL,'(1a,f8.3)')
     *        ' The biog. NO flux in 10**9 molec. cm-2 s-1 is:',
     *        1.E-9*1.E-4*NOMOLEC
         IF (NSTEP.EQ.0)
     *        WRITE(NUNMDFL,'(1a,f8.3)')
     *        ' The init. biog. NO flux in ug N m-2 hr-1 is:',
     *        NOMHR

C     -- Calling of DayCent model to calculate the the soil-biogenic
C     NO/N2O emission fluxes. The model is only called for the first
C     and last timestep, for a change in the day and for the last timestep
C     of the month, when all the daily accumulated statistics are available.
C     The actual calculation of the N fluxes is done after collecting 
C     weather information for one week which gives a one-week delay in the
C     N fluxes (and other relevant parameters). When the fluxes will be 
C     used as input to the SCM then the frequency of calling the routines 
C     must be changed. 

      ELSEIF (LDAYCENT) THEN
        IF (NSTEP.EQ.0.OR.NSTEP.EQ.NSTOP.OR.LCHDAY.OR.LEOMONTH) THEN 

          ! mz_lg_20050109+ added year and change of year
          CALL DAYCENT(NSTEP,NSTOP,IYEAR,IMON,JDAY,NDAY_L, 
     &                 LCHDAY,LCHWEEK,LCHMONTH,LCHYEAR,LEOMONTH)

          NOMOLEC=REAL(NOflux_dc,dp)*  ! NOflux_dc is in gN m-2 day-1
     &        (1./XMN)*(1./(60.*60.*24.))*ZAVO ! molec. m-2 s-1
          WRITE(*,'(1a,f6.2)')
     &    ' emiss.f: DAYCENTs NO emission flux in 1e-14 molec m-2 s-1: ',
     &     1e-14*NOMOLEC
        ENDIF

! mz_lg_20050125+ also assigning the flux when the value if not updated
        NOMOLEC=REAL(NOflux_dc,dp)*  ! NOflux_dc is in gN m-2 day-1
     &        (1./XMN)*(1./(60.*60.*24.))*ZAVO ! molec. m-2 s-1
        NOEM=REAL(NOflux_dc,dp)*1e-3*(1./(60.*60.*24.))  ! [NOEM]=kgN m-2 s-1

! mz_lg_20050406+ using DayCent's LAI in the model

        IF (LDAYCENT_LAI) THEN
	  LAI=REAL(lai_dc,dp)
          IF (NSTEP.EQ.0) THEN
	    PRINT *,'EMISS line 476: NOTE that DayCents LAI is being used!'
  	    PRINT *,'ENTER to continue'
	    READ (*,*)
	  ENDIF
        ENDIF

      ENDIF

C     LG- assigning the CO2 emission/respiration flux, which is a value
C     typical for tropical forests, according to a paper by Malhi et al, 
C     JGR 103, 31,593-31612, 1998. They reported a value for the respiration
C     of 6-8 umol CO2 m-2 s-1. The CO2 emission flux is only considered for
C     LAGS is TRUE, using the physiological model to describe CO2 assimilation

      IF (LAGS) THEN
         CO2MOLEC=(VEGFRAC(1)+BSFRAC(1)+WSFRAC(1))*(7.E-6)*ZAVO

         IF (NSTEP.EQ.0.) THEN  
            IF (VEGFRAC(1)+BSFRAC(1)+WSFRAC(1).LT.1E-10) THEN
               WRITE(*,'(1a,f8.3)')
     *              ' Location over water, the flux is corrected for land/sea mask:'
               WRITE(*,'(1a,f8.3)')
     *              ' The overland CO2 flux in 10**15 molecules cm-2 s-1 is:',
     *              1.E-4*FEMCO2*(7.E-6)*ZAVO*1E-15
            ELSE
               WRITE(*,'(1a,f8.3)')
     *              ' The CO2 flux in 10**15 molecules cm-2 s-1 is:',
     *              1.E-4*FEMCO2*CO2MOLEC*1E-15
            ENDIF
         ENDIF
      ENDIF      

C     LG- assigning the NH3 emission flux, there is as well a soil emission flux
C     as an antropogenic flux, largely controlled by the land use. The 
C     antropogenic flux has to be considered as an source of NH3 by advection
C     from the source regions, similar to NOx, CO etc

      IF (LNH3CHEM) THEN

C     LG-  September, 2001, a soil emission flux representative for the tropical
C     savanna has been assigned (see paper by Bouwman et al., 1997, 
C     Global Biogeochemical cycles, table 10). This has to be redefined in 
C     future. A possible ecosystem dependent emission factor, as presented
C     in the paper by Bouwman can be introduced. This should be consistent 
C     with other emission fluxes such as the CO2 respiration flux, linking
C     it to the NPP and soil decomposition. A flux of 3.5 g N m-2 yr-1 is 
C     recalculed to mol. NH3 m-2 s-1 by dividing it by 365*24*60*60 and then
C     by 14 to arrive at mol., the multiplication with ZAVO gives the molec.

         NH3MOLEC=(VEGFRAC(1)+BSFRAC(1)+WSFRAC(1))*(7.9E-9)*ZAVO

         IF (NSTEP.EQ.0.) THEN  
            IF (VEGFRAC(1)+BSFRAC(1)+WSFRAC(1).LT.1E-10) THEN
               WRITE(*,'(1a,f8.3)')
     *              ' Location over water, the flux is corrected for land/sea mask:'
               WRITE(*,'(1a,f8.3)')
     *              ' The overland NH3 emission flux in 10**11 molecules cm-2 s-1 is:',
     *              1.E-4*FEMNH3*(7.9E-9)*ZAVO*1E-11
            ELSE
               WRITE(*,'(1a,f8.3)')
     *              ' The NH3 soil emission flux in 10**11 molecules cm-2 s-1 is:',
     *              1.E-4*FEMNH3*NH3MOLEC*1E-11
            ENDIF
         ENDIF
      ENDIF   

C     LG- assigning the CH3CL and CHCL3 emission flux, according to Bert
C     Scheeren, estimated using a PBL budget method for the tropics

      CH3CLMOLEC=(VEGFRAC(1)+BSFRAC(1)+WSFRAC(1))*2.5e13
      CHCL3MOLEC=(VEGFRAC(1)+BSFRAC(1)+WSFRAC(1))*1.55e12

      IF (NSTEP.EQ.0.) THEN  
         IF (VEGFRAC(1)+BSFRAC(1)+WSFRAC(1).LT.1E-10) THEN
            WRITE(*,'(1a,f8.3)')
     *           ' Location over water, the flux is corrected for land/sea mask:'
            WRITE(*,'(1a,f8.3)')
     *           ' The overland CH3CL emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMCH3CL*CH3CLMOLEC*1E-11
            WRITE(*,'(1a,f8.3)')
     *           ' The overland CHCL3 emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMCHCL3*CHCL3MOLEC*1E-11
         ELSE
            WRITE(*,'(1a,f8.3)')
     *           ' The CH3CL emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMCH3CL*CH3CLMOLEC*1E-11
            WRITE(*,'(1a,f8.3)')
     *           ' The CHCL3 emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMCHCL3*CHCL3MOLEC*1E-11
         ENDIF
      ENDIF   

C     LG- end

C     LG- assigning the biogenic methanol and acetone emission flux, based on the analysis
C     by Karl et al. 2004, for tropical rainforest Costa Rica

C     LG- note that in table 1 of the manuscript by Karl et al, 2004, they propose a temperature and
C     light dependence of the emissions, that could consequently be included in vocemis.f 

      CH3OHMOLEC=VEGFRAC(1)*2.61e15  ! 0.50 mg m-2 hr-1, CH3OH (mass 32)
      ACETMOLEC=VEGFRAC(1)*1.04e15   ! 0.36 mg m-2 hr-1, CH3COCH3 (mass 58)

      IF (NSTEP.EQ.0.) THEN  
         IF (VEGFRAC(1)+BSFRAC(1)+WSFRAC(1).LT.1E-10) THEN
            WRITE(*,'(1a,f8.3)')
     *           ' Location over water/snow, the flux is corrected for the fractions:'
            WRITE(*,'(1a,f8.3)')
     *           ' The overland CH3OH emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMCH3OH*CH3OHMOLEC*1E-11
            WRITE(*,'(1a,f8.3)')
     *           ' The overland acetone emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMACET*ACETMOLEC*1E-11
         ELSE
            WRITE(*,'(1a,f8.3)')
     *           ' The CH3OH emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMCH3OH*CH3OHMOLEC*1E-11
            WRITE(*,'(1a,f8.3)')
     *           ' The acetone emission flux in 10**11 molecules cm-2 s-1 is:',
     *           1.E-4*FEMACET*ACETMOLEC*1E-11
         ENDIF
      ENDIF   

C     LG- end

C     ---
C     2: Do emission
C     

C     - - -
C     2. Set CH4-groundconcentration to default level
C     each first step of the month
C     

C     LG-   only if the switch LSETCH4 is on then the surface concentrations
C     of CH4 are set to a default value

      IF (LSETCH4) THEN

         MODROW=MOD(IROW,2)

!          ! mz_lg_20050602+ testing the calculation of the CH4 emission flux
! 	 ! as implemented in echam5/messy
!          
!          DO jl=1,NLON
!             zxtm1=1700*1e-9 ! ch4_conc_init(jl,jrow), 1700 ppbv
! 
!             ! mz_lg_20040220+ the strenght of the forcing towards the
!             !     prescribed CH4 surface layer concentrations has been reduced
!             !     by using a relaxation coefficient of 6 hrs since using the
!             !     ztmst resulted in numerical problems in the artic/antartica
!             !     regions probably due to advection problems with the L&R
!             !     advection scheme for the very small gradients. These
!             !     especially occur in those regions where there is no vertical
!             !     mixing and where the zonal mean CH4 concentrations have 
!             !     hardly any gradient.
! 
!             pxtte_ch4=0.
!             pxtm1_ch4=PM(JL,NLEV,ich4)/(ZAVO*PRHOA(JL,NLEV)/ZMAIR)
!             pxtte_ch4=pxtte_ch4-            
!      &           (pxtm1_ch4+pxtte_ch4*ptmst- 
!      &            zxtm1)/21600.
! 
!             print *,'zxtm1',pxtm1_ch4,zxtm1,pxtte_ch4
!             ! mz_lg_20030819+ calculating the surface emission flux from the
!             !     difference in the read-in monthly and zonal mean surface layer
!             !     concentrations and that calculated by the model reflecting
!             !     vertical and horizontal transport and OH oxidation. This
!             !     inferred emission flux can then be compared with previous
!             !     forward/inverse modelling inventories to check if there are
!             !     any serious flaws in the model transport/chemistry. Positive
!             !     fluxes indicate an emission flux (when the CMDL CH4 conc. >
!             !     the model simulated surface layer concentration) whereas
!             !     negative fluxes reflect a possible dry deposition flux
! 
!             ! mz_lg_20030127+ added calculation of the air density and
!             !     thickness of surface layer.
! 
!             zdensair=1e3*PRHOA(JL,NLEV)
!             zdz=64. ! assumed thickness of 64 m
! 
!             ! The recalculation from emission flux from mol CH4 mol-1 air to
!             ! molecules m-2 s-1 is *avo: molecules CH4 mol-1 air, /amd:
!             ! molecules CH4 g-1 air, /1e-3: molecules CH4 kg-1 air,
!             ! *densair [kg m-3]; molecules CH4 m-3, *zdz [m]: molecules CH4 m-2,
!             ! /ztmst [s]; molecules CH4 m-2 s-1
! 
!             ch4_emflux=                                                   
!      &            (zxtm1-                                                  
!      &            (pxtm1_ch4+ pxtte_ch4*ptmst))*                          
!      &            (zavo/(zmair*1e-3))*zdensair*(zdz/ptmst) ! molcules m-2 s-1
!             print *,'ch4_emflux',zdensair,1e-15*ch4_emflux
!          END DO
!          ! mz_lg_20050602- 

         DO JL=1,NLON

            OLDC=PM(JL,NLEV,ich4)

C     IND-PRIND
            IF (MODROW.EQ.1) THEN

C     LG- see NO emissions

C     LG- Industrial scenario

               IF (LINDUS)
     *              PM(JL,NLEV,ich4)=1.77E-6*ZAVO*PRHOA(JL,NLEV)/ZMAIR

C     LG- Pre-industrial scenario

               IF (.NOT.LINDUS)
     *              PM(JL,NLEV,ich4)=.84E-6*ZAVO*PRHOA(JL,NLEV)/ZMAIR
            ELSE

C     LG- Industrial scenario

               IF (LINDUS)
     *              PM(JL,NLEV,ich4)=1.68E-6*ZAVO*PRHOA(JL,NLEV)/ZMAIR

C     LG- Pre-industrial scenario

               IF (.NOT.LINDUS)
     *              PM(JL,NLEV,ich4)=.81E-6*ZAVO*PRHOA(JL,NLEV)/ZMAIR

            ENDIF

C     LG-       determining location

            IW=IWHERE(JL,NLEV)

            BXT(IW,IBEMIS,ich4)=BXT(IW,IBEMIS,ich4)+GRVOL(JL,NLEV)*
     *           (PM(JL,NLEV,ich4)-OLDC)
         ENDDO

C     LG- end LSETCH4

      ENDIF
      

C     - - -
C     3. NOx-emission from lightning
C     
      DO JK=1,NLEV
         DO JL=1,NLON
            IF (EMFLASH(JL).GT.0.) THEN
               IF (JK.GE.JKT(JL).AND.JK.LE.JKB(JL)) 
     *              ZCOLMSS(JL)=ZCOLMSS(JL)+GRMASS(JL,JK)
            ENDIF
         ENDDO
      ENDDO
      ZEMNOFL=0.
      DO JK=1,NLEV
         DO JL=1,NLON
            IF (EMFLASH(JL).GT.0.) THEN
               IF (JK.GE.JKT(JL).AND.JK.LE.JKB(JL)) THEN
C     addition in molec/cm3, weighted for gridmass. 

C     LG-    considering NOx emission factor

                  TOTNOFL=FEMNOX*EMFLASH(JL)*PRHOA(JL,JK)/ZCOLMSS(JL)*PTMST

C LG-             19092003, changed to deal with different treatment of the
C                 NOx species as tracers. For short timesteps, NO, NO2 are also
C                 considered as separate tracers and then the anthropogenic
C                 emissions of NOx are assumed to be in the form of NO2 

                  IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
                    PM(JL,JK,ino)=PM(JL,JK,ino)+TOTNOFL
		  ELSE
                    PM(JL,JK,inox)=PM(JL,JK,inox)+TOTNOFL
		  ENDIF

                  ZEMNOFL=ZEMNOFL+TOTNOFL*GRVOL(JL,JK)
               ENDIF
            ENDIF

C     LG-  determining location

            IW=IWHERE(JL,JK)

C     LG-  the contribution of the emission of NO through lightning
C     into the total budget is put into BXT(IW,IBEMIS,2) (no emission
C     of O3)

            BXT(IW,IBEMIS,io3)=BXT(IW,IBEMIS,io3)+ZEMNOFL

         ENDDO
      ENDDO

C     LG- two options, if the LBIOSPH=.TRUE. then the 
C     the NO emissions are in the lowest layer of the canopy whereas the
C     VOC emissions are distributed over the whole vegetation. If the switch is
C     off then the model just emits the gases in the surface layer and 
C     uses the emission input fields as they has been used in the ECHAM-3D
C     model or the biogenic NO emissions calculated in the subroutine
C     NOXEMIS.f dependent on the switch LNOEMIS

      IF (.NOT.LBIOSPH) THEN

         WRITE(NUNMDFL,*)'Emission in surface layer'
         IF (LBULKVEG.OR.LVEG_MLAY.OR.LSNOW_MLAY) ! mz_lg_20051218+
     *        WRITE(NUNMDFL,*)'(except of biogenic emis., e.g, NO, isoprene)'

C     LG-  Originally, in the L19 version the emission is in the two
C     lowest layers with 30% being emitted in layer 18 and 70%
C     in the surface layer. This has not been changed in this model 
C     version.

C     LG-  reseting the NO emission

         DO JK=1,NLEV
            DO JL=1,NLON
               XNOEMIS(JL,JK)=0.
            ENDDO
         ENDDO

         DO JK=NLEV-1,NLEV
            FFAC=0.70
            IF (JK.EQ.18) FFAC=0.30
            DO JL=1,NLON

c     EMIFAC in m2 s cm-3

               EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)

C     CO
               COEMX=FEMCO*EMIFAC*COEM(ILON,IROW)
               
C     LG-   08-2002, introduced the option to modify the CO emission fluxes

               IF (IRESET_COEMIS.EQ.1) COEMX=FEMCO*EMIFAC*XCOEMT_RESET
               
               PM(JL,JK,ico)=PM(JL,JK,ico)+COEMX

C     LG-   determining location

               IW=IWHERE(JL,JK)

               BXT(IW,IBEMIS,ico)=BXT(IW,IBEMIS,ico)+COEMX*GRVOL(JL,JK)

C     NO
C     IND-PRIND

C     LG-   different scenario's can be choosen by selecting the 
C     specific sources

C     LG-   Industrial scenario, original ECHAM emission fields

               IF ((LINDUS.AND..NOT.LNOEMIS).OR.(LINDUS.AND..NOT.LDAYCENT)) ! mz_lg_20050125+ 
     *              XNOXEMT=XNOXEM(ILON,IROW,1)+XNOXEM(ILON,IROW,2)+
     *              XNOXEM(ILON,IROW,3)

C     LG-   Industrial scenario, explicit biogenic emissions

               IF ((LINDUS.AND.LNOEMIS).OR.(LINDUS.AND.LDAYCENT)) ! mz_lg_20050125+ 
     *              XNOXEMT=XNOXEM(ILON,IROW,1)+XNOXEM(ILON,IROW,3)

C     LG-   Pre-industrial

               IF ((.NOT.LINDUS.AND..NOT.LNOEMIS).OR.(.NOT.LINDUS.AND..NOT.LDAYCENT))
     *              XNOXEMT=1.5/6.*XNOXEM(ILON,IROW,2)+XNOXEM(ILON,IROW,3)

C     LG-   Pre-industrial

               IF ((.NOT.LINDUS.AND.LNOEMIS).OR.(.NOT.LINDUS.AND.LDAYCENT))
     *              XNOXEMT=XNOXEM(ILON,IROW,3)

C     LG-   08-2002, introducing the option to use a modified total NOx emission
C     flux for IRESET_NOXEMIS=1

               IF (IRESET_NOXEMIS.EQ.1) XNOXEMT=XNOXEMT_RESET

C     LG-   end

               EMTOT=FEMNOX*EMIFAC*XNOXEMT !considering NOx emission factor

               IF (IROW.LE.12) EMTOT=0
               XNOEMIS(JL,JK)=EMTOT

C     LG-   biogenic emission is in the form of NO and not NOx!

               IF ((NLEVVEG.EQ.0.AND.LNOEMIS).OR.(NLEVVEG.EQ.0.AND.LDAYCENT).AND. ! mz_lg_20050125+
     &              IEMDD_CHEM(ino_pr).LT.1.AND.IEM_TURB(ino_pr).LT.1) THEN
                  IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
                     PM(JL,JK,ino)=PM(JL,JK,ino)+EMIFAC*NOMOLEC 
                     BXT(IW,IBEMIS,ino)=BXT(IW,IBEMIS,ino)+
     $                    EMIFAC*NOMOLEC*GRVOL(JL,JK)
                  ELSE
                     PMLOC(JL,JK,ino)=PMLOC(JL,JK,ino)+EMIFAC*NOMOLEC 
                     BG3(IW,IBEMIS,ino)=BG3(IW,IBEMIS,ino)+
     $                    EMIFAC*NOMOLEC*GRVOL(JL,JK)
                  ENDIF
               ENDIF       
               
C LG-          19092003, changed to deal with different treatment of the
C              NOx species as tracers. For short timesteps, NO, NO2 are also
C              considered as separate tracers and then the anthropogenic
C              emissions of NOx are assumed to be in the form of NO2 

               IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
                 PM(JL,JK,ino2)=PM(JL,JK,ino2)+EMTOT
	       ELSE 
                 PM(JL,JK,inox)=PM(JL,JK,inox)+EMTOT
	       ENDIF
               BXT(IW,IBEMIS,inox)=BXT(IW,IBEMIS,inox)+EMTOT*GRVOL(JL,JK)

            ENDDO
         ENDDO

 303     CONTINUE

C     LG- the hydrocarbons emission in the surface layer for the 
C     hydrocarbon chemistry scheme. The bulk hydrocarbon emission
C     flux is used here, which resembles the summed emission flux
C     in each layer within the canopy, thus assuming instant mixing of
C     between the canopy and the surface layer and not considering
C     any within-canopy interactions

         IF (LVOCEMIS) THEN
            
            DO JK=NLEV-1,NLEV
               DO JL=1,NLON
                  FFAC=0.7
                  IF (JK.EQ.NLEV) FFAC=0.3

c     EMIFAC in m2 s cm-3

                  EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)

C     LG-  recalculation from kg C m-2 s-1 to molecules C5H8 m-2 s-1,
C     there is the correction for the number of C in the isoprene
C     molecule, since the emission rate is in kg C, which is recalculated to
C     g C (1.E3), 1 g C resembles 1/12 mol C. 1/12 mol C = (1/12)*(1/5) mol 
C     isoprene and this is recalculated to molecules isoprene by multiplying
C     then with the number of avogadro  

                  IF (NLEVVEG.EQ.0.AND.IEMDD_CHEM(iisop).LT.1.AND.
     &                 IEM_TURB(iisop).LT.1) THEN
                     PM(JL,JK,iisop)=PM(JL,JK,iisop)+EMIFAC*ISOPMOLEC      

C     LG-   determining location

                     IW=IWHERE(JL,JK)

                     BXT(IW,IBEMIS,iisop)=BXT(IW,IBEMIS,iisop)+EMIFAC*
     *                    ISOPMOLEC*GRVOL(JL,JK)
                  ENDIF

C     LG-  monoterpene emissions, here the emmission factors are used to
C     distribute the total monoterpene emissions over the alpha and
C     beta species

                  IF (NLEVVEG.EQ.0.AND.IEMDD_CHEM(imaterp).LT.1.AND.
     &                 IEM_TURB(imaterp).LT.1.AND.IEMDD_CHEM(imtterp).LT.1) THEN
                     PM(JL,JK,imaterp)=PM(JL,JK,imaterp)+EMIFAC*FEMMATERP*MTERPMOLEC  
                     PM(JL,JK,imbterp)=PM(JL,JK,imbterp)+EMIFAC*FEMMBTERP*MTERPMOLEC    
                     PM(JL,JK,imtterp)=PM(JL,JK,imtterp)+EMIFAC*FEMMTTERP*MTERPMOLEC

C     LG-   determining location

                     IW=IWHERE(JL,JK)

                     BXT(IW,IBEMIS,imaterp)=BXT(IW,IBEMIS,imaterp)+EMIFAC*
     *                    FEMMATERP*MTERPMOLEC*GRVOL(JL,JK)
                     BXT(IW,IBEMIS,imbterp)=BXT(IW,IBEMIS,imbterp)+EMIFAC*
     *                    FEMMBTERP*MTERPMOLEC*GRVOL(JL,JK)
                     BXT(IW,IBEMIS,imtterp)=BXT(IW,IBEMIS,imtterp)+EMIFAC*
     *                    FEMMTTERP*MTERPMOLEC*GRVOL(JL,JK)
                  ENDIF

C     LG-  sesquiterpenes emissions

                  IF (NLEVVEG.EQ.0.AND.IEMDD_CHEM(isqterp).LT.1.AND.
     &                 IEM_TURB(isqterp).LT.1) THEN
                     PM(JL,JK,isqterp)=PM(JL,JK,isqterp)+EMIFAC*FEMOVOC*OVOCMOLEC  

C     LG-   determining location

                     IW=IWHERE(JL,JK)

                     BXT(IW,IBEMIS,isqterp)=BXT(IW,IBEMIS,isqterp)+EMIFAC*
     *                    FEMOVOC*OVOCMOLEC*GRVOL(JL,JK)
                  ENDIF

               ENDDO
            ENDDO

         ENDIF

C     LG- adding the sulfur emission (SO2/DMS), (LSULFCHEM=.TRUE) 

         IF (LSULFCHEM) THEN

            DO JK=NLEV-1,NLEV
               FFAC=0.30
               SFAC2=0.60
               SFAC1=0.
               IF (JK.EQ.18) THEN
                  FFAC=0.70
                  SFAC2=0.40
                  SFAC1=1.
               ENDIF
               DO JL=1,NLON

C     LG-  determining location

                  IW=IWHERE(JL,JK)

c     EMIFAC in m2 s cm-3

                  EMIFAC=FFAC*1.E-6*PTMST*1.E3*G*PRHOA(JL,JK)/PDP(JL,JK)

C     -- SO2 biom. burn., scale from N to S (2.3 Tg S yr-1)

                  BBNS=2.3/6.*XMN/XMS
                  SEMTOT=FEMSO2*EMIFAC*XNOXEM(ILON,IROW,2)*BBNS
                  PM(JL,JK,iso2)=PM(JL,JK,iso2)+SEMTOT
                  BXT(IW,1,iso2)=BXT(IW,1,iso2)+SEMTOT*GRVOL(JL,JK)

C     -- DMS --

                  IF (JK.EQ.NLEV) THEN
                     EMTOT=FEMDMS*EMIFAC/FFAC*DMSEM(ILON,IROW,1)
                     PM(JL,JK,idms)=PM(JL,JK,idms)+EMTOT
                     BXT(IW,1,idms)=BXT(IW,1,idms)+EMTOT*GRVOL(JL,JK)
                  ENDIF

C     -- SO2 --

                  EMTOT=FEMSO2*EMIFAC/FFAC*
     *                 (SFAC2*SO2EM(ILON,IROW,2)+SFAC1*SO2EM(ILON,IROW,1))
                  PM(JL,JK,iso2)=PM(JL,JK,iso2)+0.80*EMTOT
                  BXT(IW,1,iso2)=BXT(IW,1,iso2)+0.95*EMTOT*GRVOL(JL,JK)
                  BXT(IW,3,iso2)=BXT(IW,3,iso2)-0.15*EMTOT*GRVOL(JL,JK)
                  PM(JL,JK,iso4)=PM(JL,JK,iso4)+0.05*EMTOT
                  BXT(IW,1,iso4)=BXT(IW,1,iso4)+0.05*EMTOT*GRVOL(JL,JK)

               ENDDO
            ENDDO

C     - - -
C     volcanic
C     
            DO JL=1,NLON
               IF (IVTOP(JL,ILAT).GT.0) THEN
                  EMIFAC2=PTMST*G*AVO/XMS
                  ZEMI1=EMIFAC2*0.36
                  ZEMI2=EMIFAC2*0.36
                  ZEMI3=EMIFAC2*0.28
                  JK=IVTOP(JL,ILAT)
                  EMTOT=ZEMI1*SVOLC(ILON,ILAT)*PRHOA(JL,JK)/PDP(JL,JK)
                  PM(JL,JK,iso2)=PM(JL,JK,iso2)+EMTOT
                  BXT(IW,IBEMIS,iso2)=BXT(IW,IBEMIS,iso2)+EMTOT*GRVOL(JL,JK)
                  EMTOT16=0.5*ZEMI2*SVOLC(ILON,ILAT)*PRHOA(JL,16)/PDP(JL,16)
                  PM(JL,16,iso2)=PM(JL,16,iso2)+EMTOT16
                  BXT(IW,IBEMIS,iso2)=BXT(IW,IBEMIS,iso2)+EMTOT16*GRVOL(JL,16)
                  EMTOT15=0.5*ZEMI2*SVOLC(ILON,ILAT)*PRHOA(JL,15)/PDP(JL,15)
                  PM(JL,15,iso2)=PM(JL,15,iso2)+EMTOT15
                  BXT(IW,IBEMIS,iso2)=BXT(IW,IBEMIS,iso2)+EMTOT15*GRVOL(JL,15)
                  EMTOT11=0.5*ZEMI3*SVOLC(ILON,ILAT)*PRHOA(JL,11)/PDP(JL,11)
                  PM(JL,11,iso2)=PM(JL,11,iso2)+EMTOT11
                  BXT(IW,IBEMIS,iso2)=BXT(IW,IBEMIS,iso2)+EMTOT11*GRVOL(JL,11)
                  EMTOT10=0.5*ZEMI3*SVOLC(ILON,ILAT)*PRHOA(JL,10)/PDP(JL,10)
                  PM(JL,10,iso2)=PM(JL,10,iso2)+EMTOT10
                  BXT(IW,IBEMIS,iso2)=BXT(IW,IBEMIS,iso2)+EMTOT10*GRVOL(JL,10)
               ENDIF
            ENDDO

         ENDIF

C     LG- end sulfur emissions

C     LG- and the emission of radon being added

         DO JK=NLEV,NLEV
            DO JL=1,NLON
               FFAC=1.0

c     EMIFAC in m2 s cm-3

               EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)*
     &              MIN(1.,NSTEP/(T_SPINUP/PTMST)) ! smoothing startup

               RADEMX=FEMRAD*EMIFAC*RADEM

C     LG-  determining location

               IW=IWHERE(JL,JK)

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(irad).LT.1.AND.IEM_TURB(irad).LT.1) THEN
                  PM(JL,JK,irad)=PM(JL,JK,irad)+EMIFAC*RADEM
                  BXT(IW,IBEMIS,irad)=BXT(IW,IBEMIS,irad)+RADEMX*GRVOL(JL,JK)
               ENDIF

C     LG-  05-2001, added CO2 respiration rates 

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(ico2).LT.1.AND.IEM_TURB(ico2).LT.1) THEN
                  CO2EMX=FEMCO2*EMIFAC*CO2MOLEC
                  PM(JL,JK,ico2)=PM(JL,JK,ico2)+CO2EMX
                  BXT(IW,IBEMIS,ico2)=BXT(IW,IBEMIS,ico2)+CO2EMX*GRVOL(JL,JK)
               ENDIF

C     LG-  08-2002, added CH3CL and CHCL3 estimated emission fluxes

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(ich3cl).LT.1.AND.IEM_TURB(ich3cl).LT.1) THEN
                  CH3CLEMX=FEMCH3CL*EMIFAC*CH3CLMOLEC
                  PM(JL,JK,ich3cl)=PM(JL,JK,ich3cl)+CH3CLEMX
                  BXT(IW,IBEMIS,ich3cl)=BXT(IW,IBEMIS,ich3cl)+CH3CLEMX*GRVOL(JL,JK)
               ENDIF

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(ichcl3).LT.1.AND.IEM_TURB(ichcl3).LT.1) THEN
                  CHCL3EMX=FEMCHCL3*EMIFAC*CHCL3MOLEC
                  PM(JL,JK,ichcl3)=PM(JL,JK,ichcl3)+CHCL3EMX
                  BXT(IW,IBEMIS,ichcl3)=BXT(IW,IBEMIS,ichcl3)+CHCL3EMX*GRVOL(JL,JK)
               ENDIF

C     LG-  end

C     LG-  09-2001, added NH3 emission flux (antropogenic) 

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(inh3).LT.1.AND.IEM_TURB(inh3).LT.1) THEN
                  NH3EMX=FEMNH3*EMIFAC*NH3MOLEC
                  PM(JL,JK,inh3)=PM(JL,JK,inh3)+NH3EMX
                  BXT(IW,IBEMIS,inh3)=BXT(IW,IBEMIS,inh3)+NH3EMX*GRVOL(JL,JK)
               ENDIF

C     LG-  end

C     LG-  12-2003, added the HONO and NOX emission flux related to
C          the photodissociation of the HNO3 deposited at the surface

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(ihono).LT.1.AND.IEM_TURB(ihono).LT.1) THEN
                  HONOEMX=EMIFAC*HONO_JHNO3SEM
                  NOXEMX=EMIFAC*NOX_JHNO3SEM
                  PM(JL,JK,ihono)=PM(JL,JK,ihono)+HONOEMX
                  IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
                    PM(JL,JK,ino2)=PM(JL,JK,ino2)+NOXEMX
	          ELSE 
                    PM(JL,JK,inox)=PM(JL,JK,inox)+NOXEMX
	          ENDIF
                  BXT(IW,IBEMIS,ihono)=BXT(IW,IBEMIS,ihono)+HONOEMX*GRVOL(JL,JK)
                  BXT(IW,IBEMIS,inox)=BXT(IW,IBEMIS,inox)+NOXEMX*GRVOL(JL,JK)
               ENDIF

C     LG-  end

C     LG-  01-2006, added methanol, acetone and acetaldehyde biogenic emission flux

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(ich3oh).LT.1.AND.IEM_TURB(ich3oh).LT.1) THEN
                  CH3OHEMX=FEMCH3OH*EMIFAC*CH3OHMOLEC
                  PM(JL,JK,ich3oh)=PM(JL,JK,ich3oh)+CH3OHEMX
                  BXT(IW,IBEMIS,ich3oh)=BXT(IW,IBEMIS,ich3oh)+CH3OHEMX*GRVOL(JL,JK)
               ENDIF

               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(iacet).LT.1.AND.IEM_TURB(iacet).LT.1) THEN
                  ACETEMX=FEMACET*EMIFAC*ACETMOLEC
                  PM(JL,JK,iacet)=PM(JL,JK,iacet)+ACETEMX
                  BXT(IW,IBEMIS,iacet)=BXT(IW,IBEMIS,iacet)+ACETEMX*GRVOL(JL,JK)
               ENDIF

               ! mz_lg_20060104+ assigning acetaldehyde emissions to ald2 (higher aldehydes)
               IF (NLEVVEG.EQ.0.AND.
     &              IEMDD_CHEM(iald2).LT.1.AND.IEM_TURB(iald2).LT.1) THEN
                  ACETALDEMX=FEMACETALD*EMIFAC*ACETALDMOLEC
                  PM(JL,JK,iald2)=PM(JL,JK,iald2)+ACETALDEMX
                  BXT(IW,IBEMIS,iald2)=BXT(IW,IBEMIS,iald2)+ACETALDEMX*GRVOL(JL,JK)
               ENDIF

C     LG-  end


            ENDDO
         ENDDO

C     LG- end radon, CO2, NH3, HONO, NOx emissions

C     
C     LG- the emissions considering the biosphere
C     
      ELSEIF (LBIOSPH) THEN

         WRITE(NUNMDFL,*)'Emission within the canopy and the surface layer'

C     LG- the emission of CO is just in the surface layer since the
C     the exact locations of the sources are not well established
C     (see papers of the ABLE experiments) and also the chemical lifetime
C     is that long that only the source strenght is relevant and not 
C     the exact location within the canopy 

         DO JK=NLEV-1,NLEV
            FFAC=0.70
            IF (JK.EQ.NLEV-1) FFAC=0.30

            DO JL=1,NLON

c     EMIFAC in m2 s cm-3

               EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)

C     LG-  FEMCO is a scaling factor for the CO emission

               COEMX=FEMCO*EMIFAC*COEM(ILON,IROW)

C     LG-  08-2002, introduced the option to modify the CO emission fluxes

               IF (IRESET_COEMIS.EQ.1) COEMX=FEMCO*EMIFAC*XCOEMT_RESET

               PM(JL,JK,ico)=PM(JL,JK,ico)+COEMX

C     LG-  determining location

               IW=IWHERE(JL,JK)

               BXT(IW,IBEMIS,ico)=BXT(IW,IBEMIS,ico)+COEMX*GRVOL(JL,JK)

               IF (NSTEP.EQ.0.AND.JK.EQ.NLEV) THEN
                  COMOLEC=1E-4*FEMCO*COEM(ILON,IROW)*1E-11
                  WRITE(NUNMDFL,'(1a,f8.3)')
     *                 ' The CO flux in 10*11 molecules cm-2 s-1 is:',COMOLEC
               ENDIF      

            ENDDO
         ENDDO

         DO JL=1,NLON

C     --- emission NO and CH2O at surface:

C     NO

            DO JK=1,NLEVEL
               XNOEMIS(JL,JK)=0.
            ENDDO

C     LG- biogenic emission of NO in the lowest layer of the canopy,
C     the emitted NO is added to the NOx since this tracer is subject
C     to vertical turbulent transport, and since NO can have a relatively
C     long chemical lifetime within the canopy, vertical transport is 
C     important for NO, a possible improvement of this can be to add
C     NO to the group of long lived species (16-03-98)

            DO JK=NLEVEL,NLEVEL
               FFAC=1.0
               IF (JK.EQ.NLEVEL-1) FFAC=0.0

c     EMIFAC in m2 s cm-3

               EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)

C     LG- factor for recalculation of kg N m-2 s-1 to molecules NO m-2 s-1

               EMRECALC=1.E3*(1./XMN)*(1./1.)*ZAVO
               EMNO=CRF_YL95*FEMNO*EMIFAC*EMRECALC*NOEM
               XNOEMIS(JL,JK)=EMNO

C     LG- determining location

               IW=IWHERE(JL,JK)

C     LG- biogenic emission is in the form of NO and not NOx! For multilayer
C     approach only the biogenic emissions are considered and not other
C     sources for NO/NOx

               IF (IEMDD_CHEM(ino_pr).LT.1.AND.IEM_TURB(ino_pr).LT.1) THEN
                  IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
                     PM(JL,JK,ino)=PM(JL,JK,ino)+EMNO 
                     BXT(IW,IBEMIS,ino)=BXT(IW,IBEMIS,ino)+
     $                    EMNO*GRVOL(JL,JK)
                  ELSE
                     PMLOC(JL,JK,ino)=PMLOC(JL,JK,ino)+EMNO 
                     BG3(IW,IBEMIS,ino)=BG3(IW,IBEMIS,ino)+
     $                    EMNO*GRVOL(JL,JK)
                  ENDIF
               ENDIF
               
            ENDDO

c     -- VOC emissions
            
            IF (LVOCEMIS) THEN

               DO JK=NLEV+1,NLEVEL
                  JJ=NLEVEL+1-JK
                  FFAC=1.0

c     EMIFAC in m2 s cm-3

                  EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)

C     LG-   Isoprene, factor for recalculation of kg C m-2 s-1 to molecules 
C     C5H8 m-2 s-1
                  
                  EMRECALC=1.E3*(1./XMC)*(1./5.)*ZAVO

C     LG-   since only the vegetation emits the isoprene and not the bare soil,
C     the isoprene emission rate is corrected for the vegetation fraction
C     VEGFRAC and the wet skin fraction WSFRAC

                  EMISOP=(VEGFRAC(JL)+WSFRAC(JL))*FEMISOP*EMIFAC*EMRECALC*ISOPEM(JJ)

                  IF (IEMDD_CHEM(iisop).LT.1.AND.IEM_TURB(iisop).LT.1) THEN
                       PM(JL,JK,iisop)=PM(JL,JK,iisop)+EMISOP 
		  ENDIF

C     LG-   determining location

                  IW=IWHERE(JL,JK)

                  BXT(IW,IBEMIS,iisop)=BXT(IW,IBEMIS,iisop)+EMISOP*GRVOL(JL,JK) 

C     LG-   monoterpenes, factor for recalculation of kg C m-2 s-1 to molecules 
C     C10H18 m-2 s-1
                  
                  EMRECALC=1.E3*(1./XMC)*(1./10.)*ZAVO
                  EMMATERP=(VEGFRAC(JL)+WSFRAC(JL))*FEMMATERP*EMIFAC*
     $                 EMRECALC*MTERPEM(JJ)
                  EMMBTERP=(VEGFRAC(JL)+WSFRAC(JL))*FEMMBTERP*EMIFAC*
     $                 EMRECALC*MTERPEM(JJ)
                  EMMTTERP=(VEGFRAC(JL)+WSFRAC(JL))*FEMMTTERP*EMIFAC*
     $                 EMRECALC*MTERPEM(JJ)

                  IF (IEMDD_CHEM(imaterp).LT.1.AND.IEM_TURB(imaterp).LT.1)
     $                 PM(JL,JK,imaterp)=PM(JL,JK,imaterp)+EMMATERP 
                  IF (IEMDD_CHEM(imbterp).LT.1.AND.IEM_TURB(imbterp).LT.1)
     $                 PM(JL,JK,imbterp)=PM(JL,JK,imbterp)+EMMBTERP 
                  IF (IEMDD_CHEM(imtterp).LT.1.AND.IEM_TURB(imtterp).LT.1)
     $                 PM(JL,JK,imtterp)=PM(JL,JK,imtterp)+EMMTTERP 

C     LG-   determining location

                  IW=IWHERE(JL,JK)

                  BXT(IW,IBEMIS,imaterp)=BXT(IW,IBEMIS,imaterp)+EMMATERP*GRVOL(JL,JK) 
                  BXT(IW,IBEMIS,imbterp)=BXT(IW,IBEMIS,imbterp)+EMMBTERP*GRVOL(JL,JK) 
                  BXT(IW,IBEMIS,imtterp)=BXT(IW,IBEMIS,imtterp)+EMMTTERP*GRVOL(JL,JK) 

C     LG-   sesquiterpenes, factor for recalculation of kg C m-2 s-1 to molecules 
C     C15H24 m-2 s-1
                  
                  EMRECALC=1.E3*(1./XMC)*(1./15.)*ZAVO
                  EMSQTERP=(VEGFRAC(JL)+WSFRAC(JL))*FEMOVOC*EMIFAC*
     $                 EMRECALC*OVOCEM(JJ) !!! note that we use the OVOC emissions to mimic
                                ! the sesquiterpenes emissions

                  IF (IEMDD_CHEM(isqterp).LT.1.AND.IEM_TURB(isqterp).LT.1)
     $                 PM(JL,JK,isqterp)=PM(JL,JK,isqterp)+EMSQTERP 

C     LG-   determining location

                  IW=IWHERE(JL,JK)
                  BXT(IW,IBEMIS,isqterp)=BXT(IW,IBEMIS,isqterp)+EMSQTERP*GRVOL(JL,JK) 

               ENDDO

C     LG- end if LVOCEMIS 

            ENDIF

 5000       format (i3,5E11.2)

         ENDDO

C     LG- adding the sulfur emission (SO2/DMS),  For the 
C     biosphere mode, only the DMS and SO2 emissions through biomass
C     burning are considered, the DMS source is distributed over the whole
C     canopy, scaled using the isoprene emission (see paper Andreae, JGR, 95)

         IF (LSULFCHEM) THEN

            DO JK=NLEV-1,NLEVEL

C     LG- determining the total number of levels in order to divide
C     bulk emission rates over the total number of layers, e.g., the
C     DMS flux

               DLEV=(NLEVEL+1-(NLEV-1))

               FFAC=0.30
               IF (JK.EQ.NLEV-1) THEN
                  FFAC=0.70
               ENDIF
               DO JL=1,NLON

C     LG-  determining location

                  IW=IWHERE(JL,JK)

c     EMIFAC in m2 s cm-3, be aware of the correction for the number of levels

                  EMIFAC=(FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK))/DLEV

C     -- SO2 biom. burn., scale from N to S (2.3 Tg S yr-1)

                  BBNS=2.3/6.*XMN/XMS
                  SEMTOT=FEMSO2*EMIFAC*XNOXEM(ILON,IROW,2)*BBNS

                  IF (NSTEP.EQ.0.AND.JK.EQ.NLEV) THEN
                     SO2MOLEC=1E-4*FEMSO2*EMIFAC*XNOXEM(ILON,IROW,2)*BBNS*1E-11
                     WRITE(NUNMDFL,'(1a,f8.3)')
     *                    ' The SO2 flux in 10*11 molecules cm-2 s-1 is:',SO2MOLEC
                  ENDIF      

                  PM(JL,JK,iso2)=PM(JL,JK,iso2)+SEMTOT
                  BXT(IW,IBEMIS,iso2)=BXT(IW,IBEMIS,iso2)+SEMTOT*GRVOL(JL,JK)

C     -- DMS --

                  EMTOT=FEMDMS*EMIFAC/FFAC*DMSEM(ILON,IROW,1)
                  PM(JL,JK,idms)=PM(JL,JK,idms)+EMTOT
                  BXT(IW,IBEMIS,idms)=BXT(IW,IBEMIS,idms)+EMTOT*GRVOL(JL,JK)

               ENDDO
            ENDDO

         ENDIF

C     LG- end sulfur emissions

C     LG- and the emission of radon being added

         DO JK=NLEVEL,NLEVEL
            DO JL=1,NLON
               FFAC=1.0

c     EMIFAC in m2 s cm-3

               EMIFAC=FFAC*1.E-6*PTMST*G*1.E3*PRHOA(JL,JK)/PDP(JL,JK)

               RADEMX=FEMRAD*EMIFAC*RADEM

C     LG-  determining location

               IW=IWHERE(JL,JK)

               PM(JL,JK,irad)=PM(JL,JK,irad)+RADEMX
               BXT(IW,IBEMIS,irad)=BXT(IW,IBEMIS,irad)+RADEMX*GRVOL(JL,JK)

C     LG-  05-2001, added the CO2 respiration/emission flux

               CO2EMX=FEMCO2*EMIFAC*CO2MOLEC
               PM(JL,JK,ico2)=PM(JL,JK,ico2)+CO2EMX
               BXT(IW,IBEMIS,ico2)=BXT(IW,IBEMIS,ico2)+CO2EMX*GRVOL(JL,JK)

C     LG-  end CO2

C     LG-  09-2001, added the NH3 emission flux

               NH3EMX=FEMNH3*EMIFAC*NH3MOLEC
               PM(JL,JK,inh3)=PM(JL,JK,inh3)+NH3EMX
               BXT(IW,IBEMIS,inh3)=BXT(IW,IBEMIS,inh3)+NH3EMX*GRVOL(JL,JK)

C     LG-  end NH3

C     LG-  05-2003, added the CH3CL/CHCL3 emission flux

               CH3CLEMX=FEMCH3CL*EMIFAC*CH3CLMOLEC
               PM(JL,JK,ich3cl)=PM(JL,JK,ich3cl)+CH3CLEMX
               BXT(IW,IBEMIS,ich3cl)=BXT(IW,IBEMIS,ich3cl)+CH3CLEMX*GRVOL(JL,JK)

               CHCL3EMX=FEMCHCL3*EMIFAC*CHCL3MOLEC
               PM(JL,JK,ichcl3)=PM(JL,JK,ichcl3)+CHCL3EMX
               BXT(IW,IBEMIS,ichcl3)=BXT(IW,IBEMIS,ichcl3)+CHCL3EMX*GRVOL(JL,JK)

C     LG-  end CH3CL/CHCL3

C     LG-  05-2003, added the methanol, acetone and acetaldehyde emission flux

               CH3OHEMX=FEMCH3OH*EMIFAC*CH3OHMOLEC
               PM(JL,JK,ich3oh)=PM(JL,JK,ich3oh)+CH3OHEMX
               BXT(IW,IBEMIS,ich3oh)=BXT(IW,IBEMIS,ich3oh)+CH3OHEMX*GRVOL(JL,JK)

               ACETEMX=FEMACET*EMIFAC*ACETMOLEC
               PM(JL,JK,iacet)=PM(JL,JK,iacet)+ACETEMX
               BXT(IW,IBEMIS,iacet)=BXT(IW,IBEMIS,iacet)+ACETEMX*GRVOL(JL,JK)

               ACETALDEMX=FEMACETALD*EMIFAC*ACETALDMOLEC
               PM(JL,JK,iald2)=PM(JL,JK,iald2)+ACETALDEMX
               BXT(IW,IBEMIS,iald2)=BXT(IW,IBEMIS,iald2)+ACETALDEMX*GRVOL(JL,JK)

C     LG-  end

            ENDDO
         ENDDO

C     LG- end radon emissions

      ENDIF

C     LG- assigning the emission flux to EMIS (used for conc. estimate at zz)

      DO JL=1,NLON

C     LG-  Nitrogen oxide(s)

         IF (LNOEMIS.OR.LDAYCENT) THEN ! mz_lg_20050125+ modified
            EMIS_BS(JL,ino_pr)=NOMOLEC*1.E-4
            EMIS_VEG(JL,ino_pr)=EMIS_BS(JL,ino_pr)
            EMIS_WS(JL,ino_pr)=EMIS_BS(JL,ino_pr)
         ENDIF

C     LG-  isoprene

         IF (LVOCEMIS) THEN
            EMIS_VEG(JL,iisop)=ISOPMOLEC*1.E-4
            EMIS_WS(JL,iisop)=EMIS_VEG(JL,iisop)
         ENDIF

C     LG-  dms

         EMIS_WAT(JL,idms)=DMSEM(ILON,IROW,1)*1.E-4
         
C     LG-  radon

         EMIS_BS(JL,irad)=RADEM*1.E-4       
         EMIS_VEG(JL,irad)=EMIS_BS(JL,irad)
         EMIS_WS(JL,irad)=EMIS_BS(JL,irad)
         EMIS_WAT(JL,irad)=RADEM*1.E-4

C     LG-  summing all emission fluxes for the fractions in order to get a
C     total canopy-top/soil surface/snow surface/water surface emission
C     flux
         
         DO JT=1,NTRAC
            EMIS_TOT(JL,JT)=
     &           VEGFRAC(JL)*EMIS_VEG(JL,JT)+WSFRAC(JL)*EMIS_WS(JL,JT)+
     &           BSFRAC(JL)*EMIS_BS(JL,JT)+SNFRAC(JL)*EMIS_SN(JL,JT)+
     &           (1.-(VEGFRAC(JL)+WSFRAC(JL)+BSFRAC(JL)+SNFRAC(JL)))*
     &           EMIS_WAT(JL,JT)*MIN(1.,NSTEP/(T_SPINUP/PTMST))
         ENDDO
         
      ENDDO

C     LG- end
      
      IF (NSTEP.EQ.0) THEN
         OPEN (UNIT=NUNNOXEM,
     &        FILE='/data/ganzevl/racmo/output/NOx_emis.out',
     &        STATUS='UNKNOWN')
         WRITE(NUNNOXEM,'(1a)')'NO(x) emission rates '
      ENDIF

c     -- writing of the NO/NOx emission flux

      IF (LNOEMIS) THEN
         IF (NSTEP.EQ.0) THEN
            WRITE(NUNNOXEM,'(2a)')
     &           ' Antrop. and soil biogenic NO emission in 1E9 molec. cm-2 s-1', 
     &           ' and ng N m-2 s-1'
            WRITE(NUNNOXEM,'(a10,a15,9a14,5a20,a10)')
     &           'istep','time','emisfact-NOx',
     &           'antrop.-flx','bb-flx','NO-class','cult-index',
     &           'Tsoil[K]','soil-wetness','emisfact-NO','CRF_YL95',
     &           '[1E9-molec./cm2/s]','[ng-N/m2/s]','fert.flx [ng-N..]',
     &           'XNOXEMT [1E9...]','CRF*flx [1e9..]','f(T)'
         ENDIF

C     LG-  writing the average surface cover parameters, ofcourse only the
C     "contineous" data can be averaged and not the discrete data such as
C     as the emission classes and vegetation classes

         CRF_AVG=CRF_AVG+CRF_YL95
         CULT_AVG=CULT_AVG+CULT_INDEX
         NOMOLEC_AVG=NOMOLEC_AVG+(NOMOLEC*1.E-9*1.E-4)/CRF_YL95 ! BE AWARE THAT
                                ! THIS IS THE SOIL EMISSION FLUX, NOT CORRECTED FOR THE CRF
         NOEM_AVG=NOEM_AVG+FEMNO*1.E12*NOEM
         FERTFLX_AVG=FERTFLX_AVG+1.E12*FERTFLX
         XNOXEMT_AVG=XNOXEMT_AVG+FEMNOX*XNOXEMT*1.E-9*1.E-4
         CRFFLX_AVG=CRFFLX_AVG+NOMOLEC*1.E-9*1.E-4

         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN
            IF (NSTEP.GT.0.) THEN
               CRF_AVG=CRF_AVG/NPRINT
               CULT_AVG=CULT_AVG/NPRINT
               NOMOLEC_AVG=NOMOLEC_AVG/NPRINT
               NOEM_AVG=NOEM_AVG/NPRINT
               FERTFLX_AVG=FERTFLX_AVG/NPRINT
               XNOXEMT_AVG=XNOXEMT_AVG/NPRINT
               CRFFLX_AVG=CRFFLX_AVG/NPRINT
            ENDIF

            WRITE(NUNNOXEM,'(1x,i9.9,1x,a14,3(6x,f8.4),(12x,i2),5(6x,f8.4),
     &           5(10x,e10.4),1(5x,f5.2))')
     &           NSTEP,LDATLTIME,FEMNOX,FEMNOX*XNOXEM(ILON,IROW,1)*1.E-9*1.E-4,
     &           FEMNOX*XNOXEM(ILON,IROW,3)*1.E-9*1.E-4,INO_CLASS,CULT_AVG,
     &           SOILT(1),SOILWS(1),FEMNO,CRF_AVG,NOMOLEC_AVG,NOEM_AVG,
     &           FERTFLX_AVG,XNOXEMT_AVG,CRFFLX_AVG,FT

            CRF_AVG=0.
            CULT_AVG=0.
            NOMOLEC_AVG=0.
            NOEM_AVG=0.
            FERTFLX_AVG=0.
            XNOXEMT_AVG=0.
            CRFFLX_AVG=0.
         ENDIF
      ELSE IF (LDAYCENT) THEN
         IF (NSTEP.EQ.0) THEN
            WRITE(NUNNOXEM,'(2a)')
     &           ' Antrop. and soil biogenic NO emission in 1E9 molec. cm-2 s-1', 
     &           ' and ng N m-2 s-1'
            WRITE(NUNNOXEM,'(2a)')
     &           ' Soil biogenic NO emission calculated with DayCent model'
            WRITE(NUNNOXEM,'(2a)')
     &           ' NO CRF yet included!!'
            WRITE(NUNNOXEM,'(a10,a15,3a14,3a20)')
     &           'istep','time','emisfact-NOx',
     &           'antrop.-flx','bb-flx',
     &           '[1E9-molec./cm2/s]','[ng-N/m2/s]','XNOXEMT [1E9...]'
         ENDIF

C     LG-  writing the average surface cover parameters, ofcourse only the
C     "contineous" data can be averaged and not the discrete data such as
C     as the emission classes and vegetation classes

         NOMOLEC_AVG=NOMOLEC_AVG+(NOMOLEC*1.E-9*1.E-4)
         NOEM_AVG=NOEM_AVG+FEMNO*1.E12*NOEM
         XNOXEMT_AVG=XNOXEMT_AVG+FEMNOX*XNOXEMT*1.E-9*1.E-4

         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN
            IF (NSTEP.GT.0.) THEN
               NOMOLEC_AVG=NOMOLEC_AVG/NPRINT
               NOEM_AVG=NOEM_AVG/NPRINT
               XNOXEMT_AVG=XNOXEMT_AVG/NPRINT
            ENDIF

            WRITE(NUNNOXEM,'(1x,i9.9,1x,a14,3(6x,f8.4),
     &           3(10x,e10.4))')
     &           NSTEP,LDATLTIME,FEMNOX,FEMNOX*XNOXEM(ILON,IROW,1)*1.E-9*1.E-4,
     &           FEMNOX*XNOXEM(ILON,IROW,3)*1.E-9*1.E-4,NOMOLEC_AVG,NOEM_AVG,
     &           XNOXEMT_AVG

            NOMOLEC_AVG=0.
            NOEM_AVG=0.
            XNOXEMT_AVG=0.
         ENDIF
      ELSE
         IF (NSTEP.EQ.0) THEN
            WRITE(NUNNOXEM,'(1a)')
     &           ' NOx emission in 1E9 molec. and ng N m-2 s-1'
            WRITE(NUNNOXEM,'(a10,a15,3a20)')'istep','time','emisfact-NOx',
     &           '[1E9-molec./m2/s]','[ng-N/m2/s]'       
         ENDIF
         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0)
     &        WRITE(NUNNOXEM,'(1x,i9.9,1x,a14,3(10x,e10.4))')
     &        NSTEP,LDATLTIME,FEMNOX,FEMNOX*XNOXEMT*1E-9*1E-4,
     &        1.E12*XNOXEMT/(1.E3*(1./XMN)*(1./1.)*ZAVO)       
      ENDIF

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNNOXEM)

c     -- writing of integrated isoprene emission flux to file

      IF (LVOCEMIS.AND.EMISFACT(1).GT.0..OR.LVOCEMIS.AND.LAGRIAN) THEN

! mz_lg_20060312+ introduced the option to read a file that contains the
!        simulated VOC emission fluxes to constrain the calculations that
!        reflect exactly the same emissions for a change in meteorology and 
!        chemistry, e.g., to study the impact of a changing PBL height
    
         IF (LREAD_VOCEMIS) THEN 

           IF (NSTEP.EQ.0) THEN 
              WRITE(*,'(1a)')' emiss.f; line 1763, lread_vocemis=.TRUE.',
     &           ' Using VOC emission fluxes of previous run'
              FNAME=
     &          '/data/ganzevl/racmo/output/GABRIEL/MEGAN/3-days/INITPROF/VOCemis_nprint=1.out'
              WRITE(*,'(1a,1a)')' Reading: ',FNAME
              WRITE(*,'(1a)')' ENTER to continue '
              READ (*,*)
              OPEN (UNIT=NUNVOCEMINP,
     &           FILE=fname,
     &           STATUS='UNKNOWN')
              READ(NUNVOCEMINP,'(1a)')dummy
              READ(NUNVOCEMINP,'(1a)')dummy
              READ(NUNVOCEMINP,'(1a)')dummy
           ENDIF

C     LG-  reading the int. isoprene emission flux in ug C m-2 hr-1,

           READ(NUNVOCEMINP,'(1x,i9.9,1x,a14,9f15.3,2e15.6)') 
     &         NNSTEP,LLDATLTIME,
     &         VOCINP,VOCINP,VOCINP,
     &         VOCINP,VOCINP,
     &         VOCINP,VOCINP,
     &         VOCINP,VOCINP,
     &         ISOPEM(1),ISOPEM(2)

           WRITE(*,'(1a,i6,1x,a15,2e15.5)')
     &         ' Reading VOC emissions: ',nnstep,lldatltime,isopem(1),isopem(2)
         ENDIF

! mz_lg_20060313- end reading file produced in previous simulation

         IF (NSTEP.EQ.0) THEN 
            OPEN (UNIT=NUNVOCEM,
     &           FILE='/data/ganzevl/racmo/output/VOCemis.out',
     &           STATUS='UNKNOWN')
            WRITE(NUNVOCEM,'(1a)')'Int. VOC emissions [ug C m-2 hr-1] ',
     &           '[1E11 molec. cm-2 s-1] and [ppbv hr-1] for dZ=canopy height'
            WRITE(NUNVOCEM,'(a10,12a15)')
     &           'step','time','C5H8-emisf','C5H8','ratio','emis[1e11]',
     &           'emis[ppb/hr]','C10-emisf','C10','OVOC-emisf','OVOC',
     &           'ISOPEM(1)','ISOPEM(2)'
         ENDIF

C     LG-  writing of the int. isoprene emission flux in ug C m-2 hr-1, and the
C     ratio of the flux to the emission factor which expresses the influence
C     of the light and temperature dependent functions and of the int. 
C     isoprene emission flux in 1E11 molecules cm-2 s-1 and in ppbv hr-1
C     (to compare with tendencies). To arrive at ppbv hr-1 the emission flux
C     in molec. m-2 s-1 is multiplied with 1e-4, and divided by the canopy
C     height in cm (1.e2) yielding molec. cm-3 s-1, divided by the density
C     in g cm-3 yielding molec g-1 s-1 and then multiplied by ZMAIR [g mol-1]
C     divided by ZAVO [molec mol-1] yielding molec. C5H8 (molec. AIR)-1 s-1
C     and finally multiplied with 1e9 and 3600 to arrive at ppbv hr-1 

C     LG-  writing the average surface cover parameters, ofcourse only the
C     "contineous" data can be averaged and not the discrete data such as
C     as the emission classes and vegetation classes

         EMISFACT_AVG(1)=EMISFACT_AVG(1)+EMISFACT(1)
         EMISFACT_AVG(2)=EMISFACT_AVG(2)+EMISFACT(2)
         EMISFACT_AVG(3)=EMISFACT_AVG(3)+EMISFACT(3)
         ISOPEM_INT1_AVG=ISOPEM_INT1_AVG+FEMISOP*ISOPEM_INT*(1.E9*3600.)
         ISOPEM_INT2_AVG=ISOPEM_INT2_AVG+
     &        (FEMISOP*(ISOPEM_INT/MAX(1.E-20,DM))*1.E9*3600.)/
     &        MAX(1.E-20,EMISFACT(1))
         ISOPMOLEC1_AVG=ISOPMOLEC1_AVG+1.E-11*1.E-4*ISOPMOLEC
         ISOPMOLEC2_AVG=ISOPMOLEC2_AVG+
     &        ((1.E-4*ISOPMOLEC/(1.E2*HC))/PRHOA(1,NLEV))*(ZMAIR/ZAVO)*
     &        1.E9*3600.
         MTERPMOLEC_AVG=MTERPMOLEC_AVG+1.E-11*1.E-4*MTERPMOLEC
         OVOCMOLEC_AVG=OVOCMOLEC_AVG+1.E-11*1.E-4*OVOCMOLEC

         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN
            IF (NSTEP.GT.0.) THEN
               EMISFACT_AVG(1)=EMISFACT_AVG(1)/NPRINT
               EMISFACT_AVG(2)=EMISFACT_AVG(2)/NPRINT
               EMISFACT_AVG(3)=EMISFACT_AVG(3)/NPRINT
               ISOPEM_INT1_AVG=ISOPEM_INT1_AVG/NPRINT
               ISOPEM_INT2_AVG=ISOPEM_INT2_AVG/NPRINT
               ISOPMOLEC1_AVG=ISOPMOLEC1_AVG/NPRINT
               ISOPMOLEC2_AVG=ISOPMOLEC2_AVG/NPRINT
               MTERPMOLEC_AVG=MTERPMOLEC_AVG/NPRINT	 
               OVOCMOLEC_AVG=OVOCMOLEC_AVG/NPRINT	 
            ENDIF

            WRITE(NUNVOCEM,'(1x,i9.9,1x,a14,9f15.3,2e15.6)') 
     &           NSTEP,LDATLTIME,
     &           EMISFACT_AVG(1),ISOPEM_INT1_AVG,ISOPEM_INT2_AVG,
     &           ISOPMOLEC1_AVG,ISOPMOLEC2_AVG,
     &           EMISFACT_AVG(2),MTERPMOLEC_AVG,
     &           EMISFACT_AVG(3),OVOCMOLEC_AVG,
     &           ISOPEM(1),ISOPEM(2)

            EMISFACT_AVG(1)=0.
            EMISFACT_AVG(2)=0.
            ISOPEM_INT1_AVG=0.
            ISOPEM_INT2_AVG=0.
            ISOPMOLEC1_AVG=0.
            ISOPMOLEC2_AVG=0.
            MTERPMOLEC_AVG=0.
            OVOCMOLEC_AVG=0.

         ENDIF

c    --  writing of the radon emission flux in atoms m-2 s-1

         IF (NSTEP.EQ.0) THEN 
            OPEN (UNIT=NUNRADEM,
     &           FILE='/data/ganzevl/racmo/output/radonem.out',
     &           STATUS='UNKNOWN')
            WRITE(NUNRADEM,'(1a)')'Radon emission flux [atoms m-2 s-1] '
            WRITE(NUNRADEM,'(a10,a15,1a10)')
     &           'step','time','radon'


C LG-       and the HONO/NOx emission flux

            OPEN (UNIT=NUNHONOEM,
     &           FILE='/data/ganzevl/racmo/output/HONO-NOxem.out',
     &           STATUS='UNKNOWN')
            WRITE(NUNHONOEM,'(2a)')
     &        'HONO and NOx emission flux [1e11 molec. m-2 s-1] ',
     &        'The term N-em/DDHNO3 reflects the re-emitted fraction of ',
     &        'HNO3 deposited at the surface'
            WRITE(NUNHONOEM,'(a10,a15,7a15)')
     &           'step','time','HONOem','NOxem','VGRAT','LAI*FSLSUM',
     &           '1E11*DDHNO3','1E11*HNO3S','N-em/DDHNO3'
         ENDIF

         IF (MOD(NSTEP,NPRINT).EQ.0) 
     &        WRITE(NUNRADEM,'(1x,i9.9,1x,a14,e10.3)') 
     &        NSTEP,LDATLTIME,RADEMX

         IF (MOD(NSTEP,NPRINT).EQ.0) 
     &        WRITE(NUNHONOEM,'(1x,i9.9,1x,a14,7e15.3)') 
     &        NSTEP,LDATLTIME,1.E-11*HONO_JHNO3SEM,
     &                        1.E-11*NOX_JHNO3SEM,
     &                        VEGFRAC(1),LAI*FSLSUM,
     &                        1.E-11*DDHNO3,
     &                        1.E-11*HNO3S,
     &                        EMHNO3_ACCUM/DDHNO3_ACCUM

         IF (NSTEP.EQ.NSTOP) THEN
            CLOSE(NUNVOCEM)
            IF (LREAD_VOCEMIS) CLOSE(NUNVOCEMINP)
            CLOSE(NUNRADEM)
            CLOSE(NUNHONOEM)
         ENDIF

C     LG- end IF (LVOCEMIS)

      ENDIF

      WRITE(NUNMDFL,*)'End calculation of emission'

      RETURN
      END
