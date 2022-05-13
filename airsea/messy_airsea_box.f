! Author:
! Laurens Ganzeveld,     MPICH, 2004

    MODULE scm_airsea_box

      USE messy_airsea
      USE messy_main_constants_mem, ONLY: dp

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: messy_airsea_initialize
      PUBLIC :: messy_airsea_vdiff

    CONTAINS

      ! --------------------------------------------------------------------------

      SUBROUTINE messy_airsea_initialize(ztmst) ! ESS_lg_20150918+

	IMPLICIT NONE
	INTEGER                     :: status ! error status

! ESS_lg_20100726+ 
        REAL ZTMST

	! read CTRL_AIRSEA namelist
	CALL airsea_read_nml_ctrl(status, 99)
	IF (status /= 0) STOP

      END SUBROUTINE messy_airsea_initialize

      ! --------------------------------------------------------------------------

      SUBROUTINE messy_airsea_vdiff

	IMPLICIT NONE
	INTEGER            :: i, ii, iii
	INTEGER, PARAMETER :: nproma=3478,nasi=17,jrow=1

	REAL(dp), PARAMETER ::    vkarman=0.41  !von Karman constant [-]

	CHARACTER(LEN=*), PARAMETER :: substr='airsea_vdiff'

	LOGICAL  :: l_sea(nproma)                    ! orography: TRUE when water (no ice, no land!)
	REAL(dp) :: schmidt_air(nproma,nasi)         ! schmidt number (adimensional) in air 
	REAL(dp) :: schmidt_sea(nproma,nasi)         ! schmidt number (adimensional) in sea
	REAL(dp) :: flux_airsea(nproma,nasi)         ! flux (mol(trac) m-2 s-1)
	REAL(dp) :: densflux_airsea(nproma,nasi)     ! flux (mol(trac) kg(air) mol-1(air) m-2 s-1)
	REAL(dp) :: zflux_airsea(nproma,nasi)        ! flux (mol(trac) mol-1(air) s-1)
	REAL(dp) :: zdensair(nproma)                 ! air density!!! [kg/m3]
	REAL(dp) :: temp(nproma), sphum(nproma)      ! temperature & humidity!
	REAL(dp) :: molweight(nasi)                  ! molecolar weight
	REAL(dp) :: conc_air(nproma, nasi)           ! mixing ratio air 
	REAL(dp) :: conc_water(nproma, nasi)         ! mixing ratio water
	REAL(dp) :: delta_conc(nproma, nasi)         ! Delta_conc (Cw - Cx kh P) 
	REAL(dp) :: zdz(nproma)                      ! height of the lowest level 
	REAL(dp) :: henry_val(nproma,nasi)           ! henry value for tracer at that temperature 
	REAL(dp) :: ztmst                            ! timestep
	REAL(dp) :: cair(nproma)                     ! no lev--> we look only at surface: air concentration
                                                     ! cair => [mol(air)/m3]
	REAL(dp) :: o3_output(nasi,nproma)           ! mz_lg_20050826+ extra output               

	! mz_lg_20050825+
	REAL(dp) :: Im(nproma), C2H4(nproma), C3H6(nproma), NO3m(nproma), &
                    DMS(nproma), chlorofyll(nproma), kORG

	REAL(dp), DIMENSION(nproma) :: &
	  ! input (general)
	  sst            , & ! Sea Surface Temperature [K]
	  wind10         , & ! 10m wind speed, u component [m s-1]
	  press          , & ! surface pressure [mb]
	  zust           , & ! friction velocity [m s-1]
	  srfl           , & ! incoming net radiation [W m-2]
	  ahfs           , & ! sensible heatflux W/m**2
	  ahfl           , & ! latent heatflux W/m**2
	  gyre           , & ! Gyre index [0-1]
	  frac_coast     , & ! fraction of coastal waters [0-1]
	  z0mw           , & ! water surface roughness [m], calculated according to Charnock   
	  zref           , & ! surface layer reference height 
	  um1sl          , & ! x-component windspeed at reference height SL
	  vm1sl          , & ! y-component windspeed at reference height SL
	  geopotsl       , & ! geopotential height SL reference height
	  pblh               ! MBL height

	REAL(dp), DIMENSION(nasi,nproma) :: &
	  zoutput

	REAL(dp), DIMENSION(nproma) :: &
	  VdO3_ocean        ! mz_lg_20060218+ oceanic O3 dry deposition

	! mz_lg_20040721+ and parameters needed to calculate the ocean roughness 
	!    according to Charnock and the roughness weighted with all the surface
	!    cover fractions through the drag coeff. for surface cover fraction
	REAL(dp) :: zvisc, zustarn, zustarm, zcdragw

    ! ESS_ac_200811+ adding the header of the file
	CHARACTER*250 header,fname
	REAL(dp),DIMENSION(nproma) :: &
	  jday

    ! ESS_ac_20081111-

    ! ESS_ac-20090119+ extra parameters for reading observation files

	INTEGER  :: nobs_chl, nobs_NO3
	INTEGER, PARAMETER:: nobs_max=1e5
	REAL(dp) :: latitude, longitude, lat_obs(nproma), lon_obs(nproma),   &
                    lat_obs_chl(nproma),lon_obs_chl(nproma),chl_obs(nproma), &
                    a2,b2,c2,distance,distance_old,                          &
                    date(nobs_max),time(nobs_max),                           &
                    lat_obs_chl_salis(nobs_max),lon_obs_chl_salis(nobs_max), & 
                    salinity(nobs_max),sstC(nobs_max),                       &
                    chl_obs_salis(nproma),fcdm(nobs_max),                    &
                    oxsat(nobs_max),                                         &
                    lat_obs_NO3(nobs_max),lon_obs_NO3(nobs_max),& 
                    NO3_obs(nobs_max),NO3flag(nobs_max)
	LOGICAL  :: l_chl_RS  

    ! ESS_ac_20090119-

    !--- 0. Initialisations: -----------------------------------------------

	! -------------------------------------------------------------------------------
	! mz_lg_20040721+ start defining input parameters for the various subroutines

	l_sea(:)=.TRUE.
	l_chl_RS=.false.  ! ESS_lg_20100719+

	! mz_lg_20060217+
	VdO3_ocean(:)     = 0.045           ! O3 oceanic dry deposition velocity according to Fairall et al, 2005 
	! mz_lg_20060217-

	zdz(1:nproma) = 64.   ! layer thickness
	zref(1:nproma) =10.   ! reference height windspeed

	DO i=1,nproma
	   sst(i)    =   293.  ! 273.15+30.*COS(REAL(0.5*3.14*(nproma-i)/nproma))
	   srfl(i)   =   500.  ! REAL(i)*10. ! from 0 up to N*10 W m-2
	   temp(i)   =   293.  ! 273.15+30.*COS(REAL(0.5*3.14*(nproma-i)/nproma)) 
	   wind10(i) =  0.5+29.5*REAL(i)/nproma   ! 7.5
	   um1sl(i)  =  0.5+29.5*REAL(i)/nproma   ! 7.5    ! [m s-1]
	   vm1sl(i)  =   0.    ! [m s-1]
	   geopotsl(i)= 32.*9.8! geopotential height surface layer
	   sphum(i)  =   0.010 ! specific humidity  [g g-1]
	   ahfs(i)   =   50.   ! sensible heat flux [W/m-2]
	   ahfl(i)   =  200.   ! latent heat flux [W m-2]
	   press(i)  =  101.   ! surface pressure [mb]
	   pblh(i)   =  600.
	   ! --------------------------------------------------------------------------
	   ! mz_lg_20040721+ added the calculation of the roughness length over
	   !   water wheneveR zfrw > 0. This is required since the roughness depends
	   !   on the windspeed/friction velocity over water. The calculation requires
	   !   an iteration, which generally uses 2-3 steps to converge

	   zvisc=0.15 ! cm2 s-1
	   zustarn=wind10(i)/50. ! first-order estimate of u*
	   DO 1100 ii=1,10
               zustarm=zustarn
               z0mw(i)=0.018*(zustarm**2)/9.8  ! Charnock formula
               IF (z0mw(i) < 1.5e-5) z0mw(i) = 1.5e-5
               zcdragw=(vkarman/(LOG(zref(i)/z0mw(i))))**2
               zustarn=SQRT(zcdragw)*wind10(i)
               IF (((zustarn-zustarm)/zustarn).LT.0.0001) GOTO 1200
     1100  CONTINUE
     1200  CONTINUE

	   zust(i)=zustarn

	END DO

	! mz_lg_20050826+ intialization of tracers 

	! mz_lg_20050828+ for diagnostics and calculations, resetting oceanic
	!    NO3 concentrations over land
	NO3m(1:nproma) = 0.1764_dp ! 5 [umol NO3m l-1 ocean water], 9.8 gives Im of ~50 nmol I- l-1, 0.1764nM gives 100nM
	DMS(1:nproma)  = 0.5_dp    ! ESS_lg_20100716+, 200._dp   ! 200 nmol DMS l-1 ocean water]

	IM(1:nproma)=20.   ! prescribed fixed oceanic I- concentration [nmol l-1]
	C2H4(1:nproma)=0.5 ! prescribed fixed oceanic C2H4 concentration [nmol l-1]
	C3H6(1:nproma)=0.5 ! prescribed fixed oceanic C3H6 concentration [nmol l-1]

    ! ESS_20090119+ activate to overwrite read-in CHL observations
	chlorofyll(1:nproma)=1 ! prescribed fixed oceanic Chlorofyll concentration [ug l-1]!

	kORG          =6000. !=6000. ! organic reactivity contr. NOTE the difference units [s-1]

    ! ESS_ac_20081111+ added the reading of the input file containing the
    !   observations of the various campaings. NPROMA reflects in this case the
    !   number of datapoints in the input file

	! open file
	FNAME='input/GOMECC/expgomecc1119_new.txt'
	OPEN(unit=2,file=FNAME,form='FORMATTED',status='UNKNOWN')

	READ(2,'(1a)') header

	DO i=1,nproma

    ! ESS_lg_20100716+ note that the wind10 parameter is assigned the measured wind
    !       speed at a referenced height 18m!

            READ(2, '(2(e12.6,1x),e13.6, 6(e13.6), e13.6, 1x, e13.6)') &
               jday(i), lat_obs(i), lon_obs(i),srfl(i), sst(i), temp(i),wind10(i), zust(i), sphum(i), &
	       ahfs(i), ahfl(i)

    !        print *,'test reading',i,jday(i), lat_obs(i), lon_obs(i),srfl(i), &
    !               sst(i), temp(i),wind10(i),&
    !       	zust(i), sphum(i), ahfs(i), ahfl(i)

            sst(i)=sst(i)+273.15 ! recalculation to K
            temp(i)=temp(i)+273.15

    !	READ (*,*)

	END DO

	CLOSE(2)

    ! ESS_lg_20100719+ reading in or in-situ or remote sensing based CHL data

    ! ESS_lg_20100719+ added the reading of the GOMECC CHL observations to
    !    constrain the model, datasource: Salisbury dataset

	IF (.NOT.l_chl_RS) THEN

	  ! open file
	  FNAME='input/GOMECC/GOMECC-Salisbury-Chlorophyll-UW.txt'
	  OPEN(unit=2,file=FNAME,form='FORMATTED',status='UNKNOWN')

	  READ(2,'(1a)') header

	  i=1
    101   CONTINUE
	  READ(2,*,ERR=201,END=301) &
             time(i),lat_obs_chl_salis(i),lon_obs_chl_salis(i),&
             salinity(i),sstC(i),chl_obs_salis(i),fcdm(i),oxsat(i)

    !       print *,'reading ',&
    !           time(i),lat_obs_chl_salis(i),lon_obs_chl_salis(i),&
    !           salinity(i),sstC(i),chl_obs_salis(i),fcdm(i),oxsat(i)
    !       read (*,*)

	  i=i+1
	  GOTO 101

    201   PRINT *,'Error reading data file, check file'  ! in case of error
	  STOP

    301   NOBS_CHL=I-1                            ! IN CASE OF EOF

	  CLOSE(2)

    ! ESS_lg_20100719+ assigning the observations to the points for calculation
    !   selection of the longitude and latitude, modified

	  DO i=1,nproma
            latitude=lat_obs(i)
            longitude=lon_obs(i)
            iii=0
            distance_old=1e5
            DO ii=1,NOBS_CHL 
              b2=(latitude-lat_obs_chl_salis(ii))**2
              a2=(longitude-lon_obs_chl_salis(ii))**2
              c2=a2+b2
              distance=sqrt(c2) 
              IF (distance.LT.distance_old) THEN
		iii=ii
		distance_old=distance
              ENDIF	  
            ENDDO

            chlorofyll(i)=chl_obs_salis(iii)  

    !        print *,'chlorophyll',i,latitude,longitude,&
    !             iii,lat_obs_chl_salis(iii),lon_obs_chl_salis(iii),chl_obs_salis(iii)

    !        read (*,*)

	  ENDDO

	ENDIF

    ! ESS_lg_20100719-

    ! ESS_ac_20090119+ added the reading of the satellite CHL observations to
    !    constrain the model, datasource: NASA website

	IF (l_chl_RS) THEN

	  ! open file
	  FNAME='input/GOMECC/chlorophyll_july.txt'
	  OPEN(unit=2,file=FNAME,form='FORMATTED',status='UNKNOWN')

	  READ(2,'(1a)') header
	  READ(2,'(i3)') nobs_chl

	  DO i=1,nobs_chl
             READ(2, '(f4.1,1x,f5.1,1x,f8.6)') &
               lat_obs_chl(i), lon_obs_chl(i), chl_obs(i)
	  END DO

	  CLOSE(2)

    ! ESS_ac_20090119+ assigning the observations to the points for calculation
    !   selection of the longitude and latitude

	  DO i=1,nproma
            latitude=lat_obs(i)
            longitude=lon_obs(i)
            iii=0
            distance_old=1e5
            DO ii=1,NOBS_CHL 
              b2=(latitude-lat_obs_chl(ii))**2
              a2=(longitude-lon_obs_chl(ii))**2
              c2=a2+b2
              distance=sqrt(c2) 
              IF (distance.LT.distance_old) THEN
		iii=ii
		distance_old=distance
              ENDIF	  
            ENDDO

            chlorofyll(i)=chl_obs(iii)  

    !      print *,'chlorophyll',i,latitude,longitude,ii,lat_obs_chl(iii),lon_obs_chl(iii),iii
    !      read (*,*)

	  ENDDO

	ENDIF

    ! ESS_lg_20100719+ reading in other in-situ measurements of oceanic biogeochemical properties

	! open file
	FNAME='input/GOMECC/GOMECCBottle_NO3_15042011.txt'
	OPEN(unit=2,file=FNAME,form='FORMATTED',status='UNKNOWN')

	READ(2,'(1a)') header

	i=1
    111 CONTINUE
	READ(2,*,ERR=211,END=311) &
             time(i),lat_obs_NO3(i),lon_obs_NO3(i),&
             NO3_obs(i),NO3flag(i)

    !    print *,'reading ',&
    !         i,time(i),lat_obs_NO3(i),lon_obs_NO3(i),&
    !         NO3_obs(i),NO3flag(i)

    !    read (*,*)

	i=i+1
	GOTO 111

    211 PRINT *,'Error reading data file, check file'  ! in case of error
	  STOP

    311 NOBS_NO3=I-1                            ! IN CASE OF EOF

	CLOSE(2)

    ! ESS_lg_20100719+ assigning the observations to the points for calculation
    !   selection of the longitude and latitude, modified

	DO i=1,nproma
	  latitude=lat_obs(i)
	  longitude=lon_obs(i)
	  iii=0
	  distance_old=1e5
	  DO ii=1,NOBS_NO3 
            b2=(latitude-lat_obs_no3(ii))**2
            a2=(longitude-lon_obs_no3(ii))**2
            c2=a2+b2
            distance=sqrt(c2) 
            IF (distance.LT.distance_old.AND.NO3_obs(ii).GT.0.) THEN
	      iii=ii
	      distance_old=distance
            ENDIF	  
	  ENDDO

	  NO3m(i)=NO3_obs(iii)  ! ESS_lg_20100719+ NO3_obs is in umol kg-1 and must be recalculated
                        	!   to umol l-1

    !      print *,'NO3',i,latitude,longitude,&
    !           iii,lat_obs_NO3(iii),lon_obs_NO3(iii),NO3m(i)
    !      read (*,*)

	ENDDO

    ! ESS_lg_20100719-


    ! ESS_ac_20090119-

	WRITE(*,'(1a)')' Sucessfully read-in the data'
	WRITE(*,*)' ENTER to continue'
	READ (*,*)

    ! ESS_ac_20001111-

	!--- 1) Calculate ozone oceanic dry deposition velocities -----------------------

	CALL o3_deposition(l_sea(1:nproma),sst(1:nproma), wind10(1:nproma),   &
	   um1sl(1:nproma),vm1sl(1:nproma),geopotsl(1:nproma),                &
	   press(1:nproma), temp(1:nproma), zust(1:nproma), srfl(1:nproma),   &
	   sphum(1:nproma), ahfs(1:nproma), ahfl(1:nproma), gyre(1:nproma),   &
	   frac_coast(1:nproma), chlorofyll(1:nproma),pblh(1:nproma),         &
	   Im(1:nproma), NO3m(1:nproma), DMS(1:nproma),C2H4(1:nproma),        &
	   C3H6(1:nproma), kORG, nproma, zoutput(1:nasi,1:nproma),  &
	   VdO3_ocean(1:nproma)) ! mz_lg_20050825+ 

	! show results
	! open file
	FNAME='output/GOMECC_new.out'
	OPEN(unit=3,file=FNAME,form='FORMATTED',status='UNKNOWN')

	WRITE(3,'(16A15)') &
	  'Time','windspeed','u*','Ra','Rw','O3_Phi','SST','Henry', &
	  'NO3','Im','CHL','Aoz','Aoz_I','Aoz_DMS','Aoz_CHL','VdO3_ocean'

	WRITE(3,'(16A15)') &
        	'jday','m s-1','m s-1','s m-1','s m-1','s m-1','K','[-]',  &
		'ug l-1','nM','ug l-1','s-1','s-1','s-1','s-1','cm s-1'

	DO i=1,nproma
            WRITE(3,'(1x,f14.3,15(1x,E14.6))') &
               jday(i),wind10(i), zust(i),zoutput(2,i),zoutput(3,i), zoutput(6,i),&
		 sst(i),zoutput(17,i),zoutput(8,i),zoutput(7,i),zoutput(16,i),&
        	 zoutput(10,i),zoutput(11,i),zoutput(12,i),zoutput(13,i),&
		 100.*VdO3_ocean(i)
	END DO

	CLOSE(3)

	WRITE(*,'(1a)')'ENTER to continue'
	READ (*,*)

       END SUBROUTINE messy_airsea_vdiff

    END MODULE scm_airsea_box

    !*****************************************************************************

    PROGRAM messy_airsea

      USE messy_airsea_box, ONLY: messy_airsea_initialize, messy_airsea_vdiff

      IMPLICIT NONE

      CALL messy_airsea_initialize ! read CTRL namelist
      CALL messy_airsea_vdiff      ! calculate airsea exchanges

    END PROGRAM messy_airsea

    !*****************************************************************************
