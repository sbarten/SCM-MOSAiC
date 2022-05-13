      SUBROUTINE  READDATA(NSTEP,NSTOP,NPRINT,NYLEN,NDIN,NCBASE,
     &                     NPARAM_MAX,DTIME,OBSERVM)
C ----------------------------------------------------------------------
C  This subroutine reads in a file containing a observed parameter
C  values, the format of reading needs to be adapted to the specific
C  file structure. An example of the format of the file for which this 
C  code has been written can be found in the directory: 
C  racmo/echam/input/data/ 
C
C  Written by Laurens Ganzeveld, 26-01-99
C ----------------------------------------------------------------------

      IMPLICIT NONE
  
      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'
 
      INTEGER I,N,J,NDATA_MAX,NPARAM_MAX,NPARAM,NSTEP,NSTOP,NPRINT,
     &        NDATE,NODAY,NYMDIN,NTBASE,NDIN,NYEAR,NMONTH,NDAY,
     &        NYLEN,NCBASE,IDAT2C,JDAY_OBS,NOBS,IP,IIP,NDTG

      PARAMETER (NDATA_MAX=100000) ! ESS_lg_20100625+ increased to deal with larger files

      INTEGER NDATA,II

      REAL STARTTIME, ENDTIME, TIME(NDATA_MAX), 
     &     DTIME, OBSERV(NDATA_MAX,10), 
     &     OBSERVM(10),
     &     HOUR, MIN, SEC, DT, DTIME_INT_SEC, DTIME_OBS

      CHARACTER*100 DUMMY
      CHARACTER*100 FNAME,FNAME1,FNAME2,FILENAME,DIRNAME
      CHARACTER*15  PARNAME(0:NPARAM_MAX)

      SAVE TIME,NDATE,NPARAM,OBSERV, ! ESS_lg_20090811+ added to store NDATE,
     &     II,NDATA                  !      NPARAM, OBSERV, etc.

C ------------------PROGRAM STARTS HERE----------------------

      IF (NSTEP.EQ.0) THEN
	WRITE(*,'(1a)')
     &  ' Constraining model with observations, start reading data file'
	WRITE(*,'(1a)')
     &  ' For assigning the obs. to the model parameters see oned.f'
	WRITE(*,'(1a)')
     &   ' Give name of directory containing file(s) (default dir=case/input/observ)'
	READ(*,*)DIRNAME
	WRITE(*,'(1a)')
     &   ' Give name of file containing data (default dir=case/input/observ)'
	READ(*,*)FILENAME

	FNAME='input/observ/'
	IP=INDEX(FNAME,' ')
	FNAME1=FNAME(1:IP-1)//DIRNAME
	IP=INDEX(FNAME1,' ')
	WRITE(DUMMY,'(a1)')'/'
	FNAME2=FNAME1(1:IP-1)//DUMMY
	IP=INDEX(FNAME2,' ')
	FNAME=FNAME2(1:IP-1)//FILENAME
	PRINT *,FNAME

	OPEN(NUNDATA,FILE=FNAME,STATUS='unknown')

	READ(NUNDATA,'(a100)') DUMMY

C LG-   reading the begin and end time of the observations

	READ(NUNDATA,*)NDATE,STARTTIME,ENDTIME,NODAY,DT

C LG-   calculation of some time parameters of the observations

	NYMDIN = NDATE/100
	NTBASE = NDATE - 100*NYMDIN
	NYEAR  = NYMDIN/10000
	NMONTH =(NYMDIN-NYEAR*10000)/100
	NDAY   = NYMDIN-NYEAR*10000-NMONTH*100

C LG-   calculation of the Julian at which the observations start
C       (original parameter name for Julian day is NCBASE!)

	IF (NYLEN.EQ.365) CALL INIYL2
	NYMDIN = NDATE/100
	NTBASE = NDATE - 100*NYMDIN
C
	NYEAR   = NYMDIN/10000
	NMONTH  =(NYMDIN-NYEAR*10000)/100
	NDAY    = NYMDIN-NYEAR*10000-NMONTH*100
	IF (NYLEN.EQ.360) THEN
         JDAY_OBS = IDAT2C(NDAY,NMONTH,NYEAR)
	ELSE IF (NYLEN.EQ.365) THEN

C LG-    originally, the programm calls the julian day NCBASE
C        in the subroutine IDAT with the minimal value of zero instead of 1
C        which yields a negative julian day, this still needs to adressed
C        more specifically

         CALL IDAT(1,NMONTH,NDAY,JDAY_OBS)
	END IF

        NDTG=NDTG_ACT

C LG-   reading the header with the names of the parameters

	READ(NUNDATA,'(A100)') DUMMY
	READ(NUNDATA,'(A100)') DUMMY

	N=0
	IIP=0
 77     CONTINUE        
	IP=INDEX(DUMMY(IIP+1:100),' ')
	PARNAME(N)=DUMMY(IIP+1:IIP+IP-1)
	IF (DUMMY(IIP+1:IIP+IP+1).EQ.' ') GOTO 88
	IIP=IIP+IP
	N=N+1
	GOTO 77
 88     CONTINUE

	NPARAM=N-1

	WRITE(*,'(1a,i3,1a)')
     *   ' The read datafile contains data of: ',NPARAM,' parameters'
	DO N=1,NPARAM
         WRITE(*,'(1a)') PARNAME(N)
	ENDDO
	WRITE(*,'(1A,A10)'),' as a function of ',PARNAME(0)

	WRITE(*,'(1a,4i4)')
     *   ' The observations represent the year/month/day/Julian day ',
     *     NYEAR,NMONTH,NDAY,JDAY_OBS

	WRITE(*,'(1a,f6.2,1a,f6.2,1a,i3,1a,f6.0,1a)')
     *   ' The observations start at ',STARTTIME,' GMT and end at ',
     *     ENDTIME,' GMT after ',NODAY,' day(s), timestep: ',DT,' [s]'

	WRITE(*,'(1A)')' Press enter to continue'
	READ (*,*)

	I=1

C LG-   opening of file for checking the average observations 

	OPEN(UNIT=NUNOBS,FILE='/data/ganzevl/racmo/output/avg_obs.out',
     &     FORM='FORMATTED',STATUS='UNKNOWN')
	WRITE(NUNOBS,'(a10,a15,10(3x,a12))')
     &    'nstep','time',(PARNAME(N),N=1,NPARAM)

 100    CONTINUE

	READ(NUNDATA,*,ERR=200,END=300) 
     &      TIME(I),(OBSERV(I,N),N=1,NPARAM)

!        print *,'readdata: ',i,time(i),observ(i,1),observ(i,nparam)

        I=I+1  ! ESS_lg_20100625+ 

	GOTO 100      

 200    PRINT *,'error reading data file, check file'  ! in case of error
	STOP

 300    NDATA=I-1                            ! IN CASE OF EOF
	II=1

	CLOSE(NUNDATA,STATUS='keep')

      ENDIF
      
C LG- calculation of the time in second since the start of the 
C     simulation in order to apply the proper observations, which
C     should be defined in default format in seconds since the
C     initial time STARTTIME

C       DTIME_INT_SEC=(NTBASE-STARTTIME)*3600.+MAX(0.,NSTEP*DTIME-
C      *   MAX(0.,FLOAT(JDAY_OBS-NCBASE))*86400.)

C LG- 200301, this has been introduced to be able to read in the input
C     files that Eric Simon has produced containing the EUSTACH surface
C     observations. It stil must be checked if the same code also works
C     reading in the files that were originally being read in (NO2 conc.
C     see directories of Harvard forest data)

      DTIME_INT_SEC=NSTEP*DTIME

C LG- 200301, in case of a database with the starting time being later 
C     then the initial reference time of the simulation 

      IF (NDATE.GE.NDTG) THEN  ! ESS_lg_20100625+ GE instead of GT
        IF (NDATE.GT.NDTG_ACT) THEN ! ESS_lg_20100625+ GT instead of GE, check if this
                                    ! also works properly for a dT of 3600s or longer
          DO I=1,NDATA
            TIME(I)=TIME(I)+DTIME
          ENDDO
        ENDIF
      ENDIF  

C LG- resetting the model applied value

      DO 398 N=1,NPARAM
        OBSERVM(N)=-9999.9
  398 CONTINUE

C LG- calculation of the timestep average parameter value
           
C LG- 200301, in case of an intial time of available input data that
C     is larger then the actual time in the integration, the 
C     interpolation is not being done.

      IF (NDATE.GT.NDTG_ACT) GOTO 402 ! ESS_lg_20100625+ GT instead of GE, check if this
                                      ! also works properly for a dT of 3600s or longer

      DO 401 N=1,NPARAM

        IF (II.GT.NDATA) GOTO 400
        IF (DTIME_INT_SEC.GE.TIME(II)) THEN

          DTIME_OBS=TIME(II+1)-TIME(II)

          IF (OBSERV(II,N).LE.1e-10.OR.OBSERV(II+1,N).LE.1e-10) ! ESS_lg_20100626+ modified to avoid using negative obs.
     &      GOTO 399

          IF (II+1.LE.NDATA.AND.DTIME_OBS.GT.0.) THEN
            OBSERVM(N)=OBSERV(II,N)*
     &        (1.-((DTIME_INT_SEC-TIME(II))/DTIME_OBS))+
     &        OBSERV(II+1,N)*(1.-((TIME(II+1)-DTIME_INT_SEC)/DTIME_OBS))
          ELSE
            OBSERVM(N)=OBSERV(II,N)
          ENDIF
 399      CONTINUE
          IF (DTIME_INT_SEC+DTIME.GE.TIME(II+1).AND.N.EQ.NPARAM) II=II+1
        ENDIF

  400   CONTINUE
  401 CONTINUE
  402 CONTINUE

      IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0)
     & WRITE(NUNOBS,'(1x,i9.9,1x,a14,10(3x,f12.4))')
     &     NSTEP,LDATLTIME,(OBSERVM(N),N=1,NPARAM)

      IF(NSTEP.EQ.NSTOP) CLOSE(NUNOBS)

      END
