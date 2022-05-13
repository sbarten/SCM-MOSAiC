      SUBROUTINE  READDATA(NSTEP,NSTOP,NDIN,NTBASE,DTIME,OBSERVM)
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
 
      INTEGER I,II,J,JJ,N,MAXDATA,NDATA,NPARAM,NSTEP,NSTOP,
     &        NDATE,NODAY,NYMDIN,NTBASE,NDIN,NYEAR,NMONTH,NDAY,
     &        NOBS,IP

      PARAMETER (MAXDATA=10000,NPARAM=1)

      REAL STARTTIME, ENDTIME, TIME(MAXDATA), HOUR, MIN, SEC,
     &     DTIME, OBSERV(MAXDATA,NPARAM), OBSERVM(NPARAM),
     &     DT

      CHARACTER*50 DUMMY
      CHARACTER*40 FILENAME

      SAVE OBSERV

C ------------------PROGRAM STARTS HERE----------------------

      IF (NSTEP.EQ.0) THEN

       WRITE(*,'(1A)')
     &  ' Constraining model with observations, start reading data file'

       OPEN(NUNDATA,STATUS='unknown')

       READ(NUNDATA,'(a50)') DUMMY

C LG-  reading the begin and end time of the observations

       READ(NUNDATA,*)NDATE,STARTTIME,ENDTIME,NODAY,DT

       NYMDIN = NDATE/100
       NTBASE = NDATE - 100*NYMDIN

       NYEAR   = NYMDIN/10000
       NMONTH  =(NYMDIN-NYEAR*10000)/100
       NDAY    = NYMDIN-NYEAR*10000-NMONTH*100

C LG-  reading the header with the names of the parameters

       READ(NUNDATA,'(A50)') DUMMY
       READ(NUNDATA,'(A50)') DUMMY

       IP=INDEX(DUMMY,'   ')

       WRITE(*,'(1A,A20,1A,A9)')
     *   ' The read datafile contains data of ',
     *     DUMMY(9:IP-1),' as a function of ',DUMMY(1:8)

       WRITE(*,'(1a,3i4)')
     *   ' The observations represent the year/month/day ',
     *     NYEAR,NMONTH,NDAY

       WRITE(*,'(1a,f6.2,1a,f6.2,1a,i3,1a,f4.0,1a)')
     *   ' The observations start at ',STARTTIME,' GMT and end at ',
     *     ENDTIME,' GMT after ',NODAY,' day(s), timestep: ',DT,' [s]'

       WRITE(*,'(1A)')' Press enter to continue'
       READ (*,*)

       I=1

C LG-  opening of file for checking the average observations 

       OPEN(UNIT=NUNOBS,FILE='/data/ganzevl/racmo/output/avg_obs.out',
     *     FORM='FORMATTED',STATUS='UNKNOWN')
       WRITE(NUNOBS,'(2a8,a12)')'nstep','time','Avg_param'

 100   CONTINUE

       READ(NUNDATA,*,ERR=200,END=300) 
     &      TIME(I),(OBSERV(I,N),N=1,NPARAM)

       II=0
 101   CONTINUE

C LG-  the next statements are to deal with missing data in the 
C      input datafile, in case of missing data the parameter OBSERV
C      has been assigned a value of 9999.9999 which is then used in
C      the assignment of the observaed value to a model resolved 
C      parameter to determine each timestep if there are available data

       IF (I.GT.1.AND.(TIME(I)-(TIME(I-1)+II*DT)).GT.DT) THEN
        DO N=1,NPARAM
	 OBSERV(I+II,N)=-999.999
        ENDDO
        II=II+1
	GOTO 101
       ELSE
        IF (II.GT.0) TIME(I+II)=TIME(I)            
        GOTO 102
       ENDIF

 102   CONTINUE

       I=I+II+1
      
       GOTO 100      

 200   PRINT *,'error reading data file, check file'  ! in case of error
       STOP

 300   NDATA=I-1                            ! IN CASE OF EOF
       I=0

       CLOSE(NUNDATA,STATUS='keep')
           
      ENDIF
      
C LG- calculation of the timestep mean parameter value 
      
      IF (NDIN.EQ.NDAY.AND.NTBASE.GE.INT(STARTTIME).AND.I.EQ.0)
     &  I=1

      IF (I.GT.0) THEN

       DO 399 N=1,NPARAM
        OBSERVM(N)=0.
  399  CONTINUE

       NOBS=DTIME/DT
       JJ=0
       DO 400 N=1,NPARAM
       DO 400 J=1,NOBS
        II=(I-1)*NOBS+J

C LG-   In case of missing data, these are not considered for the averaging.
C       If there are less than half the number of available observations
C       within of time average interval then the parameter OBSERV has been 
C       assigned the value 9999.9999

        IF (OBSERV(II,N).LE.-999.999) GOTO 4000
        JJ=JJ+1
        OBSERVM(N)=OBSERVM(N)+OBSERV(II,N)
 4000   CONTINUE

        IF (J.EQ.NOBS) THEN
         IF (JJ.LT.NOBS/2) THEN
	  OBSERVM(N)=-999.999
	 ELSE
	  OBSERVM(N)=MAX(0.1,OBSERVM(N)/(JJ/NPARAM))
	 ENDIF
	ENDIF

  400  CONTINUE

       WRITE(NUNOBS,'(2(1x,i7.7),f12.4)')
     &     NSTEP,LDAYLTIME,(OBSERVM(N),N=1,NPARAM)

C LG-  in case of an integration which covers a longer time period than
C      that of the observations, the user is warned and the model continues
C      with its resolved parameter value.    

       IF (I.GT.0) I=I+1

       IF (II.GE.NDATA) THEN
        WRITE(NUNMDFL,'(1a)')
     &   ' The intergration time exceeds that of the',
     &   ' available observation range'
        WRITE(NUNMDFL,'(1a)')
     &    ' The model continues without using observed data'
        WRITE(NUNMDFL,'(1a)')
     &    ' Potential error source due to too lacking data'

C LG-   writing to screen

        WRITE(*,'(1a)')
     &    ' Potential error source due to longer intergration than',
     &    ' available data !'

        DO 402 N=1,NPARAM
         OBSERVM(N)=-0.9999
  402   ENDDO

        CLOSE(NUNOBS)

       ENDIF
       
      ENDIF

      END
