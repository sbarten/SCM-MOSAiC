      PROGRAM RandT_MECCA1_tracerlist
C ----------------------------------------------------------------------
C  This subroutine Reads in the MECCA1 species.inc file and Transforms this 
C  into a format that can easily be implemented in oned.f and help 
C  producing the IDL tracer list
C
C  to compile; ifort -extend_source RandT_MECCA1_tracerlist.f
C
C  Written by Laurens Ganzeveld, 21-03-2011
C ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NTR,IP

      CHARACTER*100 DUMMY
      CHARACTER*250 DIRNAME,FNAME,FNAME1,TRNAME

C ------------------PROGRAM STARTS HERE----------------------

      DIRNAME='/data/ganzevl/models/echam/echam5/messy/sources/v5301_messy_1.9c/messy/mbm/mecca1/boxmodel/'
      FNAME='species.inc'
      IP=INDEX(DIRNAME,' ')
      FNAME1=DIRNAME(1:IP-1)//FNAME
      PRINT *,FNAME1

      OPEN(2,FILE=FNAME1,STATUS='unknown')

C LG- reading the header with the names of the parameters
      READ(2,'(a100)') DUMMY

      NTR=1
 100  CONTINUE
      READ(2,'(a100)',ERR=200,END=300) DUMMY

! ESS_lg_20110321+ to write the tracer list as used in oned.f

      IF (NTR.EQ.1) OPEN(3,FILE='tracerlist_MECCA1_oned.out',STATUS='unknown')
      IP=INDEX(DUMMY,')')
      TRNAME=DUMMY(9:IP-1)
      WRITE(3,'(2a)')'     & ,',TRNAME

! ESS_lg_20110321+ to write the tracer list as used in IDL

      IF (NTR.EQ.1) OPEN(4,FILE='1d_chem_mecca1_trc.out',STATUS='unknown')
      IP=INDEX(DUMMY,')')
      TRNAME=DUMMY(13:IP-1)
      IP=INDEX(TRNAME,' ')
      WRITE(4,'(3a)')'1 0  0.0  1.0 4 1 ',TRNAME(1:IP-1),'                   [ppbv]'

! ESS_lg_20110321-

      NTR=NTR+1
      GOTO 100

 200  PRINT *,'error reading data file, check file'  ! in case of error
      STOP

 300  CONTINUE

      CLOSE(2,STATUS='keep')
      CLOSE(3,STATUS='keep')
      CLOSE(4,STATUS='keep')

      END
