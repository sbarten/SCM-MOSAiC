      PROGRAM DAYCENT

      implicit none

! mz_lg_20041209+ extended list of *.inc files

      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/const.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/dovars.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/jday.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/monprd.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/param.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/plot1.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/plot3.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/t0par.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/timvar.inc'
      include '/data/ganzevl/models/DayCent/DayCent4.5/DayCent/zztim.inc'


      logical lscm

      lscm = .true.

c ... Obtain startup information from user, do initializations based on
c ... answers to Modaid questions
      call detiv(lscm)  ! mz_lg_20031202+ modified

      print *,'end call detiv:',lscm

      END
     
