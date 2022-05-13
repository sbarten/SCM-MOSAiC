      SUBROUTINE DAYCENT(NSTEP,NSTOP,IMON,JDAY,NDAY_L,
     &                   LCHDAY,LCHWEEK,LCHMONTH,LEOMONTH)

c ----------------------------------------------------------------------
c     This program calculates the soil N emissions (and others) 
c     using the DAYCENT model. The original code is found in another
c     subdirectory, and linked to the SCM through a library produced
c     in the Daycent source directory, containing all the *.o files
c ----------------------------------------------------------------------
c 
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c                           DISCLAIMER
c
c        Neither the Great Plains System Research Unit - USDA (GPSR) nor
c     Colorado State University (CSU) nor any of their employees make
c     any warranty or assumes any legal liability or responsibility for
c     the accuracy, completeness, or usefulness of any information,
c     apparatus, product, or process disclosed, or represents that its
c     use would not infringe privately owned rights.  Reference to any
c     special commercial products, process, or service by tradename,
c     trademark, manufacturer, or otherwise, does not necessarily
c     constitute or imply endorsement, recommendation, or favoring by  
c     the GPSR or CSU.  The views and opinions of the authors do not
c     necessarily state or reflect those of GPSR or CSU and shall not 
c     be used for advertising or product endorsement. 

c      program main

c ... Century Soil Organic Matter Model
c ... Simulation of carbon, nitrogen, phosphorous, and sulfur cycling
c ... As of Dec. 1991, uses a 1 month time step
c ... Project - Soil Fertility in the Great Plains
c ... Modeler - Bill Parton
c ... Programmers - Vicki Kirchner, Becky McKeown, Laura Harding,
c ...               Melannie Hartman

c ... State variables and flows are grams/m2.

      implicit none

      include 'monprd.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'zztim.inc'

      include 'parchem.h'
      include 'daycent_SCM.h'

c ...              (unit 1) = plot/print file used by modaid (unformatted)
c ... <site>.100   (unit 7) = parameter values and initial values for 
c ...                         state variables; see subroutine sitein.
c ... fix.100      (unit 8) = fixed parameter values values for 
c ...                         state variables; see subroutine fixin.
c ...              (unit 9) = a file of weather data read in subroutines 
c ...                         wthini, weathr
c ... c14data     (unit 10) = a data file specifying years in which labeled
c ...                         carbon is added via plant growth and what 
c ...                         fraction of the growth is labeled.
c ... nflux.out   (unit 70) = N2/N2O fluxes computed by Trace Gas Model
c ... daily.out   (unit 80) = pet, defac, stemp, and snowpack water content
c ...                         computed by Trace Gas Model
c ... summary.out (unit 90) = tmax, tmin, prec, N2O flux, NO flux, and CH4
c ...                         computed by Trace Gas Model

c ... If you're getting floating point errors mentioned after you exit
c ... Century, uncomment the following lines, recompile, run Century
c ... in dbx with the 'catch FPE' option to find the offending code.
c ... You can also run Century outside of dbx, in which case you will
c ... get messages on your screen giving you an fpe error code (see
c ... the Floating Point Programmer's Guide, p20) and a not-very-
c ... useful-to-3rd-or-4th-generation-language-programmers location. 
c ... The error handler 'mysigal' is a Fortran callable C routine 
c ... written by Martin Fowler; it can be replaced by any user written
c ... handler or any of several library handlers, the most useful 
c ... probably being SIGFPE_ABORT.  The calls to ieee_handler won't 
c ... compile using poa's binaries.

c      external mysigal
c      ieeer=ieee_handler('set','invalid',SIGFPE_ABORT)
c      ieeer=ieee_handler('set','division',mysigal)
c      ieeer=ieee_handler('set','overflow',mysigal)
c      ieeer=ieee_handler('set','underflow',SIGFPE_ABORT)

c ... You probably won't want to uncomment the following line; inexact
c ... floating point signals occur all over the place.

c      ieeer=ieee_handler('set','inexact',mysigal)

c ... Local variables
      real month_vals(12), neg_month_vals(12)
      data month_vals /0.08, 0.17, 0.25, 0.33, 0.42, 0.50, 0.58, 0.67,
     &                 0.75, 0.83, 0.92, 1.0/
      data neg_month_vals /-0.92, -0.83, -0.75, -0.67, -0.58, -0.50,
     &                     -0.42, -0.34, -0.25, -0.17, -0.08, 0.0/

! mz_lg_20031202+ added

      integer nstep, nstop, imon, jday, nday_l
      logical lscm, lend, lchday, lchweek, lchmonth, leomonth

      lscm = .true.
      lend = .false.

! mz_lg_20031202-

! mz_lg_20031202+ only initialization being done first timestep
      if (nstep.eq.0) then

c ... Obtain startup information from user, do initializations based on
c ... answers to Modaid questions
      call detiv(lscm)  ! mz_lg_20031202+ modified

c ... Adjust parameters from crop.100 and fix.100 for weekly production
c ... if necessary. -mdh 1/95
      call adjustpar

c ... Write out starting values

      call wrtbin(time_scm)

      endif
! mz_lg_20031202- end initialization

c ... Update month

! mz_lg_20031202+ modified time-control using the SCM parameters to
!     define the month. The loop in csa_main is using timesteps of a 
!     month (dt=1/12) where blktnd and tend reflect the total number 
!     of years for which the simulation must be performed

      imonth_scm = INT(imon,i4)

      if (nstep .eq. 0) lchmonth = .true.

c ... If time is greater than the ending time for the current block,
c ... read the next block

      if ((abs(time_scm - blktnd_scm) .lt. (0.5 * dt_scm)) .and.
     &    (abs(time_scm - tend_scm)   .gt. (0.5 * dt_scm))) then
        call readblk(lscm)
      endif

c ... Perform annual tasks

! mz_lg_20030312+ modified by calling the subroutine eachyr and resetting
!     some accumulators whenever the simulation starts, NSTEP.EQ.0 and whenever
!     the simulation has run for one year 

      if (nstep .eq. 0 .or. mod(nday_l,365) .eq. 0) then   
        call eachyr
c ..... Reset accumulators for yearly trace gas output, cak - 09/23/02
        N2O_year = 0.0
        NO_year = 0.0
        CH4_year = 0.0
      endif

c ... The main driver for the model; call decomposition, growth, etc.

! mz_lg_20031202+ modified: this subroutine is being called here at a 
!     weekly frequence in constrast to the original daycent where it is
!     called once a month. The weekly frequency is applied to transfer
!     the (accumulated) parameters into simson and the other
!     subroutines.  

      if (nstep .EQ. 0. .OR. lchweek .OR. leomonth .OR. lchmonth)
     &      call simsom(lscm, lchweek, lchmonth, leomonth, 
     &                  nstep, jday)

c ... Write yearly trace gas output, cak - 09/23/02
      if ((time_scm .ge. strplt_scm) .and. (month .eq. 12)) then
        call wrtyearsum(time_scm, N2O_year, NO_year, CH4_year)
      endif

c ... Update time
c ... Add calculation to handle time drift caused by inexact floating
c ... point addition, cak - 08/23/01
c      time = time + dt
c ... Add calculation to handle incrementing the month for negative years,
c ... cak - 03/29/02

! mz_lg_20040323+ the time must be updated only for a change in the month.
!     This will change whenever dt_scm is changed to shorter timesteps.
!     The if then statement has to be modified consistently with any change
!     in dt_scm. When dt_scm will be reduced to 1 day then lchmonth should
!     be replaced by lchday
      
      if (nstep .gt. 0 .and. lchmonth) then
        if (time_scm .ge. 0.0) then
          !        time = int(time) + month_vals(month)
          ! mz_lg_20040322+ modified for the SCM timestepping
          time_scm = time_scm + dt_scm
        else
          time_scm = int(time_scm) + neg_month_vals(month)
          if (month .eq. 1) then
            time_scm = time_scm + 1.0
          endif
        endif
        if (time_scm .ge. -1.0e-07 .and. time_scm .le. 1.0e-07) then
          time_scm = 0.0
        endif
      endif ! endif (lchmonth) then
! mz_lg_20040323- 

c ... Write out values
      if ((tplt_scm - time_scm) .lt. (dt_scm * 0.5)) then
        call wrtbin(time_scm)
        tplt_scm = time_scm + dtpl_scm
      endif

! mz_lg_20031202+ modified by removing the off-line timeloop (see 
!     csa_main.f for original loop structure

      if ((tend_scm-time_scm) .lt. (dt_scm*.5)) then
         lend = .true.
      endif

! mz_lg_20031202+ at the end of the simulation, writing the output data
      if (nstep .eq. nstop .or. lend) then

      print *,'daycent.f: writing final values and closing files'

c ... Write out final values
      call wrtbin(time_scm)

c ... Close data files

c ... Close the weather file (unit=9)
      close(unit=nundcwth)
c ... Close the c14data file (unit=14)
      close(unit=nundcc14)
c ... Close the schedule file (unit=15)
      close(unit=nundcsch)
c ... Close N2/N2O flux file, nflux.out (unit=70)
      close(unit=nundcnflx)   ! mz_lg_20031213+ modified
c ... Close the daily.out file (unit=80)
      close(unit=nundcdail)
c ... Close the summary.out file (unit=90)
      close(unit=nundcsumm)

c ... Mark end of file
      endfile(unit=nundcbin)

c ... Close binary file (unit=1)
      close(unit=nundcbin)

      endif
! mz_lg_20031202-

      RETURN
      END
     
