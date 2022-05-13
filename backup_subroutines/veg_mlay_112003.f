      SUBROUTINE VEG_MLAY (NSTEP,NSTART,NSTOP,NPRINT,NLEV,NLEVEL,
     &               PTMST,DTIME,PM,PMLOC,PRHOA,PDP,PZ,PZHP,PCATYPE,
     &               SURFFLUX,XTEEMIS,XTEDRYD,XTEDIFF)

c ----------------------------------------------------------------------
c     This program calculates the influnce of within-canopy biogeochemical
c     interactions on the canopy top exchange fluxes only considering two 
c     layers within the canopy. This routine is meant to be implemented in 
c     large scale models. In this subroutine the emission within the canopy 
c     and the dry deposition process as a function of plant activity are 
c     considered. Then the flux between the lower and top canopy layer and 
c     consecutively the canopy top flux are calculated from the gradient 
c     between the vegetation layers and the surface layer concentrations 
c     and the eddy diffusivity. The diagnostic vegetation layer 
c     concentrations are then updated within the chemistry routine to 
c     correct for chemical production or destruction. The update of the 
c     surface layer concentrations is originally done with the subroutines 
c     EMISS.f and DRYDEP.f, however, this is now done within this routine. 
c
c     Laurens Ganzeveld 1999
c ----------------------------------------------------------------------
c 

      IMPLICIT NONE

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'
      INCLUDE 'comddrs.h'

C LG- adding the common block file with some used constants for the 
C     calculation of the phsysiological stomatal resistance. The compensation
C     point for CO2 is distracted from this file. paramveg.h must also be
C     included.

      INCLUDE 'paramveg.h'
      INCLUDE 'comjac.h'

      INTEGER NNFIELD
      PARAMETER (NNFIELD=6)  ! number of input fields for interpolation

      REAL EMIFACSL,EMIFAC_UNDER,PTMST,DTIME,G,QV,VX,VD,
     &     HEIGHTDD,DDFRAC,RTOT,DZ,DZ1,DZ2,DZ3,TERM1,TERM2,
     &     ALPHAPOS,ALPHANEG,CA,CB,CC,CD,PMVEG_OLD,PM_OLD,
     &     DMASS_VEG,DMASS_SL,DMASS_EM,MASS_VEG,DMASS,MASS,
     &     MASS_NEW,MASS_OLD,BCONC0,CANHEIGHT,ZCANOPY,HCMIN,
     &     RSVEG,RSWET,KHMIN,UAVG,RAH_TOT,SUMU,LADSUM,
     &     PTMST_SUB,PTMST_MIN,TS_DD,TS_EM,TS_TURB,PMNOX_OLD
     
      REAL PM(NLON,NLEVT,NTRAC),PMLOC(NLON,NLEVT,NG3X),
     &     PRHOA(NLON,NLEVT),PDP(NLON,NLEVT),
     &     EMIFACVEG(NLEVV_ML),DC(NTRAC),FLUXCTOP(NLON,NTRACT),
     &     FLUXVEG(NLON,NLEVV_ML+1,NTRAC),SURFFLUX(NLON,NTRAC),
     &     CRF(NTRACT),EMISFLUX_VEG(NLEVV_ML,NTRACT),EMISFLUX_BS(NTRACT),
     &     EMISFLUX_SN(NTRACT),EMISFLUX_WAT(NTRACT),
     &     VDDEP(NLEVV_ML,NTRACT),FLXTDIFF(NLEVT,NTRAC),PMINT(NLON,NTRACT),
     &     PZ(NLON,NLEVT),PZHP(NLON,NLEVT+1),LAD_MLAY(NLEVV_ML+1),
     &     LADRAD_MLAY(NLEVV_ML+1),ISOP_MLAY(NLEVV_ML),MTERP_MLAY(NLEVV_ML),
     &     OVOC_MLAY(NLEVV_ML),U_MLAY(NLEVV_ML),FIELDIN(NLEVV,NNFIELD),
     &     FIELDOUT(NLEVV_ML,NNFIELD),TOTWEIGHT(NLEVV_ML),RATIO(NLEVV_ML),
     &     KH(NLEVV_ML),ZMVEG(NLON,NLEVV_ML,NTRACT),ZMVEG_OLD(NLEVV_ML),
     &     ZU(NLEVV_ML+1),EM(NLEVV_ML+1),VDD(NLEVV_ML+1,NTRACT),
     &     RECALC(NLEVT),ZKH(NLEVV_ML+1),RBML(NLEVV_ML,NTRACT),
     &     CCOMP(NTR),PCATYPE(NLON),RCOX_MLAY(NLEVV_ML,NTRACT),
     &     DTIME_MIXING(NLEVV_ML)
     
      REAL TERMA(NLEVV_ML+1),TERMB(NLEVV_ML+1),TERMC(NLEVV_ML+1),
     &     TERMD(NLEVV_ML+1),TERME(NLEVV_ML+1)
 
      REAL
     &     XTEEMIS(NLEVT,NTRACT),	! emission tendency  [molec cm-3 s-1]
     &     XTEDRYD(NLEVT,NTRACT),	! dry deposition tendency  [molec cm-3 s-1]
     &     XTEDIFF(NLEVT,NTRAC)		! vertical diffusion tendency  [molec cm-3 s-1]

      INTEGER JK,JJK,JKK,JJ,JT,JJT,JL,NSTEP,NSTART,NSTOP,NPRINT,NLEV,
     &        NLEVEL,IW,NVEGLAY(NLON),NLEV1,NLEV2,NLEVPR,NSTEP_SUB,IT,
     &        NNOXSPEC

      INTEGER ISUM(NNFIELD)   ! summing up the interpolated values for ISUM=1

      DATA G/9.8/

C LG- definition of minimal canopy height, the minimum time scale at 
C     which turbulent mixing occurs and the minimal lenght of the subtimestep
C     for coupled emissions/dry deposition and turbulence calculations,
C     and the number of species forming NOX

      PARAMETER (HCMIN=0.5,PTMST_MIN=10.,NNOXSPEC=5)

      CALL RESETI(ISUM,0)       ! default setting to zero

C LG- calling of subroutine in which the compensation point (stomatal comp.
C     points) are being defined. This routine is now (10-2001) being called 
C     here in this routine but the call should be moved to a location 
C     before the "big leaf" dry deposition calculations are performed to also 
C     include the role of the compensation point in those calculations

      CALL COMP(NLEV,PRHOA,PCATYPE,CCOMP)

C LG- 

      DTIME_MIXING(1)=MAX(2.*PTMST,1.e6)	! minimal two times the timestep
      DTIME_MIXING(NLEVV_ML)=MAX(2.*PTMST,1.e6)

      IF (NSTEP.EQ.1) THEN
         WRITE(*,'(1a)')
     &      'VEG_MLAY: the applied time for nocturnal mixing events is:' 
         WRITE(*,'(1a,e10.4,1a)')
     &      'DTIME_MIXING, crown-layer: ',DTIME_MIXING(1),' sec.' 
         WRITE(*,'(1a,e10.4,1a)')
     &      'DTIME_MIXING, soil-canopy layer: ',DTIME_MIXING(NLEVV_ML),' sec.'
	 WRITE(*,'(1a)')'ENTER TO CONTINUE'
	 READ (*,*)
      ENDIF 

      WRITE(NUNMDFL,*)'Start VEG_MLAY.f'

C LG- determining the actual canopy height for the constant pressure levels
C     and the reference height of the canopy which is taken to be the 
C     displacement height + the surface roughness for heat

      CANHEIGHT=0.
      DO 408 JL=1,NLON
       DO JK=1,NLEVV_ML
	CANHEIGHT=CANHEIGHT+PDP(JL,NLEV+JK)/(PRHOA(JL,NLEV+JK)*G*1.E3)
        RECALC(NLEV+JK)=(ZMAIR*1.E9/AVO)*(1./PRHOA(1,NLEV+JK))
       ENDDO
  408 CONTINUE

      RECALC(NLEV)=(ZMAIR*1.E9/AVO)*(1./PRHOA(1,NLEV))
      UAVG=0

C LG- scaling of some of the high vertical resolution parameters to the 
C     resolution of this routine, two layers

      DO JK=1,NLEVV

C LG-  the array counter for the LAD profile and isoprene emission fluxes
C      increases with height whereas the loop of the levels is from the top
C      to the bottom layer

       JJK=NLEVV+1-JK
       FIELDIN(JK,1)=LAD(JJK)
       ISUM(1)=1
     
C LG-  recalculation of emission flux to differents units for
C      comparison and writing of data to screen/file, 1.E3 it to 
C      recalculate from kg to g, 1/XMC to recalculate to mol C and the term 
C      1/5 is to correct for the 5 C molecules. In order to recalc from mol 
C      isoprene cm-2 s-1 to molecules its multiplied with the avogadro number

C LG-  for the biogenic emission fluxes the layer average flux has to be
C      determined and not the summated since later on the fluxes are 
C      scaled up to the total source flux by multiplying it with the LAD!
C      ISUM(2/3/4)=0 !!!

       FIELDIN(JK,2)=ISOPEM(JJK)*1.E3*(1./12.)*(1./5.)*AVO

C LG-  added the monoterpene emissions, C10

       FIELDIN(JK,3)=MTERPEM(JJK)*1.E3*(1./12.)*(1./10.)*AVO

C LG-  added the sesquiterpene emissions, C15, these are mimiced using the 
C      OVOC's emissions fluxes

       FIELDIN(JK,4)=OVOCEM(JJK)*1.E3*(1./12.)*(1./15.)*AVO

C LG-  determing the canopy average windspeed with is used to scale the
C      the within-canopy average windspeed

       FIELDIN(JK,5)=U(JJK)

       IF (JK.LT.NLEVV) UAVG=UAVG+U(JK)

C LG-  determing the stomatal resistances, the interpolation is done for 
C      ozone, which is default included in the dry deposition scheme and then
C      later done for the other trace gases, using the O3 interpolated field 

       IF (LAGS) THEN
        FIELDIN(JK,6)=RCOX_AGSML(JJK,io3)
       ELSE
        FIELDIN(JK,6)=RCOX_ECHML(JJK,io3)
       ENDIF
      ENDDO

      IF (NLEVV.EQ.1) THEN
       UAVG=U(JK)
      ELSE
       UAVG=UAVG/(NLEVV-1)
      ENDIF

      CALL INTERPOLATE(CANHEIGHT,NLEV,NLEVV_ML,PZHP,NNFIELD,
     &     ISUM,FIELDIN,FIELDOUT,TOTWEIGHT)

      LADSUM=0.
      DO JK=1,NLEVV_ML
       LAD_MLAY(JK)=FIELDOUT(JK,1)

C LG-  to warn in case that the interpolation didn't work properly

       LADSUM=LADSUM+LAD_MLAY(JK)
       IF (JK.EQ.NLEVV_ML.AND.LADSUM.LT.0.99) THEN
         WRITE(*,'(1A)')'veg_mlay: problem with interpolation'
         WRITE(*,'(1A)')'The summated LAD < 0.99: STOP CALLED in veg_mlay'
         STOP
       ENDIF

       ISOP_MLAY(JK)=FIELDOUT(JK,2)
       MTERP_MLAY(JK)=FIELDOUT(JK,3)
       OVOC_MLAY(JK)=FIELDOUT(JK,4)
       U_MLAY(JK)=FIELDOUT(JK,5)
       RCOX_MLAY(JK,io3)=FIELDOUT(JK,6)

C LG-  determining the trace gas stomatal conductances by correcting
C      for differences in diffusivities.

       DO JT=1,NTR
        IF (DIFF(JT).GT.0.AND.JT.NE.irad) THEN
         RCOX_MLAY(JK,JT)=RCOX_MLAY(JK,io3)*(DIFF(JT)/DIFF(io3))
	ENDIF
       ENDDO
      ENDDO

C LG- resetting the turbulent tracer flux

      CALL RESETR (XTEDIFF,NLEVT*NTRAC,0.)

C LG- for areas with an canopy height <= HCMIN , the concentrations
C     of the canopy are assigned the concentrations of the surface
C     layer, for RERUN=1 the concentrations within the canopy are not
C     reset

      DO 101 JT=1,NTRAC
      DO 101 JL=1,NLON
      DO 101 JK=1,NLEVV_ML
        IF (HC.LE.HCMIN.OR.NSTEP.EQ.0.AND.RERUN.EQ.0) 
     &     PM(JL,NLEV+JK,JT)=PM(JL,NLEV,JT)
 101  CONTINUE

C LG- and the short-lived species

      DO 102 JT=1,NG3X
      DO 102 JL=1,NLON
      DO 102 JK=1,NLEVV_ML
        IF (HC.LE.HCMIN.OR.NSTEP.EQ.0.AND.RERUN.EQ.0) 
     &     PMLOC(JL,NLEV+JK,JT)=PMLOC(JL,NLEV,JT)
 102  CONTINUE

C LG- assigning canopy and non-vegetated surface emission fluxes

      DO JT=1,NTRAC
       DO JK=1,NLEVV_ML
        EMISFLUX_VEG(JK,JT)=0.
       ENDDO
       EMISFLUX_BS(JT)=0.
       EMISFLUX_SN(JT)=0.
       EMISFLUX_WAT(JT)=0.
      ENDDO

C LG- assigning the emission fluxes 
       
      IF (LEMIS) THEN

C LG-  emission of NO, for long timesteps it can be expected that there
C      is a subtimestep scale conversion of the NO to NO2 and that is
C      accounted for in the update of the NO/NO2 concentrations by
C      considering the reaction between NO and O3

       EMISFLUX_VEG(NLEVV_ML,ino_pr)=NOMOLEC
       EMISFLUX_BS(ino_pr)=NOMOLEC

C LG-  since the isoprene emission flux is already a grid average
C      emission flux since all the controlling parameters are grid 
C      averages, e.g. dry matter, for this routine, the 
C      grid average isoprene emission flux is scaled with the
C      relative fraction of vegetation and wet skin fraction, in 
C      order to get an overall emission flux which resembles that
C      calculated in the subroutine VOCEMIS.f 

       DO JL=1,NLON
        IF ((VEGFRAC(JL)+WSFRAC(JL)).GT.0.) THEN
         DO JK=1,NLEVV_ML
          EMISFLUX_VEG(JK,iisop)=ISOP_MLAY(JK)/(VEGFRAC(JL)+WSFRAC(JL))
          EMISFLUX_VEG(JK,imaterp)=FEMMATERP*MTERP_MLAY(JK)/
     &      (VEGFRAC(JL)+WSFRAC(JL))
          EMISFLUX_VEG(JK,imbterp)=FEMMBTERP*MTERP_MLAY(JK)/
     &      (VEGFRAC(JL)+WSFRAC(JL))
          EMISFLUX_VEG(JK,isqterp)=FEMOVOC*OVOC_MLAY(JK)/
     &      (VEGFRAC(JL)+WSFRAC(JL))
         ENDDO
        ENDIF
       ENDDO

C LG-  radon soil emissions (to check turbulent transport)
 
       EMISFLUX_VEG(NLEVV_ML,irad)=FEMRAD*RADEM

C LG-  CO2 respiration/emissions 
 
       EMISFLUX_VEG(NLEVV_ML,ico2)=FEMCO2*CO2MOLEC

C LG-  NH3 soil emission flux 
 
       EMISFLUX_VEG(NLEVV_ML,inh3)=FEMNH3*NH3MOLEC

C LG-  CH3CL canopy/soil emission flux

       EMISFLUX_VEG(NLEVV_ML,ich3cl)=CH3CLMOLEC

C LG-  CHCL3 canopy/soil emission flux

       EMISFLUX_VEG(NLEVV_ML,ichcl3)=CHCL3MOLEC

      ENDIF

C LG- calculation of the dry deposition velocities

      IF (LDRYDEP) THEN

       DO 104 JT=1,NTRACT
       DO 104 JL=1,NLON
       DO 104 JK=1,NLEVV_ML
        JJK=JK+NLEV
        QV=1.E5
        IF (RLEAF(JL,JT).EQ.0.) GOTO 104
 
C LG-   calculation of the total uptake resistance of the vegetated and
C       wet skin fraction from the specific land cover type resistances
C       and the land cover fractions, this gives the total uptake
C       velocities for these fractions, reflecting solely the within
C       canopy uptake processes. The turbulent transport is considered
C       in the calculations of the canopy top fluxes. The bulk dry deposition 
C       velocities are scaled with the LAD and radiation profile for the 
C       two layers in order to scale the "big leaf" dry deposition velocity
C       accounting for the biomass distribution and the larger removal rates
C       for the sunlit leaves

C LG-   alternative calculation of quasi-laminar boundary layer resistance
C       for canopy layers (see also ddv.f)

        ! mz_lg_20031015+ modified. The old equation used the parameter
        !     diffrb(jt) whereas the value of 180 is already reflecting the
	!     difussion properties for SO2 
        RBML(JK,JT)=DIFFRB(JT)/DIFFRB(iso2)*180.*(0.07/U_MLAY(JK))**0.5

        IF (JK.EQ.NLEVV_ML) THEN

C LG-    first order estimate of the effect of the existence of a 
C        compensation point for NO2, assuming a value for this of 0.5 ppbv

C          IF (JT.EQ.ino2_pr.AND.CCOMP(JT).GT.0.AND.NSTEP.EQ.0) THEN
C           OPEN(UNIT=2,FILE='/data/ganzevl/racmo/output/RleafNO2_Ccomp.out',
C      *    STATUS='UNKNOWN')
C           WRITE(2,'(2a)')'Rleaf NO2 as function of C, CCOMP(JT)=1 ppbv'
C           WRITE(2,'(5a15)')'CNO2','Rleaf','RstomNO2','dR','fac'
C 	  DO JJ=1,200
C            RB(JK,JT)=25.
C 	   RCOX_MLAY(JK,JT)=100.
C            RLEAF(JL,JT)=RCOX_MLAY(JK,JT)
C            IF (JJ*0.01.LT.CCOMP(JT)) THEN
C             RLEAF(JL,JT)=-RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
C      &        (-1.-(JJ*0.01-CCOMP(JT))/CCOMP(JT))**7 
C 	   ELSEIF (JJ*0.01.GT.CCOMP(JT).AND.JJ*0.01.LT.2.*CCOMP(JT)) THEN
C             RLEAF(JL,JT)=RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
C      &        (1.-(JJ*0.01-CCOMP(JT))/CCOMP(JT))**7 
C            ENDIF
C            WRITE(2,'(5E15.4)')JJ*0.01,RLEAF(JL,JT),RCOX_MLAY(JK,ino2_pr),
C      &       (1.E5-RCOX_MLAY(JK,JT)),(1.-(JJ*0.01-CCOMP(JT))/CCOMP(JT))**7
C           ENDDO
C           CLOSE(2)
C          ENDIF

C LG-    considering the potential role of the existence of a compensation
C        point for gases such as NO2, for RLEAF < 0 there is an emission
C        instead of dry deposition, for most of the gases except of CO2, 
C        some expontential increase in the uptake resistance is being
C        assumed, expressed by the term 1.E5-RSTOMX(JL,JT))*
C           &  (-1.-(PM(1,NLEV+JK,JT)-CCOMP(JT))/CCOMP(JT))**7

         IF (CCOMP(JT).GT.0.) THEN
	 
	  IF(JK.EQ.NLEVV_ML.AND.NSTEP.EQ.0.AND.CCOMP(JT).GT.0.) THEN
	    WRITE(*,'(1a,i4,1a,f10.2,1a)')
     &       ' A comp. point for tracer index: ',JT,' has been set at:',
     &         CCOMP(JT)/(PRHOA(1,NLEV)/(ZMAIR*1.E9/AVO)),' ppbv'
	    WRITE(*,'(1a)')'ENTER to continue'
	    READ(*,*)
	  ENDIF
      	 	 
          IF (JT.LE.NTRAC) THEN
           IF (PM(1,NLEV+JK,JT).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (-1.-(PM(1,NLEV+JK,JT)-CCOMP(JT))/CCOMP(JT))**7 
	   ELSEIF (PM(1,NLEV+JK,JT).GT.CCOMP(JT).AND.
     &      PM(1,NLEV+JK,JT).LT.2.*CCOMP(JT)) THEN
            RLEAF(JL,JT)=RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (1.-(PM(1,NLEV+JK,JT)-CCOMP(JT))/CCOMP(JT))**7 
           ENDIF

C LG-      short lived species
	    
	  ELSE
	   IF (PMLOC(1,NLEV+JK,JT-NTRAC).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (-1.-(PMLOC(1,NLEV+JK,JT-NTRAC)-CCOMP(JT))/CCOMP(JT))**7 
	   ELSEIF (PMLOC(1,NLEV+JK,JT-NTRAC).GT.CCOMP(JT).AND.
     &      PMLOC(1,NLEV+JK,JT-NTRAC).LT.2.*CCOMP(JT)) THEN
            RLEAF(JL,JT)=RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (1.-(PMLOC(1,NLEV+JK,JT-NTRAC)-CCOMP(JT))/CCOMP(JT))**7 
           ENDIF
          ENDIF
	 ENDIF

C LG-    19092003, overwritting the leaf resistance of HNO3

         IF (JT.EQ.ihno3) THEN
	   RCOX_MLAY(JK,JT)=1.
	   RMES(JT)=1.
	   RCUT(JT)=1.
	 ENDIF

c         IF (JT.EQ.ih2o2) THEN 
c	   RCOX_MLAY(JK,JT)=1.
c	   RMES(JT)=1.
c	   RCUT(JT)=1.
c	   print *,'veg_mlay, be carefull: RleafH2O2 set to 1 s m-1'
c	 ENDIF

C LG-    Using the multi-layer stomatal resistance, which is calculated
C        considering explicitly the extinction of radiation within the canopy.
C        Otherwise the "bulk" leaf resistance is being used

         IF (LRCO_ML) THEN
          IF (JT.LE.NTRAC) THEN
           IF (PM(1,NLEV+JK,JT).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
	   ELSEIF (PM(1,NLEV+JK,JT).GT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=
     &       (1./((AMAX1(1E-5,LAD_MLAY(JK)*LAI)/(RBML(JK,JT)+RCUT(JT)))+
     &            (AMAX1(1E-5,LAD_MLAY(JK)*LAI)/(RBML(JK,JT)+RCOX_MLAY(JK,JT)+RMES(JT)))))
           ENDIF
	  ELSE
           IF (PMLOC(1,NLEV+JK,JT-NTRAC).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
	   ELSEIF (PMLOC(1,NLEV+JK,JT-NTRAC).GT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
           ENDIF
          ENDIF

         ENDIF

C LG-    19092003, overwritting the leaf resistance of HNO3

         IF (JT.EQ.ihno3) RLEAF(JL,JT)=1.

c         IF (JT.EQ.ih2o2) THEN
c	   RLEAF(JL,JT)=1.
c	   print *,'veg_mlay, be carefull: RleafH2O2 set to 1 s m-1'
c	 ENDIF

C LG-    for the lowest canopy layer, the limiting turbulent transport from
C        the reference height to the soil surface is also considered.

         IF (RLEAF(JL,JT).GT.0.) THEN
          RSVEG=1./((1./(RSOIL(JT)+RAHCAN(JL)*PZ(JL,JK+NLEV)/CANHEIGHT))+
     &     (1./RLEAF(JL,JT)))
         ELSE
          RSVEG=1./((1./(RSOIL(JT)+RAHCAN(JL)*PZ(JL,JK+NLEV)/CANHEIGHT))+
     &     (1./((-RBML(JK,JT)+RLEAF(JL,JT))/AMAX1(1E-5,LAD_MLAY(JK)*LAI))))
         ENDIF

         RSWET=1./((1./(RWS(JT)+RAHCAN(JL)*PZ(JL,JK+NLEV)/CANHEIGHT))+
     &     (1./((RBML(JK,JT)+RWS(JT))/AMAX1(1E-5,LAD_MLAY(JK)*LAI))))

        ELSE

C LG-    considering the potential role of the existence of a compensation
C        point for gases such as NO2, for RLEAF < 0 there is an emission
C        instead of dry deposition 

         IF (CCOMP(JT).GT.0.) THEN
          IF (JT.LE.NTRAC) THEN
           IF (PM(1,NLEV+JK,JT).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (-1.-(PM(1,NLEV+JK,JT)-CCOMP(JT))/CCOMP(JT))**7 
	   ELSEIF (PM(1,NLEV+JK,JT).GT.CCOMP(JT).AND.
     &      PM(1,NLEV+JK,JT).LT.2.*CCOMP(JT)) THEN
            RLEAF(JL,JT)=RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (1.-(PM(1,NLEV+JK,JT)-CCOMP(JT))/CCOMP(JT))**7 
           ENDIF

C LG-      and the short lived species

	  ELSE
	   IF (PMLOC(1,NLEV+JK,JT-NTRAC).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (-1.-(PMLOC(1,NLEV+JK,JT-NTRAC)-CCOMP(JT))/CCOMP(JT))**7 
	   ELSEIF (PMLOC(1,NLEV+JK,JT-NTRAC).GT.CCOMP(JT).AND.
     &      PMLOC(1,NLEV+JK,JT-NTRAC).LT.2.*CCOMP(JT)) THEN
            RLEAF(JL,JT)=RCOX_MLAY(JK,JT)+(1.E5-RCOX_MLAY(JK,JT))*
     &       (1.-(PMLOC(1,NLEV+JK,JT-NTRAC)-CCOMP(JT))/CCOMP(JT))**7 
           ENDIF
          ENDIF
	 ENDIF

C LG-    19092003, overwritting the leaf resistance of HNO3

         IF (JT.EQ.ihno3) THEN
	   RCOX_MLAY(JK,JT)=1.
	   RMES(JT)=1.
	   RCUT(JT)=1.
	 ENDIF
 
C          IF (JT.EQ.ih2o2) THEN 
C 	   RCOX_MLAY(JK,JT)=1.
C 	   RMES(JT)=1.
C 	   RCUT(JT)=1.
C 	   print *,'veg_mlay, be carefull: RleafH2O2 set to 1 s m-1'
C 	 ENDIF

C LG-    Using the multi-layer stomatal resistance, which is calculated
C        considering explicitly the extinction of radiation within the canopy.
C        Otherwise the "bulk" leaf resistance is being used

         IF (LRCO_ML) THEN
          IF (JT.LE.NTRAC) THEN
           IF (PM(1,NLEV+JK,JT).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
	   ELSEIF (PM(1,NLEV+JK,JT).GT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=
     &       (1./((AMAX1(1E-5,LAD_MLAY(JK)*LAI)/(RBML(JK,JT)+RCUT(JT)))+
     &            (AMAX1(1E-5,LAD_MLAY(JK)*LAI)/(RBML(JK,JT)+RCOX_MLAY(JK,JT)+RMES(JT)))))

c            RLEAF(JL,JT)=(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
           ENDIF
	  ELSE
           IF (PMLOC(1,NLEV+JK,JT-NTRAC).LT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=-(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
	   ELSEIF (PMLOC(1,NLEV+JK,JT-NTRAC).GT.CCOMP(JT)) THEN
            RLEAF(JL,JT)=(1./((1./RCUT(JT))+(1./(RCOX_MLAY(JK,JT)+RMES(JT)))))
           ENDIF
          ENDIF
         ENDIF
         
C LG-    19092003, overwritting the leaf resistance of HNO3

         IF (JT.EQ.ihno3) RLEAF(JL,JT)=1.

c         IF (JT.EQ.ih2o2) RLEAF(JL,JT)=1.

         IF (RLEAF(JL,JT).GT.0.) THEN
          RSVEG=RLEAF(JL,JT)

c          RSVEG=(RBML(JK,JT)+RLEAF(JL,JT))/
c     &     AMAX1(1E-5,LAD_MLAY(JK)*LAI)

          RSWET=(RBML(JK,JT)+RWS(JT))/
     &     AMAX1(1E-5,LAD_MLAY(JK)*LAI)
         ELSE
          RSVEG=(-RBML(JK,JT)+RLEAF(JL,JT))/
     &     AMAX1(1E-5,LAD_MLAY(JK)*LAI)
          RSWET=(-RBML(JK,JT)+RWS(JT))/
     &     AMAX1(1E-5,LAD_MLAY(JK)*LAI)
         ENDIF

        ENDIF

        RTOT=1./
     &   ((VEGFRAC(JL)/(VEGFRAC(JL)+WSFRAC(JL)))*(1./(RSVEG))+
     &    (WSFRAC(JL)/(VEGFRAC(JL)+WSFRAC(JL)))*(1./(RSWET)))
        VDDEP(JK,JT)=MIN(10.,100./RTOT)	!max Vd of 10 cm s-1

C LG-   calculation of emission flux from leaves as a function of the
C       selected compensation point and ambient concentration, and the
C       emission velocity (-VDDEP), and the total area of stomata to
C       account for the fact that it is assumed that the emissions are 
C       through the stomata, the SAI is the stomatal area index 
C       (see Yienger and Levy, 1995), being assigned in the subroutine
C       NOXEMIS.f to the 12 discerned ecosystems for the biogenic NOx
C       emission flux calculations. The SAI multiplied with the LAI in 
C       each layer yields the total leaf emissions flux.

        IF (VDDEP(JK,JT).LT.0.) THEN
         EMISFLUX_VEG(JK,JT)=
     &     1.E4*(SAI*LAD_MLAY(JK)*LAI)*		! stomatal area index * LAI
     &     (-VDDEP(JK,JT)*(CCOMP(JT)-PM(1,NLEV+JK,JT))) !molecules m-2 s-1
         VDDEP(JK,JT)=1.E-20
        ENDIF
 104   CONTINUE

      ENDIF

C LG- the emission calculations are only done if one of the surface
C     cover types has been assigned an emission factor. Otherwise the
C     emissions are calculated within the routine emiss.f 

c LG- start-------prepare emission tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 520 JT=1,NTRAC                                     
         DO JK=1,NLEVV_ML
	   JJK=JK+NLEV
	   XTEEMIS(JJK,JT)=PM(1,NLEV+JK,JT)/PRHOA(1,JJK)
         ENDDO
	 XTEEMIS(NLEV,JT)=PM(1,NLEV,JT)/PRHOA(1,NLEV)
  520   CONTINUE
      ENDIF

C LG- calculation of the emissions within the bulk vegetation layer

      IF (LEMIS) THEN

       DO 103 JT=1,NTRAC
       DO 103 JL=1,NLON
       DO 103 JK=1,NLEVV_ML
        JJK=JK+NLEV

C LG-   no concentration update when IEM_TURB=1, this parameter
C       is default set to one for all the long lived tracers within
C       the subroutine EMISS.f. This makes the calculations of the
C       vertical diffusion being done within this subroutine by the
C       numerical procedure in which the emission process and the
C       turbulence are resolved in a coupled equation 

        IF (IEM_TURB(JT).EQ.1) GOTO 103

c EMIFAC in m2 s cm-3

        EMIFACSL=1.E-6*PTMST*G*1.E3*PRHOA(JL,NLEV)/PDP(JL,NLEV)

        EMIFACVEG(JK)=0.
        IF (HC.GT.HCMIN.AND.(VEGFRAC(JL)+WSFRAC(JL)).GT.0.)
     &   EMIFACVEG(JK)=1.E-6*PTMST/
     &      (PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))

        PM(JL,NLEV+JK,JT)=PM(JL,NLEV+JK,JT)+
     &           EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT) !EMISFLUX in molec. m-2 s-1 

C LG-   for emissions within the canopy, IW (BXT(IW,IBEMIS..) is 1
 
        BXT(1,IBEMIS,JT)=BXT(1,IBEMIS,JT)+
     &    EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT)*(VEGFRAC(JL)+WSFRAC(JL))*
     &    GRVOL(JL,JJK)

C LG-   and emissions in the surface layer

        IF (IEMDD_CHEM(JT).EQ.0.AND.IEM_TURB(JT).EQ.0) THEN
         PM(JL,NLEV,JT)=PM(JL,NLEV,JT)+
     &      SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT)+ 
     &      BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT)+ 
     &      WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT) 
  
C LG-    for emissions within the surface layer, IW (BXT(IW,IBEMIS..) is 2

         BXT(2,IBEMIS,JT)=BXT(2,IBEMIS,JT)+
     &     (SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT)+ 
     &      BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT)+ 
     &      WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT))*
     &      GRVOL(JL,NLEV)

        ENDIF

	FLUXCTOP(JL,JT)=0.
        SURFFLUX(JL,JT)=SURFFLUX(JL,JT)+
     &                  SNFRAC(JL)*EMISFLUX_SN(JT)+
     &                  BSFRAC(JL)*EMISFLUX_BS(JT)+
     &                  WTFRAC(JL)*EMISFLUX_WAT(JT)

 103   CONTINUE

C LG-  emissions of short lived species

       DO 1003 JT=1,NG3X
       DO 1003 JL=1,NLON
       DO 1003 JK=1,NLEVV_ML
        JJK=JK+NLEV

c EMIFAC in m2 s cm-3

        EMIFACSL=1.E-6*PTMST*G*1.E3*PRHOA(JL,NLEV)/PDP(JL,NLEV)

        EMIFACVEG(JK)=0.
        IF (HC.GT.HCMIN.AND.(VEGFRAC(JL)+WSFRAC(JL)).GT.0.)
     &   EMIFACVEG(JK)=1.E-6*PTMST/
     &      (PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))

C LG-   for long timesteps, not transporting the separate species of NOX
C       but solely NOX, the NO emission will already be partly chemically
C       converted to NO2, the chemical conversion is estimated using the
C       ratio of NO/NO2 of the previous timestep and the available amount
C       of O3 where it is assumed that the conversion of NO to NO2 is 
C       basically controlled by the reaction between NO and NO2

C LG-   the parameter RATIONOX is the ratio of NO and NO2, determined
C       in the chemistry subroutine from JNO2 and O3 concentration and
C       the reaction rate of the reaction between NO and O3, or from
C       the ratio of the old concentrations being update in this routine
C       to estimate the amount of NO that is being emitted in the form of 
C       NO2 at the subgridscale of the bulk vegetation model. Since it has 
C       been assumed that the NO is being transformed to NO2 by reaction 
C       with O3 this concentrations are also considered

        RATIONOX=PMLOC(JL,JJK,ino)/MAX(1.,PMLOC(JL,JJK,ino2))

C LG-   the statement MAX(1.-PM(JL,NLEV+1,io3)/(EMIFACVEG... prevents
C       the calculation of a titration of O3 larger than the old 
C       concentration (the term (PM(JL,JJK,io3)-1.) prevents the 
C       calculation of a zero concentration of O3)
 
        RATIO(JK)=1.
        
        IF (JT.EQ.ino.AND.EMISFLUX_VEG(JK,JT+NTRAC).GT.0.) THEN
         IF (LAIRCHEM.OR.LCBM4_ECH.OR.LCBM4_TM3.OR.LCHEM_MIM)
     &    RATIO(JK)=MAX(1.-(PM(JL,JJK,io3)-1.)/(EMIFACVEG(JK)*
     &     EMISFLUX_VEG(JK,JT+NTRAC)),RATIONOX/(RATIONOX+1.))

C LG-    subgrid scale chemical destruction NO, i.e., effective NO
C        emission flux

         PMLOC(JL,JJK,JT)=PMLOC(JL,JJK,JT)+
     &    RATIO(JK)*EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT+NTRAC)  

C LG-    subgrid scale chemical production NO2, i.e, NO2 emission flux

         PMLOC(JL,JJK,ino2)=PMLOC(JL,JJK,ino2)+
     &     (1.-RATIO(JK))*EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT+NTRAC)  

C LG-    updating NOX only for IEM_TURB EQ 0, otherwise the emission flux
C        of NO is assigned to NOX which is consequently being used within
C        the coupled emission/dry deposition/turbulence calculations
   
         IF (IEM_TURB(inox).EQ.0) THEN
          PMNOX_OLD=PM(JL,JJK,inox)
          PM(JL,JJK,inox)=PMLOC(JL,JJK,ino)+PMLOC(JL,JJK,ino2)+
     &     PMLOC(JL,JJK,ino3)+2*PMLOC(JL,JJK,in2o5)+PMLOC(JL,JJK,ihno4) 
          BXT(2,IBEMIS,inox)=BXT(2,IBEMIS,inox)+
     &     (PM(JL,JJK,inox)-PMNOX_OLD)*(VEGFRAC(JL)+WSFRAC(JL))*
     &      GRVOL(JL,JJK)
         ELSE
          EMISFLUX_VEG(JK,inox)=0.
          DO JJT=1,NNOXSPEC
 	   EMISFLUX_VEG(JK,inox)=EMISFLUX_VEG(JK,inox)+
     &      EMISFLUX_VEG(JK,NTRAC+JJT)
          ENDDO
	 ENDIF

C LG-    subgrid scale chemical destruction O3

         PM(JL,JJK,io3)=PM(JL,JJK,io3)-
     &     (1.-RATIO(JK))*EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT+NTRAC)  

        ELSE

         PMLOC(JL,NLEV+JK,JT)=PMLOC(JL,NLEV+JK,JT)+
     &          EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT+NTRAC) !EMISFLUX in molec. m-2 s-1 

         BG3(1,IBEMIS,JT)=BG3(1,IBEMIS,JT)+
     &     EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT+NTRAC)*(VEGFRAC(JL)+WSFRAC(JL))*
     &     GRVOL(JL,JJK)

        ENDIF

C LG-   and emissions in the surface layer

        PMLOC(JL,NLEV,JT)=PMLOC(JL,NLEV,JT)+
     &       SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT+NTRAC)+ 
     &       BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT+NTRAC)+ 
     &       WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT+NTRAC) 
  
C LG-   for emissions within the surface layer, IW (BXT(IW,IBEMIS..) is 2

        BG3(2,IBEMIS,JT)=BG3(2,IBEMIS,JT)+
     &      (SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT+NTRAC)+ 
     &       BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT+NTRAC)+ 
     &       WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT+NTRAC))*
     &       GRVOL(JL,NLEV)

 1003  CONTINUE

      ENDIF

c LG- end---------calculating emission tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN
       DO 522 JT=1,NTRAC   
        DO JK=1,NLEVV_ML
	 JJK=JK+NLEV                                     
         XTEEMIS(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &                    XTEEMIS(JJK,JT))/DTIME       
        ENDDO
        XTEEMIS(NLEV,JT)=(PM(1,NLEV,JT)/PRHOA(1,NLEV)-
     &                    XTEEMIS(NLEV,JT))/DTIME        
  522  CONTINUE
      ENDIF

c LG- end-      -------- emission  --------------------

c LG- start-------prepare dry deposition tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 524 JT=1,NTRAC  
	DO 524 JK=1,NLEVV_ML
	  JJK=JK+NLEV                                 
	  XTEDRYD(JJK,JT)=PM(1,NLEV+JK,JT)/PRHOA(1,JJK)
  524   CONTINUE
      ENDIF

C LG- calculation of the dry deposition tendency using the updated concentration
C     estimate

      IF (LDRYDEP) THEN

C LG-  updating the concentrations (only for the transported species)

       DO 105 JT=1,NTRAC
       DO 105 JL=1,NLON
       DO 105 JK=1,NLEVV_ML
        JJK=JK+NLEV
	IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
         IF (JT.EQ.inox) VDDEP(JK,JT)=1.E-20
	 IF (JT.EQ.ino3) VDDEP(JK,JT)=VDDEP(JK,ihno3)
	 IF (JT.EQ.in2o5) VDDEP(JK,JT)=VDDEP(JK,ihno3)
	ELSE

C LG-    when for the NOx deposition velocity the NO3 and N2O5 are
C        also considered then for a fair comparison with simulations
C        considering the seperate species forming NOx, the dry deposition
C        velocities of NO3 and N2O5 also need to calculated (default using
C        the HNO3 deposition velocity)

         IF (JT.EQ.inox) THEN
	  VDDEP(JK,JT)=   
     *    (VDDEP(JK,ino_pr)*PMLOC(JL,JJK,ino)+
     *     VDDEP(JK,ino2_pr)*PMLOC(JL,JJK,ino2)+
     *     VDDEP(JK,ihno3)*PMLOC(JL,JJK,ino3)+
     *     VDDEP(JK,ihno3)*2*PMLOC(JL,JJK,in2o5))
     *     /(PMLOC(JL,JJK,ino)+PMLOC(JL,JJK,ino2)+PMLOC(JL,JJK,ino3)+
     *       2*PMLOC(JL,JJK,in2o5)+PMLOC(JL,JJK,ihno4))
         ENDIF
        ENDIF

C LG-   no concentration update when IDD_TURB=1, this parameter
C       is default set to one for all the long lived tracers within
C       the subroutine DRYDEP.f. This makes the calculations of the
C       vertical diffusion being done within this subroutine by the
C       numerical procedure in which the dry deposition process and the
C       turbulence are resolved in a coupled equation 

        IF (VDDEP(JK,JT).EQ.0.OR.IDD_TURB(JT).EQ.1) GOTO 105

        QV=100./VDDEP(JK,JT)
        VX=-1./(QV*(PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))/PTMST)
        VD=100.*(PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))/PTMST*
     &     (1.-EXP(VX))
        HEIGHTDD=VD*PTMST/100.
        DDFRAC=HEIGHTDD/(PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))

        BCONC0=PM(JL,NLEV+JK,JT)
        PM(JL,NLEV+JK,JT)=(1.-DDFRAC)*PM(JL,NLEV+JK,JT)
        BXT(1,IBDDEP,JT)=BXT(1,IBDDEP,JT)+(PM(JL,NLEV+JK,JT)-BCONC0)*
     &                  (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK)
 105   CONTINUE

      ENDIF

c LG- end---------calculating dry deposition tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN
       DO 526 JT=1,NTRAC 
       DO 526 JK=1,NLEVV_ML 
         JJK=JK+NLEV                                       
         XTEDRYD(JJK,JT)=(PM(1,JK+NLEV,JT)/PRHOA(1,JJK)-
     &                    XTEDRYD(JJK,JT))/DTIME         
  526  CONTINUE
      ENDIF

c LG- end-      -------- dry deposition  --------------------

c LG- start-------prepare vertical diffusion tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 528 JT=1,NTRAC 
          DO JK=1,NLEVV_ML
	   JJK=JK+NLEV
           XTEDIFF(JJK,JT)=PM(1,NLEV+JK,JT)/PRHOA(1,JJK)
	  ENDDO
          XTEDIFF(NLEV,JT)=PM(1,NLEV,JT)/PRHOA(1,NLEV)
  528   CONTINUE
      ENDIF

      IF (LXTVDIFF) THEN

C LG-  ===numerical solver for vert. diff., dry depos. and emiss========
C      calculation of concentrations as a function of vertical transport
C      dry deposition and emissions within one equation, so without the 
C      operator splitting as it has been used in the default version. 

       DO 116 JT=1,NTRAC
       DO 116 JL=1,NLON

C LG-   no vertical diffusion of NOx for NTRAC=NTRACT since all the 
C       species are being transported separately

        IF (NTRAC.GT.ino.AND.ino.GT.1.AND.JT.EQ.inox) GOTO 116

C LG-   setting the soil surface flux to zero

        FLUXVEG(JL,NLEVV_ML+1,JT)=0.

        IF (HC.GT.HCMIN.AND.
     &     (VEGFRAC(JL)+WSFRAC(JL)).GT.0.) THEN

C LG-   calculating the eddy diffusivity from the Dz between the 
C       reference height of the surface layer and the reference height
C       of the canopy and the aerodynamic resistance RAHVEG for the canopy
C       top layer, this is already done here since the Kh can be used to
C       determine the Kh for the other canopy layers using the windspeed 
C       profile

        DZ=(PDP(JL,NLEV)/(PRHOA(JL,NLEV)*G*1.E3))/2.+
     &    (CANHEIGHT-PZ(JL,NLEV+1))
        KHMIN=(DZ**2)/DTIME_MIXING(1)	
	KH(1)=MAX(KHMIN,DZ/RAHVEG(JL))

        KH(1)=MAX(KHMIN,KH(1))*MIN(1.,NSTEP/(T_SPINUP/PTMST))
 
C LG-   starting with the lower canopy layer, going to the canopy top
C       up to the layer below the top-canopy layer

        RAH_TOT=RAHVEG(JL)
        SUMU=0.

        DO JK=NLEVV_ML,2,-1

C LG-    calculating the eddy diffusivity from the Dz between the 
C        reference height of the first and the second canopy layer. The 
C        term RAHCAN is calculated within the routine EC4_VDIFF. There 
C        are two approaches which both calculate this resistance for
C        the canopy height. Therefore the term CANHEIGHT/DZ has been 
C        introduced to correct for the difference between the vertical 
C        extent that has been represented
       
         DZ=PZ(JL,NLEV+(JK-1))-PZ(JL,NLEV+JK)
	
C LG-    definition of minimal value of Kh to study the sensitivity
C        of model resolved concentrations for small nocturnal exchange
C        fluxes. The minimal value is based on the canopy height and 
C        the number of nocturnal events during which the whole canopy is 
C        being ventilated (see paper Fitzjarrald and Moore, JGR 95,
C        16839-16850, 1990). They observed about 3 to 4 exchange events
C        implying one event each two hours [Kh]=m^2 s-1, for a canopy
C        height of 30 m and 2 hours * 3600 s=7200 s this gives a minimal
C        Kh of (30*30)/7200 ~ 0.1 m^2 s-1	
	
         KHMIN=(DZ**2)/DTIME_MIXING(JK)
         RAHCAN(JL)=MIN(RAHCAN(JL),(1./KHMIN)*CANHEIGHT)
         IF (NLEVV_ML.EQ.2) THEN

C LG-     112003, old formuluation based on using rahcan.

	  KH(JK)=MAX(0.,DZ/((MIN(RAHCAN(JL),(1./KHMIN)*CANHEIGHT))/
     &     (NLEVV_ML-1)))

C LG-     112003+ simply scaling of the surface-crown layer eddy-diffusivity 
C         with windspeed profile to estimate within canopy K profile

          KH(JK)=KH(1)*((U_MLAY(JK)+U_MLAY(1))/2.)/ 
     &                 ((U30(JL)+U_MLAY(1))/2.)

	 ELSE	
     	  KH(JK)=(U_MLAY(JK)/UAVG)*MAX(0.,DZ/
     &     ((0.5*MIN(RAHCAN(JL),(1./KHMIN)*CANHEIGHT))/(NLEVV_ML-1)))	
	 ENDIF

         KH(JK)=MAX(KHMIN,KH(JK))*MIN(1.,NSTEP/(T_SPINUP/PTMST))

        ENDDO

C LG-   definition of lenght of subtimestep used for coupled turbulence
C       dry deposition and emission calculations. The number of required 
C       subtimesteps is a function of the timescale of the dry 
C       deposition/emission/turbulent transport process calculated from 
C       the thickness of the canopy layers and the dry deposition velocity/
C       emission flux/eddy diffusivity 

C LG-   dry deposition timescale

        TS_DD=(PDP(JL,NLEV+1)/(PRHOA(JL,NLEV+1)*G*1.E3))/
     &    MAX(1.E-10,1.E-2*VDDEP(1,JT))

C LG-   emission timescale (it still needs to be adressed if you can really
C       use the concept of an emission timescale

        TS_EM=(1.E6*PM(JL,NLEV+2,JT)*(PDP(JL,NLEV+2)/
     &   (PRHOA(JL,NLEV+2)*G*1.E3)))/
     &    MAX(1.E-20,MAX(EMISFLUX_VEG(2,JT),EMISFLUX_VEG(1,JT)))

C LG-   turbulence timescale

        TS_TURB=((PDP(JL,NLEV+1)/(PRHOA(JL,NLEV+1)*G*1.E3))**2)/
     &    MAX(1.E-10,KH(1))

C LG-   the applied subtimestep is determined from the dry deposition 
C       timescale, the constant value defined in the first MIN statement 
C       determines the minimum timestep of the subtimestep scale calculations

        IF (IEM_TURB(JT).EQ.1.OR.IDD_TURB(JT).EQ.1) THEN
         PTMST_SUB=MAX(MIN(PTMST,PTMST_MIN),MIN(PTMST,0.1*MIN(TS_DD,TS_TURB))) 
	ELSE
         PTMST_SUB=PTMST
	ENDIF

        NSTEP_SUB=MAX(1,INT(PTMST/PTMST_SUB))
        PTMST_SUB=PTMST/NSTEP_SUB
 
C LG-   end

        DO JK=NLEVV_ML+1,1,-1
	 JJK=JK-1

C LG-    for IDD_TURB=1 dry deposition is considered.

         IF (IDD_TURB(JT).EQ.1.AND.JJK.GT.0) THEN
          VDD(JJK,JT)=VDDEP(JJK,JT)*MIN(1.,NSTEP/(T_SPINUP/PTMST))
	 ELSE
	  VDD(JJK,JT)=0.
         ENDIF

C LG-    for IEM_TURB=1 emission is considered.

         IF (IEM_TURB(JT).EQ.1) THEN

c EMIFAC in m2 s cm-3

          EMIFACSL=1.E-6*PTMST_SUB*G*1.E3*PRHOA(JL,NLEV)/PDP(JL,NLEV)

          EMIFACVEG(JJK)=0.
          IF (HC.GT.HCMIN.AND.(VEGFRAC(JL)+WSFRAC(JL)).GT.0.)
     &     EMIFACVEG(JJK)=1.E-6*PTMST_SUB/
     &      (PDP(JL,JJK+NLEV)/(PRHOA(JL,JJK+NLEV)*G*1.E3))*
     &       MIN(1.,NSTEP/(T_SPINUP/PTMST)) ! smoothing startup 

         ENDIF

C LG-    The updated concentrations are calculated such that the 
C        decreasing gradient occuring at a "sub" timestep scale 
C        has been considered in order to prevent the scheme to calculate
C        that large fluxes that negative concentrations are being
C        calculated. The same procedure is also applied to the 
C        calculation of decrease of the concentrations due to dry 
C        deposition. The subsript 1 refers to the top layer and 2 to
C        the lower layer (see also calculations for the top canopy layer)

         IF (JK.EQ.NLEVV_ML+1) THEN
          DZ1=(PDP(JL,NLEV+(JJK-1))/(PRHOA(JL,NLEV+(JJK-1))*G*1.E3))
          DZ2=(PDP(JL,NLEV+JJK)/(PRHOA(JL,NLEV+JJK)*G*1.E3))
          DZ3=0.
          TERMA(JK)=KH(JJK)*PTMST_SUB/(PZ(JL,NLEV+(JJK-1))*DZ2-
     &     PZ(JL,NLEV+JJK)*DZ2)
          TERMB(JK)=0.
          TERMC(JK)=1.+1.E-2*(VDD(JJK,JT)/DZ2)*PTMST_SUB+
     &      TERMA(JK)+TERMB(JK)
          EM(JK)=EMIFACVEG(JJK)*EMISFLUX_VEG(JJK,JT) !EMISFLUX in molec. m-2 s-1 
          ZMVEG_OLD(JJK)=PM(JL,NLEV+JJK,JT)
         ELSEIF(JK.GT.1) THEN
          DZ1=(PDP(JL,NLEV+(JJK-1))/(PRHOA(JL,NLEV+(JJK-1))*G*1.E3))
          DZ2=(PDP(JL,NLEV+JJK)/(PRHOA(JL,NLEV+JJK)*G*1.E3))
          DZ3=(PDP(JL,NLEV+JJK+1)/(PRHOA(JL,NLEV+JJK+1)*G*1.E3))
          TERMA(JK)=KH(JJK)*PTMST_SUB/(PZ(JL,NLEV+(JJK-1))*DZ2-
     &     PZ(JL,NLEV+JJK)*DZ2)
          TERMB(JK)=KH(JJK+1)*PTMST_SUB/(PZ(JL,NLEV+JJK)*DZ2-
     &     PZ(JL,NLEV+JJK+1)*DZ2)
          TERMC(JK)=1.+1.E-2*(VDD(JJK,JT)/DZ2)*PTMST_SUB+
     &      TERMA(JK)+TERMB(JK)
          EM(JK)=EMIFACVEG(JJK)*EMISFLUX_VEG(JJK,JT) !EMISFLUX in molec. m-2 s-1 
          ZMVEG_OLD(JJK)=PM(JL,NLEV+JJK,JT)
         ELSEIF(JK.EQ.1) THEN
          DZ1=0.
          DZ2=(PDP(JL,NLEV+JJK)/(PRHOA(JL,NLEV+JJK)*G*1.E3))
          DZ3=(PDP(JL,NLEV+JJK+1)/(PRHOA(JL,NLEV+JJK+1)*G*1.E3))
          TERMA(JK)=0.
          TERMB(JK)=KH(JJK+1)*PTMST_SUB/(PZ(JL,NLEV+JJK)*DZ2-
     &     PZ(JL,NLEV+JJK+1)*DZ2)
          TERMC(JK)=1.+TERMA(JK)+TERMB(JK)
          EM(JK)=0.
          PM_OLD=PM(JL,NLEV,JT)
         ENDIF
	 TERMD(JK)=TERMA(JK)/TERMC(JK)
	 TERME(JK)=TERMB(JK)/TERMC(JK)
        ENDDO

        DMASS=0.
	MASS=0.
        DMASS_EM=0.

C LG-   Start loop for performing emission, dry deposition, and 
C       turbulence calculations multiple times to deal with subtimestep
C       scale process timescales such as HNO3 and O3 dry deposition

        DO IT=1,NSTEP_SUB

C LG-    top canopy layer

         JK=2
         JJK=1

C LG-    budget calculations for dry deposition

         DDFRAC=(VDD(JJK,JT)*PTMST_SUB/100.)/
     &     (PDP(JL,NLEV+JJK)/(PRHOA(JL,NLEV+JJK)*G*1.E3))

         BXT(1,IBDDEP,JT)=BXT(1,IBDDEP,JT)-PM(JL,NLEV+JJK,JT)*DDFRAC*
     &               (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK+NLEV)

C LG-    for emissions within the canopy and the surface layer, 
C        IW (BXT(IW,IBEMIS..) is 1
 
         BXT(1,IBEMIS,JT)=BXT(1,IBEMIS,JT)+EM(JK)*
     &     (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK+NLEV)
         DMASS_EM=DMASS_EM+EM(JK)*(VEGFRAC(JL)+WSFRAC(JL))*
     &    GRVOL(JL,JJK+NLEV)

C LG-    end

         PM(JL,NLEV+JJK,JT)=((PM(JL,NLEV+JJK,JT)+EM(JK))/TERMC(JK)+
     &     (TERMD(JK)*PM(JL,NLEV+(JJK-1),JT))/TERMC(JK-1)+
     &     (TERME(JK)*PM(JL,NLEV+JJK+1,JT))/TERMC(JK+1))/
     &     (1.-TERMD(JK)*TERME(JK-1)-TERME(JK)*TERMD(JK+1))

C LG-    calculation of dry deposition and emission tendency which is used 
C        with the total tendency to arrive at the vertical diffusion tendency
                                     
         IF (VDD(JJK,JT).GT.0.AND.IDD_TURB(JT).EQ.1)
     &    XTEDRYD(NLEV+JJK,JT)=-1.E-2*VDD(JJK,JT)*(PM(JL,NLEV+JJK,JT)/
     &    PRHOA(JL,NLEV+JJK))/(PDP(JL,NLEV+JJK)/
     &   (PRHOA(JL,NLEV+JJK)*G*1.E3)) 
        
         IF (EM(JK).GT.0.AND.IEM_TURB(JT).EQ.1)
     &    XTEEMIS(NLEV+JJK,JT)=(EM(JK)/PRHOA(JL,NLEV+JJK))/PTMST_SUB  

C LG-    lower canopy layer

         JK=3
         JJK=2

C LG-    budget calculations for dry deposition

         DDFRAC=(VDD(JJK,JT)*PTMST_SUB/100.)/
     &     (PDP(JL,NLEV+JJK)/(PRHOA(JL,NLEV+JJK)*G*1.E3))
         BXT(1,IBDDEP,JT)=BXT(1,IBDDEP,JT)-PM(JL,NLEV+JJK,JT)*DDFRAC*
     &               (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK+NLEV)

C LG-    for emissions within the canopy and the surface layer, 
C        IW (BXT(IW,IBEMIS..) is 1
 
         BXT(1,IBEMIS,JT)=BXT(1,IBEMIS,JT)+EM(JK)*
     &     (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK+NLEV)
         DMASS_EM=DMASS_EM+EM(JK)*(VEGFRAC(JL)+WSFRAC(JL))*
     &     GRVOL(JL,JJK+NLEV)

C LG-    end

         PM(JL,NLEV+JJK,JT)=(PM(JL,NLEV+JJK,JT)+EM(JK))/TERMC(JK)+
     &    TERMD(JK)*PM(JL,NLEV+(JJK-1),JT)

C LG-    calculation of dry deposition and emission tendency which is used 
C        with the total tendency to arrive at the vertical diffusion tendency
                                     
         IF (VDD(JJK,JT).GT.0.AND.IDD_TURB(JT).EQ.1)
     &    XTEDRYD(NLEV+JJK,JT)=-1.E-2*VDD(JJK,JT)*(PM(JL,NLEV+JJK,JT)/
     &     PRHOA(JL,NLEV+JJK))/(PDP(JL,NLEV+JJK)/
     &    (PRHOA(JL,NLEV+JJK)*G*1.E3))  

         IF (EM(JK).GT.0.AND.IEM_TURB(JT).EQ.1)
     &    XTEEMIS(NLEV+JJK,JT)=(EM(JK)/PRHOA(JL,NLEV+JJK))/PTMST_SUB  

C LG-    surface layer

         JK=1
         JJK=0

C LG-    the surface layer concentrations are only updated if the canopy top
C        fluxes are not used within the vertical diffusion routine as the
C        lower boundary condition (LEM/DD_TURB=.FALSE.), else the parameter
C        SURFFLUX is used within EC4_VDIFF.f

         IF (IEMDD_CHEM(JT).EQ.0) THEN

          PM(JL,NLEV,JT)=
     &     (1.-(VEGFRAC(JL)+WSFRAC(JL)))*PM(JL,NLEV,JT)+
     &     (VEGFRAC(JL)+WSFRAC(JL))*((PM(JL,NLEV,JT)+EM(JK))/TERMC(JK)+
     &      TERME(JK)*PM(JL,NLEV+JJK+1,JT))

         ENDIF

C LG-    checking for mass conservation, this is only done whenever
C        the terms VDDEP and EM are set to zero within the vertical loop
C        in which the different terms of the numerical solver are
C        calculated! 

         IF (VDD(1,JT).EQ.0.AND.EM(3).EQ.0.AND.
     &    ABS(DMASS)/MAX(1.E-20,MASS).GT.1.E-2) THEN
          PRINT *,'VEG_MLAY.f (EM/DDTURB), no mass cons. for tracer no.: ', 
     &        JT,' and dmass: ',DMASS
          PRINT *,'Percentage mass loss (relat. to canopy mass): ',
     &        100.*(DMASS)/MAX(1.E-20,MASS)
         ENDIF

C LG-   end subtimestep loop DO IT=1,NSTEP_SUB
    
        ENDDO

C LG-   budget calculations for dry deposition, this budget is 
C       calculated as a residual term from the changes in the 
C       concentrations and the contribution by dry emissions and 
C       vertical diffusion since a straightforward calculation as 
C       in drydep.f is not possible

        MASS_OLD=PM_OLD*GRVOL(JL,NLEV)+
     &      ZMVEG_OLD(1)*(VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,NLEV+1)+
     &      ZMVEG_OLD(2)*(VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,NLEV+2)
        MASS_NEW=PM(JL,NLEV,JT)*GRVOL(JL,NLEV)+
     &      PM(JL,NLEV+1,JT)*(VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,NLEV+1)+
     &      PM(JL,NLEV+2,JT)*(VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,NLEV+2)
        DMASS=(MASS_NEW-MASS_OLD)-DMASS_EM

        DC(JT)=RECALC(NLEV)*PM(JL,NLEV,JT)-
     &     RECALC(NLEV+1)*PM(JL,NLEV+1,JT)

C LG-   calculation of the fluxes at the top of the canopy
C       layer from the dC/dt, the deposition/emission flux and the 
C       flux between the layer and the underlying layer
         
        DO JK=NLEVV_ML,1,-1
         FLUXVEG(JL,JK,JT)=	!molecules cm-2 s-1
     &    -((1.E2*(PDP(JL,NLEV+JK)/(PRHOA(JL,NLEV+JK)*G*1.E3))*
     &     (PM(JL,NLEV+JK,JT)-ZMVEG_OLD(JK)))/PTMST+
     &      VDD(JK,JT)*PM(JL,NLEV+JK,JT)-
     &      1.E-4*EM(JK+1)/MAX(1.E-20,EMIFACVEG(JK))- 
     &      FLUXVEG(JL,JK+1,JT))
        ENDDO

        IF (JT.EQ.irad) THEN
	 FLUXCTOP(JL,JT)=FLUXVEG(JL,1,JT)
     
C LG-    calculation of total surface flux considering all the land cover 
C        fractions

         SURFFLUX(JL,JT)=SURFFLUX(JL,JT)+
     &    (VEGFRAC(JL)+WSFRAC(JL))*FLUXCTOP(JL,JT)
        ELSE
         FLUXCTOP(JL,JT)=1.E-11*FLUXVEG(JL,1,JT)
     
C LG-    calculation of total surface flux considering all the land cover 
C        fractions, watch out! surfflux is in 1.e-11 molecules cm-2 s-1
C        and thus use of SURFFLUX in ec4_vdiff.f for coupled atmosphere-
C        biosphere exchange and vertical diffusion calculations must be 
C        checked for use proper units! 

         SURFFLUX(JL,JT)=1.E-11*(SURFFLUX(JL,JT)+
     &    (VEGFRAC(JL)+WSFRAC(JL))*1.E11*FLUXCTOP(JL,JT))

        ENDIF

C LG-   end IF (HC.GT.HCMIN.AND.(VEGFRAC(

        ENDIF

 116   CONTINUE

C LG-  calculating the NOx fluxes, concentrations and budget only
C      when NTRAC=NTRACT

       IF (NTRAC.GT.ino.AND.ino.GT.1) THEN

        DO 108 JL=1,NLON
       
         PM(JL,NLEV,inox)=PM(JL,NLEV,ino)+PM(JL,NLEV,ino2)+
     &    PM(JL,NLEV,ino3)+2*PM(JL,NLEV,in2o5)+PM(JL,NLEV,ihno4)

         DO JK=1,NLEVV_ML
          PM(JL,NLEV+JK,inox)=PM(JL,NLEV+JK,ino)+PM(JL,NLEV+JK,ino2)+
     &     PM(JL,NLEV+JK,ino3)+2*PM(JL,NLEV+JK,in2o5)+PM(JL,NLEV+JK,ihno4)

C LG-     determining the NOx dry deposition and emission tendency from
C         the tendencies of the species forming NOx

          XTEDRYD(NLEV+JK,inox)=
     &     XTEDRYD(NLEV+JK,ino)+XTEDRYD(NLEV+JK,ino2)+
     &     XTEDRYD(NLEV+JK,ino3)+XTEDRYD(NLEV+JK,in2o5)+
     &     XTEDRYD(NLEV+JK,ihno4)
          XTEEMIS(NLEV+JK,inox)=
     &     XTEEMIS(NLEV+JK,ino)+XTEEMIS(NLEV+JK,ino2)+
     &     XTEEMIS(NLEV+JK,ino3)+XTEEMIS(NLEV+JK,in2o5)+
     &     XTEEMIS(NLEV+JK,ihno4)

         ENDDO

         DC(inox)=RECALC(NLEV)*PM(JL,NLEV,inox)-
     &     RECALC(NLEV+1)*PM(JL,NLEV+1,inox)

         FLUXCTOP(JL,inox)=FLUXCTOP(JL,ino)+FLUXCTOP(JL,ino2)+
     &                    FLUXCTOP(JL,ino3)+FLUXCTOP(JL,in2o5)+
     &                    FLUXCTOP(JL,ihno4)

         SURFFLUX(JL,inox)=SURFFLUX(JL,ino)+SURFFLUX(JL,ino2)+
     &                    SURFFLUX(JL,ino3)+SURFFLUX(JL,in2o5)+
     &                    SURFFLUX(JL,ihno4)

 108    CONTINUE
 
      ELSE
        
C LG-   determining diagnostically the NO and NO2 flux from the
C       the concentration difference and eddy-diffusivity to interpret
C       the NOx flux in terms of the seperate species

        DO 1008 JL=1,NLON
         FLUXCTOP(JL,ino_pr)=1.E-11*
     &    (PMLOC(JL,NLEV+1,ino)-PMLOC(JL,NLEV,ino))*1.E4*KH(1)/
     &    (1.E2*(PZ(JL,NLEV)-PZ(JL,NLEV+1)))
         FLUXCTOP(JL,ino2_pr)=1.E-11*
     &    (PMLOC(JL,NLEV+1,ino2)-PMLOC(JL,NLEV,ino2))*1.E4*KH(1)/
     &    (1.E2*(PZ(JL,NLEV)-PZ(JL,NLEV+1)))
 1008   CONTINUE

      ENDIF

C LG- writing of Kh and wind profile to check profile
      
      IF (NSTEP.EQ.0) THEN
        OPEN(UNIT=NUNKHVEG,FILE='/data/ganzevl/racmo/output/Kh_veg_mlay.out',
     *    STATUS='UNKNOWN')
        WRITE(NUNKHVEG,'(2a)')'Eddy diffusivity for heat, windspeed'
        WRITE(NUNKHVEG,'(a5,2a12)')
     *   'level','KH [m2 s-1]','U [m s-1]'

        NLEV1=NLEV
        NLEV2=NLEV+NLEVV_ML
        NLEVPR=NLEV2-NLEV1+1

        WRITE (NUNKHVEG,*) NSTOP,NLEVPR,NPRINT,PTMST
        WRITE (NUNKHVEG,*) NLEV1,NLEV2
        WRITE (NUNKHVEG,*) (PZ(1,JK),JK=NLEV2,NLEV1,-1)
      ENDIF
 
      IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN

        WRITE(NUNKHVEG,*)NSTEP,JDAY,GMT,LTIME,HC,RAHVEG(1),RAHCAN(1),
     *       RAH_TOT

        DO JK=NLEV2,NLEV1,-1
         IF (JK.GT.NLEV) THEN
	  ZU(JK-NLEV+1)=U_MLAY(JK-NLEV)
          ZKH(JK-NLEV+1)=KH(JK-NLEV)
	 ELSE
	  ZU(JK-NLEV+1)=U30(1)
          ZKH(JK-NLEV+1)=KH30(1)
	 ENDIF

         WRITE(NUNKHVEG,'(I4.4,1X,3E12.3)')
     *     JK,ZKH(JK-NLEV+1),ZU(JK-NLEV+1),ZU(JK-NLEV+1)/U30(1)
        ENDDO

      ENDIF

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNKHVEG)

C LG- end writing Kh and windspeed profiles

C LG- end IF (LXTVDIFF)

      ENDIF

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 530 JT=1,NTRAC 
	  DO JK=NLEVV_ML,1,-1
	   JJK=JK+NLEV

C LG-      for considering the dry deposition flux coupled with the
C          vertical transport calculations, the total tendency is 
C          corrected with the dry deposition tendency

           IF (IDD_TURB(JT).EQ.1.AND.IEM_TURB(JT).EQ.1.) THEN
	    XTEDIFF(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &        XTEDIFF(JJK,JT))/DTIME-XTEDRYD(JJK,JT)-XTEEMIS(JJK,JT)	
           ELSEIF (IDD_TURB(JT).EQ.1.AND.IEM_TURB(JT).EQ.0) THEN   
            XTEDIFF(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &         XTEDIFF(JJK,JT))/DTIME-XTEDRYD(JJK,JT)
           ELSEIF (IEM_TURB(JT).EQ.1.AND.IDD_TURB(JT).EQ.0) THEN   
            XTEDIFF(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &         XTEDIFF(JJK,JT))/DTIME-XTEEMIS(JJK,JT)
           ELSE  
            XTEDIFF(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &         XTEDIFF(JJK,JT))/DTIME
           ENDIF
          ENDDO
          XTEDIFF(NLEV,JT)=(PM(1,NLEV,JT)/PRHOA(1,NLEV)-
     &	       XTEDIFF(NLEV,JT))/DTIME
  530   CONTINUE
      ENDIF

C LG- start calculation additional data and writing data to output
C     file

      DO 109 JL=1,NLON
       CRF(ino_pr)=0.
       CRF(ino2_pr)=0.
       CRF(inox)=0.
       IF (EMISFLUX_VEG(2,ino_pr).GT.0.OR.
     &  EMISFLUX_VEG(2,inox).GT.0.) THEN
        IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
         CRF(ino_pr)=1.E11*1.E4*FLUXCTOP(JL,ino_pr)/EMISFLUX_VEG(2,ino_pr)
         CRF(ino2_pr)=1.E11*1.E4*FLUXCTOP(JL,ino2_pr)/EMISFLUX_VEG(2,ino_pr)
	ELSE
         CRF(inox)=1.E11*1.E4*FLUXCTOP(JL,inox)/EMISFLUX_VEG(2,ino_pr)
	ENDIF
       ENDIF
       CRF(iisop)=0.
       IF (EMISFLUX_VEG(1,iisop).GT.0.)
     &   CRF(iisop)=1.E11*1.E4*FLUXCTOP(JL,iisop)/
     &     (EMISFLUX_VEG(1,iisop)+EMISFLUX_VEG(2,iisop))

C LG-  CO2

       CRF(ico2)=0.
       IF (EMISFLUX_VEG(NLEVV_ML,ico2).GT.0.)
     &   CRF(ico2)=1.E11*1.E4*FLUXCTOP(JL,ico2)/
     &     (EMISFLUX_VEG(1,ico2)+EMISFLUX_VEG(2,ico2))

C LG-  NH3

       CRF(inh3)=0.
       IF (EMISFLUX_VEG(NLEVV_ML,inh3).GT.0.)
     &   CRF(inh3)=1.E11*1.E4*FLUXCTOP(JL,inh3)/
     &     (EMISFLUX_VEG(1,inh3)+EMISFLUX_VEG(2,inh3))

109   CONTINUE

C LG- writing of canopy reduction factor and other parameters to output file

      IF (NSTEP.EQ.0) THEN
       OPEN(UNIT=NUNMVLAY,FILE='/data/ganzevl/racmo/output/veg_mlay.out',
     *    STATUS='UNKNOWN')
       WRITE(NUNMVLAY,'(1a)')
     *   'Conc., fluxes and canopy reduction factors for 2 layer veg. mode'
       WRITE(NUNMVLAY,'(1a,f8.0)')
     *   'Period for nocturnal mixing DTIME_MIXING (in canopy): ',
     *    DTIME_MIXING(NLEVV_ML)
       WRITE(NUNMVLAY,'(1a,f8.0)')
     *   'Minimum lenght subtimestep PTMST_MIN: ',PTMST_MIN
       WRITE(NUNMVLAY,'(1a)')
     *   'fluxes in 1e11 molecules cm-2 s-1, CO2 in 1e15 molecules cm-2 s-1'
       IF (LAIRCHEM) THEN
        WRITE(NUNMVLAY,'(2a10,48a12)')'nstep','time',
     *    'u*[ms-1]','RAHVEG','RAHCAN',
     *    'RBVEG(O3)','RB(1,O3)/LAI','RB(2,O3)/LAI',
     *    'Kh(1)','Kh(2)','RG','LAD(1)','LAD(2)','VdHNO3','VdNO','VdNO2(1)',
     *    'VdNO2(2)','VdNOX(1)','VdNO3','VdO3(1)','VdO3(2)','VdC5H8','VdCO2',
     *    'NOx-veg','NOx-surf','NO2-veg','NO2-surf','O3-veg','O3-surf',
     *    'OH-veg','OH-surf','DC(NOx)','DC(O3)','FLUXVEG(O3)','FLUXVEG(NO)',
     *    'FLUXCTOP(O3)','FLUX(O3)','FLUXCTOP(NO)','FLUX(NO)',
     *    'FLUXCTOP(NO2)','FLUX(NO2)','FLUXCTOP(NOx)','FLUX(NOx)',
     *    'NOMOLEC','CRF(NO)','CRF(NO2)','CRF(NOx)','FLUXCTOP(rad)',
     *    'EMFLXNOX(1)','EMFLXNOX(2)'
       ELSE
        WRITE(NUNMVLAY,'(a10,a15,76a15)')'nstep','time',
     *    'u*[ms-1]','RAHVEG','RAHCAN',
     *    'RBVEG(O3)','RB(1,O3)/LAI','RB(2,O3)/LAI',
     *    'Kh(1)','Kh(2)','RG','LAD(1)','LAD(2)','VdHNO3(1)','VdHNO3(2)',
     *    'VdNO(1)','VdNO(2)','VdNO2(1)','VdNO2(2)',
     *    'VdNOX(1)','VdNOX(2)','VdNO3','VdO3(1)','VdO3(2)','VdC5H8',
     *    'VdCO2(1)','VdCO2(2)','VdNH3(1)','VdNH3(2)',
     *    'VdH2O2(1)','VdH2O2(2)','RCO_MLAY(1)','RCO_MLAY(2)',
     *    'NOx-veg','NOx-surf','NO2-veg','NO2-surf',
     *    'O3-veg','O3-surf','C5H8-veg','C5H8-surf',
     *    'OH-veg','OH-surf','DC(NOx)','DC(O3)','FLUXVEG(O3)','FLUXVEG(NO)',
     *    'FLUXCTOP(O3)','FLUX(O3)','FLUXCTOP(NO)','FLUX(NO)',
     *    'FLUXCTOP(NO2)','FLUX(NO2)','FLUXCTOP(NOx)','FLUX(NOx)',
     *    'NOMOLEC','CRF(NO)','CRF(NO2)','CRF(NOx)',
     *    'FLUXCTOP(isop)','FLUX(isop)','ISOPMOLEC','CRF(isop)',
     *    'FLUXCTOP(rad)','EMFLXNOX(1)','EMFLXNOX(2)', 
     *    'FLUXCTOP(CO2)','FLUX(CO2)','CO2MOLEC','CRF(CO2)',
     *    'FLUXCTOP(NH3)','FLUX(NH3)','NH3MOLEC','CRF(NH3)',
     *    'FLUXCTOP(H2O2)','FLUX(H2O2)','FLUXCTOP(MHP)','FLUX(MHP)'
       ENDIF
      ENDIF

      IF (MOD(NSTEP,NPRINT).EQ.0) THEN

       IF (NTRAC.GT.ino.AND.ino.GT.1) THEN
        IF (LAIRCHEM) THEN
         WRITE(NUNMVLAY,'(1x,i9.9,1x,a14,20f12.4,21e12.4,3f12.4,3e12.4)')
     *     NSTEP,LDATLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),RBVEG(1,io3),
     *     RBML(1,io3)/AMAX1(1E-5,LAD_MLAY(1)*LAI),
     *     RBML(2,io3)/AMAX1(1E-5,LAD_MLAY(2)*LAI),KH(1),KH(2),
     *     RG,LAD_MLAY(1),LAD_MLAY(2),VDDEP(1,ihno3),VDDEP(1,ino),
     *     VDDEP(1,ino2),VDDEP(2,ino2),VDDEP(1,inox),VDDEP(1,ino3),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,iisop),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,inox),RECALC(NLEV)*PM(1,NLEV,inox),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,ino2),RECALC(NLEV)*PM(1,NLEV,ino2),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,io3),RECALC(NLEV)*PM(1,NLEV,io3), 
     *     PM(1,NLEV+1,ioh),PM(1,NLEV,ioh),DC(inox),DC(io3),
     *     FLUXVEG(1,NLEVV_ML,io3),FLUXVEG(1,NLEVV_ML,ino),
     *     FLUXCTOP(1,io3),SURFFLUX(1,io3),FLUXCTOP(1,ino),SURFFLUX(1,ino),
     *     FLUXCTOP(1,ino2),SURFFLUX(1,ino2),FLUXCTOP(1,inox),SURFFLUX(1,inox),
     *     1.e-4*1.e-11*EMISFLUX_VEG(2,ino),CRF(ino_pr),CRF(ino2_pr),
     *     CRF(inox),FLUXCTOP(1,irad),1.e-4*1.e-11*EMISFLUX_VEG(1,ino2_pr),
     *     1.e-4*1.e-11*EMISFLUX_VEG(2,ino2_pr)
        ELSE
         WRITE(NUNMVLAY,
     *      '(1x,i9.9,1x,a14,31f15.4,23e15.4,3f15.4,3e15.4,f15.4,6e15.4,f15.4,
     *        3e15.4,f15.4,4e15.4)')
     *     NSTEP,LDATLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),RBVEG(1,io3),
     *     RBML(1,io3)/AMAX1(1E-5,LAD_MLAY(1)*LAI),
     *     RBML(2,io3)/AMAX1(1E-5,LAD_MLAY(2)*LAI),KH(1),KH(2),
     *     RG,LAD_MLAY(1),LAD_MLAY(2),VDDEP(1,ihno3),VDDEP(2,ihno3),
     *     VDDEP(1,ino),VDDEP(2,ino),VDDEP(1,ino2),VDDEP(2,ino2),
     *     VDDEP(1,inox),VDDEP(2,inox),VDDEP(1,ino3),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,iisop),
     *     VDDEP(1,ico2),VDDEP(2,ico2),VDDEP(1,inh3),VDDEP(2,inh3),
     *     VDDEP(1,ih2o2),VDDEP(2,ih2o2),
     *     RCOX_MLAY(1,io3)/DIFF(io3),RCOX_MLAY(2,io3)/DIFF(io3),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,inox),RECALC(NLEV)*PM(1,NLEV,inox),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,ino2),RECALC(NLEV)*PM(1,NLEV,ino2),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,io3),RECALC(NLEV)*PM(1,NLEV,io3), 
     *     RECALC(NLEV+1)*PM(1,NLEV+1,iisop),RECALC(NLEV)*PM(1,NLEV,iisop),
     *     PM(1,NLEV+1,ioh),PM(1,NLEV,ioh),DC(inox),DC(io3),
     *     FLUXVEG(1,NLEVV_ML,io3),FLUXVEG(1,NLEVV_ML,ino_pr),
     *     FLUXCTOP(1,io3),SURFFLUX(1,io3),FLUXCTOP(1,ino_pr),
     *     SURFFLUX(1,ino_pr),FLUXCTOP(1,ino2_pr),SURFFLUX(1,ino2_pr),
     *     FLUXCTOP(1,inox),SURFFLUX(1,inox),
     *     1.e-4*1.e-11*EMISFLUX_VEG(2,ino_pr),CRF(ino_pr),CRF(ino2_pr),CRF(inox),
     *     FLUXCTOP(1,iisop),SURFFLUX(1,iisop),1.e-4*1.e-11*ISOPMOLEC,
     *     CRF(iisop),FLUXCTOP(1,irad),1.e-4*1.e-11*EMISFLUX_VEG(1,ino2),
     *     1.e-4*1.e-11*EMISFLUX_VEG(2,ino2),1.e-4*FLUXCTOP(1,ico2),
     *     1.e-4*SURFFLUX(1,ico2),1.e-4*1.e-15*CO2MOLEC,CRF(ico2),
     *     FLUXCTOP(1,inh3),SURFFLUX(1,inh3),1.e-4*1.e-11*NH3MOLEC,CRF(inh3),
     *     FLUXCTOP(1,ih2o2),SURFFLUX(1,ih2o2),
     *     FLUXCTOP(1,ich3o2h),SURFFLUX(1,ich3o2h)
        ENDIF
       ELSE
	IF (LAIRCHEM) THEN
         WRITE(NUNMVLAY,'(1x,i9.9,1x,a14,20f12.4,21e12.4,3f12.4,3e12.4)')
     *     NSTEP,LDATLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),RBVEG(1,io3),
     *     RBML(1,io3)/AMAX1(1E-5,LAD_MLAY(1)*LAI),
     *     RBML(2,io3)/AMAX1(1E-5,LAD_MLAY(2)*LAI),KH(1),KH(2),
     *     RG,LAD_MLAY(1),LAD_MLAY(2),VDDEP(1,ihno3),VDDEP(1,ino_pr),
     *     VDDEP(1,ino2_pr),VDDEP(2,ino2_pr),VDDEP(1,inox),VDDEP(1,ino3+NTRAC),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,iisop),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,inox),RECALC(NLEV)*PM(1,NLEV,inox),
     *     RECALC(NLEV+1)*PMLOC(1,NLEV+1,ino2),RECALC(NLEV)*PMLOC(1,NLEV,ino2),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,io3),RECALC(NLEV)*PM(1,NLEV,io3), 
     *     PMLOC(1,NLEV+1,ioh),PMLOC(1,NLEV,ioh),DC(inox),DC(io3),
     *     FLUXVEG(1,NLEVV_ML,io3),FLUXVEG(1,NLEVV_ML,ino_pr),
     *     FLUXCTOP(1,io3),SURFFLUX(1,io3),FLUXCTOP(1,ino_pr),
     *     SURFFLUX(1,ino_pr),FLUXCTOP(1,ino2),SURFFLUX(1,ino2_pr),
     *     FLUXCTOP(1,inox),SURFFLUX(1,inox),
     *     1.e-4*1.e-11*EMISFLUX_VEG(2,ino_pr),CRF(ino_pr),
     *     CRF(ino2_pr),CRF(inox),FLUXCTOP(1,irad),
     *     1.e-4*1.e-11*EMISFLUX_VEG(1,inox),1.e-4*1.e-11*EMISFLUX_VEG(2,inox)
        ELSE
         WRITE(NUNMVLAY,
     *      '(1x,i9.9,1x,a14,31f15.4,23e15.4,3f15.4,3e15.4,f15.4,6e15.4,f15.4,
     *        3e15.4,f15.4,4e15.4)')
     *     NSTEP,LDATLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),RBVEG(1,io3),
     *     RBML(1,io3)/AMAX1(1E-5,LAD_MLAY(1)*LAI),
     *     RBML(2,io3)/AMAX1(1E-5,LAD_MLAY(2)*LAI),KH(1),KH(2),
     *     RG,LAD_MLAY(1),LAD_MLAY(2),VDDEP(1,ihno3),VDDEP(2,ihno3),
     *     VDDEP(1,ino_pr),VDDEP(2,ino_pr),VDDEP(1,ino2_pr),VDDEP(2,ino2_pr),
     *     VDDEP(1,inox),VDDEP(2,inox),VDDEP(1,ino3+NTRAC),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,iisop),
     *     VDDEP(1,ico2),VDDEP(2,ico2),VDDEP(1,inh3),VDDEP(2,inh3),
     *     VDDEP(1,ih2o2),VDDEP(2,ih2o2),
     *     RCOX_MLAY(1,io3)/DIFF(io3),RCOX_MLAY(2,io3)/DIFF(io3),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,inox),RECALC(NLEV)*PM(1,NLEV,inox),
     *     RECALC(NLEV+1)*PMLOC(1,NLEV+1,ino2),RECALC(NLEV)*PMLOC(1,NLEV,ino2),
     *     RECALC(NLEV+1)*PM(1,NLEV+1,io3),RECALC(NLEV)*PM(1,NLEV,io3), 
     *     RECALC(NLEV+1)*PM(1,NLEV+1,iisop),RECALC(NLEV)*PM(1,NLEV,iisop),
     *     PMLOC(1,NLEV+1,ioh),PMLOC(1,NLEV,ioh),DC(inox),DC(io3),
     *     FLUXVEG(1,NLEVV_ML,io3),FLUXVEG(1,NLEVV_ML,ino_pr),
     *     FLUXCTOP(1,io3),SURFFLUX(1,io3),FLUXCTOP(1,ino_pr),
     *     SURFFLUX(1,ino_pr),FLUXCTOP(1,ino2_pr),SURFFLUX(1,ino2_pr),
     *     FLUXCTOP(1,inox),SURFFLUX(1,inox),
     *     1.e-4*1.e-11*EMISFLUX_VEG(2,ino_pr),CRF(ino_pr),
     *     CRF(ino2_pr),CRF(inox),FLUXCTOP(1,iisop),SURFFLUX(1,iisop),
     *     1.e-4*1.e-11*ISOPMOLEC,CRF(iisop),FLUXCTOP(1,irad),
     *     1.e-4*1.e-11*EMISFLUX_VEG(1,inox),1.e-4*1.e-11*EMISFLUX_VEG(2,inox),
     *     1.e-4*FLUXCTOP(1,ico2),1.e-4*SURFFLUX(1,ico2),1.e-4*1.e-15*CO2MOLEC,
     *     CRF(ico2),
     *     FLUXCTOP(1,inh3),SURFFLUX(1,inh3),1.e-4*1.e-11*NH3MOLEC,CRF(inh3),
     *     FLUXCTOP(1,ih2o2),SURFFLUX(1,ih2o2),
     *     FLUXCTOP(1,ich3o2h),SURFFLUX(1,ich3o2h)
        ENDIF

       ENDIF

      ENDIF

      IF (NSTEP.EQ.NSTOP) CLOSE(NUNMVLAY)

C LG- 

      WRITE(NUNMDFL,*)'End VEG_MLAY.f'

      RETURN
      END


