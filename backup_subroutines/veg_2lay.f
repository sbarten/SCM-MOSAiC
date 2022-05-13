      SUBROUTINE VEG_2LAY (NSTEP,NSTART,NSTOP,NPRINT,NLEV,NLEVEL,
     &               PTMST,DTIME,PM,PMLOC,PRHOA,PDP,PZ,PZHP,
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

      REAL EMIFACSL,EMIFAC_UNDER,PTMST,DTIME,G,QV,VX,VD,
     &     HEIGHTDD,DDFRAC,DZ,RTOT,TERM1,TERM2,ALPHAPOS,ALPHANEG,
     &     CA,CB,CC,CD,RECALC,PMVEG_OLD,PM_OLD,DMASS_VEG,DMASS_SL,
     &     MASS_VEG,BCONC0,CANHEIGHT,ZCANOPY,HCMIN,RSVEG
     
      REAL PM(NLON,NLEVT,NTRAC),PMLOC(NLON,NLEVT,NG3X),
     &     PRHOA(NLON,NLEVT),PDP(NLON,NLEVT),
     &     EMIFACVEG(2),DC(NTRAC),FLUXCTOP(NLON,NTRAC),
     &     FLUXVEG(NLON,NTRAC),SURFFLUX(NLON,NTRAC),CRF(NTRAC),
     &     EMISFLUX_VEG(2,NTRACT),EMISFLUX_BS(NTRACT),
     &     EMISFLUX_SN(NTRACT),EMISFLUX_WAT(NTRACT),
     &     VDDEP(2,NTRACT),FLXTDIFF(NLEVT,NTRAC),PMINT(NLON,NTRACT),
     &     PZ(NLON,NLEVT),PZHP(NLON,NLEVT+1),LAD_2LAY(2),ISOP_2LAY(2),
     &     LADRAD_2LAY(2),FIELDIN(NLEVV,2),FIELDOUT(2,2),TOTWEIGHT(2),
     &     RATIO(2),KH(2),KHMIN,DTIME_MIXING
 
      REAL
     &     XTEEMIS(NLEVT,NTRACT),	! emission tendency  [molec cm-3 s-1]
     &     XTEDRYD(NLEVT,NTRACT),	! dry deposition tendency  [molec cm-3 s-1]
     &     XTEDIFF(NLEVT,NTRAC)		! vertical diffusion tendency  [molec cm-3 s-1]

      INTEGER JK,JJK,JKK,JT,JL,NSTEP,NSTART,NSTOP,NPRINT,NLEV,NLEVEL,IW,
     &        NVEGLAY(NLON)

      DATA G/9.8/
      PARAMETER (HCMIN=0.5,DTIME_MIXING=7200.)

C LG- 

      WRITE(NUNMDFL,*)'Start VEG_2LAY.f'

C LG- determining the actual canopy height for the constant pressure levels
C     and the reference height of the canopy which is taken to be the 
C     displacement height + the surface roughness for heat

      DO 408 JL=1,NLON
        CANHEIGHT=PDP(JL,NLEV+1)/(PRHOA(JL,NLEV+1)*G*1.E3)+
     &            PDP(JL,NLEV+2)/(PRHOA(JL,NLEV+2)*G*1.E3)
        NVEGLAY(JL)=2
        IF (CANHEIGHT.LT.15.) 
     &    NVEGLAY(JL)=1
        IF (CANHEIGHT.LT.0.5)
     &    NVEGLAY(JL)=0
  408 CONTINUE

C LG- scaling of some of the high vertical resolution parameters to the 
C     resolution of this routine, two layers

      DO JK=1,NLEVV

C LG-  the array counter for the LAD profile and isoprene emission fluxes
C      increases with height whereas the loop of the levels is from the top
C      to the bottom layer

       JJK=NLEVV+1-JK
       FIELDIN(JK,1)=LAD(JJK)

c --   recalculation of emission flux to differents units for
c      comparison and writing of data to screen/file, 1.E3 it to 
c      recalculate from kg to g, 1/XMC to recalculate to mol C and the term 
C      1/5 is to correct for the 5 C molecules. In order to recalc from mol 
C      isoprene cm-2 s-1 to molecules its multiplied with the avogadro number

       FIELDIN(JK,2)=ISOPEM(JJK)*1.E3*(1./12.)*(1./5.)*AVO
      ENDDO

      CALL INTERPOLATE(CANHEIGHT,NLEV,2,PZHP,2,
     &     FIELDIN,FIELDOUT,TOTWEIGHT)

      DO JK=1,NVEGLAY(1)
       LAD_2LAY(JK)=FIELDOUT(JK,1)
       ISOP_2LAY(JK)=FIELDOUT(JK,2)
      ENDDO

C LG- the parameter LADRAD is a scaling factor which reflects both the
C     distribution of biomass and the radiation within the canopy. As a
C     a first order of estimate the ratio of the isoprene emission fluxes
C     is being used since these are a function of the LAD and radiation

      IF (ISOP_2LAY(2).GT.0.) THEN	! daytime, net radiation > 0
       LADRAD_2LAY(1)=
     &  (ISOP_2LAY(1)/ISOP_2LAY(2))/(1.+ISOP_2LAY(1)/ISOP_2LAY(2))
       LADRAD_2LAY(2)=1.-LADRAD_2LAY(1)
      ELSE	! nocturnal conditions
       LADRAD_2LAY(1)=LAD_2LAY(1)
       LADRAD_2LAY(2)=LAD_2LAY(2)
      ENDIF

C LG- resetting the turbulent tracer flux

      CALL RESETR (XTEDIFF,NLEVT*NTRAC,0.)

C LG- for areas with an canopy height <= HCMIN , the concentrations
C     of the canopy are assigned the concentrations of the surface
C     layer

      DO 101 JT=1,NTRAC
      DO 101 JL=1,NLON
      DO 101 JK=1,NVEGLAY(JL)
        IF (HC.LE.HCMIN.OR.NSTEP.EQ.0) 
     &     PM(JL,NLEV+JK,JT)=PM(JL,NLEV,JT)
 101  CONTINUE

C LG- and the short-lived species

      DO 102 JT=1,NG3X
      DO 102 JL=1,NLON
      DO 102 JK=1,NVEGLAY(JL)
        IF (HC.LE.HCMIN.OR.NSTEP.EQ.0) 
     &     PMLOC(JL,NLEV+JK,JT)=PMLOC(JL,NLEV,JT)
 102  CONTINUE

C LG- assigning canopy and non-vegetated surface emission fluxes

      DO JT=1,NTRACT
       EMISFLUX_VEG(1,JT)=0.
       EMISFLUX_VEG(2,JT)=0.
       EMISFLUX_BS(JT)=0.
       EMISFLUX_SN(JT)=0.
       EMISFLUX_WAT(JT)=0.
      ENDDO

C LG- assigning the emission fluxes 

      IF (LEMIS) THEN

       EMISFLUX_VEG(2,inoemdd)=NOMOLEC
       EMISFLUX_BS(inoemdd)=NOMOLEC

C LG-  since the isoprene emission flux is already a grid average
C      emission flux since all the controlling parameters are grid 
C      averages, e.g. dry matter, for this routine, the 
C      grid average isoprene emission flux is scaled with the
C      relative fraction of vegetation and wet skin fraction, in 
C      order to get an overall emission flux which resembles that
C      calculated in the subroutine VOCEMIS.f 

       DO JL=1,NLON
        IF ((VEGFRAC(JL)+WSFRAC(JL)).GT.0.) THEN
         EMISFLUX_VEG(1,iisop)=ISOP_2LAY(1)/(VEGFRAC(JL)+WSFRAC(JL))
         EMISFLUX_VEG(2,iisop)=ISOP_2LAY(2)/(VEGFRAC(JL)+WSFRAC(JL))
        ENDIF
       ENDDO
      ENDIF

C LG- the emission calculations are only done if one of the surface
C     cover types has been assigned an emission factor. Otherwise the
C     emissions are calculated within the routine emiss.f 

c LG- start-------prepare emission tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 520 JT=1,NTRAC                                     
         DO JK=1,NVEGLAY(1)
	   JJK=JK+NLEV
	   XTEEMIS(JJK,JT)=PM(1,NLEV+JK,JT)/PRHOA(1,JJK)
         ENDDO
	 XTEEMIS(NLEV,JT)=PM(1,NLEV,JT)/PRHOA(1,NLEV)
  520   CONTINUE
      ENDIF

C LG- calculation of the emissions within the bulk vegetation layer

      IF (LEMIS) THEN

       DO 103 JT=1,NTRACT
       DO 103 JL=1,NLON
       DO 103 JK=1,NVEGLAY(JL)
        JJK=JK+NLEV

c EMIFAC in m2 s cm-3

        EMIFACSL=1.E-6*PTMST*G*1.E3*PRHOA(JL,NLEV)/PDP(JL,NLEV)

        EMIFACVEG(JK)=0.
        IF (HC.GT.HCMIN.AND.(VEGFRAC(JL)+WSFRAC(JL)).GT.0.)
     &   EMIFACVEG(JK)=1.E-6*PTMST/
     &      (PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))

         IF (JT.LE.NTRAC) THEN
          PM(JL,NLEV+JK,JT)=PM(JL,NLEV+JK,JT)+
     &                   EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT) !EMISFLUX in molec. m-2 s-1 
         ELSE
          PMLOC(JL,NLEV+JK,JT-NTRAC)=PMLOC(JL,NLEV+JK,JT-NTRAC)+
     &                   EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT) !EMISFLUX in molec. m-2 s-1 
	 ENDIF

C LG-   for emissions within the canopy, IW (BXT(IW,IBEMIS..) is 1
 
        IF (JT.LE.NTRAC) THEN
         BXT(1,IBEMIS,JT)=BXT(1,IBEMIS,JT)+
     &     EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT)*(VEGFRAC(JL)+WSFRAC(JL))*
     &     GRVOL(JL,JJK)
        ELSE
          BG3(1,IBEMIS,JT-NTRAC)=BG3(1,IBEMIS,JT-NTRAC)+
     &     EMIFACVEG(JK)*EMISFLUX_VEG(JK,JT)*(VEGFRAC(JL)+WSFRAC(JL))*
     &     GRVOL(JL,JJK)
        ENDIF

C LG-   and emissions in the surface layer

        IF (JT.NE.inox) THEN
         IF (IEMDD_CHEM(JT).EQ.0.AND.IEMDD_TURB(JT).EQ.0) THEN
          IF (JT.LE.NTRAC) THEN
           PM(JL,NLEV,JT)=PM(JL,NLEV,JT)+
     &       SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT)+ 
     &       BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT)+ 
     &       WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT) 
  
C LG-      for emissions within the surface layer, IW (BXT(IW,IBEMIS..) is 2

           BXT(2,IBEMIS,JT)=BXT(2,IBEMIS,JT)+
     &      (SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT)+ 
     &       BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT)+ 
     &       WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT))*
     &       GRVOL(JL,NLEV)
          ELSE
           PMLOC(JL,NLEV,JT-NTRAC)=PMLOC(JL,NLEV,JT-NTRAC)+
     &       SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT)+ 
     &       BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT)+ 
     &       WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT) 
  
C LG-      for emissions within the surface layer, IW (BXT(IW,IBEMIS..) is 2

           BG3(2,IBEMIS,JT-NTRAC)=BG3(2,IBEMIS,JT-NTRAC)+
     &      (SNFRAC(JL)*EMIFACSL*EMISFLUX_SN(JT)+ 
     &       BSFRAC(JL)*EMIFACSL*EMISFLUX_BS(JT)+ 
     &       WTFRAC(JL)*EMIFACSL*EMISFLUX_WAT(JT))*
     &       GRVOL(JL,NLEV)
          ENDIF
         ENDIF

	 FLUXCTOP(JL,JT)=0.
         SURFFLUX(JL,JT)=SNFRAC(JL)*EMISFLUX_SN(JT)+
     &                   BSFRAC(JL)*EMISFLUX_BS(JT)+
     &                   WTFRAC(JL)*EMISFLUX_WAT(JT)
	ENDIF

 103   CONTINUE
 
      ENDIF

c LG- end---------calculating emission tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN
       DO 522 JT=1,NTRAC   
        DO JK=1,NVEGLAY(1)
	 JJK=JK+NLEV                                     
         XTEEMIS(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &                    XTEEMIS(JJK,JT))/DTIME       
        ENDDO
        XTEEMIS(NLEV,JT)=(PM(1,NLEV,JT)/PRHOA(1,NLEV)-
     &                       XTEEMIS(NLEV,JT))/DTIME        
  522  CONTINUE
      ENDIF

c LG- end-      -------- emission  --------------------

c LG- start-------prepare dry deposition tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 524 JT=1,NTRAC  
	DO 524 JK=1,NVEGLAY(1) 
	  JJK=JK+NLEV                                 
	  XTEDRYD(JJK,JT)=PM(1,NLEV+JK,JT)/PRHOA(1,JJK)
  524   CONTINUE
      ENDIF

C LG- calculation of the dry deposition tendency using the updated concentration
C     estimate

      IF (LDRYDEP) THEN

       DO 104 JT=1,NTRACT
       DO 104 JL=1,NLON
       DO 104 JK=1,NVEGLAY(JL)
        JJK=JK+NLEV
        QV=1.E5
        IF (RLEAF(JL,JT).EQ.0.OR.JT.EQ.inox) GOTO 104
 
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

        IF (JK.EQ.1) THEN
         RSVEG=RLEAF(JL,JT)/AMAX1(1E-5,LAD_2LAY(JK)*LAI)
        ELSE
         RSVEG=1./((1./(RSOIL(JT)))+
     &     (1./(RLEAF(JL,JT)/AMAX1(1E-5,LAD_2LAY(JK)*LAI))))
        ENDIF

        RTOT=1./MAX(1.E-20,
     &   ((VEGFRAC(JL)/(VEGFRAC(JL)+WSFRAC(JL)))*(1./(RSVEG))+
     &    (WSFRAC(JL)/(VEGFRAC(JL)+WSFRAC(JL)))*(1./(RWS(JT)))))
        VDDEP(JK,JT)=MIN(10.,100./RTOT)	!max Vd of 10 cm s-1

        IF (JT.EQ.inox) THEN
	 IF (NTRAC.EQ.NTRACT) THEN
	  VDDEP(JK,JT)=0.
	 ELSE
	  VDDEP(JK,JT)=
     *     (VDDEP(JK,inoemdd)*PMLOC(JL,JK,ino)+
     *      VDDEP(JK,ino2emdd)*PMLOC(JL,JK,ino2)
     *       )/PM(JL,JK,inox)
         ENDIF
	ENDIF

        QV=100./VDDEP(JK,JT)
        VX=-1./(QV*(PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))/PTMST)
        VD=100.*(PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))/PTMST*
     &     (1.-EXP(VX))
        HEIGHTDD=VD*PTMST/100.
        DDFRAC=HEIGHTDD/(PDP(JL,JJK)/(PRHOA(JL,JJK)*G*1.E3))

        IF (JT.LE.NTRAC) THEN
         BCONC0=PM(JL,NLEV+JK,JT)
         PM(JL,NLEV+JK,JT)=(1.-DDFRAC)*PM(JL,NLEV+JK,JT)
         BXT(1,IBDDEP,JT)=BXT(1,IBDDEP,JT)+(PM(JL,NLEV+JK,JT)-BCONC0)*
     &                  (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK)
        ELSE
         BCONC0=PMLOC(JL,NLEV+JK,JT-NTRAC)
         PMLOC(JL,NLEV+JK,JT-NTRAC)=
     &     (1.-DDFRAC)*PMLOC(JL,NLEV+JK,JT-NTRAC)
         BG3(1,IBDDEP,JT-NTRAC)=BG3(1,IBDDEP,JT-NTRAC)+
     &       (PMLOC(JL,NLEV+JK,JT-NTRAC)-BCONC0)*
     &       (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,JJK)
        ENDIF

 104   CONTINUE

      ENDIF

c LG- end---------calculating dry deposition tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN
       DO 526 JT=1,NTRAC 
       DO 526 JK=1,NVEGLAY(1) 
         JJK=JK+NLEV                                       
         XTEDRYD(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &                    XTEDRYD(JJK,JT))/DTIME         
  526  CONTINUE
      ENDIF

c LG- end-      -------- dry deposition  --------------------

C LG- updating the NOx concentrations and budget only
C     when NTRAC<NTRACT and LEMIS and LDRYDEP=.TRUE.

      IF (NTRAC.LT.NTRACT.AND.LEMIS.AND.LDRYDEP) THEN

        DO 105 JL=1,NLON
         PM(JL,NLEV,inox)=PMLOC(JL,NLEV,ino)+PMLOC(JL,NLEV,ino2)+
     &    PMLOC(JL,NLEV,ino3)+2*PMLOC(JL,NLEV,in2o5)+PMLOC(JL,NLEV,ihno4)

         DO JK=1,NVEGLAY(JL)
          PM(JL,NLEV+JK,inox)=PMLOC(JL,NLEV+JK,ino)+
     &     PMLOC(JL,NLEV+JK,ino2)+PMLOC(JL,NLEV+JK,ino3)+
     &     2*PMLOC(JL,NLEV+JK,in2o5)+PMLOC(JL,NLEV+JK,ihno4)
         ENDDO

 105    CONTINUE
 
      ENDIF

c LG- start-------prepare vertical diffusion tendency--------------------

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 528 JT=1,NTRAC 
          DO JK=1,NVEGLAY(1)
	   JJK=JK+NLEV
           XTEDIFF(JJK,JT)=PM(1,NLEV+JK,JT)/PRHOA(1,JJK)
	  ENDDO
          XTEDIFF(NLEV,JT)=PM(1,NLEV,JT)/PRHOA(1,NLEV)
  528   CONTINUE
      ENDIF

      IF (LXTVDIFF) THEN

       DO 106 JT=1,NTRAC
       DO 106 JL=1,NLON

C LG-   no vertical diffusion of NOx for NTRAC=NTRACT since all the 
C       species are being transported separately

        IF (NTRAC.EQ.NTRACT.AND.JT.EQ.inox) GOTO 106

        IF (HC.GT.HCMIN.AND.
     &     (VEGFRAC(JL)+WSFRAC(JL)).GT.0.) THEN

C LG-   starting with the lowest canopy layer

C LG-   calculating the eddy diffusivity from the Dz between the 
C       reference height of the first and the second canopy layer. The 
C       term RAHCAN is calculated within the routine EC4_VDIFF. There 
C       are two approaches which both calculate this resistance for
C       the canopy height. Therefore the term CANHEIGHT/DZ has been 
C       introduced to correct for the difference between the vertical 
C       extent that has been represented
       
        DZ=PZ(JL,NLEV+1)-PZ(JL,NLEV+2)
	
C LG-   definition of minimal value of Kh to study the sensitivity
C       of model resolved concentrations for small nocturnal exchange
C       fluxes. The minimal value is based on the canopy height and 
C       the number of nocturnal events during whcih the whole canopy is 
C       being ventilated (see paper Fitzjarrald and Moore, JGR 95,
C       16839-16850, 1990). They observed about 3 to 4 exchange events
C       implying one event each two hours [Kh]=m^2 s-1, for a canopy
C       height of 30 m and 2 hours * 3600 s=7200 s this gives a minimal
C       Kh of (30*30)/7200 ~ 0.1 m^2 s-1	
	
        KHMIN=(DZ**2)/DTIME_MIXING	
	KH(2)=MAX(KHMIN,(DZ/(RAHCAN(JL)/(CANHEIGHT/DZ))))	
        DC(JT)=1.E6*(PM(JL,NLEV+2,JT)-PM(JL,NLEV+1,JT))

C LG-   The updated concentrations are calculated such that the 
C       decreasing gradient occuring at a "sub" timestep scale 
C       has been considered in order to prevent the scheme to calculate
C       that large fluxes that negative concentrations are being calculated.
C       The same procedure is also applied to the calculation of decrease
C       of the concentrations due to dry deposition. 

        TERM1=KH(2)/((PDP(JL,NLEV+2)/(PRHOA(JL,NLEV+2)*G*1.E3))*DZ)
        TERM2=KH(2)/((PDP(JL,NLEV+1)/(PRHOA(JL,NLEV+1)*G*1.E3))*DZ)

        ALPHANEG=(-(TERM2+TERM1)+SQRT((TERM2+TERM1)**2))/2.
        ALPHAPOS=(-(TERM2+TERM1)-SQRT((TERM2+TERM1)**2))/2.
	 
        CD=TERM2/(ALPHANEG+TERM2)
	CC=TERM2/(ALPHAPOS+TERM2)
        CB=(PM(JL,NLEV+2,JT)-CC*PM(JL,NLEV+1,JT))/(CD-CC)
        CA=PM(JL,NLEV+1,JT)-CB

        PM_OLD=PM(JL,NLEV+2,JT)
        PM(JL,NLEV+2,JT)=
     &    (CC*CA*EXP(ALPHAPOS*PTMST)+CD*CB*EXP(ALPHANEG*PTMST))

        PM(JL,NLEV+1,JT)=
     &    (CA*EXP(ALPHAPOS*PTMST)+CB*EXP(ALPHANEG*PTMST))

C LG-   calculation of flux between upper and understorey from the dC/dt 
C       and the thickness of the layer (dc/dt=F/dZ => F=dZ*dc/dt)

        FLUXVEG(JL,JT)=-(PDP(JL,NLEV+2)/(PRHOA(JL,NLEV+2)*G*1.E3))*
     &     1.E6*(PM(JL,NLEV+2,JT)-PM_OLD)/PTMST

C LG-   continuing with the canopy top layer

C LG-   calculating the eddy diffusivity from the Dz between the 
C       reference height of the surface layer and the reference height
C       of the canopy and the aerodynamic resistance RAHVEG.

        DZ=PZ(JL,NLEV+1)-PZ(JL,NLEV+2)
        DZ=(PDP(JL,NLEV)/(PRHOA(JL,NLEV)*G*1.E3))/2.+
     &    (CANHEIGHT-PZ(JL,NLEV+1))
        KHMIN=(DZ**2)/DTIME_MIXING	
	KH(1)=MAX(KHMIN,(DZ/RAHVEG(JL)))

C LG-   only calculating a canopy top flux for the areas with a 
C       canopy height > HCMIN and with a vegetation and wet skin fraction
C       > 0, in order to prevent the scheme to calculate fluxes for
C       the snow covered surfaces with a canopy height > HCMIN 

        DC(JT)=1.E6*(PM(JL,NLEV+1,JT)-PM(JL,NLEV,JT))

C LG-   updating the vegetation and surface layer concentration.
C       The updated concentrations are calculated such that the 
C       decreasing gradient occuring at a "sub" timestep scale 
C       has been considered in order to prevent the scheme to calculate
C       that large fluxes that negative concentrations are being calculated.
C       The same procedure is also applied to the calculation of decrease
C       of the concentrations due to dry deposition. The subscript 1 and 2
C       in each parameter refer respectively to the surface layer and 
C       canopy parameters.

        TERM1=KH(1)/((PDP(JL,NLEV)/(PRHOA(JL,NLEV)*G*1.E3))*DZ)
        TERM2=KH(1)/((PDP(JL,NLEV+1)/(PRHOA(JL,NLEV+1)*G*1.E3))*DZ)

        ALPHANEG=(-(TERM2+TERM1)+SQRT((TERM2+TERM1)**2))/2.
        ALPHAPOS=(-(TERM2+TERM1)-SQRT((TERM2+TERM1)**2))/2.
	 
        CD=TERM2/(ALPHANEG+TERM2)
	CC=TERM2/(ALPHAPOS+TERM2)
        CB=(PM(JL,NLEV+1,JT)-CC*PM(JL,NLEV,JT))/(CD-CC)
        CA=PM(JL,NLEV,JT)-CB

        PMVEG_OLD=PM(JL,NLEV+1,JT)
	PM(JL,NLEV+1,JT)=
     &    (CC*CA*EXP(ALPHAPOS*PTMST)+CD*CB*EXP(ALPHANEG*PTMST))

C LG-   the overall grid average concentration is determined by the 
C       both the canopy top and the surface fluxes of the other land cover
C       fractions. We consider the influence of within-canopy interactions
C       on the effective surface fluxes for the wet skin and the vegetation 
C       fraction. 
C       The fact that the surface layer concentrations are calculated from
C       the fractions land cover suggest that we assume that the 
C       concentrations at the reference height are in equilibrium with the 
C       the grid average surface cover characteristics, i.e. that the 
C       surface cover fractions are rather heterogeneously distributed over
C       the whole grid square. Moreover, the concentrations at the reference 
C       heigth are assumed to be located above the blending height, which is 
C       consistent with the calculations of the local aerodynamic resistances
C       (see EC4_VDIFF.f)

C LG-   the surface layer concentrations are only updated if the canopy top
C       fluxes are not used within the vertical diffusion routine as the
C       lower boundary condition (LEMDD_TURB=.FALSE.), else the parameter
C       SURFFLUX is used within EC4_VDIFF.f

        IF (IEMDD_CHEM(JT).EQ.0.AND.IEMDD_TURB(JT).EQ.0) THEN
          PM_OLD=PM(JL,NLEV,JT)
          PM(JL,NLEV,JT)=
     &     (1.-(VEGFRAC(JL)+WSFRAC(JL)))*PM(JL,NLEV,JT)+
     &     (VEGFRAC(JL)+WSFRAC(JL))*(CA*EXP(ALPHAPOS*PTMST)+
     &                               CB*EXP(ALPHANEG*PTMST))

C LG-     checking for mass conservation from the differences in the
C         concentrations for the vegetation layer and the surface layer.
C         The total mass of the vegetation compartment is corrected
C         for the vegetation and wet skin fraction. This is also done
C         since the mass increase of the surface layer is assumed to be
C         the total grid mass increase.

          DMASS_VEG=(PM(JL,NLEV+1,JT)-PMVEG_OLD)*
     &       (VEGFRAC(JL)+WSFRAC(JL))*GRVOL(JL,NLEV+1)
     	  DMASS_SL=(PM(JL,NLEV,JT)-PM_OLD)*GRVOL(JL,NLEV)
	  MASS_VEG=PM(JL,NLEV+1,JT)*(VEGFRAC(JL)+WSFRAC(JL))*
     &             GRVOL(JL,NLEV+1)

          IF ((ABS(DMASS_VEG+DMASS_SL))/MAX(1.E-20,MASS_VEG).GT.
     &         1.E-2) THEN
           PRINT *,'VEG_2LAY.f, no mass cons. for tracer no.: ', 
     &        JT,' and dmass: ',DMASS_VEG+DMASS_SL
           PRINT *,'Percentage mass loss (relat. to canopy mass): ',
     &         100.*(DMASS_VEG+DMASS_SL)/MAX(1.E-20,MASS_VEG)
          ENDIF

         ENDIF

C LG-    calculation of canopy top flux from the dC/dt and the thickness
C        of the layer (dc/dt=F/dZ => F=dZ*dc/dt)

         FLUXCTOP(JL,JT)=-(PDP(JL,NLEV+1)/(PRHOA(JL,NLEV+1)*G*1.E3))*
     &     1.E6*(PM(JL,NLEV+1,JT)-PMVEG_OLD)/PTMST

C LG-    calculation of total surface flux considering all the land cover 
C        fractions

         SURFFLUX(JL,JT)=SURFFLUX(JL,JT)+
     &      (VEGFRAC(JL)+WSFRAC(JL))*FLUXCTOP(JL,JT)

        ENDIF

 106   CONTINUE
 
C LG-  calculating the NOx fluxes, concentrations and budget only
C      when NTRAC=NTRACT

       IF (NTRAC.EQ.NTRACT) THEN

        DO 107 JL=1,NLON
       
         PM(JL,NLEV,inox)=PM(JL,NLEV,ino)+PM(JL,NLEV,ino2)+
     &    PM(JL,NLEV,ino3)+PM(JL,NLEV,in2o5)+PM(JL,NLEV,ihno4)

         DO JK=1,NVEGLAY(JL)
          PM(JL,NLEV+JK,inox)=PM(JL,NLEV+JK,ino)+PM(JL,NLEV+JK,ino2)+
     &     PM(JL,NLEV+JK,ino3)+PM(JL,NLEV+JK,in2o5)+PM(JL,NLEV+JK,ihno4)
         ENDDO

         FLUXCTOP(JL,inox)=FLUXCTOP(JL,ino)+FLUXCTOP(JL,ino2)+
     &                    FLUXCTOP(JL,ino3)+2.*FLUXCTOP(JL,in2o5)+
     &                    FLUXCTOP(JL,ihno4)

         SURFFLUX(JL,inox)=SURFFLUX(JL,ino)+SURFFLUX(JL,ino2)+
     &                    SURFFLUX(JL,ino3)+2.*SURFFLUX(JL,in2o5)+
     &                    SURFFLUX(JL,ihno4)

 107    CONTINUE
 
      ENDIF

C LG- end IF (LXTVDIFF)

      ENDIF

      IF (NSTEP.GE.NSTART.AND.NSTEP.LE.NSTOP) THEN 
        DO 530 JT=1,NTRAC 
	  DO JK=1,NVEGLAY(1)
	   JJK=JK+NLEV
           XTEDIFF(JJK,JT)=(PM(1,NLEV+JK,JT)/PRHOA(1,JJK)-
     &         XTEDIFF(JJK,JT))/DTIME
          ENDDO
          XTEDIFF(NLEV,JT)=(PM(1,NLEV,JT)/PRHOA(1,NLEV)-
     &	       XTEDIFF(NLEV,JT))/DTIME
  530   CONTINUE
      ENDIF


      IF (LEMIS.AND.LXTVDIFF) THEN
 
       DO 108 JL=1,NLON
        CRF(ino)=0.
        CRF(ino2)=0.
        CRF(inox)=0.
        IF (EMISFLUX_VEG(2,inoemdd).GT.0.) THEN
         IF (NTRAC.EQ.NTRACT) THEN
          CRF(ino)=FLUXCTOP(JL,ino)/EMISFLUX_VEG(2,inoemdd)
          CRF(ino2)=FLUXCTOP(JL,ino2)/EMISFLUX_VEG(2,inoemdd)
	 ENDIF
         CRF(inox)=FLUXCTOP(JL,inox)/EMISFLUX_VEG(2,inoemdd)
        ENDIF
        CRF(iisop)=0.
        IF (EMISFLUX_VEG(1,iisop).GT.0.)
     &   CRF(iisop)=FLUXCTOP(JL,iisop)/
     &     (EMISFLUX_VEG(1,iisop)+EMISFLUX_VEG(2,iisop))

 108   CONTINUE

C LG-  writing of canopy reduction factor and other parameters to output file

       IF (NSTEP.EQ.0) THEN
        OPEN(UNIT=NUN2VLAY,FILE='/data/ganzevl/racmo/output/veg_2lay.out',
     *    STATUS='UNKNOWN')
        WRITE(NUN1VLAY,'(1a)')
     *   'Conc., fluxes and canopy reduction factors for 2 layer veg. mode'
        IF (LAIRCHEM) THEN
         WRITE(NUN2VLAY,'(2a10,34a12)')'nstep','time',
     *    'u*[ms-1]','RAHVEG','RAHCAN','Kh(1)','Kh(2)',
     *    'RG','LAD(1)','LAD(2)',
     *    'VdNO','VdNO2','VdO3(1)','VdO3(2)','VdHNO3',
     *    'NOx-veg','NOx-surf','NO2-veg','NO2-surf','O3-veg','O3-surf',
     *    'OH-veg','OH-surf','DC(NOx)','FLUXVEG(O3)','FLUXVEG(no)',
     *    'FLUXCTOP(NO)','FLUXCTOP(O3)','FLUXCTOP(NOx)','FLUX(NOx)',
     *    'NOMOLEC','RATIO(1)','RATIO(2)','CRF(NO)','CRF(NO2)','CRF(NOx)'
    
        ELSE
         WRITE(NUN2VLAY,'(2a10,40a15)')'nstep','time',
     *    'u*[ms-1]','RAHVEG','RAHCAN','Kh(1)','Kh(2)',
     *    'RG','LAD(1)','LAD(2)',
     *    'VdNO','VdNO2','VdO3(1)','VdO3(2)','VdHNO3',
     *    'NOx-veg','NOx-surf','NO2-veg','NO2-surf',
     *    'O3-veg','O3-surf','C5H8-veg','C5H8-surf',
     *    'OH-veg','OH-surf','DC(NOx)','FLUXVEG(O3)','FLUXVEG(NO)',
     *    'FLUXCTOP(O3)','FLUXCTOP(NO)','FLUXCTOP(NOx)','FLUX(NOx)',
     *    'NOMOLEC','RATIO(1)','RATIO(2)','CRF(NO)','CRF(NO2)','CRF(NOx)',
     *    'FLUXCTOP(isop)','FLUX(isop)','ISOPMOLEC','CRF(isop)'
        ENDIF
       ENDIF

       RECALC=(ZMAIR*1.E9/AVO)*(1./PRHOA(1,NLEV))

       IF (MOD(NSTEP,NPRINT).EQ.0) THEN

        IF (NTRAC.EQ.NTRACT) THEN
         IF (LAIRCHEM) THEN
          WRITE(NUN2VLAY,'(1x,i9.9,1x,a9,13f12.4,16e12.4,5f12.4)')
     *     NSTEP,LHRLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),KH(1),KH(2),
     *     RG,LAD_2LAY(1),LAD_2LAY(2),VDDEP(1,ino),VDDEP(1,ino2),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,ihno3),
     *     RECALC*PM(1,NLEV+1,inox),RECALC*PM(1,NLEV,inox),
     *     RECALC*PM(1,NLEV+1,ino2),RECALC*PM(1,NLEV,ino2),
     *     RECALC*PM(1,NLEV+1,io3),RECALC*PM(1,NLEV,io3), 
     *     PM(1,NLEV+1,ioh),PM(1,NLEV,ioh),DC(inox),
     *     FLUXVEG(1,io3),FLUXVEG(1,ino),FLUXCTOP(1,io3),FLUXCTOP(1,ino),
     *     FLUXCTOP(1,inox),SURFFLUX(1,inox),EMISFLUX_VEG(2,ino),
     *     RATIO(1),RATIO(2),CRF(ino),CRF(ino2),CRF(inox)
         ELSE
          WRITE(NUN2VLAY,
     *      '(1x,i9.9,1x,a9,13f15.4,18e15.4,5f15.4,3e15.4,f15.4)')
     *     NSTEP,LHRLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),KH(1),KH(2),
     *     RG,LAD_2LAY(1),LAD_2LAY(2),VDDEP(1,ino),VDDEP(1,ino2),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,ihno3),
     *     RECALC*PM(1,NLEV+1,inox),RECALC*PM(1,NLEV,inox),
     *     RECALC*PM(1,NLEV+1,ino2),RECALC*PM(1,NLEV,ino2),
     *     RECALC*PM(1,NLEV+1,io3),RECALC*PM(1,NLEV,io3), 
     *     RECALC*PM(1,NLEV+1,iisop),RECALC*PM(1,NLEV,iisop),
     *     PM(1,NLEV+1,ioh),PM(1,NLEV,ioh),DC(inox),
     *     FLUXVEG(1,io3),FLUXVEG(1,ino),FLUXCTOP(1,io3),FLUXCTOP(1,ino),
     *     FLUXCTOP(1,inox),SURFFLUX(1,inox),EMISFLUX_VEG(2,ino),
     *     RATIO(1),RATIO(2),CRF(ino),CRF(ino2),CRF(inox),
     *     FLUXCTOP(1,iisop),SURFFLUX(1,iisop),ISOPMOLEC,
     *     CRF(iisop)
         ENDIF
        ELSE
	 IF (LAIRCHEM) THEN
          WRITE(NUN2VLAY,'(1x,i9.9,1x,a9,13f12.4,16e12.4,5f12.4)')
     *     NSTEP,LHRLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),KH(1),KH(2),
     *     RG,LAD_2LAY(1),LAD_2LAY(2),VDDEP(1,inoemdd),VDDEP(1,ino2emdd),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,ihno3),
     *     RECALC*PM(1,NLEV+1,inox),RECALC*PM(1,NLEV,inox),
     *     RECALC*PMLOC(1,NLEV+1,ino2),RECALC*PMLOC(1,NLEV,ino2),
     *     RECALC*PM(1,NLEV+1,io3),RECALC*PM(1,NLEV,io3), 
     *     PMLOC(1,NLEV+1,ioh),PMLOC(1,NLEV,ioh),DC(inox),
     *     FLUXVEG(1,io3),FLUXVEG(1,inoemdd),FLUXCTOP(1,io3),
     *     FLUXCTOP(1,inoemdd),FLUXCTOP(1,inox),SURFFLUX(1,inox),
     *     EMISFLUX_VEG(2,inoemdd),RATIO(1),RATIO(2),CRF(inoemdd),
     *     CRF(ino2emdd),CRF(inox)
         ELSE
          WRITE(NUN2VLAY,
     *      '(1x,i9.9,1x,a9,13f15.4,18e15.4,5f15.4,3e15.4,f15.4)')
     *     NSTEP,LHRLTIME,USTAR(1),RAHVEG(1),RAHCAN(1),KH(1),KH(2),
     *     RG,LAD_2LAY(1),LAD_2LAY(2),VDDEP(1,inoemdd),VDDEP(1,ino2emdd),
     *     VDDEP(1,io3),VDDEP(2,io3),VDDEP(1,ihno3),
     *     RECALC*PM(1,NLEV+1,inox),RECALC*PM(1,NLEV,inox),
     *     RECALC*PMLOC(1,NLEV+1,ino2),RECALC*PMLOC(1,NLEV,ino2),
     *     RECALC*PM(1,NLEV+1,io3),RECALC*PM(1,NLEV,io3), 
     *     RECALC*PM(1,NLEV+1,iisop),RECALC*PM(1,NLEV,iisop),
     *     PMLOC(1,NLEV+1,ioh),PMLOC(1,NLEV,ioh),DC(inox),
     *     FLUXVEG(1,io3),FLUXVEG(1,inoemdd),FLUXCTOP(1,io3),
     *     FLUXCTOP(1,inoemdd),FLUXCTOP(1,inox),SURFFLUX(1,inox),
     *     EMISFLUX_VEG(2,inoemdd),RATIO(1),RATIO(2),CRF(inoemdd),
     *     CRF(ino2emdd),CRF(inox),FLUXCTOP(1,iisop),SURFFLUX(1,iisop),
     *     ISOPMOLEC,CRF(iisop)
         ENDIF
        
	ENDIF

       ENDIF

       IF (NSTEP.EQ.NSTOP) CLOSE(NUN2VLAY)

C LG- end IF (LEMIS.AND.LXTVDIFF)

      ENDIF

C LG- 

      WRITE(NUNMDFL,*)'End VEG_2LAY.f'

      RETURN
      END


