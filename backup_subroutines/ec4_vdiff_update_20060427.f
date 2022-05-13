      SUBROUTINE EC4_VDIFF ( KIDIA,KFDIA,KLON,KLP2,KTDIA,KLEV,KLEVM1
     *     , KLEVP1

C     LG- added the number of soil layers, required for AGS.f

     *     , KLEVS

C     LG- some chemistry parameters

     *     , KLEVEL , KLEVELM1, KLEVELP1 
     *     , KTRAC  , NPRINT  

C     LG- end

     *     , CONACC
     *     , APZERO
     *     , NSTEP,NSTART,NSTOP,TWODT,EPS
     *     , LVDIFF,LSURF
C-----------------------------------------------------------------------
C     - INPUT 2D .
     *     , PACLCM,   PAPHM1,   PAPM1,   PGEOM1,  PQM1,    PTKEM
     *     , PTKEM1M,  PTM1,     PUM1,    PVM1,    PXM1
      
C     LG- update 03-2000, GL scheme 

     *     , PTVM1   
      
C     LG- adding the tracers

     *     , PXTM1,    PXTMVEG, PXTMSN     
      
C     LG- and the soil type, required for the subroutine AGS.f

     *     , PSOTYPE

C     LG- end     
      
C     - INPUT 1D .
     *     , LALAND
     *     , PAHFLM,   PAHFSM,   PAZ0M,   PDEW2M,  PEVAPM,  PFORESTM
     *     , PSEAICE,  PSNM1M,   PSRFL,   PTEMP2M, PTSM1M,  PT2MAXM
     *     , PT2MINM,  PUSTAR3M, PUSTRM,  PU10M,   PVDISM,  PVSTRM
     *     , PV10M,    PWIMAXM,  PWIND10M,PWLM1M,  PWSM1M
     *     , PWSMXM,   PVLTM
     *     , PBLHM,    PTKEVIM
     *     , PTSM,     PWSM
C     - INPUT 1D LAM SPECIFIC
     *     , PSINLAT
C     - OUTPUT 2D .
     *     , PTKE,     PTKEM1,   DTKEDT
     *     , PKE
C     - OUTPUT 1D .
     *     , PAHFL,    PAHFS,    PAZ0,    PCVGHL,  PCVS,    PCVW
     *     , PDEW2,    PDHFQS,   PDHFQW,  PDHFT,   PEVAP
     *     , PQHFL,    PRSFL,    PTEMP2,  PTHFL,   PT2MAX,  PT2MIN
     *     , PUSTAR3,  PUSTR,    PU10,    PVDIS,   PVSTR,   PV10
     *     , PWIMAX,   PWLMX,    PWIND10, PXHFL
     *     , PBLH,     PTKEVI
C     - INPUT/OUTPUT 2D .
     *     , PVOL,     PVOM,     PQTE,     PTTE,     PXTE

C     LG- adding the tracers

     *     , PXTTE,PXTEEMIS,PXTEDRYD,PXTETRNS,PXTESTRL
     *     , PXTEDIFF,PXTEDIFF_VEG

C     LG- end

C     - INPUT/OUTPUT 1D .
     *     , PVGRAT 
      
C     LG- pressure differences within canopy, surface temperature, 
C     aerodynamic resistance for reference height zz, emissivity and
C     surface flux

     *     , DAPM1, DAPHM1, PTS, RAHZZ, PEMTERM, SURFFLUX, PTD3
     *     , UOBS , VOBS,   WOBS,TOBS,  QOBS,    XOBS)     

C     LG- end
C     
C**** *VDIFF* - DOES THE VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
C     
C     
C     SUBJECT.
C     --------
C     
C     THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE FOUR
C     PROGNOSTIC VARIABLES U,V,T AND Q DUE TO THE VERTICAL EXCHANGE BY
C     TURBULENT (= NON-MOIST CONVECTIVE) PROCESSES. THESE TENDENCIES ARE
C     OBTAINED AS THE DIFFERENCE BETWEEN THE RESULT OF AN IMPLICIT
C     TIME-STEP STARTING FROM VALUES AT T-1 AND THESE T-1 VALUES. ALL
C     THE DIAGNOSTIC COMPUTATIONS (EXCHANGE COEFFICIENTS, ...) ARE DONE
C     FROM THE T-1 VALUES. AS A BY-PRODUCT THE ROUGHNESS LENGTH OVER SEA
C     IS UPDATED ACCORDINGLY TO THE *CHARNOCK FORMULA. HEAT AND MOISTURE
C     SURFACE FLUXES AND THEIR DERIVATIVES AGAINST TS, WS AND WL
C     (THE LATTER WILL BE LATER WEIGHTED WITH THE SNOW FACTOR IN
C     *VDIFF*), LATER TO BE USED FOR SOIL PROCESSES TREATMENT, ARE ALSO
C     COMPUTED AS WELL AS A STABILITY VALUE TO BE USED AS A DIAGNOSTIC
C     OF THE DEPTH OF THE WELL MIXED LAYER IN CONVECTIVE COMPUTATIONS.
C     
C**   INTERFACE.
C     ----------
C     
C     *VDIFF* IS CALLED FROM *PHYSC*.
C     
C     INPUT ARGUMENTS.
C     ----- ----------
C     
C     - 3D
C     PXTM1    : TRACER VARIABLES (T-DT)
C     - 2D
C     PACLCM   : CLOUD COVER (OLD VALUE)
C     PAPHM1   : HALF LEVEL PRESSURE (T-DT)
C     PAPM1    : FULL LEVEL PRESSURE (T-DT)
C     PGEOM1   : GEOPOTENTIAL ABOVE SURFACE (T-DT)
C     PQM1     : HUMIDITY (T-DT)
C     PTKEM    : TURBULENT KINETIC ENERGY
C     PTKEM1M  : TURBULENT KINETIC ENERGY (T-DT)
C     PTM1     : TEMPERATURE (T-DT)
C     PUM1     : ZONAL WIND (T-DT)
C     PVM1     : MERIDIONAL WIND (T-DT)
C     PXM1     : CLOUD WATER (T-DT)
C     - 1D
C     LALAND   : LAND-SEA FLAG
C     PAHFLM   : SURFACE LATENT HEAT FLUX (OLD VALUE)
C     PAHFSM   : SURFACE SENSIBLE HEAT FLUX (OLD VALUE)
C     PAZ0M    : ROUGHNESS LENGTH (OLD VALUE)
C     PDEW2M   : DEW POINT TEMPERATURE AT 2 METER (ACCUMULATED, OLD VALUE)
C     PEVAPM   : SURFACE EVAPORATION (ACCUMULATED, OLD VALUE)
C     PFORESTM : FOREST COVERAGE
C     PSEAICE  : SEA ICE COVER (NEW VALUE)
C     PSNM1M   : SNOW DEPTH (T-DT)
C     PSRFL    : NET SOLAR RADIATIVE FLUX AT THE SURFACE
C     PTEMP2M  : TEMPERATURE AT 2 METER (ACCUMULATED, OLD VALUE)
C     PTSM1M   : SURFACE TEMPERATURE (T-DT)

C     LG- update 03-2000

C     #ifdef EVM
C     PTVM1    : VIRTUAL TEMPERATURE AT T-DT
C     #endif

C     LG- end

C     PT2MAXM  : MAXIMUM TEMP. AT 2 M. BETWEEN OUTPUT INTERVALS (OLD VALUE)
C     PT2MINM  : MINIMUN TEMP. AT 2 M. BETWEEN OUTPUT INTERVALS (OLD VALUE)
C     PUSTAR3M : TKE FOR OCEAN MIXED LAYER (ACCUMULATED, OLD VALUE)
C     PUSTRM   : U-STRESS (ACCUMULATED, OLD VALUE)
C     PVSTRM   : V-STRESS (ACCUMULATED, OLD VALUE)
C     PU10M    : U-WIND AT 10 METER (ACCUMULATED, OLD VALUE)
C     PV10M    : V-WIND AT 10 METER (ACCUMULATED, OLD VALUE)
C     PWIND10M : WIND SPEED AT 10 METER (ACCUMULATED, OLD VALUE)
C     PWIMAXM  : MAXIMUM WINDSPEED AT 10 M. BETW. OUTP. INTERV. (OLD VALUE)
C     PVDISM   : BOUNDARY LAYER DISSIPATION (ACCUMULATED, OLD VALUE)
C     PWLM1M   : SKIN RESERVOIR CONTENT (T-DT)
C     PWSM1M   : SURFACE SOIL WETNESS (T-DT)
C     PWSMXM   : FIELD CAPACITY OF SOIL
C     PVLTM    : LEAF AREA INDEX
C     PBLHM    : PLANETARY BOUNDARY LAYER HEIGHT (HOLTSLAG ET AL.)
C     PTKEVIM  : VERTICALLY INTEGRATED TURBULENT KINETIC ENERGY
C     EVM
C     PTSM     : SURFACE TEMPERATURE
C     PWSM     : SURFACE SOIL WETNESS
C     
C     OUTPUT ARGUMENTS.
C     ------ ----------
C     
C     - 2D
C     PTKE     : TURBULENT KINETIC ENERGY (T+DT)
C     PTKEM1   : TURBULENT KINETIC ENERGY (FILTERED)
C     DTKEDT   : TKE-TENDENCY
C     PKE      : KINETIC ENERGY LOSS FOR EACH LAYER
C     - 1D
C     KTROPO   : TROPOPAUSE INDEX
C     PAHFL    : SURFACE LATENT HEAT FLUX (NEW VALUE)
C     PAHFS    : SURFACE SENSIBLE HEAT FLUX (NEW VALUE)
C     PAZ0     : ROUGHNESS LENGTH (NEW VALUE)
C     PCVGHL   : RATIO OF MOISTURE FLUXES
C     PCVS     : SNOW COVER FRACTION
C     PCVW     : WET SKIN FRACTION
C     PDEW2    : DEW POINT TEMPERATURE AT 2 METER (ACCUMULATED, NEW VALUE)
C     PDHFQS   : DERIVITAVE OF MOISTURE FLUX OVER SNOW WITH RESP. TO SNOW D
C     PDHFQW   : DERIVITAVE OF MOISTURE FLUX WITH RESPECT TO SKIN RESERVOIR
C     PDHFT    : DERIVITAVE OF SENSIBLE HEAT FLUX WITH RESP. TO SURF. TEMP.
C     PEVAP    : SURFACE EVAPORATION (ACCUMULATED, NEW VALUE)
C     PQHFL    : MOISTURE FLUX AT THE SURFACE
C     PRSFL    : LARGE SCALE RAIN FLUX AT THE SURFACE
C     PTHFL    : SENSIBLE HEAT FLUX AT THE SURFACE
C     PXHFL    : LIQUID WATER FLUX AT THE SURFACE
C     PBLH     : PLANETARY BOUNDARY LAYER HEIGHT (HOLTSLAG ET AL.)
C     PTKEVI   : VERTICALLY INTEGRATED TURBULENT KINETIC ENERGY
C     PTEMP2   : TEMPERATURE AT 2 METER (ACCUMULATED, NEW VALUE)
C     PT2MAX   : MAXIMUM TEMP. AT 2 M. BETWEEN OUTPUT INTERVALS (NEW VALUE)
C     PT2MIN   : MINIMUN TEMP. AT 2 M. BETWEEN OUTPUT INTERVALS (NEW VALUE)
C     PUSTAR3  : TKE FOR OCEAN MIXED LAYER (ACCUMULATED, NEW VALUE)
C     PUSTR    : U-STRESS (ACCUMULATED, NEW VALUE)
C     PVSTR    : V-STRESS (ACCUMULATED, NEW VALUE)
C     PU10     : U-WIND AT 10 METER (ACCUMULATED, NEW VALUE)
C     PV10     : V-WIND AT 10 METER (ACCUMULATED, NEW VALUE)
C     PVDIS    : BOUNDARY LAYER DISSIPATION (ACCUMULATED, NEW VALUE)
C     PWIMAX   : MAXIMUM WINDSPEED AT 10 M. BETW. OUTP. INTERV. (NEW VALUE)
C     PWLMX    : MAXIMUM SKIN RESERVOIR CONTENT
C     PWIND10  : WIND SPEED AT 10 METER (ACCUMULATED, NEW VALUE)
C     
C     INPUT/OUTPUT ARGUMENTS.
C     -----------------------
C     
C     - 3D
C     PXTTE    : TENDENCIES OF TRACER VARIABLES
C     - 2D
C     PVOL     : TENDENCY OF MERIDIONAL WIND
C     PVOM     : TENDENCY OF ZONAL WIND
C     PQTE     : TENDENCY OF HUMIDITY
C     PTTE     : TENDENCY OF TEMPERATURE
C     PXTE     : TENDENCY OF CLOUD WATER
C     - 1D
C     PVGRAT   : VEGETATION RATIO
C     
C     
C     METHOD.
C     -------
C     
C     FIRST AN AUXIALIARY VARIABLE CP(Q)T+GZ IS CREATED ON WHICH
C     THE VERTICAL DIFFUSION PROCESS WILL WORK LIKE ON U,V AND Q. THEN
C     ALONG THE VERTICAL AND AT THE SURFACE, EXCHANGE COEFFICIENTS (WITH
C     THE DIMENSION OF A PRESSURE THICKNESS) ARE COMPUTED FOR MOMENTUM
C     AND FOR HEAT (SENSIBLE PLUS LATENT). THE LETTERS M AND H ARE USED
C     TO DISTINGUISH THEM. THE DIFFUSIONCOEFFICENTS DEPEND ON THE
C     TURBULENT KINETIC ENERGY (TKE) CALCULATED BY AN ADDITIONAL
C     PROGNOSTIC EQUATION, WHICH CONSIDERS ADVEKTION OF TKE.
C     IN THE SECOND PART OF THE ROUTINE THE IMPLICIT LINEAR
C     SYSTEMS FOR U,V FIRST AND T,Q SECOND ARE SOLVED BY A *GAUSSIAN
C     ELIMINATION BACK-SUBSTITUTION METHOD. FOR T AND Q THE LOWER
C     BOUNDARY CONDITION DEPENDS ON THE SURFACE STATE.
C     FOR TKE THE LOWER BOUNDARY CONDITION DEPENDS ON THE SQUARE OF
C     THE FRICTIONAL VELOCITY.
C     OVER LAND, TWO DIFFERENT REGIMES OF EVAPORATION PREVAIL:
C     A STOMATAL RESISTANCE DEPENDENT ONE OVER THE VEGETATED PART
C     AND A SOIL RELATIVE HUMIDITY DEPENDENT ONE OVER THE
C     BARE SOIL PART OF THE GRID MESH.
C     POTENTIAL EVAPORATION TAKES PLACE OVER THE SEA, THE SNOW
C     COVERED PART AND THE LIQUID WATER COVERED PART OF THE
C     GRID MESH AS WELL AS IN CASE OF DEW DEPOSITION.
C     FINALLY ONE RETURNS TO THE VARIABLE TEMPERATURE TO COMPUTE
C     ITS TENDENCY AND THE LATER IS MODIFIED BY THE DISSIPATION'S EFFECT
C     (ONE ASSUMES NO STORAGE IN THE TURBULENT KINETIC ENERGY RANGE) AND
C     THE EFFECT OF MOISTURE DIFFUSION ON CP. Z0 IS UPDATED AND THE
C     SURFACE FLUXES OF T AND Q AND THEIR DERIVATIVES ARE PREPARED AND
C     STORED LIKE THE DIFFERENCE BETWEEN THE IMPLICITELY OBTAINED
C     CP(Q)T+GZ AND CP(Q)T AT THE SURFACE.
C     
C     EXTERNALS.
C     ----------
C     ....
C     
C     REFERENCE.
C     ----------
C     
C     SEE VERTICAL DIFFUSION'S PART OF THE MODEL'S DOCUMENTATION
C     FOR DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
C     
C     AUTHOR.
C     -------
C     U. SCHLESE     DKRZ-HAMBURG  FEB-93
C     MODIFIED     E. ROECKNER  - 1994
C     
C     BASED  ON  ORIGINAL ECMWF VERSION BY J.F. GELEYN  - 1982
C     MODIFIED BY C.B. BLONDIN - 1986
C     H. FEICHTER  - 1991
C     S. BRINKOP   - 1992
C     M. CLAUSSEN  - 1993
C     
C     MODIFICATIONS.
C     --------------
C     
C     JHC*CALL PARAM
C     JHC*CALL COMCTL
C     JHC*CALL COMGAU
C     JHC*CALL COMDIAZ
C     JHC*CALL COMPSW
C     JHC*CALL COMPH2
C     JHC*CALL COMCON
C     JHC*CALL COMSDS
C     JHC*CALL COMVEG
C     JHC*CALL COMTRC
C     JHC*CALL YOTLUC

C     LG- 
      
      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'
      INCLUDE 'comddrs.h'
      
C     LG- end     

      INCLUDE 'comcon.h'

C     LG- it turned out that something went wrong reading input files with
C     surface flux data due the fact that the common block comio.h was not 
C     included in this routine

      INCLUDE 'comio.h'

C     LG- end

      INCLUDE 'comph2.h'
      INCLUDE 'paramsoil.h'
      INCLUDE 'comveg.h'
      INCLUDE 'yotluc.h'
C     ---------------------------------
      LOGICAL LO,LO1
C     DIR$ VFUNCTION EXPHF, SQRTHF
C     *vdir vectfunc (sqrthf)
C     DIR$ VFUNCTION LOGHF
C     ----------------
C     
      LOGICAL LALAND(KLP2) ,LVDIFF,LSURF

      REAL PACLCM(KLP2,KLEV),  PAPHM1(KLP2,KLEVP1), PAPM1(KLP2,KLEV)
     *     ,PGEOM1(KLP2,KLEV),  PQM1(KLP2,KLEV),     PTKEM(KLP2,KLEV)
     *     ,PTKEM1M(KLP2,KLEV), PTM1(KLP2,KLEV),     PUM1(KLP2,KLEV)
     *     ,PVM1(KLP2,KLEV),    PXM1(KLP2,KLEV)

C     LG- update 03-2000, GL scheme

     *     ,PTVM1(KLP2,KLEV)

C     LG- adding the tracers

      REAL PXTM1(KLP2,NLEVT,KTRAC), PXTMVEG(KLP2,NLEVV_ML,KTRAC),
     *     PXTMSN(KLP2,NLEVS_ML,KTRAC)

C     LG- end

      REAL PAHFLM(KLP2),       PAHFSM(KLP2),        PAZ0M(KLP2)
     *     ,PDEW2M(KLP2),       PEVAPM(KLP2),        PFORESTM(KLP2)
     *     ,PSEAICE(KLP2),      PSNM1M(KLP2),        PSRFL(KLP2)
     *     ,PTEMP2M(KLP2),      PTSM1M(KLP2),        PT2MAXM(KLP2)
     *     ,PT2MINM(KLP2),      PUSTAR3M(KLP2),      PUSTRM(KLP2)
     *     ,PU10M(KLP2),        PVDISM(KLP2),        PVSTRM(KLP2)
     *     ,PV10M(KLP2),        PWIMAXM(KLP2),       PWIND10M(KLP2)
     *     ,PWLM1M(KLP2),       PWSM1M(KLP2),        PWSMXM(KLP2)
     *     ,PVLTM(KLP2) ,       PSINLAT(KLP2)
     *     ,PBLHM(KLP2),        PTKEVIM(KLP2)
     *     ,PTSM(KLP2),         PWSM(KLP2)
      REAL PTKE(KLP2,KLEV),    PTKEM1(KLP2,KLEV),   DTKEDT(KLP2,KLEV)
     *     ,PKE(KLP2,KLEV)
      REAL PAHFL(KLP2),        PAHFS(KLP2),         PAZ0(KLP2)
     *     ,PCVGHL(KLP2),       PCVS(KLP2),          PCVW(KLP2)
     *     ,PDEW2(KLP2),        PDHFQS(KLP2)
     *     ,PDHFQW(KLP2),       PDHFT(KLP2),         PEVAP(KLP2)
     *     ,PQHFL(KLP2),        PRSFL(KLP2),         PTEMP2(KLP2)
     *     ,PTHFL(KLP2),        PT2MAX(KLP2),        PT2MIN(KLP2)
     *     ,PUSTAR3(KLP2),      PUSTR(KLP2),         PU10(KLP2)
     *     ,PVDIS(KLP2),        PVSTR(KLP2),         PV10(KLP2)
     *     ,PWIMAX(KLP2),       PWLMX(KLP2),         PWIND10(KLP2)
     *     ,PXHFL(KLP2)
     *     ,PBLH(KLP2),         PTKEVI(KLP2)

C     LG- adding the tracers

      REAL PXTTE(KLON,NLEVT,NTRAC)
     *     ,PXTEEMIS(NLEVT,NTRACT) ! emission tendency  [molec cm-3 s-1]
     *     ,PXTEDRYD(NLEVT,NTRACT) ! dry deposition tendency  [molec cm-3 s-1]
     *     ,PXTETRNS(NLEVT,NTRACT) ! neg. conc. tendency  (see TRANSP.f)
     *     ,PXTESTRL(NLEVT,NTRACT) ! strat. proc. tendency (see STRFILL.f)
     *     ,PXTEDIFF(NLEVT,NTRAC) ! vertical diffusion tendency [molec cm-3 s-1] 
     *     ,PXTEDIFF_VEG(NLEVT,NTRAC) ! vertical diffusion tendency (see bulkveg.f) 

C     LG- end

      REAL PVOL(KLP2,KLEV),    PVOM(KLP2,KLEV),     PQTE(KLP2,KLEV)
     *    ,PTTE(KLP2,KLEV),    PXTE(KLP2,KLEV)
      REAL PVGRAT(KLP2)

! mz_lg_20060427+
CGEERT
      REAL PMFU(KLP2,KLEV), PLPARCEL(KLP2)
! mz_lg_20060427-

C     LG- declaration of surface temperature, the soil temperature and emissivity

      REAL PTS(KLP2),	PTD3(KLP2),	PEMTERM (KLP2,KLEVP1)

C     LG- end

C     
C     TEMPORARY ARRAYS
C     
      INCLUDE 'paramh.h'
      INCLUDE 'paramv100.h'

C     LG- adding the tracers

      REAL ZXTDIF(KLON,NLEVT,KTRAC), ZXTEMS(KLON,KTRAC),
     *     EM(KLON,NLEVT,NTRACT),   EM_TOT(KLON,NTRACT),
     *     DD(KLON,NLEVT,NTRACT),   DD_TOT(KLON,NTRACT),
     *     SURFFLUX(KLON,NTRACT)

C     LG- end

      REAL ZCFM(JPHR,MLEV),    ZDIS(JPHR,MLEV)
     *     ,ZCFH(JPHR,MLEV),    ZCPTGZ(JPHR,MLEV),   ZEBSM(JPHR,MLEV)
     *     ,ZUDIF(JPHR,MLEV),   ZVDIF(JPHR,MLEV)
     *     ,ZWET(JPHR),         ZQS(JPHR),           ZDQS(JPHR)
     *     ,ZCPTS(JPHR),        ZTVS(JPHR),          ZRI(JPHR)
     *     ,ZUCF(JPHR),         ZSCF(JPHR),          ZCFNC(JPHR)
     *     ,ZCDN(JPHR),         ZTCOE(JPHR),         ZWLMXI(JPHR)
      REAL ZCR(JPHR),          ZRS0(JPHR),          ZHUM(JPHR)
     *     ,ZCSAT(JPHR),        ZCAIR(JPHR),         ZTDIF(JPHR,MLEV)
     *     ,ZQDIF(JPHR,MLEV),   ZEBSH(JPHR,MLEV),    ZVIDIS(JPHR)
     *     ,Z1MXTM1(JPHR)
      REAL ZBN(JPHR),          ZBM(JPHR),           ZBH(JPHR)
     *     ,ZCHN(JPHR),  ZCH(JPHR),  ZUSTAR(JPHR), ZWST(JPHR)
     *     ,ZRICLS(JPHR), ZTESS(JPHR), ZCFNCH(JPHR)
     *     ,ZHSOIL(JPHR)
     *     ,ZHDYN(JPHR),        ZTETA1(JPHR,MLEV),   ZLTETA1(JPHR,MLEV)
     *     ,ZTVIR1(JPHR,MLEV),  ZHH(JPHR,MLEVM1),    ZQSS(JPHR,MLEV)
     *     ,ZXDIF(JPHR,MLEV),   ZEDIF(JPHR,MLEV),    ZTKEVN(JPHR,MLEV)
      REAL ZQSSM(JPHR,MLEVM1),ZTMITTE(JPHR,MLEVM1),ZTVIRMIT(JPHR,MLEVM1)
     *     ,ZFAXEN(JPHR,MLEVM1),ZFAXE(JPHR,MLEV),    ZCCOVER(JPHR,MLEVM1)
     *     ,ZLWCMIT(JPHR,MLEVM1),ZTEMIT(JPHR,MLEVM1),ZQMIT(JPHR,MLEVM1)
     *     ,ZCDUM(JPHR,MLEV)
      REAL ZSENKF (JPHR),       ZLATKF (JPHR),      ZUSTAR1(JPHR)
     *     ,ZCDH(JPHR)   ,       ZCDM   (JPHR)
     *     ,ZPBLH  (JPHR),       ZOBUKL (JPHR),      ZBUOYPR(JPHR)

      REAL zxtvn(JPHR,klev,ktrac) ! op_ck_20031001

C     LG- update 03-2000, GL scheme

C     GL...
      REAL  ZGBUOY(JPHR,MLEV), ZGSHEAR(JPHR,MLEV), ZGRI(JPHR,MLEV)
     *     ,ZMIXDWH(JPHR, MLEVP1), ZMIXUPH(JPHR, MLEVP1)          
     *     ,ZMIXDWM(JPHR, MLEVP1), ZMIXUPM(JPHR, MLEVP1)
     *     ,ZMIXQUADH(JPHR,MLEVP1),ZMIXQUADM(JPHR,MLEVP1)
     *     ,ZGMONIN(JPHR), ZGWST(JPHR)
     *     ,PTM1PR(JPHR,MLEV)  ! mz_lg_20060427+

! mz_lg_20060427+
      REAL 
     *    DZ_PARCEL(JPHR), TKE_PARCEL_UP(JPHR), TH_PARCEL_UP(JPHR), 
     *    TKE_CBL(JPHR)

      REAL  ZVARQTMF(JPHR,MLEV) ,ZQLMF(JPHR,MLEV), ZPMFU(JPHR,MLEV)

! mz_lg_20060427-
C     GL...

      INTEGER IHPBL(JPHR),     IHPBLC(JPHR),        IHPBLD(JPHR)
C     EVM
      REAL ZDPH   (JPHR,MLEVP1),ZBETA(JPHR)
     *     ,ZCFHSRF(JPHR)       ,WRIH(JPHR)

C     LG- declaration of extra parameters

      REAL ZXTM1(JPHR,NLEVT,NTRAC),ZCFTR(JPHR,NLEVT),
     *     ZTCOETR(JPHR,NTRAC),ZKH(JPHR,NLEVT),HPBL(JPHR),
     *     ZLM(NLEVV),ZIW(NLEVV),ZKH_NEW(JPHR,NLEVT),ZU(NLEVT),
     *     ZTETAS(JPHR),ZTETASOIL(JPHR),ZRIB(JPHR,NLEVT),
     *     ZQSA(JPHR),ZXTMCO2(JPHR,NLEVT),ZDIR,DQWS_DQAVG

      REAL FREQ,GUST,ZSTAB,ZU_RATIO,X,ZLEAFWIDTH,DISPSNOW

      INTEGER NLEVMIX,IZPBLH(JPHR),NPRINT,N,NSTEP_GUST,ISETZ0,
     *     TIME_SUNRISE,TIME_MIXING, ! time of sunrise and onset of mixing
     *     DHR,DMIN

C     LG- LGUST is a switch to consider the role of gusts on the turbulent exchange
C     by simply scaling the calculated K values. LMXL is a switch to calculate 
C     the eddy-diffusivity applying the mixing lenght theory with the mixing 
C     lenght being estimated within the canopy according to Goudriaan, 
C     Crop Micrometeorology. LSENS_ZTMST is the switch to define some parameters
C     to test the sensitivity of the vertical exchange calculations for the
C     timestep

      LOGICAL LGUST, LMXL, LSENS_ZTMST

      SAVE FREQ,ISETZ0,TIME_SUNRISE,TIME_MIXING

C     LG- added some parameters for the routine AGS.f, May 2001

      INTEGER KLEVS		! the number of soil layers

      REAL PSOTYPE(JPHR)	! the soil type
      REAL PCATYPE(JPHR)	! the C3/C4 vegetation type
      REAL PWSAM1M(JPHR,KLEVS)	! the soil water content of KLEVS layers
      REAL PWET_AGS(JPHR)	! the physiological stomatal resistance,
                                ! corrected for soil water waterstress
      REAL PRS0_AGS(JPHR)	! the physiological stomatal resistance,
                                ! not corrected 
      REAL PWET(JPHR)		! the original ECHAM stomatal resistance
      REAL PRS0_AGSML(JPHR,NLEVV) ! the physiological stomatal resistance,
                                ! of the multilayer model

C     LG- parameters required for the calculation of rho, height of half pressure
C     levels and the halp level pressure within the canopy

      REAL ZRHOA(JPHR,NLEVT),ZT(JPHR,NLEVT),ZP(JPHR,NLEVT),
     *     ZQ(JPHR,NLEVT),ZDP(JPHR,NLEVT),
     *     Z(JPHR,NLEVT),ZZDZ(JPHR,NLEVT),ZHP(JPHR,NLEVT+1),
     *     ZAPHM1(JPHR,NLEVT+1),ZAPM1(JPHR,NLEVT),
     *     DP(NLEVT),DHP(NLEVT+1),DDP(NLEVT+1),DDP_Z(NLEVT),
     *     DDP_ZHP(NLEVT+1)

C     LG- pressure differences between levels within the canopy. 

      REAL DAPM1(JPHR,NLEVT),DAPHM1(JPHR,NLEVT+1)

      REAL TVIRT

C     LG- end

      REAL CSB,CEPS,CLV,CR,CEMISIR,DQDT,FRACRAD_SOIL,
     &     LABDA,ZRHOCHU,TERM1,TERM2,SOILFLUX

C     LG- added the observed windspeed

      REAL UOBS(JPHR,KLEV),VOBS(JPHR,KLEV),UVOBS(JPHR,KLEV),
     &     WOBS(JPHR,KLEV),TOBS(JPHR,KLEV),QOBS(JPHR,KLEV),
     &     XOBS(JPHR,KLEV)

C     LG- added the parameters to calculate explicitly the water resistance

      REAL REYN(JPHR),WEIGHT(JPHR),AZ0S(JPHR),AVGWS(JPHR)
      REAL CHENRY(JPHR,NTRACT),CENTHAL(NTRACT),
     &     CKL(JPHR,NTRACT),ENHANC(JPHR,NTRACT),
     &     RWATER(JPHR,NTRACT)

C     LG- 20031118+ added some PBL parameters (see also PBLH.out)

      REAL ZPBLHMAX,ZPBLHOLD,ZPBLH_RESID,ZRGOLD
      INTEGER IRESET

C     LG- end

C======================================================================
C     DEFINE NAMELIST *NAMECH4* (ONLY SUITABLE FOR SINGLE-COLUMN)
C     ------ -------- -------- ----- -------- --- -------------- 
C     (OPPORTUNITY TO INITIALIZE LOCAL SWITCHES)
C     
      LOGICAL LINIT
      NAMELIST/NAMECH4/LSRFFLUX_OFF,LSRFFLUX_READ,LTS_OBS,LDISS_OFF
     +     ,CREFLUX,LAZ0FIX,FEPDU,LLSCALE_MOD

      LOGICAL LSRFFLUX_OFF,LSRFFLUX_READ,LTS_OBS,LDISS_OFF,LAZ0FIX
     +     ,LLSCALE_MOD

      DATA LINIT/.TRUE./
      DATA LSRFFLUX_OFF,LDISS_OFF/.FALSE.,.FALSE./
      DATA LSRFFLUX_READ,LTS_OBS/.FALSE.,.FALSE./
      DATA LAZ0FIX/.FALSE./
      DATA LLSCALE_MOD/.TRUE./
      DATA CREFLUX/0./
      DATA FEPDU/0.01/
C     

! mz_lg_20060427+
CGEERT 
      DATA ZENTR /0.3/
      DATA ZPRANDTL /1.00/
      DATA LLSCALE_MOD /.TRUE./
      DATA L_UPDATED   /.FALSE./
      DATA L_STABLE    /.FALSE./
      DATA LK_MFU /.FALSE./
CGEERT
! mz_lg_20060427-

C     LG- declaration of parameters used in the dry deposition scheme
C     as it has been incorporated in the ECHAM-3D model. The scheme calculates
C     the dry deposition velocities of O3, NO, NO2, HNO3, SO2 and SO4 according 
C     to the "big leaf" approach" (see Ganzeveld and Lelieveld, JGR, 1995 and 
C     1997 (in press, 23-07-1997). The deposition velocities are used finally in
C     the routine DRYDEP.f in PRECHEM.f to calculate the dry deposition flux!
C     For SO4 the aerosol dry deposition velocities is calculated using an
C     bimodal mass size distribution over sea and a unimodal distribution
C     over land. The resistances for SO4 are defined, however they are not 
C     relevant since the deposition velocity is controlled by other processes
C     The resistances are only defined to avoid numerical problems within the
C     loop of the tracer dry deposition velocity calculations.

*     C VDIFF
*     I VDIFF.136

C     ---- Parameterization of Vd, scheme Ganzeveld

C     Dimensioning output arrays
C     assigning array number to component different to sequence by 
C     Roelofs !!

C     -- Trace gases of the dry deposition scheme
C     O3, HNO3, NO, NO2, SO2, SO4

C     -- GOOSE, extension with the Wesely scheme, 1989 (Atmospheric Environ.)
C     in which the resistances of trace gases, which have not been observed
C     at all or not that often, are estimated from the reactivity relatively 
C     to ozone and the dissolvation relatively to SO2. The deposition process
C     of these two trace gases is relatively well known and used as a 
C     reference for the estimation.
C     
C     H2O2, CH3CHO, HCHO, CH3O2H, CH3C(O)O2H, HCOOH,
C     NH3, CH3COO2(NO2), HNO2 (see table 2, Wesely, 1989, Atmos. Environ.)

C     LG- end a extension to CO2, May 2001

C     -- GOOSE, deposition parameters
      REAL PSIM(NLON),PSIMZZ(NLON),PSIH(NLON),PSIH_OLD(NLON),
     *     PSIHZZ(NLON),RAHZZ(NLON),UMZZ(NLON),UMZZN(NLON),
     *     ZREF,ZRISURF,ZSURF,ZMONIN,ZOVERL,ZA,ZB,ZC,ZBLEND,ZRICRIT,
     *     KMKH,ZETA,CM(NLON),RBWAT(NLON,NTR),
     *     RCOX(NLON,NTR),RSVEG(NLON,NTR),VDLAND(NLON,NTR),VDWAT(NLON,NTR),
     *     RCO(NLON),CANHT(NLON),ZCDNVEG(NLON),
     *     ZCDNSLSN(NLON),VDSO4WAT(NLON),
     *     USTSLSN(NLON),CMVEG(NLON),CMSLSN(NLON),
     *     STHETA(NLON),VDSO4SLSN(NLON)

      REAL TSM(NLON),PGLACM(NLON)

C     LG- SOILPH is the soil pH 

      REAL SOILPH(96,48,5),DOC(64,32),AZ0MLOC,AZ0MSNOW,AZ0MSOIL,Z0MOLD

      SAVE SOILPH,DOC,Z0MOLD

C     LG- more declarations

*     I VDIFF.200
      REAL ZGJCONS(NLON),ZGJVEFF(NLON)

C     LG- setting some switches

      LMXL=.FALSE.
      LSENS_ZTMST=.FALSE.
      
C     LG- it turned out that the calculation of the vegetation temperature, 
C     as it has been included in the model (up to 08-2002) can only be done
C     for dry vegetation. To mimic this the interception capacity is set to
C     an extremly large value

      IF (LTVEG) CWLMAX=2.E20

C     LG- end     

C     LG- resetting of parameters TIME_SUNRISE and TIME_MIXING to zero,
C     to study the time between the sunrise and the onset of the mixing

      IF (NSTEP.EQ.NRESUM)  THEN
         TIME_SUNRISE=0
         TIME_MIXING=0

C     LG- and the old PBL height

         ZPBLHOLD=0.

      ENDIF 

C     =======================================================================
C     LG- calling of routine in which the surface resistances for use in the
C     dry deposition scheme are being calculated

      IF (NSTEP.EQ.0) CALL VD_RS(NSTEP) ! mz_lg_20060429+ multiple call of Vd_rs

C     LG- end
C     =======================================================================

      WRITE(NUNMDFL,*)'Start EC4_VDIFF.f'

C     LG- end

      IF (MLEV.LT.KLEV) STOP ' INCREASE JPNLEV IN PARAMV100 '
C     
C     READ NAMELIST *NAMECH4*
C     ---- -------- --------
C     
      IF (LINIT) THEN
         NUNECH=19
         NUNNAM= 7
         OPEN (UNIT=NUNECH ,FORM='FORMATTED',STATUS='OLD',FILE='namech4')
         READ  (NUNECH,NAMECH4)
         WRITE (NUNNAM,NAMECH4)
         CLOSE (NUNECH)
         WRITE (6,NAMECH4)
         IF (LSRFFLUX_OFF) LSRFFLUX_READ = .FALSE.
         IF (.NOT.LSRFFLUX_READ) LTS_OBS = .FALSE.
         DO JUN=1,2
            IF (JUN.EQ.1) IUNOUT=6
            IF (JUN.EQ.2) IUNOUT=NUNNAM
            IF (LSRFFLUX_OFF) THEN
               WRITE (IUNOUT,'(A)') 
     +              ' *** SURFACE HEAT FLUXES IN *EC_VDIFF* TURNED OFF BY HAND *'
            END IF
            IF (LSRFFLUX_READ) THEN
               IF (LTS_OBS) THEN
                  WRITE (IUNOUT,'(A)') 
     +                 ' *** OBSERVED SURFACE TEMPERATURES USED ***'
               ELSE
                  WRITE (IUNOUT,'(A)') 
     +                 ' *** SURFACE TEMPERATURES CALCULATED ***'
               END IF
               IF (CREFLUX.NE.0.0) THEN
                  WRITE (IUNOUT,'(A,F12.5,A)') 
     +                 ' *** FRACTION ',CREFLUX,' MOVED FROM LATF TO SENF ***'
               END IF
            END IF
            IF (LDISS_OFF) THEN
               WRITE (IUNOUT,'(A)') 
     +              ' *** DISSIPATION IN *EC_VDIFF* TURNED OFF BY HAND ***'
            END IF
            IF (LAZ0FIX) THEN
               WRITE (IUNOUT,'(A)')
     +              ' *** SEA SURFACE ROUGHNESS LENGTH FIXED *** '
            END IF
            IF (FEPDU.NE.0.01) THEN
               WRITE (IUNOUT,'(A,E8.2,A)')
     +              ' *** MINIMUM VALUE FOR SHEAR SET TO ',FEPDU,' M/S ***' 
            END IF
            IF (LLSCALE_MOD) THEN
               WRITE (IUNOUT,'(2A)')
     +              ' *** MODIFIED LENGTH SCALE FORMULATION ACCORDING '
     +              ,' TO THE 1995-PAPER BY BRINKOP AND ROECKNER USED '
            ENDIF

C     LG- 
            IF (LXTVDIFF) THEN
               WRITE (IUNOUT,'(A,E8.2,A)')
     +              ' *** TRACER VERTICAL TRANSPORT ***' 
            END IF

C     LG- end

         ENDDO
         LINIT = .FALSE.
      END IF
C     ------------------------------------------------
C     
C     *    PHYSICAL CONSTANTS.
C     -------- ----------
C     
C     *ZLAM* IS THE ASYMPTOTIC MIXING LENGTH FOR MOMENTUM EXCHANGE,
C     *ZKAP* IS THE VON KARMAN CONSTANT, *ZB*, *ZC* AND *ZD* ARE SOME
C     CONSTANTS FOR THE FORMULAE ABOUT STABILITY DEPENDENCY RESPECTIVELY
C     NEAR THE NEUTRAL CASE, IN THE UNSTABLE CASE AND IN THE STABLE
C     CASE AND *ZCHAR* IS THE CONSTANT OF THE *CHARNOCK FORMULA.
C     *ZQWSSAT* AND *ZQSNCR* ARE THE INVERSES OF CRITICAL VALUES FOR
C     SOIL WATER AND SNOW DEPTH THAT ARE USED IN THE COMPUTATION OF THE
C     EVAPOTRANSPIRATION'S EFFICIENCY.
C     
C     
C     
      ZLAM=CLAM
      ZKAP=CKAP
      ZB=CB
      ZC=CC
      ZD=CD

C     LG- it turned out that using ZD=5 yields a critical Richardson number of
C     of 0.13 which is too low compared top the value which has generally been
C     used, 0.25. Therefore ZD has been reset to 8./3. so that Ricrit becomes
C     0.25 (see ECHAM3 manual, 3.3.3.17)

      ZD=8./3.

!      IF (LFRICHSTAB_ECMWF) ZD = 1.  ! mz_lg_20060427+

      ZCHAR=CCHAR
      ZVA=CVA
      ZVB=CVB
      ZVC=CVC

      ZVBC=CVBC
      ZVK=CVK
      ZVKC=CVKC
      ZVABC=CVABC
      ZVRAD=CVRAD
      ZWLMAX=CWLMAX
      ZQSNCR=CQSNCR
      ZDA1=15.0
      ZUSTF=3.75

! mz_lg_20060427+
      ZDA1=ZUSTF**2                  ! dissipation length scale
      ZNEUTR = 1./(ZUSTF**0.5)       ! scaling neutral lengthscale at surf
! mz_lg_20060427-

      ZWSTF=0.2
      ZTKEMIN=1.E-4
      ZTMELT=TMELT
      ZRVRD=VTMPC1+1.
      ZRDRV=1./ZRVRD

C     LG- for tropical forest, see paper by Sellers et al., 1989

      IF (LSIB89) THEN

         IF (NSTEP.EQ.0) THEN
            print *,'Modified parameters for calculations Rstom (Sellers et 1989)'
            print *,'The parameters are representative for tropical rainforest!'
            print *,'ENTER TO CONTINUE'
            read (*,*)
         ENDIF

C     LG-   07-2002, in case of tropical forest, some of the initial parameters
C     are modified and an alternative method to calculate the exchange is
C     is applied, Reinder Ronda, personal communications, July 2002

         IF (OLSONVEG.EQ.33.OR.OLSONVEG.EQ.29) THEN 
            ZVA=2335.9
            ZVB=0.0145
            ZVC=153.5
            ZVBC=ZVB*ZVC
            ZVKC=ZVK*ZVC
            ZVABC=(ZVA+ZVBC)/ZVC
         ENDIF

      ENDIF 

C     LG- end
C     
      ZCPD=CPD
      ZRD=RD
      ZKAPPA=ZRD/ZCPD
      ZC3LES=C3LES
      ZC3IES=C3IES
      ZC4LES=C4LES
      ZC4IES=C4IES
C     
C     *      PARAMETERS FOR BOUNDARY LAYER DIAGNOSTICS
C     ---------- --- -------- ----- -----------
C     
      ZHUV=10.*G
      ZHTQ=2.*G
      ZEPHUM=5.E-2
      ZRHOS=RHOH2O*1.025
C     
C     
C     *    SECURITY PARAMETERS.
C     --------------------
C     
C     ZEPDU2 IS A MINIMUM SQUARED WIND INCREMENT TO AVOID DIVIDING BY
C     ZERO IN THE *RICHARDSON NUMBER'S CALCULATION AND ZEPZZO IS A
C     MINIMUM ROUGHNESS LENGTH.
C     
C     ZEPDU2=0.1
      ZEPDU2=FEPDU**2
      ZEPZZO=1.5E-05
      ZEPZ0O=2.
      ZEPCOR=5.E-05
C     
C     ZEPSW IS THE MINIMUM RELATIVE HUMIDITY OF THE GROUND,
C     ZEPSR IS A MINIMUM VALUE FOR THE RADIATION IN THE
C     VISIBLE PART OF THE SPECTRUM USED TO COMPUTE THE
C     CANOPY RESISTANCE.
C     
      ZEPSW=1.E-3
      ZEPSR=1.E-10
C     
C     ZEPEVAP IS THE MINIMUM ATMOSPHERIC DEMAND
C     
      ZEPEVAP=1.E-10
C     
C     
C     ZEPSEC IS A MINIMUM VALUE FOR THE DRAG COEFFICIENT
C     
      ZEPSEC=1.E-2
C     
C     *    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C     
      ZTMST=TWODT

C     LG- determining of half the timestep

      DTIME=0.5*TWODT
      
C     LG- end

      IF (NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZDIAGT=CONACC*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
C     
      ZTPFAC1=CVDIFTS
      ZTPFAC2=1./ZTPFAC1
      ZTPFAC3=1.-ZTPFAC2
      ZTPFAC4=1.+ZTPFAC3
C     
      ZZZLAM=30.
      ZCONS2=0.5*ZKAP/G
      ZCONS3=ZLAM
      ZCONS5=3.*ZB*ZC*G**2
      ZCONS6=1./3.
      ZCONS8=2.*ZB
      ZCONS9=3.*ZB
      ZCONS10=1./CPD
      ZCONS11=3.*ZB*ZC
      ZCONS12=ZTPFAC1*ZTMST*G/RD
      ZCONS13=1./ZTMST
      ZCONS14=ZCHAR*RD/(G**2*ZTMST)
      ZCONS15=1./(G*ZTMST)
      ZCONS16=CPD*VTMPC2
      ZCONS18=ZTPFAC1*ZTMST*G**2
      ZCONS17=1./ZKAP**2
      ZCONS25=ZCONS2/ZCONS3
C     
      ZPLMAX=0.75
      ZPLMIN=0.35
C     
      ZCHNEU=.3
C     
C     CONSTANT FOR FREE CONVECTION
C     EXPONENT FOR THE INTERPOLATION BETWEEN FREE CONVECTION
C     AND NEUTRAL OVER SEA
C     
      ZFREEC=0.0016
      ZGAM=1.25
      Z1DGAM=1./ZGAM
C     
C     
C     NEUTRAL STABILITY FUNCTIONS (MELLOR/YAMADA, 1982)
C     -------------------------------------------------------
C     
      ZH1= 2.22
      ZH2= 0.22
      ZM1= 1.24
      ZM2= 2.37
      ZM4= 3.69
      ZSHN=ZH1*ZH2*SQRT(2.)
      ZSMN=ZSHN*ZM1*ZM2/ZM4
C     
      ITOP=1
      ITOPP1=ITOP+1


! mz_lg_20060427+

      IF (LLSCALE_MOD) THEN

!       GEERT 
!       DO SOME TIME FILTERING TO FILTER OUT 2DT MODE LEAP FROG SCHEME
!
        ZTIMFIL = MIN (1. , ZTMST / 600.)        ! RELAXATION IN ABOUT 10 MINUTES
        DO JK=KTDIA+1,KLEV-1
          DO JL=KIDIA,KFDIA 
            ZPMFU(JL,JK) = ZTIMFIL * PMFU (JL,JK)  + (1 - ZTIMFIL)*ZPMFU(JL,JK)
          ENDDO
        ENDDO 

      ENDIF
! mz_lg_20060427-

C     
C     
C     ------------------------------------------------------------------
C     
C     *         1.     LOCATE AND POSITION SPACE.
C     ------ --- -------- ------
C     
 100  CONTINUE

C     LG- 

      IF (LCHEM) THEN

C     LG- Originally a file with the LAI data and land cover mask
C     is read (AREA) but this is not done anymore in this model
C     since for the LAI the NDVI data and the Olson database is used
C     and the landmask is also available at the highest resolution.
C     Also the file containing the local surface roughness fields
C     at ECHAM resolution are not read anymore.

C     Reading of input file with soil pH data,
C     SOILPH(JL,JR,1) - soil pH <5.5
C     SOILPH(JL,JR,2) - soil 5.5 <pH <7.3
C     SOILPH(JL,JR,3) - soil 7.3< pH <8.5
C     SOILPH(JL,JR,4) - soil 8.5 <pH
C     SOILPH(JL,JR,5) - soil 4 < pH <8.5

         IF (NSTEP.EQ.0) THEN

            DO 96 JT=1,5
               READ (125,699) ((SOILPH(JL,JR,JT),JL=1,96),JR=1,47,2),
     *              ((SOILPH(JL,JR,JT),JL=1,96),JR=48,2,-2)
 96         CONTINUE

 699        FORMAT (32f5.3/32f5.3/32F5.3)
            
            CLOSE(125)

C     LG-     reading in DOC(JL,JR) - dissolved organic matter in the oceans, 
C     used in the calculation of the ozone surface resistance

            DO 697 JR=1,31,2
               READ (126,*) (DOC(JL,JR),JL=1,64)
 697        CONTINUE                
            DO 698 JR=32,2,-2
               READ (126,*) (DOC(JL,JR),JL=1,64)
 698        CONTINUE
            
         ENDIF

C     LG-
C     05-2001, the call of the subroutine to calculate the vertical
C     profiles of radiation has originally being done in the 
C     chemistry routine PRECHEM1.f. However to calculate the vertical
C     profiles of the stomatal resistance according to the ECHAM
C     stomatal resistance and that calculated by the physiologial model
C     by Ronda et al. [2001], the radiation profiles are already 
C     required within EC4_VDIFF before the call of the PRECHEM1.f

C     Call subroutine to calculate the distribution of direct
C     and diffuse PAR within the canopy as a function of the 
C     sun angle, LAI and PAR at the top of the canopy. The fraction
C     of sunlit leaves (FSUN), shaded leaves (FSHADE), the flux density
C     of PAR on sunlit leaves (QSUN) and on shaded leaves (QSHADE) are 
C     determined. The subroutine is only called if LBIOSPH=.TRUE. in 
C     order to resolve the vertical profiles of radiation required for
C     the scaling of the photodissocation rates and the multilayer
C     dry deposition calculations, and if LVOCEMIS=.TRUE., in order to 
C     calculate the isoprene/monoterpene emission rates which are
C     sensitive for the vertical profiles of radiation. There is the 
C     opportunity to study the impact of changing the vertical resolution
C     of radiation profiles by putting the switch LRAD_1LAY to TRUE. for
C     which the subroutine CALCRAD.f is called (see PRECHEM1.f) 

         IF (LBIOSPH.OR.LVOCEMIS.OR.LRCO_ML)
     &        CALL CALCPROF(NSTEP,NSTOP)

C     LG-   end

      ENDIF

C     LG- assigning of the local surface roughness to the
C     parameter AZ0MLOC

      AZ0MLOC=MAX(1.E-3,Z0M)
      AZ0MSNOW=MAX(1.E-3,Z0MSNOW)
      AZ0MSOIL=MAX(1.E-3,Z0MSOIL)
      
C     LG- Calling of subroutine in which the surface roughness is being defined as 
C     a function of the wind direction. This done since measurements often show
C     a direction dependent roughness. To ensure a fair comparison this 
C     direction dependency can therefor optionally being introduced by setting
C     the switch LZ0MWDIR

      IF (LZ0MWDIR) THEN

         IF (NSTEP.EQ.0) THEN
            WRITE(*,'(1a)')
     &           ' The surface roughness is wind direction dependent since LZ0MWDIR=.TRUE.!'
            WRITE(*,'(1a)')
     &           ' For more details see the subroutine ASSIGN_Z0MWDIR'
            WRITE (*,'(1a)')' Enter to continue'
            READ (*,*)
         ENDIF

         ZDIR  = -90. - 180./API*ATAN2(PVM1(JPHR,KLEV),PUM1(JPHR,KLEV))
         IF (ZDIR.LE.0.) ZDIR = ZDIR + 360.
         CALL ASSIGN_Z0MWDIR(ZDIR,Z0M,Z0MWDIR)

C        LG- analysis of the simulated K values applying the nudging technique
C        including the winddirection dependent roughness that fast changes in
C        roughness result in large increases in the turbulence. It seems more
C        appropriate therefore to smoothen the roughness changes by also
C        applying some relaxation based on the new and old value   
         
         IF (NSTEP.EQ.0) PAZ0M(JL)=Z0M
	 AZ0MLOC = MAX(1.E-3,PAZ0M(JL) - (PAZ0M(JL) -Z0MWDIR) / (1800./ZTMST)) 
      ENDIF

C     LG- initialization of canopy temperature

      IF (NSTEP.EQ.0) TVEGM(JPHR)=PTSM(JPHR)

C     LG- end

C     ---------------------------------------------------------------
C     
C     *         2.     NEW THERMODYNAMIC VARIABLE AND BOUNDARY CONDITIONS.
C     --- ------------- -------- --- -------- -----------
C     
 200  CONTINUE
C     
C     *         2.04    DETERMINE PRESSURE THICKNESS ZDPH
C     
      CALL COPYRE(PAPM1(KIDIA,1) , ZDPH(KIDIA,1) , KLON)
      DO 204 JK=1,KLEVM1
         CALL DIFFRE(PAPM1(KIDIA,JK+1),PAPM1(KIDIA,JK)
     +        ,ZDPH(KIDIA,JK+1),KLON)
 204  CONTINUE
      CALL DIFFRE(PAPHM1(KIDIA,KLEVP1),PAPM1(KIDIA,KLEV)
     +     ,ZDPH(KIDIA,KLEVP1),KFDIA-KIDIA+1)

! mz_lg_20060427+
      IF (LLSCALE_MOD) THEN
!
!       GEERT 
!       IF L_UPDATED COMPUTE WITH UPDATE TENDENCY FOR PTM1 
!
        RFUZ = ZTPFAC1
        DO JK=KTDIA,KLEV
          DO JL=KIDIA,KFDIA
            PTM1PR(JL,JK)=PTM1(JL,JK)
            IF (L_UPDATED) 
     *        PTM1PR(JL,JK)=PTM1(JL,JK) + RFUZ*ZTMST*PTTE(JL,JK)
          ENDDO
        ENDDO

      ENDIF
! mz_lg_20060427+


C     
C     *         2.1     REPLACE T BY CP(Q)*T+GZ IN THE ATMOSPHERE.
C     
 210  CONTINUE
      DO 212 JK=KTDIA,KLEV
         DO 211 JL=KIDIA,KFDIA
            ZCPTGZ(JL,JK)=PGEOM1(JL,JK)+PTM1(JL,JK)*CPD
     *           *(1.+VTMPC2*PQM1(JL,JK))
            ZTETA1(JL,JK)=PTM1(JL,JK)*(100000./PAPM1(JL,JK))**ZKAPPA
            ZTVIR1(JL,JK)=ZTETA1(JL,JK)*(1.+VTMPC1*PQM1(JL,JK)-PXM1(JL,JK))
C     
            LO=PTM1(JL,JK).GE.ZTMELT
            ZFAXE(JL,JK)=CVMGT(ALV,ALS,LO)
            ZBET=ZFAXE(JL,JK)/ZCPD
            ZUSUS1=ZBET*ZTETA1(JL,JK)/PTM1(JL,JK)*PXM1(JL,JK)
            ZLTETA1(JL,JK)=ZTETA1(JL,JK)-ZUSUS1
C     
            IT=PTM1(JL,JK)*1000.
            ZES=TLUCUA(IT)/PAPM1(JL,JK)
            ZES=MIN(ZES,0.5)
            ZCOR=1./(1.-VTMPC1*ZES)
            ZQSS(JL,JK)=ZES*ZCOR
 211     CONTINUE
 212  CONTINUE
C     
      DO 214 JK=KTDIA,KLEVM1
         DO 215 JL=KIDIA,KFDIA
            ZHH(JL,JK)=(PGEOM1(JL,JK)-PGEOM1(JL,JK+1))
            ZSDEP1=(PAPHM1(JL,JK)-PAPHM1(JL,JK+1))/(PAPHM1(JL,JK)-
     1           PAPHM1(JL,JK+2))
            ZSDEP2=(PAPHM1(JL,JK+1)-PAPHM1(JL,JK+2))/(PAPHM1(JL,JK)-
     1           PAPHM1(JL,JK+2))
C     
            ZQSSM(JL,JK)=ZSDEP1*ZQSS(JL,JK)+ZSDEP2*ZQSS(JL,JK+1)
            ZTMITTE(JL,JK)=ZSDEP1*PTM1(JL,JK)+ZSDEP2*PTM1(JL,JK+1)
            ZTVIRMIT(JL,JK)=ZSDEP1*ZTVIR1(JL,JK)+ZSDEP2*ZTVIR1(JL,JK+1)
            ZFAXEN(JL,JK)=ZSDEP1*ZFAXE(JL,JK)+ZSDEP2*ZFAXE(JL,JK+1)
            ZLWCMIT(JL,JK)=ZSDEP1*PXM1(JL,JK)+ZSDEP2*PXM1(JL,JK+1)
            ZQMIT(JL,JK)=ZSDEP1*PQM1(JL,JK)+ZSDEP2*PQM1(JL,JK+1)
            ZTEMIT(JL,JK)=ZSDEP1*ZTETA1(JL,JK)+ZSDEP2*ZTETA1(JL,JK+1)
            ZCCOVER(JL,JK)=PACLCM(JL,JK)*ZSDEP1+PACLCM(JL,JK+1)*ZSDEP2
 215     CONTINUE
 214  CONTINUE

C     -----------------------------------------------------------------------------
C     LG- calculation of the within canopy exchange coefficients which
C     are applied for the calculation of the vertical turbulent exchange
C     of the tracers. There is the opportunity to perform  a complete mixing 
C     of the mass for a specific number of layers for a specific frequency 
C     which can be defined. This complete mixing is introduced in order to 
C     simulate the effect of gusts on the net exchange fluxes between the 
C     canopy and the surface layer.
C     -----------------------------------------------------------------------------

C     LG- definition of pressure levels within the canopy. The pressure
C     at the HALF!! pressure level of the surface layer and the surface 
C     pressure has been used and scaling is performed using the difference
C     in height to perform the extrapolation (dP=rho*G*dZ, and rho=p/RT
C     some manipulation yields dP=P*G/RT*dZ 

      DO JK=1,KLEV
         DO JL=1,NLON

C     LG-  determining the height of the half pressure levels from 
C     the thickness of the layer and the reference height of the full 
C     pressure level

            ZDP(JL,JK)=PAPHM1(JL,JK+1)-PAPHM1(JL,JK)
            ZT(JL,JK)=PTM1(JL,JK)+PTTE(JL,JK)*DTIME
            ZQ(JL,JK)=PQM1(JL,JK)+PQTE(JL,JK)*DTIME
            TVIRT=ZT(JL,JK)*(1.+0.607717*ZQ(JL,JK))
            ZRHOA(JL,JK)=PAPM1(JL,JK)*ZMAIR*1E-6/(TVIRT*ZGASC)

C     LG-  Calculation of depth of layers and altitude,
C     the canopy height is accounted for by increasing the
C     model's altitude with the displacement height, which is about
C     2/3 the canopy height. 

            ZZDZ(JL,JK)=ZDP(JL,JK)/(ZRHOA(JL,JK)*1.E3*G)
            Z(JL,JK)=PGEOM1(JL,JK)/G+HC
            ZHP(JL,JK)=Z(JL,JK)+0.5*ZZDZ(JL,JK)

            ZAPHM1(JL,JK)=PAPHM1(JL,JK)
            ZAPHM1(JL,JK+1)=PAPHM1(JL,JK+1)
            ZAPM1(JL,JK)=PAPM1(JL,JK)

         ENDDO
      ENDDO

C     LG- defining different physical parameters for the vegetation/snow layers,
C     e.g. pressure, mass, volume, density.

      IF (LBIOSPH.OR.LBULKVEG.OR.LVEG_MLAY.OR.LSNOW_MLAY) THEN

         IF (LBIOSPH) THEN
            KLEVELS=NLEVT
         ELSEIF (LBULKVEG) THEN
            KLEVELS=KLEV+1
         ELSEIF (LVEG_MLAY) THEN
            KLEVELS=KLEV+NLEVV_ML
         ELSEIF (LSNOW_MLAY) THEN
            KLEVELS=KLEV+NLEVS_ML
         ENDIF

         DO JK=KLEV+1,KLEVELS
            DO JL=1,NLON  

C     LG-   for lagragian simulations over land the reference height is 
C     corrected for the actual number of layers as a function of the
C     canopy height. This done by correction the pressure difference
C     for the half and full pressure levels (DAPM1 and DAPHM1) with 
C     the reduced number of vegetation layers and the specific height
C     of each layer expressed in its pressure difference. The corrected
C     pressure differences DP and DHP are finally used to determine
C     the canopy pressure levels and reference heights
               
               IF (LAGRIAN.AND.ILSTYPE.GT.1.AND.
     &              JK-((KLEVELS-KLEV)-NLEVVEG).LE.KLEV) THEN
                  JJK=KLEV
               ELSE
                  JJK=JK
               ENDIF

               IF (JK.EQ.KLEV+1) THEN
                  JJ=0
                  DDP(JK-1)=0.
                  DDP_Z(JK)=0.
                  DDP_ZHP(JK)=0.
               ENDIF

               ZT(JL,JK)=ZT(JL,KLEV)
               ZQ(JL,JK)=ZQ(JL,KLEV)

               IF (JJK.EQ.KLEV) THEN
                  JJ=JJ+1
                  DDP(JK)=DDP(JK-1)+(DAPHM1(JL,JK+1)-DAPHM1(JL,JK))
                  DP(JK)=0.
                  DDP_Z(JK)=0.
                  DHP(JK)=0.
                  DHP(JK+1)=0.
                  DDP_ZHP(JK)=0.
               ELSE

C     LG-    in order to ensure that the pressure correction for 
C     a reduced number of actual canopy layers becomes larger than
C     the pressure difference between the surface and the surface 
C     layer for a canopy height less than the thickness of the 
C     layers, these canopy is incorporated in this pressure difference
C     DDP

                  DDP(JK)=DDP(JK-1)
                  DHP(JK)=DAPHM1(JL,JK)-DDP(JK)
                  DHP(JK+1)=DAPHM1(JL,JK+1)-DDP(JK)
                  DP(JK)=DAPM1(JL,JK)-DDP(JK)

C     LG-    first estimate of height at half pressure level

                  ZHP(JL,JK)=ZHP(JL,KLEV)-
     *                 DHP(JK)/(PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
                  ZHP(JL,JK+1)=ZHP(JL,KLEV)-
     *                 DHP(JK+1)/(PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))

C     LG-    determining the correction of the term DDP in order to
C     prevent the calculation of a negative height

                  IF (ZHP(JL,JK+1).LT.0.) THEN
                     DDP_ZHP(JK+1)=ZHP(JL,JK+1)*(PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
                     IF (.NOT.LBULKVEG) THEN
                        DDP_Z(JK)=DDP_ZHP(JK+1)
                        DDP_ZHP(JK)=DDP_ZHP(JK+1)
                     ELSEIF (LBULKVEG) THEN
                        DDP_Z(JK)=DDP_ZHP(JK+1)/2.
                     ENDIF
                  ENDIF

                  DDP(JK)=DDP(JK-1)
                  DP(JK)=DAPM1(JL,JK)-DDP(JK)
                  DHP(JK)=DAPHM1(JL,JK)-DDP(JK)
                  DHP(JK+1)=DAPHM1(JL,JK+1)-DDP(JK)
               ENDIF

            ENDDO
         ENDDO

         DO JK=KLEV+1,KLEVELS
            DO JL=1,NLON

C     LG-   pressure and reference height at full pressure level

               ZAPM1(JL,JK)=PAPM1(JL,KLEV)+DP(JK)+DDP_ZHP(JK)
               Z(JL,JK)=Z(JL,KLEV)-(DP(JK)+DDP_Z(JK))/
     *              (PAPM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
               HGHT(JK)=Z(JL,JK)

C     LG-   pressure and reference height at half pressure level, for both
C     the upper and lower boundary of the layer

               ZAPHM1(JL,JK)=PAPHM1(JL,KLEV)+DHP(JK)+DDP_ZHP(JK)
               ZAPHM1(JL,JK+1)=PAPHM1(JL,KLEV)+DHP(JK+1)+DDP_ZHP(JK+1)
               ZHP(JL,JK)=ZHP(JL,KLEV)-
     *              (DHP(JK)+DDP_ZHP(JK))/(PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
               ZHP(JL,JK+1)=ZHP(JL,KLEV)-(DHP(JK+1)+DDP_ZHP(JK+1))/
     *              (PAPHM1(JL,KLEV)*G/(RD*ZT(JL,JK)))
               HGHTHP(JK)=ZHP(JL,JK)
               HGHTHP(JK+1)=ZHP(JL,JK+1)

C     LG-   density

               TVIRT=ZT(JL,JK)*(1.+0.607717*ZQ(JL,JK))
               ZRHOA(JL,JK)=ZAPM1(JL,JK)*ZMAIR*1E-6/(TVIRT*ZGASC)

            ENDDO
         ENDDO

C     LG- end IF LBIOSPH

      ENDIF
C     
C     COMPUTE FRACTIONAL SURFACE COVERAGES
C     
      DO 213 JL=KIDIA,KFDIA
         PCVS(JL)=MIN(1.,PSNM1M(JL)*ZQSNCR)
         
C     LG- end

C         zwlmax=1e10 ! mz_lg_20060225+ to simulate zero wet skin fraction

         PWLMX(JL)=ZWLMAX*((1.-PVGRAT(JL))+PVGRAT(JL)*PVLTM(JL))
         ZWLMXI(JL)=1./PWLMX(JL)
         PCVW(JL)=PWLM1M(JL)*ZWLMXI(JL)

C     LG- writing the fractions of surface cover to an output file and
C     assigning these fractions to the parameters put in COMMON block
C     INP (see comchem.h)

         VEGFRAC(JL)=(1.-PCVS(JL))*(1.-PCVW(JL))*PVGRAT(JL)
         SNFRAC(JL)=PCVS(JL)
         WSFRAC(JL)=(1.-PCVS(JL))*PCVW(JL)
         BSFRAC(JL)=0.
         IF (LALAND(JL))
     &        BSFRAC(JL)=(1.-PCVS(JL))*(1.-PVGRAT(JL))*(1.-PCVW(JL))

C     LG- determining the water/ocean surface cover being the residual

         WTFRAC(JL)=(1.-(VEGFRAC(JL)+WSFRAC(JL)+SNFRAC(JL)+BSFRAC(JL)))

C     LG- and assigning the forest fraction

         FSTFRAC(JL)=PFORESTM(JL)

         IF (NSTEP.EQ.0) THEN 
            OPEN (UNIT=NUNFRAC,FILE='/data/ganzevl/racmo/output/surffrac.out',
     *           STATUS='UNKNOWN')
            WRITE(NUNFRAC,'(1a)')'Fractions of surface cover and some add. param.'
            WRITE(NUNFRAC,'(a10,a15,9a12)')
     *           'step','time','snow','wet-skin','veg.','soil','water',
     *           'forest','ws+veg','wlmx','wlm1m'
         ENDIF
         
         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0)
     *        WRITE(NUNFRAC,'(1x,i9.9,1x,a14,7f12.2,7e12.5)')
     *        NSTEP,LDATLTIME,SNFRAC(JL),WSFRAC(JL),VEGFRAC(JL),
     *        BSFRAC(JL),WTFRAC(JL),FSTFRAC(JL),WSFRAC(JL)+VEGFRAC(JL),
     *        PWLMX(JL),PWLM1M(JL)

         IF (NSTEP.EQ.NSTOP) CLOSE(NUNFRAC)

C     LG- end

 213  CONTINUE

C     
C     *         2.2     SATURATION PARAMETERS,
C     *                 RELATIVE HUMIDITY OVER THE BARE LAND PART
C     *                 AND VIRTUAL TEMPERATURE AT THE SURFACE.
C     
 220  CONTINUE
C     
      DO 221 JL=KIDIA,KFDIA
         LO=(PTSM1M(JL)-TMELT).GT.0.

C     LG- replacing the surface temperature by the skin temperature

         IF (LTVEG)
     &        LO=(TVEGM(JL)-TMELT).GT.0.
         
         ZCVM3=CVMGT(C3LES,C3IES,LO)
         ZCVM4=CVMGT(C4LES,C4IES,LO)
         ZCVM5=CVMGT(C5LES,C5IES,LO)
         IT=PTSM1M(JL)*1000.

C     LG- replacing the surface temperature by the skin temperature

         IF (LTVEG)
     &        IT=(1.-(VEGFRAC(JL)+WSFRAC(JL)))*PTSM1M(JL)*1000.+
     &        (VEGFRAC(JL)+WSFRAC(JL))*TVEGM(JL)*1000.

         ZES=TLUCUA(IT)/PAPHM1(JL,KLEVP1)
         ZCOR=1./(1.-VTMPC1*ZES)
         ZQS(JL)=ZES*ZCOR
         ZDQS(JL)=ZQS(JL)*ZCVM5*ZCOR*(1./(PTSM1M(JL)-ZCVM4))**2

C     LG- replacing the surface temperature by the skin temperature

         IF (LTVEG)
     &        ZDQS(JL)=(1.-(VEGFRAC(JL)+WSFRAC(JL)))*
     &        ZQS(JL)*ZCVM5*ZCOR*(1./(PTSM1M(JL)-ZCVM4))**2+
     &        (VEGFRAC(JL)+WSFRAC(JL))*
     &        ZQS(JL)*ZCVM5*ZCOR*(1./(TVEGM(JL)-ZCVM4))**2
         
         PWSM1M(JL)=MIN(PWSM1M(JL),PWSMXM(JL))

C     LG- assigning the value of the old soil moisture level PWSM1, to the 
C     parameter SOILWS which is required in the calculations of the NO
C     soil biogenic emissions (see noxemis.f)

         SOILWS(JL)=PWSM1M(JL)

C     LG- end

         ZWSTOP=MIN(0.1,PWSMXM(JL))
         ZWSLEV=PWSMXM(JL)-ZWSTOP
         IF(PWSM1M(JL).GT.ZWSLEV) THEN
            ZHUM(JL)=0.5*(1.-COS((PWSM1M(JL)-ZWSLEV)*API/ZWSTOP))
         ELSE
            ZHUM(JL)=0.
         ENDIF
         ZHSOIL(JL)=PCVS(JL)+(1.-PCVS(JL))*(PCVW(JL)+(1.-PCVW(JL))
     >        *ZHUM(JL))
         ZHSOIL(JL)=CVMGT(ZHSOIL(JL),1.,LALAND(JL))

C     LG- update 03-2000, GL scheme

C     GL.  CORRECTION IN MOIST CONDITIONS

         IF (LLSCALE_MOD) THEN
            ZHSOIL(JL)= ZHSOIL(JL) + PACLCM(JL,KLEV)*(1.-ZHSOIL(JL))
         ENDIF

         LO=PQM1(JL,KLEV).GT.ZQS(JL)
         ZHSOIL(JL)=CVMGT(1.,ZHSOIL(JL),LO)
         ZTESS(JL)=PTSM1M(JL)*(1.E5/PAPHM1(JL,KLEVP1))**ZKAPPA

C     LG- replacing the surface temperature by the skin temperature

         IF (LTVEG)
     &        ZTESS(JL)=((1.-(VEGFRAC(JL)+WSFRAC(JL)))*PTSM1M(JL)+
     &        (VEGFRAC(JL)+WSFRAC(JL))*TVEGM(JL))*
     &        (1.E5/PAPHM1(JL,KLEVP1))**ZKAPPA

C     LG- end

         ZTVS(JL)=ZTESS(JL)*(1.+VTMPC1*ZHSOIL(JL)*ZQS(JL))

C     GL.  NO MOIST CORRECTION STABILITY
         IF (LLSCALE_MOD) THEN
            ZTVS(JL)=ZTESS(JL)*(1.+VTMPC1*PQM1(JL,KLEV))
         ENDIF

 221  CONTINUE

C     
C     *         2.3     DEFINITION OF THE STOMATAL RESISTANCE
C     
 230  CONTINUE
C     
      DO 231 JL=KIDIA,KFDIA
         ZWCRIT=ZPLMAX*PWSMXM(JL)
         ZWPWP=ZPLMIN*PWSMXM(JL)
         ZQWEVAP=1./(ZWCRIT-ZWPWP)
         ZSOIL=MAX(ZEPSW,MIN(1.,(PWSM1M(JL)-ZWPWP)*ZQWEVAP))
         ZSRFL=MAX(ZEPSR,PSRFL(JL)*ZVRAD)
         ZABCS=(ZVA+ZVBC)/(ZVC*ZSRFL)
         ZVKLT=ZVK*PVLTM(JL)
         ZVXPKLT=EXP(ZVKLT)
         ZVXMKLT=EXP(-ZVKLT)
         ZLN1=LOG((ZABCS*ZVXPKLT+1.)/(ZABCS+1.))
         ZLN2=LOG((ZABCS+ZVXMKLT)/(ZABCS+1.))
         ZRSI=(ZVB*ZLN1/ZVABC-ZLN2)/ZVKC
         ZRS0(JL)=1./ZRSI       ! this is the resistance, uncorrected for water stress
         
C     LG- for net radiation of zero, ZRO(JL) is set to large value, since in the
C     the original Sellers scheme for RG=0, rstom is about 600 s m-1 (LAI=1)

         IF (ZSRFL.LE.1E-10) ZRS0(JL)=1.E5

C     LG- testing the impact of using a constant FWS factor for tropical forest
C     pasture to get a better agreement between the simulated and observed
C     surface energy balances

         IF (LFWS_RESET) THEN

           IF (NSTEP.EQ.0) THEN
             print *,
     &        'ec4_vdiff: The initial soil moisture stress value is: ec4_vdiff',ZSOIL
             print *,'ENTER to continue'
             read (*,*)
           ENDIF
    
           ZSOIL=FWS_RESET
           IF (NSTEP.EQ.0) 
     &        WRITE(*,'(1a,f4.1)')' EC_VDIFF; FWS SET TO CONSTANT VALUE OF: ',
     &          FWS_RESET 
         
         ENDIF
C        LG- end 
         
         ZWET(JL)=ZRS0(JL)/ZSOIL

C     LG- 
C     Stomatal resistance computation
*     I VDIFF.467
         FWS(JL)=ZSOIL
C     recalculated rstom for a LAI of 1 (new ln1 and ln2 calculation)
         RCO(JL)=1./((ZVB*(ALOG((ZABCS*2.71**0.9+1.)/(ZABCS+1.)))
     *        /ZVABC-(ALOG((ZABCS+2.71**(-0.9))/(ZABCS+1.))))/ZVKC)

C     LG- for net radiation of zero, ZRO(JL) is set to large value, since in the
C     the original Sellers scheme for RG=0, rstom is about 600 s m-1 (LAI=1)

         IF (ZSRFL.LE.1E-10) RCO(JL)=1.E5

C     LG- end

         IF (LCHEM) THEN

C     GOOSE, calculation of mesophyll resistance of NO/NO2 from stomatal resistance O3
            RCOX(JL,io3)=DIFF(io3)*RCO(JL)

C     LG-    in the "big leaf" dry deposition description NO and NO2 dry deposition
C     are being considered with for NO2 and mesophyllic resistance that yields
C     a stomatal conductance being about 2/3 that of ozone. The NO stomatal
C     conductance is very small. These approach is based on observations for
C     which it is not quite clear if a sink in the form of stomatal uptake is
C     present or that a chemical sink is responsible for the removal. 

            IF (LBULKVEG.OR.LVEG_MLAY) THEN
               RMES(ino_pr)=1.E5 !5.*RCOX(JL,io3)
               RMES(ino2_pr)=RMES(io3) !0.5*RCOX(JL,io3)
            ELSE
               RMES(ino_pr)=5.*RCOX(JL,io3)
               RMES(ino2_pr)=0.5*RCOX(JL,io3)
            ENDIF

         ENDIF

C     LG- end

         LO=PQM1(JL,KLEV).GT.ZQS(JL)
         ZWET(JL)=CVMGT(0.,ZWET(JL),LO)

C     LG- assigning the canopy stomatal resistance for comparison with the
C     physiological stomatal resistance

         PWET(JL)=ZWET(JL)

 231  CONTINUE

C     =========================================================================
C     LG- introduction of physiological model for stomatal exchange by Ronda
C     et al., 2001, the model needs the definition of the CO2 concentration
C     and is therefore only called whenever the switch LCHEM=true

      IF (LAGS) THEN

C     LG-   some additional parameters need to be calculated for use within the
C     subroutine. Originally the subroutine is implemented in the surface
C     scheme of the ECMWF model, and consequently some of the calculations
C     of the routine EC4VB_VDIFF. 

         DO 3313 JL=KIDIA,KFDIA
            ZQSA(JL)=PQM1(JL,KLEV)*(1.-ZCAIR(JL))+ZCSAT(JL)*ZQS(JL)
            ZQSA(JL)=MAX(1.E-12,ZQSA(JL))

C     LG-     assigning the C3/C4 vegetation type

            PCATYPE(JL)=PC3C4TYPE

C     LG-     the scheme uses the soil water of all the distinguished soil layers.
C     However, in EC4_VDIFF this information is not available and therefore
C     the bucket soil water content of EC4_VDIFF is divided over the number
C     of soil layers and subsequently integrated in AGS.f and used within
C     the calculations of the stomatal conductance

            DO JK=1,KLEVS
               PWSAM1M(JL,JK)=PWSM1M(JL)/KLEVS
            ENDDO

C     LG-     recalculating the CO2 concentration to ppmv, for use within the 
C     routine AGS.f

            DO JK=1,KLEVEL
               IF (ZXTMCO2(JL,JK).LT.1E-20) THEN

C     LG- initialization of in-canopy CO2 concentrations for 
C     LBIOSPH=.TRUE.

                  ZXTMCO2(JL,JK)=PXTM1(JL,KLEV,ico2)/(1.E3*ZRHOA(JL,KLEV)/
     &                 (ZMAIR*1.E9/AVO)) 
               ELSE
                  ZXTMCO2(JL,JK)=PXTM1(JL,JK,ico2)/(1.E3*ZRHOA(JL,JK)/
     &                 (ZMAIR*1.E9/AVO))
               ENDIF
            ENDDO

C     LG-     for the bulk and multi-layer vegetation models

            IF (LBULKVEG.OR.LVEG_MLAY) THEN
               DO JK=1,NLEVV
                  ZXTMCO2(JL,JK+KLEV)=PXTMVEG(JL,JK,ico2)/
     &                 (1.E3*ZRHOA(JL,JK)/(ZMAIR*1.E9/AVO))
                  IF (ZXTMCO2(JL,JK+KLEV).EQ.0.) ZXTMCO2(JL,JK+KLEV)=350.E3
               ENDDO
            ENDIF

C     LG-     end

 3313    CONTINUE

         IF (LALAND(JPHR))
     *        CALL AGS(KIDIA,KFDIA,KLON,KLEV,KLEVS,
     *        PSOTYPE,PCATYPE,PQM1, ! PCATYPE = C3/C4 vegetation type
     *        PTSM1M, PWSAM1M,  ! C3, PCATYPE=1, C4, PCATYPE=2
     *        PSRFL , PVLTM, PAPHM1,

C     LG-   added the permanent wilting point and some other soil water parameters
C     relevant to the calculation of the soil water stress effect, moreover
C     the CO2 concentration is added (in ppmv)

     *        ZWPWP,  ZQWEVAP, ZXTMCO2,

C     LG-   end

     *        ZQS   , ZQSA , PRS0_AGS  ,PWET_AGS, PRS0_AGSML)	

C     LG-   end

      ENDIF

C     ==========================================================================

C     LG- calculation of multi-layer ECHAM and AGS stomatal resistance for use 
C     in subroutine VEG_MLAY.f, the interpolation from the NLEVV layers to 
C     the number of layers actually being used in the VEG_MLAY.f is 
C     done within that subroutine

C LG- 092004- testing of the next parameterization has revealed that it 
C     initially (without modification of adding (1-FSL) term in front of
C     RVD(JK) term), resulted in much smaller deposition velocities through
C     providing much larger leaf stomatal uptake resistances. With the
C     modification there is already a significant increase in uptake but it
C     still smaller compared to using the SiB86 bulk parameterization. 


      IF (LCHEM.AND.LBIOSPH.OR.LBULKVEG.OR.LVEG_MLAY) THEN
         DO JL=KIDIA,KFDIA
            DO JK=1,NLEVV
               DO  JT=1,NTR

C                 LG- for net radiation of zero, ZRO(JL) is set to large value, since in the
C                 the original Sellers scheme for RG=0, rstom is about 600 s m-1 (LAI=1)

                  IF (ZSRFL.GT.1E-10) THEN
                    ZABCS=(ZVA+ZVBC)/(ZVC*RBVD*FSL(JK)+(1-FSL(JK))*RVD(JK))
                    RCOX_ECHML(JK,JT)=DIFF(JT)/
     *                 ((ZVB*(ALOG((ZABCS*2.71**0.9+1.)/(ZABCS+1.)))
     *                 /ZVABC-(ALOG((ZABCS+2.71**(-0.9))/(ZABCS+1.))))/ZVKC)
                    RCOX_ECHML(JK,JT)=RCOX_ECHML(JK,JT)/FWS(JL)
                  ELSE 
		    RCOX_ECHML(JK,JT)=1.E5
                  ENDIF

                  ! mz_lg_20060225+ modified for lagrangian simulations
                  IF (LAGS.AND.LAI.GT.1.E-10) THEN
		    RCOX_AGSML(JK,JT)=DIFF(JT)*PRS0_AGSML(JL,JK)
		  ELSE
		    RCOX_AGSML(JK,JT)=1.E5
		  ENDIF
                  ! mz_lg_20060225-
               ENDDO
            ENDDO
	 ENDDO
      ENDIF

C     LG- end

C***  
      IF (LVDIFF) THEN
C***  
C     
C     ------------------------------------------------------------------
C     
C     *         3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
C     ----------- -- --- -------- -------------
C     
C     THE SURFACE LAYER IS NOW COMPUTED BEFORE THE OTHER LEVELS
C     
 300     CONTINUE
C     
C     *       3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
C     *                 RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
C     *                 AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
C     *                 COMMON PART OF THE DRAG COEFFICIENTS.
C     
 310     CONTINUE
C     

C     LG- in case of LXTVDIFF=.TRUE., and LBIOSPH=.TRUE., 
C     calculation of some within canopy exchange parameters. This is
C     already done here in order to make it possible to correct the
C     surface layer windspeed with some value determined using a 
C     function which mimics the gustiness of the exchange. In order to
C     have an optimal consistency in all parameters, the surface layer
C     windspeed is corrected.

         IF (LXTVDIFF) THEN
C     
C     --  Definition of the frequence of the gusts which is used
C     for a total mixing of mass over a specific domain in the vertical
C     this is only possible for intergrations with timestep less than 
C     five minutes, which is the typical frequency of the occurence of
C     these gusts
C     
            IF (NSTEP.EQ.0 ) THEN
               LGUST=.FALSE.
               IF (LBIOSPH.AND.ZTMST.LE.300.) THEN
                  N=0
                  LGUST=.TRUE.
                  GUST=1
                  PRINT *,'Frequency of occurence of gusts [seconds] ?'
                  PRINT *,'(for not considering gusts enter value > 1e3)'
                  READ *,FREQ
                  IF (FREQ.GT.1.E3) LGUST=.FALSE.
               ELSE
                  FREQ=1.E20
               ENDIF
            ENDIF

C     LG-  determining the shape of the function being used to mimic gustiness
C     (see page 228 of Crop Micrometeorology by Goudriaan), The model 
C     resolved average surface layer windspeed is adjusted using the 
C     function value. The corrected windspeed is then applied to determine
C     the vertical profiles of the tracer exchange coefficients within
C     the canopy. The values of the function applied by Goudriaan are 
C     derived from the figure normalising the time by the period of the 
C     gust (thus 100 seconds) and using a polynomial fit to reproduce the
C     increasing and decreasing branch

            IF (LGUST) THEN
               NSTEP_GUST=INT(FREQ/DTIME)
               X=N*DTIME/FREQ
               IF (N.LT.0.5*NSTEP_GUST) THEN
                  GUST=0.25242-0.29707*X+8.0591*X**2
               ELSE
                  GUST=9.1071-17.642*X+8.7946*X**2
               ENDIF       
               IF (N.EQ.NSTEP_GUST) THEN 
                  N=1
               ELSE
                  N=N+1
               ENDIF
            ENDIF

C     LG-  definition of layer being the upper boundary of the domain
C     where the impact of gusts are considered. KPBLHE is the top of 
C     the PBL, be careful with definition of KLEV/NLEV!!!
            
            NLEVMIX=MIN(KPBLHE(JL),KPBLHE(JL))

C     LG-  only for the daytime, gusts are being considered

            IF (LGUST.AND.RG.GT.0.) THEN
               WRITE(NUNMDFL,*)
     *              'Considering gustiness from the level: ',NLEVMIX,' up to ',KLEVEL 
            ELSE
               GUST=1.
            ENDIF

         ENDIF

         DO 311 JL=KIDIA,KFDIA

C     LG- modifying the surface layer windspeed for LGUST=.TRUE.

            IF (LGUST) THEN
               PUM1(JL,KLEV)=GUST*PUM1(JL,KLEV)
               PVM1(JL,KLEV)=GUST*PVM1(JL,KLEV)
            ENDIF

C     LG- endif

            ZDU2=MAX(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)
            ZQMITTE=(PQM1(JL,KLEV)+ZQS(JL)*ZHSOIL(JL))/2.
            ZQTMIT=PXM1(JL,KLEV)*0.5+ZQMITTE
            ZTMIT=(PTM1(JL,KLEV)+PTSM1M(JL))/2.

C     LG- replacing the surface temperature by the skin temperature

            IF (LTVEG)
     &           ZTMIT=(PTM1(JL,KLEV)+
     &           ((1.-(VEGFRAC(JL)+WSFRAC(JL)))*PTSM1M(JL)+
     &           (VEGFRAC(JL)+WSFRAC(JL))*TVEGM(JL)))/2.

C     LG- end

            ZTEMITTE=(ZTETA1(JL,KLEV)+ZTESS(JL))/2.
            ZVIRMITTE=(ZTVIR1(JL,KLEV)+ZTVS(JL))/2.
            ZQSMIT=(ZQSS(JL,KLEV)+ZQS(JL))/2.
            ZQLWI1=PQM1(JL,KLEV)+PXM1(JL,KLEV)
            ZQLWI2=ZQS(JL)*ZHSOIL(JL)

C           LG- update 03-2000, GL scheme

            IF (LLSCALE_MOD) ZQLWI2 = ZQLWI1

C           LG- end

            ZFUX=ZFAXE(JL,KLEV)/(ZCPD*ZTMIT)
            ZFOX=ZFAXE(JL,KLEV)/(ZRD*ZTMIT)
            ZMULT1=1.+VTMPC1*ZQTMIT
            ZMULT2=ZFUX*ZMULT1-ZRVRD
            ZMULT3=ZRDRV*ZFOX*ZQSMIT/(1.+ZRDRV*ZFOX*ZFUX*ZQSMIT)
            ZMULT5=ZMULT1-ZMULT2*ZMULT3
            ZMULT4=ZFUX*ZMULT5-1.
C     
            ZDUS1=PACLCM(JL,KLEV)*ZMULT5+(1.-PACLCM(JL,KLEV))*ZMULT1
            ZDUS2=PACLCM(JL,KLEV)*ZMULT4+(1.-PACLCM(JL,KLEV))*VTMPC1
            ZTELDIF=ZLTETA1(JL,KLEV)-ZTESS(JL)
            ZQDDIF=ZQLWI1-ZQLWI2

            ZBUOY=ZDUS1*ZTELDIF+ZDUS2*ZTEMITTE*ZQDDIF     

! mz_lg_20060427+ further extension of GL mixing lenght concept

            IF (LLSCALE_MOD) THEN

CGEERT  DRY BUOY. 

c        ZTESS(JL)=PTSM1M(JL)*(1.E5/PAPHM1(JL,KLEVP1))**ZKAPPA
c        ZTVS(JL)=ZTESS(JL)*(1.+VTMPC1*0.8*ZQS(JL))
c        ZTVNLEV = ZTETA1(JL,KLEV)*(1.+VTMPC1*PQM1(JL,KLEV))
c        ZBUOY = ZTVNLEV - ZTVS(JL)
     
              ZGBUOY(JL,KLEV) = G/ZVIRMITTE*ZBUOY*(G/PGEOM1(JL,KLEV)) 

c        dtkebuoy(klev)   = ZGBUOY(JL,klev)
c        zmixlength(klev) = PWSM1M(JL)
c        dtkeshear(klev)  = ZQLWI2*1e3 
c        dtkediss(klev)   = ZHSOIL(JL) 
c        zrichard(klev)   = PCVW(JL)
c        dtkenld(klev)    = PCVS(JL)
         
              DZHELP = PGEOM1(JL,KLEV)/G
              ZRI(JL)=PGEOM1(JL,KLEV)*ZBUOY/(ZVIRMITTE*
     &                 (ZDU2 + 1e-6*PTKEM1M(JL,KLEV)*DZHELP**2 ))
              ZGRI(JL,KLEV) = ZRI(JL)
            ELSE         
              ZRI(JL)=PGEOM1(JL,KLEV)*ZBUOY/(ZVIRMITTE*ZDU2)
            ENDIF
! mz_lg_20060427-

            ZRICLS(JL)=ZRI(JL)

C     LG- assigning the time of sunrise and mixing, indicated by the first time 
C     that ZRI < 0 
            
            IF (PSRFL(JL).GT.0.1.AND.TIME_SUNRISE.EQ.0) THEN
               TIME_SUNRISE=LTIME ! assigning the sunrise time
            ELSEIF (PSRFL(JL).LT.0.1.AND.TIME_SUNRISE.GT.0) THEN
               TIME_SUNRISE=0   ! setting to zero again in the evening
               TIME_MIXING=0
            ENDIF

            IF (ZRI(JL).LT.-0.001.AND.TIME_MIXING.EQ.0) THEN
               TIME_MIXING=LTIME
               IF (TIME_SUNRISE.GT.0) THEN
                  DHR=INT(FLOAT(TIME_MIXING-TIME_SUNRISE)/100.)
                  IF (INT(FLOAT(TIME_MIXING-TIME_SUNRISE)-DHR*100).GT.60) THEN 
                     DMIN=INT((FLOAT(TIME_MIXING-TIME_SUNRISE)-DHR*100)*(60./120.))
                  ELSE
                     DMIN=INT(FLOAT(TIME_MIXING-TIME_SUNRISE)-DHR*100)
                  ENDIF
                  WRITE(*,'(1a,i5,1a)')
     &                 ' The time of the onset of the mixing (Ri < -0.001) is: ',TIME_MIXING,
     &                 ' [hr:min]'
                  WRITE(*,'(1a,i3,1a,i3,1a)')
     &                 ' The time between sunrise and the onset of the mixing is: ',
     &                 DHR,' hr(s) and: ',DMIN,' minutes'
               ENDIF
            ENDIF           	

C     LG- the local surface roughness has explicitly been derived from the
C     roughness for all the ecosystem classes as defined in the Olson database. 
C     The grid average surface roughness is then calculated from the fractions 
C     of land cover PCVS, PCVW, PVGRAT and an assumed roughness of 0.005 m for
C     bare snow and snow/ice covered surfaces (see calculation of local
C     drag coefficient over vegetated and bare soil and snow/ice covered
C     surfaces in the calculation of the "local" aerodynamic resistance for the
C     dry deposition velocity) and the vegetation roughness. This vegetation
C     roughness has also been used for the wet skin fraction. This is all done 
C     if the switch LAZ0LOC is on, otherwise the model uses the RACMO surface
C     roughness. Applying this local roughness does ofcourse only make sense
C     for an over land grid (LALAND(JL)=.TRUE.), otherwise, the roughness
C     is calculated according to the Charnock formula.

C     For the recalculation of grid surface roughness out of drag coefficients
C     see DKRZ (MPI, Hamburg) report nr 135, Claussen et al., pg 7, see also
C     this report for the definition of ZBLEND and its value

            ZBLEND=100.
            
C     LG- for a vegetation + wet skin fraction of about zero and a forest
C     fraction > 0 the surface cover properties reflect a site with 
C     snow cover on the forest floor and therefore the roughness must be
C     corrected since the rouhgness is also affected by the branches and
C     trunk

            IF (PCVS(1).GT.1.E-10.AND.FSTFRAC(1).GT.0.) THEN
               IF (NSTEP.EQ.0) THEN
                  WRITE(*,'(1a)')
     &                 ' Snow cover fraction > 0., and Forest fraction > 0'
                  WRITE(*,'(1a)')
     &                 ' The snow surface roughness is corrected for the height of the',
     &                 ' leafless canopy'
               ENDIF

C     LG-   call subroutine for calculation of the surface roughness and 
C     displacement heigth from the LAI and canopy height according
C     to Raupach, "Simplified expressions for vegetation roughness.."
C     in Boundary Layer Meteorology, 71, 211-216, 1994

C     LG-   a minimal LAI of 0.01 has been selected since this gives a z0
C     of about 0.5 m for a canopy height of 15 m 

               CALL CALCZ0D(0.1,HC,AZ0MSNOW,DISPSNOW)

C     LG-   correction for the forest fraction, with an assumed roughness of 
C     the non-forested fraction of about 0.005 m

               AZ0MSNOW=ZBLEND/(EXP(SQRT(1./(
     *              FSTFRAC(JL)/
     *              ((ALOG(ZBLEND/AZ0MSNOW))**2.)+
     *              (1.-FSTFRAC(JL))/
     *              ((ALOG(ZBLEND/0.005))**2.)))))

               IF (.NOT.LAGRIAN.AND.NSTEP.EQ.0.OR.LAGRIAN.AND.
     *              MOD(NSTEP,NPRINT).EQ.0)
     *              WRITE(*,'(1a,f6.4,1a)')
     *              ' The snow fraction surface roughness is : ',AZ0MSNOW,' [m]'

            ENDIF

            IF(LAZ0LOC.AND.LALAND(JL)) THEN
               PAZ0M(JL)=ZBLEND/(EXP(SQRT(1./(
     *              PCVS(JL)/   !snow fraction
     *              ((ALOG(ZBLEND/AZ0MSNOW))**2.)+
     *              (1.-PCVS(JL))*PCVW(JL)/ ! wet skin fraction
     *              ((ALOG(ZBLEND/AZ0MLOC))**2.)+
     *              (1.-PCVS(JL))*(1.-PCVW(JL))*(1.-PVGRAT(JL))/ !bare soil
     *              ((ALOG(ZBLEND/AZ0MSOIL))**2.)+
     *              (1.-PCVS(JL))*(1.-PCVW(JL))*PVGRAT(JL)/ !vegetation fraction
     *              ((ALOG(ZBLEND/AZ0MLOC))**2.)))))

C     LG-  for a vegetation + wet skin fraction less than 1, the surface 
C     roughness will be less than the selected value for initialization
C     therefore the option exists to re-define the z0 

               IF (NSTEP.EQ.0.AND..NOT.LAGRIAN.AND.((1.-PCVS(JL))*PCVW(JL)+
     &              (1.-PCVS(JL))*(1.-PCVW(JL))*PVGRAT(JL)).LT.1.) THEN
                  WRITE(*,'(1a)')
     &                 ' Vegetation + wet skin fraction < 1., this will decrease the selected z0'
                  WRITE(*,'(1a,f7.4,1a)')
     &                 ' Type 1 to use the initial surface roughness of :',AZ0MLOC,' [m]'
                  ISETZ0=0
                  READ(*,*)ISETZ0
                  WRITE(*,'(1a,f7.4,1a)')
     &                 ' The surface roughness being used : ',AZ0MLOC,' [m]'
                  WRITE (*,'(1a)')' Enter to continue'
                  READ (*,*)
               ENDIF

            ENDIF

            IF (ISETZ0.EQ.1) PAZ0M(JL)=AZ0MLOC

C     LG- end

            ZCDN(JL)=(ZKAP/LOG(1.+PGEOM1(JL,KLEV)/(G*PAZ0M(JL))))**2
            Z0H=PAZ0M(JL)*EXP(2.-86.276*PAZ0M(JL)**0.375)
            
C     LG- definition of surface roughness for heat over land, being 1/10 the 
C     roughness of momentum for vegetated areas and resembling
C     the surface roughness for momentum for other land cover 

            IF (LALAND(JL)) Z0H=ZBLEND/(EXP(SQRT(1./(
     *           PCVS(JL)/
     *           ((ALOG(ZBLEND/0.005))**2.)+
     *           (1.-PCVS(JL))*PCVW(JL)/
     *           ((ALOG(ZBLEND/(0.1*AZ0MLOC)))**2.)+
     *           (1.-PCVS(JL))*(1.-PCVW(JL))*(1.-PVGRAT(JL))/
     *           ((ALOG(ZBLEND/0.005))**2.)+
     *           (1.-PCVS(JL))*(1.-PCVW(JL))*PVGRAT(JL)/
     *           ((ALOG(ZBLEND/(0.1*AZ0MLOC)))**2.)))))

            IF (ISETZ0.EQ.1) Z0H=0.1*AZ0MLOC

C     LG- end

            ZALO=LOG(1.+PGEOM1(JL,KLEV)/(G*PAZ0M(JL)))
            ZALOH=LOG(1.+PGEOM1(JL,KLEV)/(G*Z0H))
            ZCHN(JL)=ZKAP**2/(ZALO*ZALOH)
            ZUCF(JL)=1./(1.+ZCONS11*ZCDN(JL)*SQRT(ABS(ZRI(JL))*(1.
     &           +PGEOM1(JL,KLEV)/(G*PAZ0M(JL)))))
            ZSCF(JL)=SQRT(1.+ZD*ABS(ZRI(JL)))
            ZCONS=ZCONS12*PAPHM1(JL,KLEVP1)/
     &           (PTM1(JL,KLEV)*(1.+VTMPC1*PQM1(JL,KLEV)-PXM1(JL,KLEV)))
            ZCFNC(JL)=ZCONS*SQRT(ZDU2)*ZCDN(JL)
            ZCFNCH(JL)=ZCONS*SQRT(ZDU2)*ZCHN(JL)
            ZDTHV=MAX(0.,(ZTVS(JL)-ZTVIR1(JL,KLEV)))
            ZWST(JL)=ZDTHV*SQRT(ZDU2)/ZVIRMITTE
            ZCR(JL)=(ZFREEC/(ZCHN(JL)*SQRT(ZDU2)))*ABS(ZBUOY)**ZCONS6
C     
C     CANOPY RESISTANCE
C     

C     LG- assigning the physiological stomatal resistance to the parameter
C     ZWET for LAGS=.TRUE.

            IF (LAGS) THEN
               ZWET(JL)=PWET_AGS(JL)/ZCONS
               ZRS0(JL)=PRS0_AGS(JL)/ZCONS
            ELSE

C     LG-  using the original ECHAM resistance

               ZWET(JL)=ZWET(JL)/ZCONS
               ZRS0(JL)=ZRS0(JL)/ZCONS
            ENDIF

 311     CONTINUE

C     
C     *    3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
C     *         BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
C     
 320     CONTINUE
C     

         DO 321 JL=KIDIA,KFDIA
            IF(ZRI(JL).GE.0.) THEN
               ZCFM(JL,KLEV)=ZCFNC(JL)/(1.+ZCONS8*ZRI(JL)/ZSCF(JL))
               IF(LALAND(JL)) THEN

C     LG-    the constant 3 is reduced to 2 in ZCONS9, to study the sensitivity 
C     of the PBL growth for very stable conditions, based on a presentation
C     by Anton Beljaars, mentioning this change for the ECMWF model, 
C     which has problems in representing the surface energy balance and
C     PBL structure over snow and ice covered surfaces (march 1998, and
C     see also the book of the 1998 workshop, "Clear and Cloudy Boundary
C     Layers" eds. Holtslag and Duynkerke, 1998)

                  ZCFH(JL,KLEV)=ZCFNC(JL)/(1.+(2./3.)*ZCONS9*ZRI(JL)*ZSCF(JL))
                  ZCH(JL)=ZCFH(JL,KLEV)/ZCFNC(JL)*ZCDN(JL)
               ELSE
                  ZCFH(JL,KLEV)=ZCFNCH(JL)/(1.+(2./3.)*ZCONS9*ZRI(JL)*ZSCF(JL))
                  ZCH(JL)=ZCFH(JL,KLEV)/ZCFNCH(JL)*ZCHN(JL)
               END IF
            ELSE
               ZCFM(JL,KLEV)=ZCFNC(JL)*(1.-ZCONS8*ZRI(JL)*ZUCF(JL))
               IF(LALAND(JL)) THEN
                  ZCFH(JL,KLEV)=ZCFNC(JL)*(1.-ZCONS9*ZRI(JL)*ZUCF(JL))

! mz_lg_20060427+ 
                  IF (LLSCALE_MOD) THEN

CGEERT              DO CORRECTION FOR LOW SHEAR CASES
C
C                   SQUARED SURFACE CONVECTIVE SCALE WSSURF
C                   W_C^2 = G/THETA * DTHETAV * DZ

                    WSSURF = 1.*MAX(0.,-PGEOM1(JL,KLEV)*ZGBUOY(JL,KLEV)/(ZVIRMITTE))

C                   CORRECTION
                    ZDU2=MAX(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)
                    ZCFH(JL,KLEV) = ZCFH(JL,KLEV)* SQRT(1. + WSSURF / ZDU2)
                  ENDIF
! mz_lg_20060427-

                  ZCH(JL)=ZCFH(JL,KLEV)/ZCFNC(JL)*ZCDN(JL)
               ELSE
C     
C     *             SPECIAL FREE CONVECTION LIMIT OVER SEA
C     
                  ZCFH(JL,KLEV)=ZCFNCH(JL)*((1.+ZCR(JL)**ZGAM)**Z1DGAM)
                  ZCH(JL)=ZCFH(JL,KLEV)/ZCFNCH(JL)*ZCHN(JL)
               END IF
            END IF

C     LG- 20030707+, added the assignment of the surface layer tracer eddy-
C     diffusivity

            ZCFTR(JL,KLEV)=ZCFH(JL,KLEV)

C     LG- 
*     /-----    parameters for dry deposition
*     I VDIFF.536
            ZGJCONS(JL)=ZCONS

C     LG- assigning of value of windspeed at reference height surface layer

            U30(JL)=SQRT(AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+
     *           PVM1(JL,KLEV)**2))

C     LG- calculation of parameter which is used to compute the 
C     aerodynamic resistance. This aerodynamic resistance is used 
C     in the routine "drydep.f" to calculate the deposition velocities 
C     of species which are not included in the Ganzeveld et al. 
C     "big leaf" dry deposition model

*     I VDIFF.566
            PCDV(JL)=ZCFH(JL,KLEV)/ZGJCONS(JL)
            PCDV(JL)=MAX(PCDV(JL),1.E-10)

C     LG- assigning the surface Richardson number which is used in the calculation
C     of the 2 m concentration after the calculation of the bulk dry
C     deposition velocity

            ZRISURF=ZRI(JL)

            CM(JL)=ZCDN(JL)*ZCFM(JL,KLEV)/ZCFNC(JL)
            USTAR(JL)=SQRT(CM(JL))*U30(JL)

C     LG- 
*     I VDIFF.567
C     Richardson number
C     computation of stability correction term, 
C     see paper Williams and Hicks et al. (LGP313)

C     LG- this is the original method used in the dry deposition scheme

            ZSURF=(PGEOM1(JL,KLEV)/G)+HC
            IF (ZRI(JL).GE.0) THEN
               ZMONIN=(USTAR(JL)*((ZTVIR1(JL,KLEV)+ZTVS(JL))/2.)*
     *              SQRT(AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)))/
     *              (ZKAP*G*(ZTVIR1(JL,KLEV)-ZTVS(JL)))
               ZOVERL=ZSURF/ZMONIN
               PSIH_OLD(JL)=-5.*ZOVERL
            ELSE
               ZMONIN=ZSURF/ZRI(JL)
               ZOVERL=ZSURF/ZMONIN
               PSIH_OLD(JL)=EXP(0.589+0.390*ALOG(-ZOVERL)-
     *              0.09*(ALOG(-ZOVERL))**2)
            ENDIF

C     LG- end

C     LG- these are the integrated stability correction functions over the range 
C     z0-z (surface roughness and reference height) !! The stability
C     correction functions are taken from Stull (page 383-385) (08-11-98)
C     and are slightly different from those by Williams and Hicks et al., 
C     (LGP313) which were originaly being used in the dry deposition scheme.
C     KMKH is the ratio of the eddy-diffusivity of momentum and heat. The sign
C     is opposite to that mentioned in Stull, however in the original scheme
C     there is minus sign in front of the integrated stability correction 
C     function is whereas Stull has a positive sign there!!!! 

            ZSURF=(PGEOM1(JL,KLEV)/G)+HC
            DZ=ZSURF-(AZ0MLOC+HC)
            ZRICRIT=2./(3.*ZD)  ! see ECHAM manual
            KMKH=0.74
            ALPHA=4.7

C     LG- calculation of bulk Richardson number according to equation 5.6.3
C     in Stull, page 177 for comparison with the ECHAM/RACMO calculated Ri

            ZRIB(JL,KLEV)=(G*(ZTVIR1(JL,KLEV)-ZTVS(JL))*DZ)/
     &           (((ZTVIR1(JL,KLEV)+ZTVS(JL))/2.)*
     &           AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))

            IF (ZRI(JL).GT.0) THEN

C     LG-   for stable conditions the Monin-Obukhov lenght is calculated from the 
C     Richardson number using the relationship Ri=z/L*(KMKH+4.7*Z/L)/
C     (1+4.7*z/L)^2 (see Stull page 183, figure 5.23), the parameter ZOVERL,
C     representing z/L, is calculated such that it resembles the assymptotic
C     behaviour for very stable conditions, with Ri approaching the value
C     of the critical Richardson number with increasing z/L.
               
               ZA=(ALPHA**2.)*MIN(ZRICRIT,ZRI(JL))-ALPHA
               ZB=2.*ALPHA*MIN(ZRICRIT,ZRI(JL))-KMKH
               ZC=MIN(ZRICRIT,ZRI(JL))
               ZOVERL=(-ZB-SQRT(ZB**2.-4.*ZA*ZC))/(2.*ZA)

C     an alternative method to calculate the Monin-Obukhov lenght 
C     directly, for stable conditions, applying the formula given by
C     Stull, 9.7.5k, page 386

               ZMONIN=(USTAR(JL)*((ZTVIR1(JL,KLEV)+ZTVS(JL))/2.)*
     *              SQRT(AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)))/
     *              (ZKAP*G*(ZTVIR1(JL,KLEV)-ZTVS(JL)))

C     LG-   calculation of Monin-Obukhov lenght according to formula 5.7c,
C     page 181 of Stull, with the sensible heat flux in K m s-1 

               IF (ABS(PAHFSM(JL)).LT.0.) ! set LT.0 to GT.0 for using equation
     *              ZMONIN=(((ZTVIR1(JL,KLEV)+ZTVS(JL))/2.)*USTAR(JL)**3.)/
     *              (ZKAP*G*(PAHFSM(JL)/(ZDIAGT*ZCPD*1.E3*ZRHOA(JL,KLEV))))

C     LG-   calculation of stability correction functions for a reference
C     height in the surface layer ZREFUM, which basically resembles the
C     height of the observed wind speed. This is done to be able to
C     calculate the modeled wind speed at this reference height using
C     the equation u=(u*/k)*(ln((zrefum-disp)/z0)+psim(z/L))

               ZOVERL=(ZREFUM+HC-DISP)/ZMONIN
               PSIMZZ(JL)=-ALPHA*ZOVERL ! momentum
               ZOVERL=ZSURF/ZMONIN
               PSIM(JL)=-ALPHA*ZOVERL ! momentum
               PSIH(JL)=-ALPHA*ZOVERL ! heat
            ELSE
               ZMONIN=ZSURF/ZRI(JL)
               ZOVERL=(ZREFUM+HC-DISP)/ZMONIN

C     LG-   see for the proper value of the exponent figure 9.9 of Stull and also
C     page 366 of "A review of flux-profile relationships" by Dyer, BLM 7,
C     363-372, 1974.

               XZSURF=(1.-9.*(ZOVERL))**(0.5) 
               XZREF=1.

               PSIMZZ(JL)=      ! momentum
     *              (2.*ALOG((1.+XZSURF)/2.)+ALOG((1.+XZSURF**2.)/2.)- 
     *              2.*ATAN(XZSURF))- ! primitive function value for z
     *              (2.*ALOG((1.+XZREF)/2.)+ALOG((1.+XZREF**2.)/2.)-
     *              2.*ATAN(XZREF)) ! primitive function value for zz

               ZOVERL=ZSURF/ZMONIN
               XZSURF=(1.-9.*(ZOVERL))**(0.5) 
               PSIM(JL)=	! momentum
     *              (2.*ALOG((1.+XZSURF)/2.)+ALOG((1.+XZSURF**2.)/2.)- 
     *              2.*ATAN(XZSURF))- ! primitive function value for z
     *              (2.*ALOG((1.+XZREF)/2.)+ALOG((1.+XZREF**2.)/2.)-
     *              2.*ATAN(XZREF)) ! primitive function value for zz

               ZMONIN=ZSURF/ZRI(JL)
               ZOVERL=ZSURF/ZMONIN
               XZSURF=KMKH*(1.-9.*(ZOVERL))**(0.5) 
               XZREF=KMKH
               PSIH(JL)=	! heat
     *              (2.*ALOG((1.+XZSURF)/2.)+ALOG((1.+XZSURF**2.)/2.)- 
     *              2.*ATAN(XZSURF))- ! primitive function value for z
     *              (2.*ALOG((1.+XZREF)/2.)+ALOG((1.+XZREF**2.)/2.)-
     *              2.*ATAN(XZREF)) ! primitive function value for zz
            ENDIF 

C     LG- computation of wind speed at reference height ZREFUM above
C     surface for comparison with observed wind speeds.

            UMZZ(JL)=U30(JL)
            IF (LUMZZ.AND.ZREFUM+HC.GT.DISP) THEN
               UMZZ(JL)=U30(JL)*(ALOG((ZREFUM+HC-DISP)/PAZ0M(JL))-PSIMZZ(JL))/
     *              (ALOG((ZSURF+HC)/PAZ0M(JL))-PSIM(JL))

C     LG-  added the estimated windspeed at zz without using the stability 
C     correction functions 

               UMZZN(JL)=U30(JL)*(ALOG((ZREFUM+HC-DISP)/PAZ0M(JL)))/
     *              (ALOG((ZSURF+HC)/PAZ0M(JL)))
            ENDIF

C     Computation of rah (s m-1), z0 for momentum!!, GEOM1/G=altitude,
            
            RAH(JL)=AMAX1(1.,(1./(USTAR(JL)*ZKAP))*
     *           (ALOG((PGEOM1(JL,KLEV)/G)/PAZ0M(JL))-PSIH(JL)))

C     -- The minimal z0 for vegetation is taken from the paper of Wieringa
C     0.02 m is the average of the reported range for short grass/moss

            IF (AZ0MLOC.GT.0.OR.LALAND(JL)) THEN
               ZCDNVEG(JL)=(ZKAP/ALOG(1.+PGEOM1(JL,KLEV)/
     *              (G*AMAX1(0.02,AZ0MLOC))))**2
               ZCDNSLSN(JL)=(ZKAP/ALOG(1.+PGEOM1(JL,KLEV)/
     *              (G*0.005)))**2
               CMVEG(JL)=ZCDNVEG(JL)*ZCFM(JL,KLEV)/ZCFNC(JL)
               CMSLSN(JL)=ZCDNSLSN(JL)*ZCFM(JL,KLEV)/ZCFNC(JL)
               USTVEG(JL)=SQRT(CMVEG(JL))*SQRT(AMAX1(ZEPDU2,
     *              PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))
               USTSLSN(JL)=SQRT(CMSLSN(JL))*SQRT(AMAX1(ZEPDU2,
     *              PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))
               RAHVEG(JL)=AMAX1(1.,(1./(USTVEG(JL)*ZKAP))*
     *              (ALOG((PGEOM1(JL,KLEV)/G)/
     *              AMAX1(0.02,AZ0MLOC))-PSIH(JL)))
               RAHSLSN(JL)=AMAX1(1.,(1./(USTSLSN(JL)*ZKAP))*
     *              (ALOG((PGEOM1(JL,KLEV)/G)/
     *              0.005)-PSIH(JL)))
            ELSE
               CMVEG(JL)=CM(JL)
               CMSLSN(JL)=CM(JL)
               USTVEG(JL)=USTAR(JL)
               USTSLSN(JL)=USTAR(JL)
               RAHVEG(JL)=RAH(JL)
               RAHSLSN(JL)=RAH(JL)
            ENDIF

C     Computation of STHETA from the calculated RAH (see DDIM)

            IF (PSRFL(JL).GT.10) THEN
               STHETA(JL)=AMIN1(0.5,SQRT(AMAX1(1.E-5,1./(RAHVEG(JL)/9.+
     *              SQRT(AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))))))
            ELSE
               STHETA(JL)=AMIN1(0.5,SQRT(AMAX1(1.E-5,1./(RAHVEG(JL)/4.+
     *              SQRT(AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))))))
            ENDIF

C     LG- end

C     ==========================================================================
C     LG- 07-2002 including a correcting for the water vapor pressure deficit 
C     to check if this plays a significant role in increasing the afternoon
C     stomatal resistance for tropical forest

C     LG- calculating the surface layer RH

            ES=0.611*EXP(19.59*(PTM1(JL,KLEV)-273.)/(PTM1(JL,KLEV)))
            E=0.611*EXP(19.59*(1.-273./AMIN1(PTM1(JL,KLEV),(TMELT-ZFRAC*ZCVM4)
     *           /(1.-ZFRAC))))
            RH(JL)=E/ES

C     LG- calculation of quasi-laminar boundary layer resistance

            RBVEG(JL,1)=(2./(USTVEG(JL)*0.40))

C     LG- calculation of the water vapour pressure of the air and within
C     the leaf

            TA=PTM1(JL,KLEV)-273.15
            VPAIR = 6.108*10**(7.5*TA/(237.3+TA))*(RH(JL)*100.)/100.
            VPLEAF = 6.108*10**(7.5*TA/(237.3+TA)) 
            CALL FVPD(VPLEAF,VPAIR,RBVEG(JL,1),RCO(JL),FV)
            IF (.NOT.LFVPD) THEN
               FV=1.
            ELSE
               ZWET(JL)=ZWET(JL)/FV
            ENDIF

C     LG- 07-2002, it seems that the quasi-laminar boundary layer resistance is
C     not included in the calculation of the latent heat flux. It still must
C     checked if this parameter is implicitly included in the Sellers 1986
C     SiB parameterization. This term could be important for forest canopies
C     with a low canopy stomatal resistance.

C     ZWET(JL)=(ZWET(JL)*ZCONS+RBVEG(JL,1))/ZCONS
C     print *,'ec4_vdiff, Rb added to ZWET!',zwet(jl)*fv*zcons,zwet(jl)*zcons,
C     &       rbveg(jl,1)

C     LG- end
C     ============================================================================

            ZCDUM(JL,KLEV)=ZCFM(JL,KLEV)
C     EVM------------------------------------------------------
            IF (LSRFFLUX_READ) THEN
               ZCFHSRF(JL)=ZCFH(JL,KLEV)
               WRIH   (JL)=ZRI (JL)
            ENDIF
C     EVM------------------------------------------------------
C     --------------------------------------------------------
C     IF *LSRFFLUX_OFF* = .TRUE. TURN OFF THE SURFACE FLUXES
C     FOR SENSIBLE AND LATENT HEAT
C     
            IF (LSRFFLUX_OFF) THEN
C     
               CALL RESETR (ZCFH (KIDIA,KLEV),KFDIA-KIDIA+1,0.)
               CALL RESETR (ZCFHSRF (KIDIA),KFDIA-KIDIA+1,0.)
C     
            END IF

C     EVM------------------------------------------------------
C     
C     INTERPOLATIONFUNCTIONS FOR DIAGNOSTICS
C     
            ZBN(JL)=ZKAP/SQRT(ZCDN(JL))
            ZBM(JL)=MAX(ZEPSEC,SQRT(ZCFM(JL,KLEV)*ZCDN(JL)*
     &           ZCONS17/ZCFNC(JL)))
            ZBH(JL)=MAX(ZEPSEC,ZCH(JL)/ZBM(JL)*ZCONS17)
            ZBM(JL)=1./ZBM(JL)
            ZBH(JL)=1./ZBH(JL)
 321     CONTINUE

C     LG- in case of LXTVDIFF=.TRUE., and LBIOSPH=.TRUE., 
C     calculation of some within canopy exchange parameters

         IF (LXTVDIFF.AND.LBIOSPH.OR.LXTVDIFF.AND.LBULKVEG.OR.
     &        LXTVDIFF.AND.LVEG_MLAY) THEN

C     LG-  calling of subroutine WNDPROF in which the vertical profile of the
C     the windspeed within the canopy is calculated, these calculations
C     are default done for the maximum amount of vegetation layers 

            CALL WNDPROF
            
C     LG-  assigning the windspeed at the reference height 

            U(NLEVV+1)=U30(JPHR)

C     LG-  calculation of the within-canopy windspeed at half the canopy height, 
C     which is being used in the 2+layer vegetation routine to calculate 
C     the within-canopy eddy diffusivity and the aerodynamic resistance. 

            DO JK=1,NLEVV
               IF((JK*HC/NLEVV)-0.5*(HC/NLEVV).LT.HC/2.) 
     &              UCAN(JPHR)=U(JK)
            ENDDO

C     LG-  calculation of Kh in order to apply these surface layer Eddy
C     diffusivity to determine the within-canopy ZCFNCH values based
C     on the scaling of the assigned within-canopy Kh profiles with
C     the surface layer Kh. The eddy diffusivity is determined from the
C     half the thickness of the surface layer (since the aerodynamic 
C     resistance is calculated over the height interval up to the
C     reference height of the surface layer) and the aerodynamic resistance
C     of the surface layer which resembles the reciprocal value of the
C     drag coefficient times the windspeed (see also calculation of 
C     canopy temperature)

            DO JL=KIDIA,KFDIA
               ZKH(JL,KLEV)=(PGEOM1(JL,KLEV)/G)/
     &              (1./((ZCFH(JL,KLEV)/ZCONS)*U30(JL)))
            ENDDO

            DO JL=KIDIA,KFDIA
               DO JK=KLEV+1,KLEVEL
                  JJ=KLEVEL+1-JK
                  JJJ=JJ+NLEVV-NLEVVEG ! for lagragian simulations, LAD profile

                  IF (JK.LT.KLEVEL) THEN
                     ZZDZ(JL,JK)=Z(JL,JK)-Z(JL,JK+1)
                  ELSE 
                     ZZDZ(JL,JK)=Z(JL,JK)
                  ENDIF
                  
C     LG-   simply scaling of the surface layer eddy-diffusivity with
C     windspeed profile to estimate within canopy K profile

                  ZKH(JL,JK)=ZKH(JL,KLEV)*(MIN(U30(JL),U(JJ)))/ 
     &                 ((U30(JL)+U(NLEVV))/2.)

C     LG-   calculating the eddy-diffusivity applying the mixing lenght
C     theory with the mixing lenght being estimated within the canopy
C     according to Goudriaan, Crop Micrometeorology

                  IF (LMXL) THEN
                     ZLEAFWIDTH=0.05
                     ZLM(JJ)=
     &                    MIN(ZZDZ(JL,JK),((4.*ZLEAFWIDTH)/
     &                    (API*MAX(1.E-10,LAD(JJJ))*LAI))**0.5)
                     ZIW(JJ)=0.5
                     ZKH_NEW(JL,JK)=ZLM(JJ)*ZIW(JJ)*U(JJ)       
                  ENDIF

C     LG-   the parameter ZCONS is calculated from surface layer properties
C     of the physical parameters also since these parameters are not yet
C     explicitly resolved within the canopy.

                  ZCONS=ZCONS12*ZAPHM1(JL,JK+1)/
     &                 (ZT(JL,JK)*(1.+VTMPC1*ZQ(JL,JK)-PXM1(JL,KLEV)))

                  ZDU2=MAX(ZEPDU2,(U(JJ+1)-U(JJ))**2)
                  ZCFNCH(JL)=ZCONS*SQRT(ZDU2)*ZCHN(JL)

C     LG-   calculating the eddy-diffusivity applying the mixing lenght
C     theory

                  IF (LMXL) THEN
                     ZDUDZ=(U(JJ+1)-U(JJ))/(Z(JL,JK-1)-Z(JL,JK))
                     ZCFNCH(JL)=ZCONS*ZLM(JJ)*ZDUDZ
                  ENDIF

C     LG-   assigning the bulk Richardson number for the vegetation layers
C     in order to study the impact of the within stability on the
C     the vertical turbulent exchange. Default, the models surface
C     layer stability is being used to define the stability within 
C     the canopy, thus implicitly assuming that there is good coupling 
C     between the surface layer and the canopy

                  ZRIB(JL,JK)=ZRI(JL)

C     LG-   set HC < 50 m, then observed (Kruijt et al., 1999) Ri numbers 
C     for tropical rainforest are being used

                  IF (HC.GT.50.) THEN
                     IF (LTIME.GT.0.AND.LTIME.LT.400) THEN
                        IF (Z(JL,JK).LT.9.) ZRIB(JL,JK)=-150.
                        IF (Z(JL,JK).GE.9.AND.Z(JL,JK).LT.17.) ZRIB(JL,JK)=108.
                        IF (Z(JL,JK).GE.17.) ZRIB(JL,JK)=1.57
                     ELSEIF (LTIME.GT.400.AND.LTIME.LT.800) THEN
                        IF (Z(JL,JK).LT.9.) ZRIB(JL,JK)=-543.
                        IF (Z(JL,JK).GE.9.AND.Z(JL,JK).LT.17.) ZRIB(JL,JK)=3.1
                        IF (Z(JL,JK).GE.17.) ZRIB(JL,JK)=14.7
                     ELSEIF (LTIME.GT.800.AND.LTIME.LT.1200) THEN 
                        IF (Z(JL,JK).LT.9.) ZRIB(JL,JK)=-44.
                        IF (Z(JL,JK).GE.9.AND.Z(JL,JK).LT.17.) ZRIB(JL,JK)=-430.
                        IF (Z(JL,JK).GE.17.) ZRIB(JL,JK)=31.3
                     ELSEIF (LTIME.GT.1200.AND.LTIME.LT.1600) THEN
                        IF (Z(JL,JK).LT.9.) ZRIB(JL,JK)=575.
                        IF (Z(JL,JK).GE.9.AND.Z(JL,JK).LT.17.) ZRIB(JL,JK)=540
                        IF (Z(JL,JK).GE.17.) ZRIB(JL,JK)=17.8
                     ELSEIF (LTIME.GT.1600.AND.LTIME.LT.2000) THEN
                        IF (Z(JL,JK).LT.9.) ZRIB(JL,JK)=75.5
                        IF (Z(JL,JK).GE.9.AND.Z(JL,JK).LT.17.) ZRIB(JL,JK)=241.8
                        IF (Z(JL,JK).GE.17.) ZRIB(JL,JK)=14.5
                     ELSEIF (LTIME.GT.2000.AND.LTIME.LT.2400) THEN
                        IF (Z(JL,JK).LT.9.) ZRIB(JL,JK)=-13.7
                        IF (Z(JL,JK).GE.9.AND.Z(JL,JK).LT.17.) ZRIB(JL,JK)=19.1
                        IF (Z(JL,JK).GE.17.) ZRIB(JL,JK)=0.57
                     ENDIF
                  ENDIF

                  IF(ZRIB(JL,JK).GE.0.) THEN
                     ZSCF(JL)=SQRT(1.+ZD*ABS(ZRIB(JL,JK)))

C     LG-    for the correction term in front of ZCONS9 (2./3.), see the comments
C     added to the calculation of the stability correction functions of 
C     heat

                     ZCFTR(JL,JK)=ZCFNCH(JL)/(1.+(3./3.)*ZCONS9*ZRIB(JL,JK)*ZSCF(JL))
                  ELSE
                     ZUCF(JL)=1./(1.+ZCONS11*ZCDN(JL)*SQRT(ABS(ZRIB(JL,JK))*(1.
     &                    +PGEOM1(JL,KLEV)/(G*PAZ0M(JL)))))
                     ZCFTR(JL,JK)=ZCFNCH(JL)*(1.-ZCONS9*ZRIB(JL,JK)*ZUCF(JL))
                  ENDIF

               ENDDO
            ENDDO

         ENDIF

C     LG- end IF (LXTVDIFF.AND.LBIOSPH)
         
 3200    CONTINUE

C     -----------------------------------------------------------------
C     
C     
C     INITIALIZE SURFACE EMISSION FOR TRACERS
C     

C     LG- originally out-commented, however reintroduced here for the chemistry

         IF(KTRAC.GT.0) THEN
            DO 3211 JT=1,KTRAC
               DO 3210 JL=KIDIA,KFDIA
                  ZXTEMS(JL,JT)=0.
 3210          CONTINUE
 3211       CONTINUE

C     LG- the surface emission and dry deposition are not calculated
C     since these processes are already accounted for in the chemistry
C     routine

C     SURFACE EMISSIONS AND DRY DEPOSITION
C     
C     IF(LEMIS) THEN
C     DO 3212 JL=KIDIA,KFDIA
C     Z1MXTM1(JL)=
C     *      PAPM1(JL,KLEV)/(PTM1(JL,KLEV)*RD*(1.+VTMPC1*PQM1(JL,KLEV)))
C     3212     CONTINUE
C     
C     CALL XTEMISS(KLON,   KLEV,     IROW,   CVDIFTS,  DTIME, ITASK,
C     *                 PXTM1,  ZXTEMS,   Z1MXTM1,
C     *                 LALAND, PFORESTM, PSNM1M)
C     ENDIF

         ENDIF

C---------------------------------------------------------------------
C     
C     *       3.3       EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY COEFFICIENT
C     
 330     CONTINUE
C     
         DO 331 JL=KIDIA,KFDIA

            ZWET(JL)=PCVS(JL)+(1.-PCVS(JL))*(PCVW(JL)+(1.-PCVW(JL))/
     *           (1.+ZCFH(JL,KLEV)*ZWET(JL)))
            ZWET(JL)=CVMGT(ZWET(JL),1.,LALAND(JL))
C     OBC
            LO=(ZHUM(JL)*ZQS(JL)).LE.PQM1(JL,KLEV)
            ZCSAT(JL)=PCVS(JL)+
     *           (1.-PCVS(JL))*(PCVW(JL)+(1.-PCVW(JL))*CVMGT(0.,ZHUM(JL),LO))
            ZCAIR(JL)=PCVS(JL)+
     *           (1.-PCVS(JL))*(PCVW(JL)+(1.-PCVW(JL))*CVMGT(0.,1.,LO))
            ZCSAT(JL)=CVMGT(ZCSAT(JL),1.,LALAND(JL))
            ZCAIR(JL)=CVMGT(ZCAIR(JL),1.,LALAND(JL))
            LO=PQM1(JL,KLEV).GT.ZQS(JL)
            ZCSAT(JL)=CVMGT(1.,ZCSAT(JL),LO)
            ZCAIR(JL)=CVMGT(1.,ZCAIR(JL),LO)
            ZCSAT(JL)=PVGRAT(JL)*ZWET(JL)+(1.-PVGRAT(JL))*ZCSAT(JL)
            ZCAIR(JL)=PVGRAT(JL)*ZWET(JL)+(1.-PVGRAT(JL))*ZCAIR(JL)
C     EVM------------------------------------------------------
            IF (LSRFFLUX_READ) THEN
               IF (NSTEP.GT.0) THEN
                  ZCSAT(JL)=ZBETA(JL)
                  ZCAIR(JL)=ZBETA(JL)
               END IF
            END IF
C     EVM------------------------------------------------------
            ZCPTS(JL)=PTSM1M(JL)*CPD*(1.+VTMPC2*
     *           (ZCSAT(JL)*ZQS(JL)+(1.-ZCAIR(JL))*PQM1(JL,KLEV)))


C     LG- replacing the surface temperature by the skin temperature

            IF (LTVEG)
     &           ZCPTS(JL)=(1.-(VEGFRAC(JL)+WSFRAC(JL)))*
     &           PTSM1M(JL)*CPD*(1.+VTMPC2*
     &           (ZCSAT(JL)*ZQS(JL)+(1.-ZCAIR(JL))*PQM1(JL,KLEV)))+
     &           (VEGFRAC(JL)+WSFRAC(JL))*
     &           TVEGM(JL)*CPD*(1.+VTMPC2*
     &           (ZCSAT(JL)*ZQS(JL)+(1.-ZCAIR(JL))*PQM1(JL,KLEV)))

C     LG========================================================================
C     calculation of canopy temperature from old temperature, 
C     energy balance components, soil temperature and coupling
C     parameter for daytime stable and nocturnal free convection conditions

            IF (LTVEG) THEN
               CSB=5.67E-8      ! Stefan-Boltzmann constant (see Stull page 259)
               CEPS=0.622       ! ratio of gas constant for air and water vapor (Stull 275)
               CLV=2.5E6        ! latent heat of vaporization (see Stull 641)
               CR=287.04        ! gas constant for dry air (see Stull 640)
               CEMISIR=1.       ! infrared emissivity (see Stull 259)

C     LG-  definition of coupling coefficient  (ECMWF model and Bosveld),
C     see thesis of v/d Hurk, page 113 for definition

               LABDA=9.5        ! coupling coefficient for Tveg > Tsoil [W m-2 K-1]
               IF (TVEGM(JL).LT.PTD3(JL)) LABDA=20. ! and for Tveg < Tsoil

C     LG-  determining the fraction of shortwave radiation that reaches the
C     the soil surface, 1 - this fraction is the amount of net radation
C     being intercepted by the canopy and thus being used to warm up the 
C     canopy

               IF (PSRFL(JL).GT.0.) THEN       
                  FRACRAD_SOIL=(FSL(1)*RBVD+RVD(1))/PSRFL(JL) 
		  ! the fraction of the net radiation reaching the soil
		  ! FSL(1) is the fraction of sunlit leaves in the lowest
		  ! canopy layer discerned for the radiation calculations
               ELSE
                  FRACRAD_SOIL=0.
               ENDIF

               IT=TVEGM(JL)*1000.
               ZES=TLUCUA(IT)/PAPHM1(JL,KLEVP1)
               ZCOR=1./(1.-VTMPC1*ZES)
               ZQS(JL)=ZES*ZCOR

C     LG-  Clausius-Clapeyron equation [g g-1 K-1]

               DQDT=(CEPS*CLV*ZQS(JL))/(CR*TVEGM(JL)**2.) ! [kg K-1] 
               ZRHO=PAPHM1(JL,KLEVP1)/ZRD/
     &              (PTM1(JL,KLEV)*(1.+VTMPC1*PQM1(JL,KLEV))) !density [kg m-3]
               ZCONS=ZCONS12*PAPHM1(JL,KLEVP1)/ ! see calculation ZCFH and ZWET
     &              (PTM1(JL,KLEV)*(1.+VTMPC1*PQM1(JL,KLEV)-PXM1(JL,KLEV)))

C     LG-  this term is the density muliplied with the drag coefficient and the
C     windspeed. It is divided by ZCONS since this term has initially already
C     being added to the term ZCFH. The term ZRHOCHU actually reflects the
C     aerodynamic resistance multiplied with the density, the stomatal 
C     resistance for the evapotranspiration is represented by the terms 
C     ZCAIR and ZSCAT

               ZRHOCHU=ZRHO*ZCFH(JL,KLEV)/ZCONS ! [kg m-2 s-1]

C     LG-  manipulation of the equation 4.15 of the thesis by v/d Hurk, page 113
C     to derive an equation for the new canopy temperature yields two
C     terms, one in which the old canopy temperature and some terms, that
C     are independent of the canopy temperature, are being found [W m-2], and 
C     a second term in which all the terms that are a function of the new
C     canopy temperature are being found [W m-2 K-1], dividing term1 by term2
C     yields the new canopy temperature. 

               TERM1=
     &              PSRFL(JL)*(1.-FRACRAD_SOIL)          ! intercepted radiation
     &              +CSB*PEMTERM(JL,KLEVP1)*PTSM(JL)**4. ! incoming longwave radiation
     &              +CSB*CEMISIR*PTD3(JL)**4.	
     &              +3.*CSB*CEMISIR*TVEGM(JL)**4.        ! outgoing longwave radiation
     &              +ZRHOCHU*ZCPTGZ(JL,KLEV)             ! sensible heatflux
     &              +ZRHOCHU*CLV*ZCAIR(JL)*PQM1(JL,KLEV) ! latent heatflux
     &              -ZRHOCHU*CLV*ZCSAT(JL)*ZQS(JL)       ! latent heatflux
     &              +ZRHOCHU*CLV*ZCSAT(JL)*DQDT*TVEGM(JL)! latent heatflux
     &              +LABDA*PTD3(JL) ! [W m-2]	         ! conductivity, soil heatflux

               TERM2=
     &              4.*CSB*CEMISIR*TVEGM(JL)**3.         ! longwave radiation
     &              +ZRHOCHU*CPD                         ! sensible heatflux 
     &              +ZRHOCHU*CLV*ZCSAT(JL)*DQDT          ! latent heatflux
     &              +LABDA      ! [W m-2 K-1]	         ! conductivity, soil heatflux

               TVEG(JL)=TERM1/TERM2 ! [W m-2 / W m-2 K-1 => K]
            ELSE
               TVEG(JL)=PTSM(JL) ! assigning the surface temperature 
            ENDIF

C     PTD3(JL)=TVEGM(JL)+2.
C     IF (PSRFL(JL).GT.0.) PTD3(JL)=TVEGM(JL)-2.
C     SOILFLUX=LABDA*(TVEGM(JL)-PTD3(JL))
C     
C     print *,'soilflux',SOILFLUX
C     TVEGM(JL)=TVEG(JL)
C     IT=TVEGM(JL)*1000.
C     ZES=TLUCUA(IT)/PAPHM1(JL,KLEVP1)
C     ZCOR=1./(1.-VTMPC1*ZES)
C     ZQS(JL)=ZES*ZCOR
C     
C     SOILFLUX=PSRFL(JL)*(1.-FRACRAD_SOIL)
C     C      &        +CSB*PEMTERM(JL,KLEVP1)*PTSM(JL)**4
C     C      &        +CSB*CEMISIR*PTD3(JL)**4
C     C      &        -CSB*CEMISIR*TVEGM(JL)**4
C     &        +CSB*PEMTERM(JL,KLEVP1)*TVEGM(JL)**4.
C     &        -CSB*CEMISIR*TVEGM(JL)**4
C     &        +ZRHOCHU*(ZCPTGZ(JL,KLEV)-CPD*TVEGM(JL))
C     &        +ZRHOCHU*CLV*(ZCAIR(JL)*PQM1(JL,KLEV)-ZCSAT(JL)*ZQS(JL))
C     
C     PRINT *,ptd3(jl)
C     c      PTD3(JL)=TVEGM(JL)-SOILFLUX/LABDA
C     PRINT *,'new',
C     &        tvegm(jl),ptd3(jl),soilflux,
C     &        PSRFL(JL)*(1.-FRACRAD_SOIL),
C     &        CSB*PEMTERM(JL,KLEVP1)*PTM1(JL,KLEV)**4,
C     &        CSB*CEMISIR*PTD3(JL)**4,
C     &        -CSB*CEMISIR*TVEGM(JL)**4,
C     &        ZRHOCHU*(ZCPTGZ(JL,KLEV)-CPD*TVEGM(JL)),
C     &        ZRHOCHU*CLV*(ZCAIR(JL)*PQM1(JL,KLEV)-ZCSAT(JL)*ZQS(JL))
            
C     LG- calculation of bulk Richardson number for the canopy according to 
C     equation 5.6.3 in Stull, page 177, moisture and pressure corrections
C     are ignored based on the thin layer and the small changes in moisture.
C     The calculations are only performed when the vertical diffusion is
C     being calculated. 
            
            IF (LXTVDIFF.AND.LTVEG) THEN
               RIBVEG(JL)=(G*(TVEG(JL)-PTD3(JL))*(DISP+Z0H))/
     &              (((TVEG(JL)+PTD3(JL))/2.)*AMAX1(ZEPDU2,UCAN(JL))) ! see comment below
            ELSE
               RIBVEG(JL)=0.
            ENDIF
            ZRIB(JL,KLEV+1)=RIBVEG(JL)

C     LG- end

C     LG- calculation of within-canopy aerodynamic resistance from the 
C     within canopy windspeed and the temperature gradient between the 
C     reference height DISP+z0H and the soil temperature

            IF(RIBVEG(JL).GE.0.) THEN
               ZSCF(JL)=SQRT(1.+ZD*ABS(RIBVEG(JL)))

C     LG-   for the correction term in front of ZCONS9 (2./3.), see the comments
C     added to the calculation of the stability correction functions of 
C     heat

               ZSTAB=1./(1.+(3./3.)*ZCONS9*RIBVEG(JL)*ZSCF(JL))
            ELSE
               ZUCF(JL)=1./(1.+ZCONS11*ZCDN(JL)*SQRT(ABS(RIBVEG(JL))*(1.
     &              +PGEOM1(JL,KLEV)/(G*PAZ0M(JL)))))
               ZSTAB=(1.-ZCONS9*RIBVEG(JL)*ZUCF(JL))
            ENDIF

            ZKHN=(PGEOM1(JL,KLEV)/G)/(AMAX1(1.,(1./(USTVEG(JL)*ZKAP))*
     &           (ALOG((PGEOM1(JL,KLEV)/G)/AMAX1(0.02,AZ0MLOC)))))

C     LG===========================================================================

C     LG- in case of a stable surface layer it is assumed that the within-canopy
C     resistance is determined by the within-canopy stability regime. This
C     is basically occuring during the night when free convective conditions 
C     will prevail within the canopy. For an unstable surface layer, the 
C     within-canopy exchange is being controlled by non-local mixing events
C     and as a first estimate of the mixing a neutral eddy diffusivity is 
C     being applied, scaled with the average within canopy windspeed

            IF (LXTVDIFF.AND.KLEVEL.GT.KLEV) THEN
               IF (ZRISURF.GT.0.) THEN
                  RAHCAN(JL)=HC/(ZKHN*ZSTAB*(UCAN(JL)/U30(JL)))
               ELSE
                  RAHCAN(JL)=HC/(ZKHN*(UCAN(JL)/U30(JL)))
               ENDIF
            ELSE

C     LG-  Computation of in-canopy aerodynamic resistance of the original
C     dry deposition code (see Erisman & Van Pul, LGP209)

               RAHCAN(JL)=AMAX1(1.E-5,14.*LAI*HC/USTVEG(JL))

C     LG-  formulation of RAHCAN using the neutral K value and scaling with the
C     windspeed, it seems that the parameterization by Erisman and Van Pul
C     results in a much too large resistance (July 2000)

               IF (LXTVDIFF.AND.LBULKVEG.OR.LXTVDIFF.AND.LVEG_MLAY) 
     &              RAHCAN(JL)=HC/(ZKHN*(UCAN(JL)/U30(JL)))
            ENDIF

C     LG- end 

C     -- GOOSE, calculation of the integrated Vd from the mass size 
C     distribution for rural continental aerosol as defined by Warnecke,
C     observations by Mehlmann (1986), three integrations are available
C     rural continental aerosol over land, rural continental aerosol over water
C     and the marine aerosol over water. The empirical relationships are 
C     polynominals (3th degree) determined by the curve fit procedure in Kaleida
C     graph and express the relationship between u* and Vd SO4

C     -- rural continental Vd SO4
            
            VDSO4SLSN(JL)=AMAX1(1E-5,0.08-0.57*USTSLSN(JL)+
     *           1.66*USTSLSN(JL)**2-0.36*USTSLSN(JL)**3)
            
C     -- Marine Vd SO4 for marine sulfate distribution

            VDSO4WAT(JL)=AMAX1(1E-5,0.12-0.24*USTAR(JL)+
     *           5.25*USTAR(JL)**2-1.84*USTAR(JL)**3)

C     -- Marine Vd SO4 for rural continental SO4 mass size distribution

c     VDSO4WAT(JL)=AMAX1(1E-5,0.07+0.1*USTAR(JL)+
c     *    2.1*USTAR(JL)**2-0.20*USTAR(JL)**3)

C     Temperature correction term for HNO3 above snow, ice and soil
C     (Wesely, 1989), recomputation from K to 0C
C     The surface temperature of the previous timestep has been used
C     based on the very strong zonal pattern in tsurf as extracted from
C     the subroutine SURF 
C     -- GOOSE, parameterization based on Dasch and Cadle observations

C     LG- originally the parameter TSM (which is the old surface temperature)
C     has been used in the ECHAM-3D model but has however been replaced
C     here by PTS!

            RSNOWHNO3SO2(JL)=AMIN1(AMAX1(10.,10.**(-0.09*(PTS(JL)-273.)
     *           +2.4)),1.E+5)

C     Computation of relative humidity out T and Tdew (2m) (Monteith)

C     LG- determining the 2 m temperature being the average of the maximum and 
C     minimum 2 m temperature

            ZT2=(PT2MAXM(JL)+PT2MINM(JL))/2.

C     LG- end

            ES=0.611*EXP(19.59*(ZT2-273.)/(ZT2))
            E=0.611*EXP(19.59*(1.-273./AMIN1(ZT2,(TMELT-ZFRAC*ZCVM4)
     *           /(1.-ZFRAC))))
            RH(JL)=E/ES

C     Computation of rsoil for SO2 from soil pH two different rh classes
C     The rsoil values for each class are determined out of regression and
C     the average value of pH of the class. Concerning rh, a wet and dry class
C     is used with the threshold value of 60%

            RSOILSO2(JL)=AMAX1(50.,
     *           SOILPH(ILON,IROW,1)*115.+
     *           SOILPH(ILON,IROW,2)*65.+
     *           SOILPH(ILON,IROW,3)*25.+
     *           SOILPH(ILON,IROW,4)*25.+
     *           SOILPH(ILON,IROW,5)*70.+
     *           AMIN1(AMAX1(0.,1000.*EXP(269.-PTS(JL))),
     *           1.E+5))
            IF (RH(JL).LT.0.6.AND.RH(JL).GT.0.4) THEN
               RSOILSO2(JL)=AMAX1(50.,(RSOILSO2(JL)*3.41-85.)+
     *              AMIN1(AMAX1(0.,1000.*EXP(269.-PTS(JL))),
     *              1.E+5))
            ENDIF
            IF (RH(JL).LE.0.4) THEN
               RSOILSO2(JL)=AMAX1(50.,AMIN1(1000.,(RSOILSO2(JL)*3.41-85.+
     *              ((0.4-RH(JL))/0.4)*1.E+3)+
     *              AMIN1(AMAX1(0.,1000.*EXP(269.-PTS(JL))),
     *              1.E+5)))
            ENDIF

C     LG- end

 331     CONTINUE

C     LG- start calculation of "big-leaf" dry deposition velocities
         
*     I VDIFF.1266
C     Loop for tracers

         DO 702 JT=1,NTR

c     --   Only if a value for DIFF (diffusivity) is defined, the program 
c     calculates a surface resistance

            IF (DIFF(JT).LE.0.OR.JT.EQ.irad) GOTO 702

            DO 701 JL=KIDIA,KFDIA

               VDLG(JL,0)=0.
               VDLG(JL,JT)=1.E-10

C     Boundary layer resistance, Van Pul, see pag 82, k=0.4, 
C     ln(z0/zs)=2 and diffusivity correction, diffusivity for 
C     other components still need to be defined

               RBVEG(JL,JT)=(2./(USTVEG(JL)*0.40))*DIFFRB(JT)
               RBSLSN(JL,JT)=(1./(USTSLSN(JL)*0.40))*DIFFRB(JT)
               RB(JL,JT)=(1./(USTAR(JL)*0.40))*DIFFRB(JT)

C     Computation of stomatal resistance for component X, correction 
C     based on differences in molecular diffusion of H2O and X (see Forgro)

C     LG- 07-2002 added an estimate of the vapor pressure deficit effect

               RCOX(JL,JT)=DIFF(JT)*RCO(JL)/FV
               
C     Water stress correction fws

               RCOX(JL,JT)=RCOX(JL,JT)/FWS(JL)

C     LG-      assigning the stomatal resistance to the parameter RSTOMX 
C     which has been used in the subroutine BULKVEG.f to calculate 
C     the emission velocity of the leaf stomata (isoprene), and in 
C     other routines to calculate first-order estimates of the 
C     compensation point effect

               RSTOMX(JL,JT)=RCOX(JL,JT)

C     LG-      end

               RLEAF(JL,JT)=(1./((1./RCUT(JT))+(1./(RCOX(JL,JT)+RMES(JT)))))

C                IF (JT.EQ.ih2o2) THEN
C 	        RLEAF(JL,JT)=1.
C 	        print *,'ec4_vdiff, be carefull: RleafH2O2 set to 1 s m-1'
C 	       ENDIF

C     LG-      added the stomatal resistance of the physiological model, for 
C     use in the other subroutines, May 2001

               IF (LAGS) THEN

C     LG-       the stomatal resistance of the AGS model, PWET_AGS is already upscaled
C     to the canopy scale. Therefore to calculate the trace gas leaf resistance
C     this parameter is multiplied with the LAI to arrive at an estimate of the
C     "bulk" leaf resistance

                  RSTOMX_AGS(JL,JT)=DIFF(JT)*PWET_AGS(JL)*LAI
                  RLEAF(JL,JT)=
     &                 (1./((1./RCUT(JT))+(1./(RSTOMX_AGS(JL,JT)+RMES(JT)))))
               ENDIF

C     -- GOOSE, LAI derived from Olson! ECHAM's stomatal resistance is
C     corrected for light extinction, page 1355, Sellers, 1986

               RSVEG(JL,JT)=(1./((1./(RAHCAN(JL)+
     *              (1./(USTVEG(JL)*0.40))*DIFFRB(JT)+RSOIL(JT)))+
     *              (1./(RLEAF(JL,JT)/AMAX1(1E-5,LAI)))))

C     It can be assumed that HNO3 deposition only depends on ra
C     so that surface resistance = 0, however definition of minimal 
C     surface resistance in order to avoid extreme Vd values.

               IF (JT.EQ.ihno3.OR.JT.EQ.iso4) 
     *              RSVEG(JL,JT)=AMAX1(10.,RSVEG(JL,JT))

C     Computation of deposition velocity Vd (m s-1) above land
C     Vd for surface with vegetation

               VDVEG(JL,JT)=(1./(RAHVEG(JL)+RBVEG(JL,JT)
     *              +RSVEG(JL,JT)))        

C     DDIM sulfate deposition calculation according to Wesely

               IF(STHETA(JL).GT.0.175.AND.PSRFL(JL).GT.250) THEN
                  RSVEG(JL,iso4)=1./(0.01*USTVEG(JL))
                  VDVEG(JL,iso4)=1./(RAHVEG(JL)+RSVEG(JL,JT))
               ELSE
                  RSVEG(JL,iso4)=1./(0.002*USTVEG(JL))
                  VDVEG(JL,iso4)=1./(RAHVEG(JL)+RSVEG(JL,JT))
               ENDIF

C     Vd for surface with bare soil, soil boundary layer resistance
C     still needs to be defined 

               VDSOIL(JL,JT)=(1./(RAHSLSN(JL)+RBSLSN(JL,JT)
     *              +RSOIL(JT)))
               IF (JT.EQ.iso2) VDSOIL(JL,JT)=(1./(RAHSLSN(JL)+RBSLSN(JL,JT)
     *              +RSOILSO2(JL)))

C     For sulfate over snow and soil the integrated parameters are used
C     the terms are already in cm s-1 and thus recalculated to m s-1

               IF (JT.EQ.iso4) VDSOIL(JL,JT)=VDSO4SLSN(JL)/100.

C     Wet skin fraction

               VDWS(JL,JT)=(1./(RAHVEG(JL)+RBVEG(JL,JT)+RWS(JT)))

C     GOOSE it is assumed that the wet skin fraction mainly consists of
C     vegetation

               IF (JT.EQ.iso4) VDWS(JL,JT)=VDVEG(JL,JT)

C     Snow coverage

               IF (JT.EQ.ihno3.OR.JT.EQ.iso2) 
     *              RSNOW(JT)=RSNOWHNO3SO2(JL)

               VDSN(JL,JT)=(1./(RAHSLSN(JL)+RBSLSN(JL,JT)
     *              +RSNOW(JT)))

               IF (JT.EQ.iso4) VDSN(JL,JT)=VDSO4SLSN(JL)/100.

C     Total areal weighted deposition velocity

               VDLAND(JL,JT)=PCVS(JL)*VDSN(JL,JT)+
     *              (1.-PCVS(JL))*(1.-PCVW(JL))*PVGRAT(JL)*VDVEG(JL,JT)+
     *              (1.-PCVS(JL))*(1.-PVGRAT(JL))*(1.-PCVW(JL))*VDSOIL(JL,JT)+
     *              (1.-PCVS(JL))*PCVW(JL)*VDWS(JL,JT)

               IF (LRWATER.AND.WTFRAC(JL).GT.0.) THEN

C     LG-        20030129, added the calculation of the O3 water resistance as a 
C     function of parameters like the henry coefficient, turbulent mixing 
C     efficiency, etc. This is a code developed in echam4 and based on
C     various studies by Liss and Merlivat, Schwartz, Joffre and offers
C     the opportunity to describe more mechanistically the uptake of
C     gases by water     

                  CENTHAL(io3)=2560. ! units?
                  CHENRY(JL,io3)=3.56 ! units?

                  IF (PTS(JL).GT.225.) THEN
                     CHENRY(JL,io3)=0.08314*PTS(JL)*0.0409*1./(3.56*2.72**
     *                    (-2560.*(1./PTS(JL)-1./298.)))
                  ELSE
                     CHENRY(JL,io3)=0.4 
                  ENDIF

C     See Liss and Merlivat, the derived relationship is for a windspeed 
C     at a height of 10 m !!!, the extra z0 term in the first ALOG is 
C     introduced in order to avoid negative windspeeds in above water 
C     surfaces bordering land with extreme large surface roughness 
C     (> 10 m gives negative ln)

                  AVGWS(JL)=SQRT(AMAX1(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))
     *                 *ALOG((10.+PAZ0M(JL))/PAZ0M(JL))/
     *                 ALOG((PGEOM1(JL,KLEV)/G)/PAZ0M(JL))

                  IF (ZSPEED.GT.0.) AVGWS(JL)=ZSPEED

                  IF (AVGWS(JL).LE.3.6) THEN
                     CKL(JL,io3)=0.17*AVGWS(JL)
                  ENDIF
                  IF (AVGWS(JL).GT.3.6.AND.AVGWS(JL).LE.13.) THEN
                     CKL(JL,io3)=2.85*AVGWS(JL)-9.65
                  ENDIF
                  IF (AVGWS(JL).GT.13.) THEN
                     CKL(JL,io3)=5.9*AVGWS(JL)-49.3
                  ENDIF

C     LG-         Recalculation to cm s-1

                  CKL(JL,io3)=CKL(JL,io3)/3600.

C     LG-         Temperature correction kl (see Liss 1986, and Duce et al.)
C     CKL(JL,1)=CKL(JL,1)*600.**(0.5)*(1858-108.7*(PTS(JL)-273)+
C     *        2.77*((PTS(JL)-273)**2.)-0.026*((PTS(JL)-273)**3.))**(-0.5)

C     LG-         Computation of enhancement deposition velocity O3 by surfactants
C     the DOC is in mol C l-1, recalculated to mg surfactant l-1 
C     (M=12, g to mg 1e+3, 12/10 see Hunter and Liss, 
C     (0.036/0.0004)/6, see paper McKay et al.)

                  IF (JT.EQ.io3) DOC(ILONT21,IROWT21)=DOC(ILONT21,IROWT21)*1.001 ! just to check the
		                  ! sensitivity slowly increasing the initial DOC value

                  ENHANC(JL,io3)=MAX(1.,DOC(ILONT21,IROWT21)*12.E+3*(12./10.)*
     *                 (0.036/0.0004)/6.)

C     LG-         rwat in s m-1

                  IF (JT.EQ.io3) THEN
                     RWATER(JL,JT)=100.*(1./(0.25*4.*5.)+
     *                    (1./(CHENRY(JL,JT)*CKL(JL,JT)*ENHANC(JL,JT))))


C LG- 102004- new formula based on the paper by Chang et al., AE2004, where 
C     a hybrid version of the Liss and Merlivat and Wanninkhof's parameterization 
C     is proposed to reproduce the observed deposition velocities for low wind speeds
C     and which is not captured by the Liss and Merlivat parameterization

C LG- check: what are the terms (0.25*4.*5)? and why is CHENRY multiplied with ENHANC (should be ENHANC/HENRY)?
C                    RWATER(JL,JT)=100.*(1./(0.25*4.*5.)+
C     *                    (1./(CHENRY(JL,JT)*CKL(JL,JT)*ENHANC(JL,JT)+         ! Liss and Merlivat, pkw, Eq 12
C     *                         CHENRY(JL,JT)*(CONCIM*KIMO3+CONCDMS*KDMSO3))))  ! Wanninkhof, q, Eq 12

!                      print *,'ec4_vdiff, 2880: assigning rwater O3 to rwat(io3)!!'

                     RWAT(JT)=RWATER(JL,JT)
                  ENDIF

C     LG-         Boundary layer resistance over water (Joffre, 1988, LGP125)

                  REYN(JL)=USTAR(JL)*PAZ0M(JL)/(0.15*10.**(-4.))
                  REYN(JL)=AMIN1(AMAX1(REYN(JL),10.**(-4.)),10.**(4.))

C     LG-         Weighting factor

                  IF (REYN(JL).GT.0.1) THEN
                     WEIGHT(JL)=1./(1.+20.*(ALOG(REYN(JL))/
     *                    ALOG(10.)+0.92)**6.)
                  ELSE
                     WEIGHT(JL)=1.
                  ENDIF
                  AZ0S(JL)=(1-WEIGHT(JL))*(7.4*PAZ0M(JL)*
     *                 2.72**(-7.3*0.4*1.*REYN(JL)**0.25*
     *                 (0.15/0.25)**0.5))+
     *                 WEIGHT(JL)*30.*(0.15*10.**(-4.)/USTAR(JL))*
     *                 2.72**(-13.6*0.4*1.*(0.15/0.25)**(2./3.))
                  AZ0S(JL)=AMAX1(1.5*10.**(-6.),AZ0S(JL))
                  AZ0S(JL)=AMIN1(100.*PAZ0M(JL),AZ0S(JL))
                  RBWAT(JL,JT)=1./(USTAR(JL)*0.4)*ALOG(PAZ0M(JL)/AZ0S(JL))*
     *                 DIFFRB(JT)

C     LG-         20030129- end calculation explicit water resistances

                  IF (JT.EQ.io3) THEN

                    IF (NSTEP.EQ.0) THEN
                      OPEN(UNIT=NUNVDWAT,FILE='/data/ganzevl/racmo/output/Vdwater.out',
     &                   STATUS='UNKNOWN')
                      WRITE(NUNVDWAT,'(1a)')'Oceanic O3 dry deposition parameters'     
                      WRITE(NUNVDWAT,'(a10,a15,5a12)')'istep','time','AVGWS',
     &                       'RwaterO3','RbwaterO3','DOC','Enhanc.'
                    ENDIF

                    IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN  
                       WRITE(NUNVDWAT,'(1x,i9.9,1x,a14,5e12.4)')
     &                      NSTEP,LDATLTIME,AVGWS(JL),RWATER(JL,JT),RBWAT(JL,JT),
     &                      DOC(ILONT21,IROWT21),ENHANC(JL,JT)
                    ENDIF

                    IF (NSTEP.EQ.NSTOP) CLOSE(NUNVDWAT)

                  ENDIF ! ENDIF (JT.EQ.IO3)

               ENDIF            ! end if (LRWATER)

               IF (PSEAICE(JL).GT.0) THEN
                  VDWAT(JL,JT)=PSEAICE(JL)*(1./(RAHSLSN(JL)
     *                 +RBSLSN(JL,JT)+RSNOW(JT)))+
     *                 (1.-PSEAICE(JL))*(1./(RAH(JL)+RB(JL,JT)+RWS(JT)))
                  IF (JT.EQ.ihno3.OR.JT.EQ.iso2) VDWAT(JL,JT)=
     *                 PSEAICE(JL)*(1./(RAHSLSN(JL)
     *                 +RBSLSN(JL,JT)+RSNOWHNO3SO2(JL)))+
     *                 (1.-PSEAICE(JL))*(1./(RAH(JL)+RB(JL,JT)+RWAT(JT)))
                  IF (JT.EQ.iso4) VDWAT(JL,JT)=PSEAICE(JL)*VDSO4SLSN(JL)/100.+
     *                 (1.-PSEAICE(JL))*VDSO4WAT(JL)/100.
               ELSE
                  VDWAT(JL,JT)=(1./(RAH(JL)+RB(JL,JT)+RWAT(JT)))
                  IF(JT.EQ.iso4) VDWAT(JL,JT)=VDSO4WAT(JL)/100.
               ENDIF

C     --       Computation of Vd for water/land, land grids, without
C     snow/ice cover

               IF(LALAND(JL)) THEN
                  VDLG(JL,JT)=100.*VDLAND(JL,JT)
                  IF (VEGFRAC(JL)+WSFRAC(JL).GT.0.999) THEN
                     RSURF(JL,JT)=(1./(0.01*VDLG(JL,JT)))-
     &                    (RAHVEG(JL)+RBVEG(JL,JT))
                  ELSE
                     RSURF(JL,JT)=1E5
                  ENDIF
               ELSE
                  VDLG(JL,JT)=100.*VDWAT(JL,JT)
               ENDIF

C     --       In the case of ice/snow cover ECHAM's landmask is used since
C     otherwise the snow/ice cover does not completely covers land
C     where it should be covered !!!

               IF(LALAND(JL)) THEN
                  IF(PCVS(JL).GT.0.) THEN
                     VDLG(JL,JT)=100.*VDLAND(JL,JT)
                  ENDIF
               ELSE
                  IF(PSEAICE(JL).GT.0.) THEN
                     VDLG(JL,JT)=100.*VDWAT(JL,JT)
                  ENDIF
               ENDIF

               VDLG(JL,JT)=AMAX1(1.E-5,VDLG(JL,JT))

 701        CONTINUE
 702     CONTINUE

C     LG- calculation of integrated stability correction function for the interval
C     zz (2 m) to z, this calculation is done by some manipulations of the
C     flux-profile relationships, which yield for an integral from zz (0 m) up
C     to z the stability correction function PSIH as it has been used in the
C     calculation of Ra, used to calculate the overall Vd, see Stull (page
C     383-385). This results in the calculation of the RAH for the specific 
c     height interval which is then used in the chemistry routine prechem.f
c     to calculate the concentration at height zz from the total emission and
c     dry deposition fluxes and the chemistry for these concentrations at zz

         DO 703 JL=KIDIA,KFDIA

C     LG-  defining the minimal height for which surface layer concentrations
C     are calculated being the surface roughness for the grid square PAZ0M

            ZREF=AMAX1(PAZ0M(JL),ZREFXTM)

C     LG-  stable conditions

            IF (ZRISURF.GE.0) THEN
               PSIHZZ(JL)=-4.7*(ZSURF-ZREFXTM)/ZMONIN

C     LG-  unstable conditions

            ELSE
               XZSURF=KMKH*(1.-9.*(ZSURF/ZMONIN))**(0.5)
               XZREF=KMKH*(1.-9.*(ZREF/ZMONIN))**(0.5)
               PSIHZZ(JL)=
     *              (2.*ALOG((1.+XZSURF)/2.)+ALOG((1.+XZSURF**2.)/2.)- 
     *              2.*ATAN(XZSURF))- ! primitive function value for z
     *              (2.*ALOG((1.+XZREF)/2.)+ALOG((1.+XZREF**2.)/2.)-
     *              2.*ATAN(XZREF)) ! primitive function value for zz
            ENDIF
            RAHZZ(JL)=(1./(USTAR(JL)*ZKAP))*
     *           AMAX1(1.,(ALOG(ZSURF/ZREF)-PSIHZZ(JL)))

 703     CONTINUE
         
C     END scheme Ganzeveld

C     LG- 

         IF (LCHEM) THEN

! mz_lg_20031203+ added a routine in which the required input parameters
!     to constrain the DayCent soil biogeochemical model are being
!     calculated/assigned. For more details see comchem.h and the DayCent 
!     routines

            IF (LDAYCENT) 
     &          CALL GET_WTH(NSTEP,DTIME,T2M,RAINTM,SNOWTM,
     &                       PSRFL,WS10M)

! mz_lg_20031203-

            IF (LEM_TURB) THEN

C     LG-   setting the parameter IEM/DD_TURB to 1 for those gases for which the
C     canopy top flux is being considered within this routine. This switch
C     is also used within the routine BULKVEG.f to determine for which
C     gases the surface layer concentration must be updated as a result of
C     vertical turbulent exchange calculated within that specific routine

               IEM_TURB(iisop)=1

               IF (ino_pr.LE.KTRAC) THEN
                  IEM_TURB(ino_pr)=0
               ELSE
                  IEM_TURB(inox)=0
               ENDIF
            ENDIF
            
C     LG-  calling of subroutine PRECHEM which is originally called in the 3-D model
C     version at the end of PHECHAM.f, however, since some of the surface
C     fluxes are considered within the vertical diffusion calculations to study
C     the impact of splitting these processes, that part of PRECHEM.f containing
C     the subroutines to calculate the atmosphere-biosphere/surface exchange 
C     fluxes is now called here before the loop of trace gas turbulent exchange
C     calculations

            CALL PRECHEM1
     *           (PGEOM1,PAPM1,PAPHM1,DAPM1,DAPHM1,PTM1,PTTE,PQM1,PQTE,PXTM1,
     *           PXTTE,PXTMVEG,PXTMSN,PXTEEMIS,PXTEDRYD,PXTETRNS,PXTESTRL,
     *           PXTEDIFF_VEG,SURFFLUX,KLEV,KLEVEL,KLEVELM1,KLEVELP1,NSTEP,
     *           DTIME,NSTOP,TWODT,NPRINT,PGLACM,PSEAICE,PSNM1M,PCATYPE)

C     LG-  end IF (LCHEM)

         ENDIF
         

C     LG- adding the tracers, in order to be able to calculate the vertical
C     diffusion tendency, this tendency is assigned here the value of
C     tracer tendency calculated up to this location within the code. 
C     Since there is already a vertical diffusion tendency being calculated
C     within the subroutine bulkveg.f, this vertical diffusion tendency 
C     must be corrected for this tendency

         DO JT=1,KTRAC
            DO JK=1,KLEVEL
               DO JL=KIDIA,KFDIA 
                  PXTEDIFF(JK,JT) = PXTTE(JL,JK,JT)
               ENDDO
            ENDDO
         ENDDO

         IF (LCHEM) THEN

C     LG- there is the opportunity to consider the emission and dry deposition
C     fluxes within this routine instead of resolving these within the 
C     subroutines emiss.f and drydep.f. The emission fluxes are accounted for 
C     by assigning the emission and dry deposition fluxes in kg m-2 s-1 
C     (see ECHAM manual) for the proper model levels. It should be noted that
C     in the emission and dry deposition routines, and also the routine 
C     bulkveg.f for LBULKVEG=.TRUE., then the same emission and dry 
C     deposition fluxes should be incorporated within the IF (LEM/DD_TURB) 
C     statement to avoid calculating twice the emissions/dry deposition fluxes 

            DO JT=1,KTRAC
               DO JL=KIDIA,KFDIA
                  DO JK=1,KLEVEL
                     EM(JL,JK,JT)=0.
                     DD(JL,JK,JT)=0.
                  ENDDO
               ENDDO
            ENDDO

            IF (LEM_TURB) THEN

               ZAVO=6.022045E23     
               XMC=12.
               
               DO JL=KIDIA,KFDIA
                  EM_TOT(JL,iisop)=0.
                  IF (ino_pr.LE.KTRAC) THEN
                     EM_TOT(JL,ino_pr)=0
                  ELSE
                     EM_TOT(JL,inox)=0
                  ENDIF
                  DO JK=1,KLEVEL

                     IF (.NOT.LBIOSPH.AND.JK.EQ.KLEVEL) THEN

C     LG-    for the bulk vegetation mode the surface fluxes, consisting
C     of the canopy top flux and the emission fluxes of the other
C     surfaces, is used as the lower boundary condition within the 
C     vertical diffusion routine 

                        IF (LBULKVEG.OR.LVEG_MLAY) THEN

C     The fluxes are generally defined in molecules m-2 s-1,
C     and are divided by the density in g m-3 to yield an concentration
C     increase in molecules g-1 m s-1. EMIS is then multiplied within 
C     the vertical loop with ZTMST to yield molec g-1 m , and 
C     G*ZDQP (1/DPD) in m s-2 Pa-1. To get from m s-2 Pa-1 to m-1
C     this term must be multiplied with the density in kPa s2 m-2, 
C     which therefore cancels out in the formulas except of the term 
C     1E3, The term CVDIFTS (1/ZTPFAC2) is introduced since in the 
C     final calculations, the emission flux is multiplied with 
C     ZTPFAC2

                           DO JT=1,NTRAC

C     LG-      when both IEM_TURB(JT) and IDD_TURB(JT) eq 1 then the
C     coupling between the emissions, dry deposition and turbulence
C     is already made within the routine VEG_MLAY.f, SURRFLUX is
C     generally in 1.e-11 molecules cm-2 s-1. Therefor the term must
C     multiplied with 1e4 to get the flux in molecules m-2 s-1

                              IF (IEM_TURB(JT).EQ.1.AND.LDRYDEP.AND.IDD_TURB(JT).EQ.0) THEN
                                 EM(JL,KLEVEL,JT)=1.E11*1.E4*SURFFLUX(JL,JT)*CVDIFTS/1.E3
                              ENDIF
                           ENDDO

                        ELSE

C     LG-     the isoprene emission flux is already a grid average
C     emission flux since the applied input parameters for the
C     isoprene emission calculations, dry matter, LAI etc, are
C     also averaged over the whole grid in the preprocessing.
C     Therefore, no weighting using the fraction of surface cover
C     is required in these calculations        

                           EM(JL,KLEVEL,iisop)=ISOPMOLEC*CVDIFTS/1.E3

c     EM(JL,KLEVEL,ino_pr)=NOMOLEC*CVDIFTS/1.E3

                        ENDIF
                     ELSEIF (LBIOSPH.AND.JK.GT.KLEV) THEN
                        JJ=KLEVEL+1-JK+(NLEVV-NLEVVEG)

C     LG-    ISOPEM is in kg C m-2 s-1 and thus needs to be recalculated to
C     molecules m-2 s-1, and the emission scaling factor also needs
C     to be considered

                        IEM_TURB(iisop)=1
                        EMRECALC=1.E3*(1./XMC)*(1./5.)*ZAVO
                        EM(JL,JK,iisop)=
     &                       FEMISOP*EMRECALC*ISOPEM(JJ)*CVDIFTS/1.E3

C     IF (JK.EQ.KLEVEL) THEN
C     IEM_TURB(ino_pr)=1
C     EM(JL,JK,ino_pr)=NOMOLEC*CVDIFTS/1.E3
C     ENDIF

C     LG-    determining the parameter DD for HNO3, due to the large
C     dry deposition fluxes, numerical problems are introduced
C     and therefore the dry deposition and vertical diffusion are
C     also coupled for HNO3

                        DD(JL,JK,ihno3)=0.

C     LG-    calculation of total intergrated emission flux in molec. g-1 s-1

                        EM_TOT(JL,iisop)=EM_TOT(JL,iisop)+
     &                       EM(JL,JK,iisop)/CVDIFTS

C     EM_TOT(JL,ino_pr)=EM_TOT(JL,ino_pr)+
C     &      EM(JL,JK,ino_pr)/CVDIFTS

                     ENDIF

                  ENDDO
               ENDDO

C     LG-  writing emission fluxes to output file for checking the applied
C     emission rates, be carefull with the interpretation of these emission
C     rates for the different modes (big leaf, bulk vegetation and 
C     multilayer) due the fact that the fractions of surface cover are
C     differently considered in these approaches!

               IF (NSTEP.EQ.0) THEN
                  OPEN(UNIT=NUNEMVD,FILE='/data/ganzevl/racmo/output/emis_vd_turb.out',
     &                 STATUS='UNKNOWN')
                  WRITE(NUNEMVD,'(2a)')'Emission fluxes used in the routine ',
     &                 'ec4_vdiff.f in [molec. g-1 s-1]'     
                  WRITE(NUNEMVD,'(a10,a15,1a12)')'istep','time','Isoprene'
               ENDIF

               IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN  
                  IF (LBIOSPH) THEN
                     WRITE(NUNEMVD,'(1x,i9.9,1x,a14,1e12.4)')
     &                    NSTEP,LDATLTIME,EM_TOT(JPHR,iisop)
                  ELSE      
                     WRITE(NUNEMVD,'(1x,i9.9,1x,a14,1e12.4)')
     &                    NSTEP,LDATLTIME,EM(JPHR,KLEVEL,iisop)/CVDIFTS
                  ENDIF
               ENDIF
               
               IF (NSTEP.EQ.NSTOP) CLOSE(NUNEMVD)

C     LG-  end IF (EM_TURB) 

            ENDIF

            ZDIR  = -90. - 180./API*ATAN2(PVM1(JPHR,KLEV),PUM1(JPHR,KLEV))
            IF (ZDIR.LE.0.) ZDIR = ZDIR + 360.

            UVOBS(JPHR,KLEV)=-9999.99
            IF (NSTEP.GT.0.AND.UOBS(JPHR,KLEV).GT.-9999.AND.
     *           VOBS(JPHR,KLEV).GT.-9999.) 
     *           UVOBS(JPHR,KLEV)=SQRT(UOBS(JPHR,KLEV)**2+VOBS(JPHR,KLEV)**2)

C     LG- writing of RI, stability parameter, and some other turbulent
C     exchange parameters

            IF (NSTEP.EQ.NRESUM) THEN
               OPEN(UNIT=NUNRI,FILE='/data/ganzevl/racmo/output/micromet.out',
     *              STATUS='UNKNOWN')
               WRITE(NUNRI,'(5a)')
     *              'Surface roughness,drag coefficent (CM), windspeed z and zz, ',
     *              'friction velocity, Ri number, Monin-Obukhov lenght, Z/L, ',
     *              'PSIM,PSIH and CFH, surface layer, canopy top and soil temp. ',
     *              'and specific and relative humidity in surface layer.',
     *              'The stomatal resistance for the physiological model is added'
               WRITE(NUNRI,'(a10,a15,37a14)')
     *              'step','time','RG','z0m','CM','UM+VM','UMZZ','UMZZN','U*',
     *              'U*-veg','DIR','Rah','Rah_veg','Rah-can','Ri','Ri-bulk-sl',
     *              'HF(old)','L','z/L','PSIM','PSIH','PSIH_OLD','CFH','Tair',
     *              'Can.-temp','soiltemp.','Ribulk-can','q-sl','RH-sl','Rstom','FWS',
     *              'Rstom-leaf','Rstom-AGS','wsfrac','U+V-obs','T-obs',
     *              'Q-obs','X-obs','FVPD'
            ENDIF

            IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) 
     *           WRITE(NUNRI,'(1x,i9.9,1x,a14,37e14.5)')
     *           NSTEP,LDATLTIME,
     *           RG,PAZ0M(JPHR),CM(JPHR),SQRT(AMAX1(ZEPDU2,PUM1(JPHR,KLEV)**2+
     *           PVM1(JPHR,KLEV)**2)),UMZZ(JPHR),UMZZN(JPHR),USTAR(JPHR),
     *           USTVEG(JPHR),ZDIR,RAH(JPHR),RAHVEG(JPHR),RAHCAN(JPHR),ZRI(JPHR),
     *           ZRIB(JPHR,KLEV),PAHFSM(JPHR),ZMONIN,ZOVERL,PSIM(JPHR),PSIH(JPHR),
     *           PSIH_OLD(JPHR),ZCFH(JPHR,KLEV),PTM1(JPHR,KLEV),TVEG(JPHR),
     *           PTD3(JPHR),RIBVEG(JPHR),PQM1(JPHR,KLEV),RH(JPHR),
     *           PWET(JPHR)/FV,FWS(JPHR),RCO(JPHR)/(FWS(JPHR)*FV),
     *           PWET_AGS(JPHR)/FV,WSFRAC(JPHR),UVOBS(JPHR,KLEV),
     *           TOBS(JPHR,KLEV),QOBS(JPHR,KLEV),XOBS(JPHR,KLEV),FV      
            
            IF (NSTEP.EQ.NSTOP) CLOSE(NUNRI)

C     LG-  end IF (LCHEM) 

         ENDIF

C     LG- end

C     
C     *       3.4       COMPUTATION OF THE PBL EXTENSION.
C     
 340     CONTINUE
C     
         DO 341 JL=KIDIA,KFDIA
            ZDU2=MAX(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)
C     JHC
C     
C     ------------------------------------------------
C     ECHAM CORIOLIS-FUNCTION CHANGED INTO HIRLAM SINLAT ARRAY
C     ------------------------------------------------
C     
            ZCOR=MAX(ABS(2.*OMEGA*PSINLAT(JL)),ZEPCOR)
C     JHC
C     ZCOR=MAX(ABS(CORIOL(IROW)),ZEPCOR)
            LO=PAZ0M(JL).GT.ZEPZ0O
            ZCDN2M=CVMGT((ZKAP/LOG(1.+PGEOM1(JL,KLEV)/(G*ZEPZ0O)))**2,
     &           ZCDN(JL),LO)
            ZCDNR=ZCDN2M/ZCDN(JL)
            ZCFM2M=CVMGT(ZCFNC(JL)*ZCDNR*(1.-ZCONS8*ZRI(JL))/(1.+ZCONS11*
     &           ZCDN2M*SQRT(ABS(ZRI(JL))*(1.+PGEOM1(JL,KLEV)/(G*ZEPZ0O)))),
     &           ZCFM(JL,KLEV)*ZCDNR,LO.AND.ZRI(JL).LT.0.)
            ZUSTAR(JL)=SQRT(ZCFM2M*SQRT(ZDU2)*PTM1(JL,KLEV)*
     &           (1.+VTMPC1*PQM1(JL,KLEV)-PXM1(JL,KLEV))
     &           /(ZCONS12*PAPHM1(JL,KLEVP1)))
            ZHDYN(JL)=MIN(PGEOM1(JL,1)/G,ZCHNEU*ZUSTAR(JL)/ZCOR)
C     
            IHPBLC(JL)=KLEV
            IHPBLD(JL)=KLEV
 341     CONTINUE
C     
         DO 343 JK=KLEVM1,1,-1
            DO 342 JL=KIDIA,KFDIA
               ZDS=ZCPTGZ(JL,JK)-ZCPTGZ(JL,KLEV)
               ZDZ=PGEOM1(JL,JK)/G-ZHDYN(JL)
               IHPBLC(JL)=ICVMGT(JK,IHPBLC(JL),IHPBLC(JL).EQ.KLEV.AND.ZDS.GT.0.)
               IHPBLD(JL)=ICVMGT(JK,IHPBLD(JL),IHPBLD(JL).EQ.KLEV.AND.ZDZ.GE.0.)
 342        CONTINUE
 343     CONTINUE
C     
C     CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
C     TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
C     
         DO 344 JL=KIDIA,KFDIA

! mz_lg_20060427+ 
            IF (LLSCALE_MOD) THEN

CGEERT        A BETTER MONIN OB. LENGTH SCALE
              ZDU2=MAX(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)
              ZGWST(JL) = (ZCPTS(JL) - ZCPTGZ(JL,KLEV))
              DCPTMIN = 0.1*ZCPD
              IF (ABS(ZGWST(JL)).LE.DCPTMIN) THEN
                ZGWST(JL) = DCPTMIN 
              ENDIF 
       
              ZGWST(JL) = ZGWST(JL)*SQRT(ZDU2)/ZCPTS(JL)
              ZMONOB=(ZUSTAR(JL)**3)/(ZKAP*G*ZGWST(JL)*ZCH(JL))
        
c        zmixlength(klev) = ZMONOB
c        dtkebuoy(klev)   = ZUSTAR(JL)
c        dtkeshear(klev)  = ZGWST(JL)
c        dtkediss(klev)   = ZCH(JL)
c        zrichard(klev)   = ZGRI(JL,KLEV)

            ENDIF
! mz_lg_20060427-

            IHPBL(JL)=MIN0(IHPBLC(JL),IHPBLD(JL),KLEV-3)
            ZGHABL=MIN(50000.,PGEOM1(JL,IHPBL(JL)))
            IF(ZWST(JL).GT.ZEPSR) THEN
               ZCONVS=(ZWST(JL)*ZCH(JL)*ZGHABL)**ZCONS6

               IF (.NOT.LLSCALE_MOD) ! mz_lg_20060427+
     &            ZMONOB=(ZUSTAR(JL)**3)/(ZKAP*G*ZWST(JL)*ZCH(JL))

               ZSTABF=(PGEOM1(JL,KLEV)/(G*ZMONOB))**(ZCONS6*2.)
               ZSTABF=MIN(ZUSTF*3.,ZSTABF)
            ELSE
               ZCONVS=0.
               ZSTABF=0.
            ENDIF
            ZTKEVN(JL,KLEV)=(ZUSTF+ZSTABF)*(ZUSTAR(JL)**2)
     >           +ZWSTF*(ZCONVS**2)
            ZTKEVN(JL,KLEV)=MAX(ZTKEMIN,ZTKEVN(JL,KLEV))

! mz_lg_20060427+ 
            IF (LLSCALE_MOD) THEN

CGEERT testing simple BL formulation TKE

              ZTKEVN(JL,KLEV)=ZUSTF*ZUSTAR(JL)**2
     &           + ZWSTF*ZCONVS**2

c        ZTKEVN(JL,KLEV)=ZUSTF*ZUSTAR(JL)**2
c     &      * (MAX(0.1, 1 - ZGRI(JK,KLEV) ))**0.5

              ZTKEVN(JL,KLEV)=MAX(ZTKEMIN,ZTKEVN(JL,KLEV))
            ENDIF
! mz_lg_20060427-

C     LG- update 03-2000, GL scheme

C     GL..

            IF (LLSCALE_MOD) THEN
               ZUO = 0.2
               ZGMONIN(JL) = ZMONOB  ! mz_lg_20060427+ modified, was
	                             ! zgmonin(jl) = 300*( 1. -  EXP( -(ZUSTAR(JL)/ZUO)**2) )
            ENDIF
C     GL.

C     LG- end

 344     CONTINUE
C     
         IF(NSTEP.EQ.NSTART) THEN
            DO 345 JL=KIDIA,KFDIA
               PTKEM1M(JL,KLEV)=ZTKEVN(JL,KLEV)
               PTKEM(JL,KLEV)=ZTKEVN(JL,KLEV)
 345        CONTINUE
         ENDIF
C     *       3.5       VERTICAL LOOP.
C     
 350     CONTINUE
C***  
         DO 372 JK=KTDIA,KLEVM1
C***  
C     
C     *       3.6       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
C     *                 BUOYANCY, RICHARDSON NUMBER, MIXING LENGTHS.
C     
 360        CONTINUE
C     
C     MODIFIED RICHARDSON NUMBER IN VERTICAL LOOP
C     
C     DIR$ IVDEP
            DO 361 JL=KIDIA,KFDIA
               ZQTMIT=ZLWCMIT(JL,JK)+ZQMIT(JL,JK)
               ZFUX=ZFAXEN(JL,JK)/(ZCPD*ZTMITTE(JL,JK))
               ZFOX=ZFAXEN(JL,JK)/(ZRD*ZTMITTE(JL,JK))
               ZMULT1=1.+VTMPC1*ZQTMIT
               ZMULT2=ZFUX*ZMULT1-ZRVRD
               ZMULT3=ZRDRV*ZFOX*ZQSSM(JL,JK)/
     1              (1.+ZRDRV*ZFUX*ZFOX*ZQSSM(JL,JK))
               ZMULT5=ZMULT1-ZMULT2*ZMULT3
               ZMULT4=ZFUX*ZMULT5-1.
C     
               ZDUS1=ZCCOVER(JL,JK)*ZMULT5+(1.-ZCCOVER(JL,JK))*ZMULT1
               ZDUS2=ZCCOVER(JL,JK)*ZMULT4+(1.-ZCCOVER(JL,JK))*VTMPC1
               ZTELDIF=(ZLTETA1(JL,JK)-ZLTETA1(JL,JK+1))/ZHH(JL,JK)*G
               ZDQTOT=(PQM1(JL,JK)+PXM1(JL,JK))-(PQM1(JL,JK+1)+PXM1(JL,JK+1))
               ZQDDIF=ZDQTOT/ZHH(JL,JK)*G
               ZBUOY=(ZTELDIF*ZDUS1+ZTEMIT(JL,JK)*ZDUS2*
     1              ZQDDIF)*G/ZTVIRMIT(JL,JK)
               ZDIVV=(PUM1(JL,JK)-PUM1(JL,JK+1))**2
               ZDIVV1=(PVM1(JL,JK)-PVM1(JL,JK+1))**2
               ZSHEAR=MAX(ZEPDU2,(ZDIVV+ZDIVV1))*(G/ZHH(JL,JK))**2
               ZRI(JL)=ZBUOY/ZSHEAR

C     LG- update 03-2000, GL scheme

C     GL...
               ZGBUOY(JL,JK) = ZBUOY
               ZGSHEAR(JL,JK) = ZSHEAR 
               ZGRI(JL,JK) = ZRI(JL)

! mz_lg_20060427+ 
               IF (LLSCALE_MOD) THEN

CGEERT    CUIJPERS COEFFICIENTS USED
        	 ZMULT1=1.+VTMPC1*ZQTMIT
        	 ZMULT2=1.-ZQTMIT + ZRVRD*ZQSSM(JL,JK)*
     *                      (1.+ ZRDRV*ZFOX)
        	 ZMULT5=ZMULT2/
     *               (1.+ZRDRV*ZFUX*ZFOX*ZQSSM(JL,JK))
        	 ZMULT4=ZFUX*ZMULT5-1.
        	 ZXB = ZCCOVER(JL,JK)
        	 ZDUS1=ZXB*ZMULT5+(1.-ZXB)*ZMULT1
        	 ZDUS2=ZXB*ZMULT4+(1.-ZXB)*VTMPC1

        	 ZTELDIF=(ZLTETA1(JL,JK)-ZLTETA1(JL,JK+1))/ZHH(JL,JK)*G
        	 ZTEDIF = (ZTETA1(JL,JK)-ZTETA1(JL,JK+1))/ZHH(JL,JK)*G
        	 ZTVIRDIF = (ZTVIR1(JL,JK)-ZTVIR1(JL,JK+1))/ZHH(JL,JK)*G


        	 ZDQTOT=(PQM1(JL,JK)+PXM1(JL,JK))-(PQM1(JL,JK+1)+PXM1(JL,JK+1))
        	 ZQDDIF=ZDQTOT/ZHH(JL,JK)*G
        	 ZBUOY=( ZTELDIF*ZDUS1 +
     *               ZTEMIT(JL,JK)*ZDUS2*ZQDDIF ) * G/ZTVIRMIT(JL,JK)

        	 ZBUOY_MOIST =( ZTELDIF*ZMULT5 +
     *             ZTEMIT(JL,JK)*ZMULT4*ZQDDIF ) * G/ZTVIRMIT(JL,JK)

        	 ZBUOY_DRY = ZTVIRDIF * G/ZTVIRMIT(JL,JK)

        	 ZSDDIF = (ZCPTGZ(JL,JK)-ZCPTGZ(JL,JK+1))/ZHH(JL,JK)*G
        	 CPM = CPD * (1.+VTMPC2*PQM1(JL,JK))
        	 ZBUOY_SD = G / (CPM * ZTMITTE(JL,JK)) * ZSDDIF * 
     *                    (1. +  VTMPC1 * ZQTMIT) +
     *                     1.*  G * VTMPC1 * ZQDDIF   

        	 ZDIVV=(PUM1(JL,JK)-PUM1(JL,JK+1))**2
        	 ZDIVV1=(PVM1(JL,JK)-PVM1(JL,JK+1))**2
        	 ZSHEAR=MAX(ZEPDU2,(ZDIVV+ZDIVV1))*(G/ZHH(JL,JK))**2

        	 ZRI(JL)=ZBUOY/ZSHEAR
CGEERT
C                to prevent becoming infinit in case of low shear
        	 ZRI(JL) = ZBUOY/( ZSHEAR + 1e-6*PTKEM1M(JL,JK))
        	 ZGBUOY(JL,JK) = ZBUOY
        	 ZGSHEAR(JL,JK) = ZSHEAR 
        	 ZGRI(JL,JK) = ZRI(JL)
              ENDIF

! mz_lg_20060501+ assigning the Richardson number to ZRIB for diagnostics
               ZRIB(JL,JK) = ZRI(JL)

! mz_lg_20060427-

 361        CONTINUE
C***  
 372     CONTINUE
C***  
         
         IF (LLSCALE_MOD) THEN
C     GL...
C     
C     CGL...
C     
C     COMPUTE A SIMPLE QUADRATIC LENGTH SCALE
C     THIS IS A PRELIMINARY CODE, WHICH GIVES A "QUADRATIC" LENGTHSCALE
C     ACCORDING TO 1/ZMIXQUAD = 1/ZMIXUP + 1/ZMIXDW
C     ZMIXUP: LENGTHSCALE FROM BOTTOM (SAY GROUND)
C     LOWER BOUNDARY: GROUND OR LOWEST POSITION UNSTABLE LAYER
C     ZMIXDW: LENGTHSCALE FROM TOP (SAY INVERSION)
C     UPPER BOUNDARY: "INVERSION"
C     BOTH ZMIXUP AND ZMIXDW INCREASE FOR UNSTABLE AND NEUTRAL 
C     CONDITIONS, AND DECREASE FOR STABLE CONDITIONS.
C     THIS IS GOVERNED BY ZBUOYTHRESS; SMALL ZBUOYTHRESS ==> 
C     STRONG DAMPING STABLE CONDITIONS (INVERSION), AND 
C     INCREASE FOR UNSTABLE CONDITIONS. 
C     IN NEUTRAL CONDITIONS THE LENGTHSCALE ZMIXUP ~ ZMIXQUAD ~ 0.4*Z 
C     NEAR THE SURFACE
C     ZPRANDT: GOVERNS DIFFERENCE BETWEEN LENGTHSCALE FOR MOMENTUM
C     AND HEAT.
C     

! mz_lg_20060427+ 
      
	DO JL=KIDIA,KFDIA   
          TH_PARCEL_UP(JL) = PTM1(JL,KLEV)*(1.E5/PAPM1(JL,KLEV))**ZKAPPA + 0.2      
          TKE_PARCEL_UP(JL) = 0.1 
          DZ_PARCEL(JL) = 0.
          TKE_CBL(JL) = 0.

          DO JK = KLEVM1, KTDIA+1, -1          ! VERTICAL LOOP BOTTOM-UP

            RHO = 0.5*(PAPHM1(JL,JK+1)+PAPHM1(JL,JK+2))
     *             /(ZRD*(PTVM1(Jl,JK+1)))
            DZ = (PAPHM1(JL,JK+2) - PAPHM1(JL,JK+1))/(G*RHO)

            TH_ENV = PTM1(JL,JK)*(1.E5/PAPM1(JL,JK))**ZKAPPA  ! improve interpolation !!!
            ZEPSILON = DZ/500.
            TH_PARCEL_UP(JL) =  TH_PARCEL_UP(JL) +            ! simple entrainment 
     *        ZEPSILON*(TH_ENV - TH_PARCEL_UP(JL))          
            DBUOY = G / TH_ENV * (TH_PARCEL_UP(JL) - TH_ENV)  
            IF (TKE_PARCEL_UP(JL) .GT. -DBUOY*DZ) THEN
              DZ_PARCEL(JL) = DZ_PARCEL(JL) + DZ
              TKE_PARCEL_UP(JL) = TKE_PARCEL_UP(JL)  + DBUOY*DZ
              TKE_CBL(JL) =  TKE_CBL(JL) + DZ*0.5*
     *          (PTKEM1M(JL,JK) + PTKEM1M(JL,MIN(JK+1,KLEVM1)))
            ELSE
              ZFAC = TKE_PARCEL_UP(JL) / (-DBUOY*DZ) 
              DZ_PARCEL(JL) = DZ_PARCEL(JL) + ZFAC*DZ
              PLPARCEL(JL) = DZ_PARCEL(JL)
              TKE_PARCEL_UP(JL) = TKE_PARCEL_UP(JL) + DBUOY*DZ
              TKE_CBL(JL) =  TKE_CBL(JL) + ZFAC*DZ*0.5*
     *          (PTKEM1M(JL,JK) + PTKEM1M(JL,MIN(JK+1,KLEVM1)))
              TKE_CBL(JL) =  TKE_CBL(JL)/DZ_PARCEL(JL)
              GOTO 999
            ENDIF

          ENDDO    ! JK LOOP
999     ENDDO      ! JL LOOP


c-----------------------------------------------------------------------

	ZBUOYTHRESS = 0.2E-4            ! TUNING PARAMETERS !
	ZPRANDT =  0.6 ! 0.6 
	ZALPHA = 0.
	PI2 = 2. * ATAN(1.)
	ZFINEUTR = ZNEUTR *ZKAP 
	ZFIAS =  ZNEUTR *ZKAP * 3. / ZPRANDT         !   1.666

c	print *, 'zfias momentum = ', ZFIAS*ZPRANDT

	ZBRI = 4.0     ! 5
	IF (L_STABLE) ZBRI = 8.
!      
!     COMPUTE BOTTOM UP LENGTH SCALE
!
	DO JL=KIDIA,KFDIA                   
          ZMIXUPH(JL, KLEV) = 0.
          ZMIXUPM(JL, KLEV) = 0.
	ENDDO

	DO JK = KLEVM1, KTDIA+1, -1          ! VERTICAL LOOP BOTTOM-UP
           DO JL=KIDIA,KFDIA   

            RHO = 0.5*(PAPHM1(JL,JK+1)+PAPHM1(JL,JK+2))
     *             /(ZRD*(PTVM1(Jl,JK+1)))
            DZ = (PAPHM1(JL,JK+2) - PAPHM1(JL,JK+1))/(G*RHO)

            ZBUOY_FUL = 0.5*(ZGBUOY(JL,JK) + ZGBUOY(JL,MIN(JK+1,KLEVM1))) 
            X = ZBUOY_FUL / ZBUOYTHRESS 
            ZRI_FUL = 0.5*(ZGRI(JL,JK) + ZGRI(JL,MIN(JK+1,KLEV))) 
cgeert test 
c            IF (L_STABLE) ZRI_FUL = ZRI_FUL - 0.05

            XARM = ZBRI * PI2 * ZNEUTR * ZKAP  / (ZFIAS*ZPRANDT - ZFINEUTR)
            XARH = ZBRI * PI2 * ZNEUTR * ZKAP  / (ZFIAS - ZFINEUTR)

            XM = XARM * ZRI_FUL
            XH = XARH * ZRI_FUL

            IF (X.GT.0) THEN
              DZM = ZFINEUTR*DZ - (ZFIAS*ZPRANDT - ZFINEUTR)/PI2
     *                          *DZ*XM
              DZH = ZFINEUTR*DZ - (ZFIAS - ZFINEUTR)/PI2
     *                          *DZ*XH
            ELSE
              DZM = ZFINEUTR*DZ - (ZFIAS*ZPRANDT - ZFINEUTR)/PI2
     *                          *DZ*ATAN(XM)
              DZH = ZFINEUTR*DZ - (ZFIAS - ZFINEUTR)/PI2
     *                          *DZ*ATAN(XH)
            ENDIF

!           BOTTOM UP LENGTHSCALE

c            print *,'zmixuph',jk,x,dzh,ZFINEUTR*DZ,zfineutr,dz,(ZFIAS -ZFINEUTR) / PI2
c     *                          *DZ*XH,zmixuph(jl,jk+1),ZMIXUPH (JL, JK + 1) + DZH, 
c     *                 ZGBUOY(JL,JK),ZGBUOY(JL,JK+1)
c            read (*,*)

            ZMIXUPH(JL, JK) =  ZMIXUPH (JL, JK + 1) + DZH
            ZMIXUPH(JL, JK) = MAX(ZMIXUPH(JL, JK) , 1.e-20) ! mz_lg_20060427+ was zero
            ZMIXUPM(JL, JK) =  ZMIXUPM (JL, JK + 1) + DZM
            ZMIXUPM(JL, JK) = MAX(ZMIXUPM(JL, JK) , 1.e-20) ! mz_lg_20060427+ was zero

          ENDDO
	ENDDO

!
!     COMPUTE TOP DOWN LENGTH SCALE 
!
	DO JL=KIDIA,KFDIA                   
          ZMIXDWH(JL, KTDIA) = 0.
          ZMIXDWM(JL, KTDIA) = 0.
	ENDDO

	DO JK= KTDIA + 1, KLEVM1            ! VERTICAL LOOP TOP, DOWN
          DO JL=KIDIA,KFDIA     

            RHO = 0.5*(PAPHM1(JL,JK+1)+PAPHM1(JL,JK))
     *             /(ZRD*(PTVM1(Jl,JK)))
            DZ = (PAPHM1(JL,JK+1) - PAPHM1(JL,JK))/(G*RHO)

            ZBUOY_FUL = 0.5*(ZGBUOY(JL,JK) + ZGBUOY(JL,JK-1)) 
            X = ZBUOY_FUL / ZBUOYTHRESS

            ZRI_FUL = 0.5*(ZGRI(JL,JK) + ZGRI(JL,JK-1)) 
cgeert test 
c            IF (L_STABLE) ZRI_FUL = ZRI_FUL - 0.05

            XARM = ZBRI * PI2 * ZNEUTR * ZKAP  / (ZFIAS*ZPRANDT - ZFINEUTR)
            XARH = ZBRI * PI2 * ZNEUTR * ZKAP  / (ZFIAS - ZFINEUTR)

            XM = XARM * ZRI_FUL
            XH = XARH * ZRI_FUL

            IF (X.GT.0) THEN
              DZM = ZFINEUTR*DZ - (ZFIAS*ZPRANDT - ZFINEUTR) / PI2
     *                          *DZ*XM
              DZH = ZFINEUTR*DZ - (ZFIAS -ZFINEUTR) / PI2
     *                          *DZ*XH
            ELSE
              DZM = ZFINEUTR*DZ - (ZFIAS*ZPRANDT - ZFINEUTR) / PI2
     *                          *DZ*ATAN(XM)
              DZH = ZFINEUTR*DZ - (ZFIAS - ZFINEUTR) / PI2
     *                          *DZ*ATAN(XH)
            ENDIF

!           TOP DOWN LENGTHSCALE 

c            print *,'zmixdwh',jk,x,dzh,ZFINEUTR*DZ,zfineutr,dz,(ZFIAS -ZFINEUTR) / PI2
c     *                          *DZ*XH

            ZZZLAM = 75.
            Z2GEOMF=PGEOM1(JL,JK)+PGEOM1(JL,JK+1)
            ZHEIGHT = Z2GEOMF / (G*2.) 
c            ZMIX_HELP = ZZZLAM / (1 + ZZZLAM / (ZNEUTR *ZKAP * ZHEIGHT))
            ZMIX_HELP = ZZZLAM * EXP ( - ZHEIGHT / 500.)   ! Equation 18, Lenderink and Holtslag, 2004.
            IF ( .NOT. L_STABLE) ZMIX_HELP = 1.e-20  ! mz_lg_20060427+ was zero

            ZMIXDWH(JL,JK) =  ZMIXDWH (JL,JK - 1)  +  DZH             
            ZMIXDWH(JL,JK) = MAX (ZMIXDWH(JL,JK),ZMIX_HELP) 

c            print *,'zmixdwh',jk,zmixdwh(jl,jk),ZMIXDWH (JL,JK - 1),ZMIXUPH(JL,JK),dzh,zmix_help

            ZMIXDWM(JL,JK) =  ZMIXDWM (JL,JK - 1)  +  DZM
            ZMIXDWM(JL,JK) = MAX (ZMIXDWM(JL,JK),ZMIX_HELP)    

!           COMPOSED LENGTH SCALE OF TOP-DOWN LS AND BOTTOM-UP LS

            ZAVE = 1
            ZLENGTHRH = 1. / ZMIXDWH(JL,JK)**ZAVE 
     *              + 1. / ZMIXUPH(JL,JK)**ZAVE
            ZMIXQUADH(JL,JK) = 1. / ZLENGTHRH**(1./ZAVE)
            ZLENGTHRM = 1. / ZMIXDWM(JL,JK)**ZAVE 
     *              + 1. / ZMIXUPM(JL,JK)**ZAVE
            ZMIXQUADM(JL,JK) = 1. / ZLENGTHRM**(1./ZAVE)

           ENDDO
	ENDDO                                  ! END VERTICAL LOOP

C     END COMPUTATION FREE TURBULENCE LENGTH SCALE 
C     LENGTHSCALES ARE IN ZMIXQUADH (HEAT) AND ZMIXQUADM (MOMENTUM)
C          
! mz_lg_20060427-

      ENDIF                  ! (LLSCALE_MOD)     

C     BEGIN NEW LOOP
C     
C***  
         DO JK=KTDIA,KLEVM1
C***  
            DO JL=KIDIA,KFDIA

C     LG- end

C     
C     ASYMPTOTIC MIXING LENGTH FOR MOMENTUM AND
C     HEAT (ZLAM) ABOVE THE PBL AS A FUNCTION OF HEIGHT
C     ACCORDING TO HOLTSLAG AND BOVILLE (1992), J. CLIMATE.
C     
               ZHEXP=EXP(1.-PGEOM1(JL,JK)/PGEOM1(JL,IHPBL(JL)))
               ZLAM=ZZZLAM+(ZCONS3-ZZZLAM)*ZHEXP
               IF(JK.GE.IHPBL(JL)) THEN
                  ZCONS23=ZCONS25
               ELSE
                  ZCONS23=ZCONS2/ZLAM
               END IF
C     
C     MIXING LENGTH (BLACKADAR) + STABILITY DEPENDENT FUNCTION
C     
               Z2GEOMF=PGEOM1(JL,JK)+PGEOM1(JL,JK+1)
               ZZ2GEO=ZCONS2*Z2GEOMF
               ZMIX=ZZ2GEO/(1.+ZCONS23*Z2GEOMF)

C     LG- update 03-2000, GL scheme

C     GL...
               IF (LLSCALE_MOD) THEN

! mz_lg_20060427+
CGEERT
C         BRINKOP AND ROECKNER (1995) FORMULATION FOR THE 
C         MIXING LENGTH SCALES
C         ZNBRUNT = BRUNT VAISALA FREQUENCY (SQUARED)
C         ZENTR  = ENTRAINMENT PARAMETER
C         ZMAILBEN = CONSTANT MAILHOT AND BENOIT (1982) 
   
          ZMAILBEN= 1.00
c          ZENTR = 0.20
          ZENTRH  = ZENTR
c          ZENTRM  = ZENTR * MAX(1.,MIN(1.+ZGRI(JL,JK)*4.,10.))
c          ZRIHELP = MAX (ZGRI(JL,JK),0.0)
c          ZHELP = EXP ( - ZRIHELP *4)  
c     *          + 4. * ZRIHELP * EXP (-ZRIHELP/4.)
c     *          + 10 * (1-EXP(-ZRIHELP/10))

          ZHELP =  MAX(1.,MIN(1.+ZGRI(JL,JK)*2.,100.)) 
          ZENTRM  = ZENTR * MAX(1.,MIN(ZHELP,4.))
c          ZENTRM = ZENTR

          ZZZLAM = 75.0
          Z2GEOMF=PGEOM1(JL,JK)+PGEOM1(JL,JK+1)
          ZHEIGHT = Z2GEOMF / (G*2.) 

          ! limited to half the neutral value
          ZMIX2 = ZZZLAM / (1 + ZZZLAM / ( 0.5 * ZNEUTR *ZKAP * ZHEIGHT))

          ZMIXH =  MAX(ZMIXQUADH(JL,JK),ZMIX2)         
          ZMIXM =  MAX(ZMIXQUADM(JL,JK),ZMIX2)

          RPOW = 2                                          ! NO PHYS.
					                    ! ADDED FOR CONTINUITY BETWEEN UNSTABLE AND STABLE
          ZMIXH =  (ZMIX2**RPOW + ZMIXQUADH(JL,JK)**RPOW)**(1./RPOW) 
          ZMIXM =  (ZMIX2**RPOW + ZMIXQUADM(JL,JK)**RPOW)**(1./RPOW)

          ZNBRUNT = ZGBUOY(JL,JK)
          ZTKESQ  = SQRT(MAX(PTKEM1M(JL,JK),ZTKEMIN))

          IF (ZNBRUNT.GT.0) THEN
            ZMIXCH = 1. +  ZMIXH*SQRT(ZNBRUNT)/(ZENTRH*ZTKESQ)     
            ZMIXCM = 1. +  ZMIXM*SQRT(ZNBRUNT)/(ZENTRM*ZTKESQ)     
          ELSE
            ZMIXCH = 1.
            ZMIXCM = 1.
          ENDIF        

          ZSM=ZMAILBEN/ZMIXCM
          ZSH=ZMAILBEN/ZMIXCH

c new      
           RPOW =  - 2
           IF (ZNBRUNT.GT.0) THEN
             ZMIXH_STABLE = ZENTRH * ZTKESQ / SQRT(ZNBRUNT)
             ZMIXCH = (ZMIXH_STABLE**(RPOW) + ZMIXH**(RPOW))**(1./RPOW) 
             ZMIXCH = ZMIXH / ZMIXCH 
             ZMIXM_STABLE = ZENTRM * ZTKESQ / SQRT(ZNBRUNT)
             ZMIXCM = (ZMIXM_STABLE**(RPOW) + ZMIXM**(RPOW))**(1./RPOW) 
             ZMIXCM = ZMIXM / ZMIXCM 
           ELSE
             ZMIXH_STABLE = 0.
             ZMIXCH = 1.
             ZMIXCM = 1.
           ENDIF        

          ZSM=ZMAILBEN/ZMIXCM
          ZSH=ZMAILBEN/ZMIXCH

c          DTKENLD(JK) = ZFAC_INH 
CGEERT
! mz_lg_20060427-

               ELSE             ! ECHAM4 lengthscale 

C     LG- end

C     
C     STABILITY FUNCTIONS (LOUIS, 1979)
C     
                  ZALH2=ZMIX*ZMIX
                  ZUCF(JL)=1./(1.+ZCONS5*ZALH2*SQRT(ABS(ZRI(JL))*(((PGEOM1(JL,JK)
     *                 /PGEOM1(JL,JK+1))**ZCONS6-1.)/(PGEOM1(JL,JK)
     *                 -PGEOM1(JL,JK+1)))**3/(PGEOM1(JL,JK+1))))
                  IF(ZRI(JL).LT.0.) THEN
                     ZSH=ZSHN*(1.-ZCONS9*ZRI(JL)*ZUCF(JL))
                     ZSM=ZSMN*(1.-ZCONS8*ZRI(JL)*ZUCF(JL))
                  ELSE
                     ZSH=ZSHN/(1.+ZCONS8*ZRI(JL)*SQRT(1.+ZRI(JL)))
                     ZSM=ZSMN/(1.+ZCONS8*ZRI(JL)/SQRT(1.+ZD*ZRI(JL)))
                  END IF
C     GL
C     
C     ENDIF (LLSCALE_MOD)
               ENDIF
C     GL
C     
C     *       3.7       DIMENSIONLESS COEFFICIENTS MULTIPLIED BY PRESSURE
C     *                 THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE.
C     
 370           CONTINUE
C     
C     

C     LG- update 03-2000, GL scheme

C     GL...      ZZB=ZSHEAR*ZMIX*ZSM-ZBUOY*ZMIX*ZSH
C     GL...      ZDISL=ZDA1*ZMIX/ZTMST

               ZTKESQ=SQRT(MAX(PTKEM1M(JL,JK),ZTKEMIN)) ! mz_lg_20060427+
               ZZB=ZGSHEAR(JL,JK)*ZMIXM*ZSM-ZGBUOY(JL,JK)*ZMIXH*ZSH
               ZDISL=ZDA1*ZMIXM/ZTMST

C     LG- end

               ZKTEST=1.+(ZZB*ZTMST+SQRT(PTKEM1M(JL,JK))*2.)/ZDISL
               IF (ZKTEST.LE.1.) THEN
                  ZTKEVN(JL,JK)=ZTKEMIN
               ELSE
                  ZTKEVN(JL,JK)=MAX(ZTKEMIN,(ZDISL*(SQRT(ZKTEST)-1.))**2)
               ENDIF
               IF(NSTEP.EQ.NSTART) THEN
                  PTKEM1M(JL,JK)=ZTKEVN(JL,JK)
                  PTKEM(JL,JK)=ZTKEVN(JL,JK)
               END IF

C     LG- update 03-2000, GL scheme

C     GL...      ZZZM=ZMIX*ZSM*ZTKESQ
C     GL...      ZZZH=ZMIX*ZSH*ZTKESQ
               ZZZM=ZMIXM*ZSM*ZTKESQ 
               ZZZH=ZMIXH*ZSH*ZTKESQ

C     LG- end

C     LG- end
C     
               ZALF=PAPHM1(JL,JK+1)/(ZTVIRMIT(JL,JK)*ZHH(JL,JK)*ZRD)
               ZCFM(JL,JK)=ZZZM*ZCONS18*ZALF

C     LG- In some simulations with a updated surface cover description,
C     it seems that there are initialization problems resulting in an
C     extreme unstable energy balance which gives unrealistic changes
C     in the soil temperatures. This can be reduced by applying some
C     spin-up time (defined in namchem) that has already been used for
C     some of the chemistry calculations. It is also used to reduce the
C     initial fluxes of heat and moisture (07-2000) 

               ZZZH=ZZZH*MIN(1.,NSTEP/(T_SPINUP/ZTMST))

               ZCFH(JL,JK)=ZZZH*ZCONS18*ZALF

C     LG- assigning the tracer exchange coeff. 

               IF (LXTVDIFF) THEN

C     LG-  defining Kz vertical profiles for night and daytime in order
C     to compare model simulations for different timestep with a 
C     comparable vertical diffusion so that the influence of
C     changing timesteps on other processes can be studied

                  IF (LSENS_ZTMST) THEN
                     IF (RG.LT.1.E-5) THEN
                        IF (JK.LT.17) THEN 
                           ZZZH=1.E-5*MIN(1.,NSTEP/(T_SPINUP/ZTMST))
                        ELSE
                           ZZZH=0.1*MIN(1.,NSTEP/(T_SPINUP/ZTMST))
                        ENDIF
                     ELSE
                        IF (JK.LT.14) THEN
                           ZZZH=1.E-5*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ELSEIF(JK.EQ.14) THEN
                           ZZZH=0.1*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ELSEIF(JK.EQ.15) THEN
                           ZZZH=1.*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ELSEIF(JK.EQ.16) THEN
                           ZZZH=10.*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ELSEIF(JK.EQ.17) THEN
                           ZZZH=100.*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ELSEIF(JK.EQ.18) THEN
                           ZZZH=10.*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ELSEIF(JK.EQ.19) THEN
                           ZZZH=1.*MIN(1.,NSTEP/(T_SPINUP/ZTMST))*(RG/600.)
                        ENDIF
                     ENDIF
                  ENDIF         

C     LG-  ZCFTR is the exchange coefficient of the tracers and
C     is assumed to resemble that of heat and moisture. 

                  ZCFTR(JL,JK)=ZZZH*ZCONS18*ZALF

               ENDIF

C     LG- assigning the value of ZZZH to the eddy-diffusivity for
C     heat and the windspeed, being written away to an outputfile

               ZKH(JL,JK)=ZZZH

C     LG- assigning of value of eddy diffusivity between surface layer and
C     second layer
               
               IF (JK.EQ.KLEV-1) KH30(JL)=ZKH(JL,JK)

               ZU(JK)=SQRT(PUM1(JL,JK)**2+PVM1(JL,JK)**2)

C     LG- end

               ZCDUM(JL,JK)=ZCFM(JL,JK)/ZTKESQ*SQRT(ZTKEVN(JL,JK))

C     LG- update 03-2000, GL scheme

C     GL...
            ENDDO               ! JL-LOOP 
C***  
         ENDDO                  ! JK-LOOP

C***  
C     GL...  361 CONTINUE
C     GL...***
C     GL...  372 CONTINUE
C     GL...***

C     LG- end

C     
C     *       3.8        DIFFUSION IMPLICIT COMPUTATIONS FOR TKE
C     --------- -------- ------------ --- ---
C     
         DO 380 JK=KTDIA,KLEV
            DO 381 JL=KIDIA,KFDIA
               ZEDIF(JL,JK)=ZTPFAC2*ZTKEVN(JL,JK)
 381        CONTINUE
 380     CONTINUE
C     
         DO 385 JL=KIDIA,KFDIA
            ZTCOE(JL)=(ZCDUM(JL,ITOP)+ZCDUM(JL,ITOPP1))*0.5
            ZQDP=1./(PAPM1(JL,ITOPP1)-PAPM1(JL,ITOP))
            ZDISC=1./(1.+(ZCDUM(JL,ITOP)+ZCDUM(JL,ITOPP1))*0.5*ZQDP)
            ZEBSM(JL,ITOP)=ZDISC*(ZCDUM(JL,ITOP)+ZCDUM(JL,ITOPP1))*0.5*ZQDP
            ZEDIF(JL,ITOP)=ZDISC*ZEDIF(JL,ITOP)
 385     CONTINUE
C     
         DO 386 JK=ITOPP1,KLEV-2
            DO 387 JL=KIDIA,KFDIA
               ZQDP=1./(PAPM1(JL,JK+1)-PAPM1(JL,JK))
               ZFAC=ZTCOE(JL)*ZQDP
               ZTCOE(JL)=(ZCDUM(JL,JK+1)+ZCDUM(JL,JK))*0.5
               ZDISC=1./(1.+ZFAC*(1.-ZEBSM(JL,JK-1))+(ZCDUM(JL,JK+1)+
     1              ZCDUM(JL,JK))*0.5*ZQDP)
               ZEBSM(JL,JK)=ZDISC*(ZCDUM(JL,JK+1)+ZCDUM(JL,JK))*0.5*ZQDP
               ZEDIF(JL,JK)=ZDISC*(ZEDIF(JL,JK)+ZFAC*ZEDIF(JL,JK-1))
 387        CONTINUE
 386     CONTINUE
C     
         DO 390 JL=KIDIA,KFDIA
            ZQDP=1./(PAPM1(JL,KLEV)-PAPM1(JL,KLEVM1))
            ZFAC=ZTCOE(JL)*ZQDP
            ZTCOE(JL)=(ZCDUM(JL,KLEV)+ZCDUM(JL,KLEVM1))*0.5
            ZDISC=1./(1.+ZFAC*(1.-ZEBSM(JL,KLEV-2))+(ZCDUM(JL,KLEV)+
     1           ZCDUM(JL,KLEVM1))*0.5*ZQDP)
            ZEDIF(JL,KLEVM1)=ZDISC*((ZCDUM(JL,KLEV)+ZCDUM(JL,KLEVM1))*0.5
     1           *ZQDP*ZEDIF(JL,KLEV)+ZEDIF(JL,KLEVM1)+ZFAC*ZEDIF(JL,KLEV-2))
 390     CONTINUE
C     
         DO 392 JK=KLEV-2,ITOP,-1
            DO 393 JL=KIDIA,KFDIA
               ZEDIF(JL,JK)=ZEDIF(JL,JK)+ZEBSM(JL,JK)*ZEDIF(JL,JK+1)
 393        CONTINUE
 392     CONTINUE
C     
C     --------------------------------------------------------------------
C     
C     *    TIME INTEGRATION OF TURBULENT KINETIC ENERGY AND CHECK
C     ---- ----------- -- --------- ------- ------ --- -----
C     
         DO 394 JK=ITOP,KLEV
            ZTEST=0.
            DO 395 JL=KIDIA,KFDIA
               PTKE(JL,JK)=ZEDIF(JL,JK)+ZTPFAC3*ZTKEVN(JL,JK)
               ZTEST=ZTEST+CVMGM(1.,0.,PTKE(JL,JK))
C     EVM
               DTKEDT(JL,JK)=(PTKE(JL,JK)-PTKEM1M(JL,JK))/ZTMST
 395        CONTINUE
            IF(ZTEST.NE.0.) STOP ' STOP IN EC4_VDIFF: TKE IS NEGATIVE '
 394     CONTINUE
C     
C     *    TIME FILTER FOR TURBULENT KINETIC ENERGY
C     ---- ------ --- --------- ------- ------
C     
         IF(NSTEP.NE.NSTART) THEN
            ZEPS=EPS
         ELSE
            ZEPS=0.
         ENDIF
         DO 397 JK=KTDIA,KLEV
            DO 396 JL=KIDIA,KFDIA
               PTKEM1(JL,JK)=PTKEM(JL,JK)
     *              +ZEPS*(PTKEM1M(JL,JK)-2.*PTKEM(JL,JK)+PTKE(JL,JK))
 396        CONTINUE
 397     CONTINUE
C     EVM
C     *    VERTICALLY INTEGRATED TURBULENT KINETIC ENERGY
C     ---------- ---------- --------- ------- ------
C     
         DO 1391 JL=KIDIA,KFDIA
            PTKEVI(JL)=0.
 1391    CONTINUE
         DO 1393 JK=KTDIA,KLEV
            DO 1394 JL=KIDIA,KFDIA
               PTKEVI(JL)=PTKEVI(JL)
     *              + (PAPHM1(JL,JK+1)-PAPHM1(JL,JK))*PTKEM1M(JL,JK)
 1394       CONTINUE
 1393    CONTINUE
         DO 1397 JL=KIDIA,KFDIA
            PTKEVI(JL)=PTKEVIM(JL)+PTKEVI(JL)*ZDIAGT/G
 1397    CONTINUE
C     EVM  
C     
C     -----------------------------------------------------------------
C     
C     
C     *         4.     DIFFUSION IMPLICIT COMPUTATIONS FOR MOMENTUM.
C     --------- -------- ------------ --- ---------
C     
 400     CONTINUE

C     
C     
C     *         4.1     SETTING OF RIGHT HAND SIDES.
C     
 410     CONTINUE
         DO 412 JK=ITOP,KLEV
            DO 411 JL=KIDIA,KFDIA
               ZUDIF(JL,JK)=ZTPFAC2*PUM1(JL,JK)
               ZVDIF(JL,JK)=ZTPFAC2*PVM1(JL,JK)
 411        CONTINUE
 412     CONTINUE
C     
C     *         4.2     TOP LAYER ELIMINATION.
C     
 420     CONTINUE
C     
         DO 421 JL=KIDIA,KFDIA
            ZQDP=1./(PAPHM1(JL,ITOPP1)-PAPHM1(JL,ITOP))
            ZDISC=1./(1.+ZCFM(JL,ITOP)*ZQDP)
            ZEBSM(JL,ITOP)=ZDISC*(ZCFM(JL,ITOP)*ZQDP)
            ZUDIF(JL,ITOP)=ZDISC*ZUDIF(JL,ITOP)
            ZVDIF(JL,ITOP)=ZDISC*ZVDIF(JL,ITOP)
 421     CONTINUE
C     
C     *         4.3     ELIMINATION FOR MIDDLE LAYERS.
C     
 430     CONTINUE
C     
         DO 432 JK=ITOPP1,KLEVM1
            DO 431 JL=KIDIA,KFDIA
               ZQDP=1./(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
               ZFAC=ZCFM(JL,JK-1)*ZQDP
               ZDISC=1./(1.+ZFAC*(1.-ZEBSM(JL,JK-1))+ZCFM(JL,JK)*ZQDP)
               ZEBSM(JL,JK)=ZDISC*(ZCFM(JL,JK)*ZQDP)
               ZUDIF(JL,JK)=ZDISC*(ZUDIF(JL,JK)+ZFAC*ZUDIF(JL,JK-1))
               ZVDIF(JL,JK)=ZDISC*(ZVDIF(JL,JK)+ZFAC*ZVDIF(JL,JK-1))
 431        CONTINUE
 432     CONTINUE
C     
C     *         4.4     BOTTOM LAYER ELIMINATION.
C     
 440     CONTINUE
C     
         DO 441 JL=KIDIA,KFDIA
            ZQDP=1./(PAPHM1(JL,KLEVP1)-PAPHM1(JL,KLEV))
            ZFAC=ZCFM(JL,KLEVM1)*ZQDP
            ZTCOE(JL)=ZCFM(JL,KLEV)
            ZDISC=1./(1.+ZFAC*(1.-ZEBSM(JL,KLEVM1))+ZCFM(JL,KLEV)*ZQDP)
            ZUDIF(JL,KLEV)=ZDISC*(ZUDIF(JL,KLEV)+ZFAC*ZUDIF(JL,KLEVM1))
            ZVDIF(JL,KLEV)=ZDISC*(ZVDIF(JL,KLEV)+ZFAC*ZVDIF(JL,KLEVM1))
 441     CONTINUE
C     
C     *         4.5     BACK-SUBSTITUTION.
C     
 450     CONTINUE
C     
         DO 452 JK=KLEVM1,ITOP,-1
            DO 451 JL=KIDIA,KFDIA
               ZUDIF(JL,JK)=ZUDIF(JL,JK)+ZEBSM(JL,JK)*ZUDIF(JL,JK+1)
               ZVDIF(JL,JK)=ZVDIF(JL,JK)+ZEBSM(JL,JK)*ZVDIF(JL,JK+1)
 451        CONTINUE
 452     CONTINUE
C     
C     *         4.6     INCREMENTATION OF U AND V TENDENCIES AND STORAGE OF
C     *                 THE DISSIPATION.
C     
 460     CONTINUE
C     
         DO 461 JL=KIDIA,KFDIA
            ZVIDIS(JL)=0.
 461     CONTINUE
C***  
         DO 471 JK=ITOP,KLEV
C***  
            DO 462 JL=KIDIA,KFDIA
               ZDUDT=(ZUDIF(JL,JK)-ZTPFAC2*PUM1(JL,JK))*ZCONS13
               PVOM(JL,JK)=PVOM(JL,JK)+ZDUDT
               ZDVDT=(ZVDIF(JL,JK)-ZTPFAC2*PVM1(JL,JK))*ZCONS13
               PVOL(JL,JK)=PVOL(JL,JK)+ZDVDT
               ZDIS(JL,JK)=0.5*((ZTPFAC2*PUM1(JL,JK)-ZUDIF(JL,JK))*(ZTPFAC4*
     *              PUM1(JL,JK)+ZUDIF(JL,JK))+(ZTPFAC2*PVM1(JL,JK)-
     *              ZVDIF(JL,JK))*(ZTPFAC4*PVM1(JL,JK)+ZVDIF(JL,JK)))
               PKE(JL,JK)=ZDIS(JL,JK)*ZCONS13
               ZVIDIS(JL)=ZVIDIS(JL)+ZDIS(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
 462        CONTINUE
C     
C***  
 471     CONTINUE
C     
C     --------------------------------------------------------
C     IF *LDISS_OFF* = .TRUE. TURN OFF THE DISSIPATION OF
C     KINETIC ENERGY
C     
         IF (LDISS_OFF) THEN
C     
            DO 478 JK=ITOP,KLEV
               CALL RESETR (ZDIS (KIDIA,JK),KFDIA-KIDIA+1,0.)
               CALL RESETR (PKE  (KIDIA,JK),KFDIA-KIDIA+1,0.)
 478        CONTINUE
            CALL RESETR (ZVIDIS (KIDIA),KFDIA-KIDIA+1,0.)
C     
         END IF
C***  
C     
C     *         4.8     UPDATING OF Z0 FOR OPEN SEA.
C     
 480     CONTINUE
C     
         DO  481 JL=KIDIA,KFDIA
            PAZ0(JL)=CVMGT(PAZ0M(JL),MAX(ZCONS14*ZTCOE(JL)
     *           *SQRT(ZUDIF(JL,KLEV)**2+ZVDIF(JL,KLEV)**2)
     *           *PTM1(JL,KLEV)*(1.+VTMPC1*PQM1(JL,KLEV)-PXM1(JL,KLEV))
C     EVM *        /PAPHM1(JL,KLEVP1),ZEPZZO),LALAND(JL))
     *           /PAPHM1(JL,KLEVP1),ZEPZZO),LALAND(JL).OR.LAZ0FIX)
            LO=(.NOT.LALAND(JL)).AND.(PSEAICE(JL).GT.0.5)
            PAZ0(JL)=CVMGT(CZ0ICE,PAZ0(JL),LO)
            ZTAUX=ZCONS15*ZTCOE(JL)*ZUDIF(JL,KLEV)
            ZTAUY=ZCONS15*ZTCOE(JL)*ZVDIF(JL,KLEV)
            PUSTR(JL)=PUSTRM(JL)+ZDIAGT*ZTAUX
            PVSTR(JL)=PVSTRM(JL)+ZDIAGT*ZTAUY
            ZTAU=SQRT(ZTAUX**2+ZTAUY**2)
C     PUSTAR3(JL)=PUSTAR3M(JL)+ZDIAGT*SQRT(ZTAU/ZRHOS)**3
C     
C     EVM  951107 COMPUTE USTAR ON THE BASIS OF DENSITY OF AIR
C     
            ZRHOL=PAPHM1(JL,KLEVP1)/(ZRD*PTM1(JL,KLEV)*
     +           (1.+VTMPC1*PQM1(JL,KLEV)))
            PUSTAR3(JL)=PUSTAR3M(JL)+ZDIAGT*SQRT(ZTAU/ZRHOL)**3
            PVDIS(JL)=PVDISM(JL)+ZDIAGT*ZCONS15*ZVIDIS(JL)
 481     CONTINUE



C     JHC
C     ZBUD=BUDW(IROW)
C     DVDISZ(IROW)=ZDIAGT*ZCONS15*ZBUD*SSUM(KLON,ZVIDIS,1)
C     JHC
C     
C     ------------------------------------------------------------------
C     
C     *         5.     DIFFUSION IMPLICIT COMPUTATIONS FOR HEAT (S.+L.).
C     --------- -------- ------------ --- ---- --------
C     
 500     CONTINUE
         DO 502 JK=1,KLEV
            DO 501 JL=KIDIA,KFDIA
               ZTDIF(JL,JK)=0.
               ZQDIF(JL,JK)=0.
               ZXDIF(JL,JK)=0.
 501        CONTINUE
 502     CONTINUE
C     
C     *         5.1     SETTING OF RIGHT HAND SIDES.
C     
 510     CONTINUE
         DO 512 JK=ITOP,KLEV
            DO 511 JL=KIDIA,KFDIA
               ZTDIF(JL,JK)=ZTPFAC2*ZCPTGZ(JL,JK)
               ZQDIF(JL,JK)=ZTPFAC2*PQM1(JL,JK)
               ZXDIF(JL,JK)=ZTPFAC2*PXM1(JL,JK)
 511        CONTINUE
 512     CONTINUE
C     
C     *         5.2     TOP LAYER ELIMINATION.
C     
 520     CONTINUE
C     
         DO 521 JL=KIDIA,KFDIA
            ZQDP=1./(PAPHM1(JL,ITOPP1)-PAPHM1(JL,ITOP))
            ZDISC=1./(1.+ZCFH(JL,ITOP)*ZQDP)
            ZEBSH(JL,ITOP)=ZDISC*(ZCFH(JL,ITOP)*ZQDP)
            ZTDIF(JL,ITOP)=ZDISC*ZTDIF(JL,ITOP)
            ZQDIF(JL,ITOP)=ZDISC*ZQDIF(JL,ITOP)
            ZXDIF(JL,ITOP)=ZDISC*ZXDIF(JL,ITOP)
 521     CONTINUE
C     
C     *         5.3     ELIMINATION FOR MIDDLE LAYERS.
C     
 530     CONTINUE
C     
         DO 532 JK=ITOPP1,KLEVM1
            DO 531 JL=KIDIA,KFDIA
               ZQDP=1./(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
               ZFAC=ZCFH(JL,JK-1)*ZQDP
               ZDISC=1./(1.+ZFAC*(1.-ZEBSH(JL,JK-1))+ZCFH(JL,JK)*ZQDP)
               ZEBSH(JL,JK)=ZDISC*(ZCFH(JL,JK)*ZQDP)
               ZTDIF(JL,JK)=ZDISC*(ZTDIF(JL,JK)+ZFAC*ZTDIF(JL,JK-1))
               ZQDIF(JL,JK)=ZDISC*(ZQDIF(JL,JK)+ZFAC*ZQDIF(JL,JK-1))
               ZXDIF(JL,JK)=ZDISC*(ZXDIF(JL,JK)+ZFAC*ZXDIF(JL,JK-1))
 531        CONTINUE
 532     CONTINUE
C     
C     *         5.4     BOTTOM LAYER ELIMINATION.
C     
 540     CONTINUE
C     
         DO 541 JL=KIDIA,KFDIA
            ZQDP=1./(PAPHM1(JL,KLEVP1)-PAPHM1(JL,KLEV))
            ZFAC=ZCFH(JL,KLEVM1)*ZQDP
            ZTCOE(JL)=ZCFH(JL,KLEV)
            ZDISC=1./(1.+ZFAC*(1.-ZEBSH(JL,KLEVM1))+ZCFH(JL,KLEV)*ZQDP)
            ZDISQ=1./(1.+ZFAC*(1.-ZEBSH(JL,KLEVM1))+ZCAIR(JL)*ZCFH(JL,KLEV)
     *           *ZQDP)
            ZTDIF(JL,KLEV)=ZDISC*(ZTDIF(JL,KLEV)+(ZCFH(JL,KLEV)*ZQDP)
     *           *ZTPFAC2*ZCPTS(JL)+ZFAC*ZTDIF(JL,KLEVM1))
            ZQDIF(JL,KLEV)=ZDISQ*(ZQDIF(JL,KLEV)+(ZCSAT(JL)*ZCFH(JL,KLEV)*
     *           ZQDP)*ZTPFAC2*ZQS(JL)+ZFAC*ZQDIF(JL,KLEVM1))
            ZXDIF(JL,KLEV)=ZDISC*(ZXDIF(JL,KLEV)+ZFAC*ZXDIF(JL,KLEVM1))
 541     CONTINUE
C     *         5.5     BACK-SUBSTITUTION.
C     
 550     CONTINUE
C     
         DO 552 JK=KLEVM1,ITOP,-1
            DO 551 JL=KIDIA,KFDIA
               ZTDIF(JL,JK)=ZTDIF(JL,JK)+ZEBSH(JL,JK)*ZTDIF(JL,JK+1)
               ZQDIF(JL,JK)=ZQDIF(JL,JK)+ZEBSH(JL,JK)*ZQDIF(JL,JK+1)
               ZXDIF(JL,JK)=ZXDIF(JL,JK)+ZEBSH(JL,JK)*ZXDIF(JL,JK+1)
 551        CONTINUE
 552     CONTINUE
C     
C     *         5.6     INCREMENTATION OF T AND Q TENDENCIES.
C     
 560     CONTINUE
C     
C***  
         DO 571 JK=ITOP,KLEV
C***  
            DO 561 JL=KIDIA,KFDIA
               ZQDIF(JL,JK)=ZQDIF(JL,JK)+ZTPFAC3*PQM1(JL,JK)
               ZDQDT=(ZQDIF(JL,JK)-PQM1(JL,JK))*ZCONS13
               PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
               ZTDIF(JL,JK)=ZTDIF(JL,JK)+ZTPFAC3*ZCPTGZ(JL,JK)
               ZDTDT=((ZTDIF(JL,JK)+ZDIS(JL,JK)-PGEOM1(JL,JK))
     *              /(CPD*(1.+VTMPC2*ZQDIF(JL,JK)))-PTM1(JL,JK))*ZCONS13
               PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
               ZXDIF(JL,JK)=ZXDIF(JL,JK)+ZTPFAC3*PXM1(JL,JK)
               ZDXMDT=(ZXDIF(JL,JK)-PXM1(JL,JK))*ZCONS13
               PXTE(JL,JK)=PXTE(JL,JK)+ZDXMDT
 561        CONTINUE
C     
C***  
 571     CONTINUE
         
C     -------------------------------------------------------------------
C     LG- Vertical transport of tracers, the code has been rewritten in
C     order to account for the vertical exchange within the canopy and
C     between the canopy and the surface layer/PBL if the integrated 
C     biosphere model is used. The vertical transport of tracers is 
C     not calculated for the first timestep since for the calculation of 
C     the tracer exchange coefficients some parameters are required which
C     are initialized in the chemistry scheme which is called after 
C     EC4_VDIFF.f  
C     -------------------------------------------------------------------

         IF (LXTVDIFF) THEN

            DO JK=1,KLEVEL
               DO JL=1,NLON 

                  DO JT=1,KTRAC
                     ZXTM1(JL,JK,JT)=PXTM1(JL,JK,JT)

C     LG-     originally the whole NOx group is transported as a group, however, 
C     since for the biosphere model, the NO emission is treated seperately
C     this tracer is transported itself and is therefore not considered
C     in the vertical exchange of the NOx group

c     IF (JT.EQ.inox) ZXTM1(JL,JK,JT)=PXTM1(JL,JK,JT)-PXTM1(JL,JK,ino)-
c     *     PXTM1(JL,JK,ino2)-PXTM1(JL,JK,ino3)-PXTM1(JL,JK,in2o5)-
c     *     PXTM1(JL,JK,ihno4)

                  ENDDO

               ENDDO
            ENDDO

C     LG- writing of Kh and profile to check profile
            
            IF (NSTEP.EQ.0) THEN
               OPEN(UNIT=NUNKH,FILE='/data/ganzevl/racmo/output/Kh.out',
     *              STATUS='UNKNOWN')
               WRITE(NUNKH,'(2a)')'Eddy diffusivity for heat, windspeed, ',
     *              'bulk Richardson numbers, and the gust function'
               WRITE(NUNKH,'(a5,4a12)')
     *              'level','KH [m2 s-1]','U [m s-1]','Rib','ZCFTR/DT'

               NLEV1=1
               NLEV2=NLEVATM
               NLEVPR=NLEV2-NLEV1+1

               WRITE (NUNKH,*) NSTOP,NLEVPR,NPRINT,ZTMST
               WRITE (NUNKH,*) NLEV1,NLEV2
               WRITE (NUNKH,*) (HGHT(JK),JK=NLEV2,NLEV1,-1)
            ENDIF
            
            IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN

               WRITE(NUNKH,'(a14)') LDATLTIME
               WRITE(NUNKH,*)NSTEP,JDAY,GMT,LTIME,HC,GUST

               DO JK=NLEVATM,1,-1
                  ZU_RATIO=0.
                  IF (ZU(NLEVATM).GT.0.) THEN
                     ZU_RATIO=ZU(JK)/ZU(NLEVATM)
                  ENDIF
                  IF (JK.EQ.NLEVATM) THEN
                     ZU(JK)=SQRT(PUM1(1,JK)**2+PVM1(1,JK)**2)
                  ENDIF
                  WRITE(NUNKH,'(I4.4,1X,5E12.3)')
     &                 JK,ZKH(JPHR,JK),ZU(JK),ZU_RATIO,
     &                 ZRIB(JPHR,JK),ZCFTR(JPHR,JK)/ZTMST
               ENDDO

            ENDIF

            IF (NSTEP.EQ.NSTOP) CLOSE(NUNKH)
            
C     LG- end writing Kh and windspeed profiles
            
C     ---------------------------------------------------------------------
C     LG- Vertical tracer transport is extended with the exchange within the 
C     the canopy and between the canopy and the surface layer/PBL. The 
C     explicit treatment of the turbulent exchange within the canopy 
C     considering the local and local contribution in the overall exchange
C     still needs to be incorporated (18-08-97). The exchange can be 
C     determined by some representative K profiles for some
C     specific conditions such as the nocturnal and daytime K profiles
C     of the ABLE experiments
C     ------------------------------------------------------------------
C     
C     *         5.     DIFFUSION IMPLICIT COMPUTATIONS FOR TRACERS
C     --------- -------- ------------ --- ---- --------
C     

 599        CONTINUE

            DO 604 JT=1,KTRAC
               DO 603 JK=1,KLEVEL
                  DO 602 JL=KIDIA,KFDIA
                     ZXTDIF(JL,JK,JT)=0.
 602              CONTINUE
 603           CONTINUE
 604        CONTINUE

            IF (LCHEM) THEN     !only when LCHEM=.TRUE.
C     
C     *         5.1     SETTING OF RIGHT HAND SIDES.
C     
 610           CONTINUE
               DO 613 JT=1,KTRAC
                  DO 612 JK=ITOP,KLEVEL
                     DO 611 JL=KIDIA,KFDIA
                       ! op_ck_20031001+
                       zxtvn(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
C                        ZXTDIF(JL,JK,JT)=ZTPFAC2*ZXTM1(JL,JK,JT)
                       zxtdif(jl,jk,jt)=ztpfac2*zxtvn(jl,jk,jt)
                       ! op_ck_20031001-
 611                 CONTINUE
 612              CONTINUE
 613           CONTINUE

C     
C     *         5.2     TOP LAYER ELIMINATION.
C     
 620           CONTINUE
C     
               DO 622 JT=1,KTRAC
                  DO 621 JL=KIDIA,KFDIA
                     ZTCOETR(JL,JT)=ZCFTR(JL,ITOP)
                     ZQDP=1./(ZAPHM1(JL,ITOPP1)-ZAPHM1(JL,ITOP))
                     ZDISC=1./(1.+ZCFTR(JL,ITOP)*ZQDP)
                     ZEBSH(JL,ITOP)=ZDISC*(ZCFTR(JL,ITOP)*ZQDP)
                     ZXTDIF(JL,ITOP,JT)=ZDISC*ZXTDIF(JL,ITOP,JT)
 621              CONTINUE
 622           CONTINUE
C     
C     *         5.3     ELIMINATION FOR MIDDLE LAYERS.
C     
 630           CONTINUE
C     

               DO 633 JT=1,KTRAC
                  DO 632 JK=ITOPP1,KLEVELM1
                     DO 631 JL=KIDIA,KFDIA
                        ZQDP=1./(ZAPHM1(JL,JK+1)-ZAPHM1(JL,JK))
                        ZFAC=ZTCOETR(JL,JT)*ZQDP
                        ZTCOETR(JL,JT)=ZCFTR(JL,JK)
                        ZDISC=1./(1.+ZFAC*(1.-ZEBSH(JL,JK-1))+ZCFTR(JL,JK)*ZQDP)
                        ZEBSH(JL,JK)=ZDISC*(ZCFTR(JL,JK)*ZQDP)
                        ZXTDIF(JL,JK,JT)=ZDISC*(ZXTDIF(JL,JK,JT)+

C     LG- modification by considering the emission flux within the routine
C     of vertical turbulent transport to study the impact of the splitting on
C     the resolved concentrations

     *                       ZTMST*G*ZQDP*EM(JL,JK,JT)*ZTPFAC2       
     *                       +ZFAC*ZXTDIF(JL,JK-1,JT))

C     LG- determining the emission tendency, which is determined seperately and
C     distracted from the vertical diffusion tendency in the subroutine
C     write_oned.f, to calculate the actual vertical diffusion flux. The
C     ABS term is added since the effective canopy top flux EM/DD can also
C     be negative for cases with a stronger sink than a source canopy layers,
C     since it reflects both the net emission and dry deposition fluxes.

                        IF (LBIOSPH.AND.ABS(EM(JL,KLEVEL,JT)).GT.0.) 
     *                       PXTEEMIS(JK,JT)=G*ZQDP*EM(JL,JK,JT)*ZTPFAC2

C     LG- end

 631                 CONTINUE
 632              CONTINUE
 633           CONTINUE

C     
C     *         5.4     BOTTOM LAYER ELIMINATION.
C     
 640           CONTINUE
C     

               DO 642 JT=1,KTRAC
                  DO 641 JL=KIDIA,KFDIA
                     ZQDP=1./(ZAPHM1(JL,KLEVELP1)-ZAPHM1(JL,KLEVEL))
                     ZFAC=ZTCOETR(JL,JT)*ZQDP
                     ZTCOETR(JL,JT)=ZCFTR(JL,KLEVEL)
                     ZDISXT=1./(1.+ZFAC*(1.-ZEBSH(JL,KLEVELM1)))
                     ZXTDIF(JL,KLEVEL,JT)=ZDISXT*(ZXTDIF(JL,KLEVEL,JT)+

C     LG-    the parameter ZXTEMIS(JL,JT) has been replaced by the 
C     parameter EM/DD(JL,JK,JT)

     *                    ZTMST*G*ZQDP*EM(JL,KLEVEL,JT)*ZTPFAC2
     *                    +ZFAC*ZXTDIF(JL,KLEVELM1,JT))

C     LG- determining the emission tendency, which is determined seperately and
C     distracted from the vertical diffusion tendency in the subroutine
C     write_oned.f, to calculate the actual vertical diffusion flux

                     IF (LBIOSPH.AND.ABS(EM(JL,KLEVEL,JT)).GT.0.) 
     *                    PXTEEMIS(KLEVEL,JT)=G*ZQDP*EM(JL,KLEVEL,JT)*ZTPFAC2

C     LG- end

 641              CONTINUE
 642           CONTINUE

C     *         5.5     BACK-SUBSTITUTION.
C     
 650           CONTINUE
C     

               DO 653 JT=1,KTRAC
                  DO 652 JK=KLEVELM1,ITOP,-1
                     DO 651 JL=KIDIA,KFDIA
                        ZXTDIF(JL,JK,JT)=ZXTDIF(JL,JK,JT)+ZEBSH(JL,JK)*ZXTDIF(JL,JK+1,JT)
 651                 CONTINUE
 652              CONTINUE
 653           CONTINUE

               DO 677 JT=1,KTRAC

C     LG- no vertical diffusion of NOx since all the species are being transported
C     separately for short timesteps (ino.LT.NTRAC)

C     ==============================================================================
C     LG- uncommented this statement removes the vertical turbulent mixing
C     in the atmospheric column

c     GOTO 677

C     ===============================================================================

                  IF (JT.EQ.inox.AND.KTRAC.GT.ino.AND.ino.GT.1) GOTO 677

C     LG- determining the mass conservation within the vertical diffusion routine,
C     initialisation of the accumulators 

                  TOTMASS_OLD=0.
                  TOTMASS_NEW=0.

                  DO 675 JK=ITOP,KLEVEL
                     DO 673 JL=KIDIA,KFDIA

! op_ck_20031001+

!                         ZXTDIF(JL,JK,JT)=ZXTDIF(JL,JK,JT)+ZTPFAC3*ZXTM1(JL,JK,JT)
!                         ZDXTDT=(ZXTDIF(JL,JK,JT)-ZXTM1(JL,JK,JT))*ZCONS13

                        zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt) +                
     &                             ztpfac3*zxtvn(jl,jk,jt) 
                        zdxtdt=(zxtdif(jl,jk,jt)-zxtvn(jl,jk,jt))*zcons13 
! op_ck_20031001-

C     LG- determining the mass conservation within the vertical diffusion routine,
C     the original mass of the column is ZXTM1*GRMASS intergrated over the
C     the column whereas the new mass of the column consists of the ZXTM1
C     plus the sum of the tracer tendency which should be zero, intergrated
C     over the whole column. Testing this routine on mass conservation that
C     is sometimes some problem with mass conservation, especially short after
C     the initialisation of the vertical profiles within the canopy, e.g.
C     for isoprene due to the strong gradients that exists initially. For
C     other trace gases with a smoother profile, e.g. ozone it turns out
C     that there is mass conservation within this routine (dM <1e-6 %) 
C     
                        TOTMASS_OLD=TOTMASS_OLD+(ZXTM1(JL,JK,JT)+PXTTE(JL,JK,JT)*ZTMST)*
     *                       GRMASS(JL,JK)

C     LG- end

                        PXTTE(JL,JK,JT)=PXTTE(JL,JK,JT)+ZDXTDT

C     LG- calculating the new mass accounting for the vertical diffusion tendency

                        TOTMASS_NEW=TOTMASS_NEW+(ZXTM1(JL,JK,JT)+PXTTE(JL,JK,JT)*ZTMST)*
     *                       GRMASS(JL,JK)
                        IF (JK.EQ.KLEVEL) THEN
                           DMASS=100.*ABS((TOTMASS_NEW-TOTMASS_OLD)/TOTMASS_NEW)

C     LG-  when the emission/dry deposition flux is considered within this 
C     routine then there is obviously an increase/decrease of mass within
C     this routine!

                           IF (DMASS.GT.0.5.AND.IEM_TURB(JT).EQ.0) THEN
                              WRITE(*,'(1a,i3)')
     *                             ' More than 0.5% mass loss in EC4_VDIFF.f for tracer no.: ',JT
                              WRITE(*,'(1a,f10.4)')' Mass loss in [%]: ',DMASS
                           ENDIF
                        ENDIF

C     LG- end

 673                 CONTINUE
 675              CONTINUE
 677           CONTINUE
            ENDIF

C     LG- end tracer vertical turbulent transport

C     LG- end IF (LCHEM)

         ENDIF

C     
C     *         5.8     STORAGE OF THE SURFACE HEAT (S.+L.) AND MOISTURE
C     *         FLUXES AND THEIR DERIVATIVES
C     *         AGAINST SURFACE VARIABLES
C     
 580     CONTINUE
C     
         DO 581 JL=KIDIA,KFDIA
            ZCOEFF=ZCONS15*ZTCOE(JL)
C     
C     OBC  Multiplied by zqs on both sides:
            LO=(ZQS(JL)*ZHUM(JL)).LE.PQM1(JL,KLEV)
            ZHUM(JL)=CVMGT(0.,ZHUM(JL),LO)
            ZCA=CVMGT(0.,1.,LO)
            ZHUM(JL)=CVMGT(ZHUM(JL),1.,LALAND(JL))
            ZCA=CVMGT(ZCA,1.,LALAND(JL))
            LO=PQM1(JL,KLEV).GT.ZQS(JL)
            ZHUM(JL)=CVMGT(1.,ZHUM(JL),LO)
            ZCA=CVMGT(1.,ZCA,LO)
            ZHUM(JL)=(1.-PCVS(JL))*(1.-PCVW(JL))*ZHUM(JL)
            ZCA=(1.-PCVS(JL))*(1.-PCVW(JL))*ZCA
C     
            ZQNLEV=ZQDIF(JL,KLEV)-ZTPFAC3*PQM1(JL,KLEV)
            ZZQS=ZTPFAC2*ZQS(JL)
            PQHFL(JL)=ZCOEFF*(ZCAIR(JL)*ZQNLEV-ZCSAT(JL)*ZZQS)
C     
            ZTNLEV=ZTDIF(JL,KLEV)-ZTPFAC3*ZCPTGZ(JL,KLEV)
            ZZCPTS=ZTPFAC2*ZCPTS(JL)
            PTHFL(JL)=ZCOEFF*(ZTNLEV-ZZCPTS)
            PDHFT(JL)=-ZCONS16*PQHFL(JL)

C     LG- replacing the surface temperature by the skin temperature

            IF (LTVEG) THEN
               PTHFL(JL)=PTHFL(JL)+(1.-(VEGFRAC(JL)+WSFRAC(JL)))*
     &              PTSM1M(JL)*PDHFT(JL)+
     &              (VEGFRAC(JL)+WSFRAC(JL))*TVEGM(JL)*PDHFT(JL)
            ELSE
               PTHFL(JL)=PTHFL(JL)+PTSM1M(JL)*PDHFT(JL)
            ENDIF

C     
            ZXNLEV=ZXDIF(JL,KLEV)-ZTPFAC3*PXM1(JL,KLEV)
            ZXHFL=ZCOEFF*ZXNLEV
            PRSFL(JL)=MAX(0.,ZXHFL)
            PXHFL(JL)=MIN(ZXHFL,0.)
C     
            PAHFS(JL)=PAHFSM(JL)+ZDIAGT*PTHFL(JL)
            PEVAP(JL)=PEVAPM(JL)+ZDIAGW*(PQHFL(JL)+PXHFL(JL))
C     

            PDHFQW(JL)=ZCOEFF*ZWLMXI(JL)*(1.-PCVS(JL))*(ZQNLEV-ZZQS)

C     LG- introducing a correction of the wet skin flux to study the sensitivity
C     of the simulated wet skin evaporation for the wet skin evaporation
C     rate using a wet skin specific moisture gradient compared to using
C     the grid average values (see spreadsheets for the polynomial fit of
C     the dq-ws/dqavg dependence on the wet skin fraction

C     DQWS_DQAVG=37.915*WSFRAC(JL)**6-95.487*WSFRAC(JL)**5+
C     &           90.636*WSFRAC(JL)**4-39.428*WSFRAC(JL)**3+
C     &           7.7496*WSFRAC(JL)**2-0.455*WSFRAC(JL)+
C     &           0.0678
C     
C     PDHFQW(JL)=ZCOEFF*ZWLMXI(JL)*(1.-PCVS(JL))*(ZQNLEV-ZZQS)*DQWS_DQAVG
C     print *,'ec4_vdiff, LE, line 4278',wsfrac(jl),
C     & ZCOEFF*ZWLMXI(JL)*(1.-PCVS(JL))*(ZQNLEV-ZZQS),
C     & pdhfqw,DQWS_DQAVG

C     LG- end  
C     
            PDHFQS(JL)=ZCOEFF*(ZQNLEV-ZZQS)
C     
            ZWET(JL)=ZWET(JL)-PCVS(JL)-(1.-PCVS(JL))*PCVW(JL)
            ZQHFLV=PVGRAT(JL)*ZWET(JL)*(ZQNLEV-ZZQS)
            ZQHFLB=(1.-PVGRAT(JL))*(ZCA*ZQNLEV-ZHUM(JL)*ZZQS)
            PCVGHL(JL)=CVMGT(ABS(ZQHFLV)/MAX(ZEPEVAP,ABS(ZQHFLV)+
     *           ABS(ZQHFLB)),1.,LALAND(JL))
 581     CONTINUE
C     
C     JCH  DHFSZ(IROW) =ZDIAGT*ZBUD*SSUM(KLON,PTHFL,1)
C     JCH  DEVAPZ(IROW)=ZDIAGW*ZBUD*SSUM(KLON,PQHFL,1)
C     
         DO 582 JL=KIDIA,KFDIA
            PQHFL(JL)=PQHFL(JL)-PCVS(JL)*PDHFQS(JL)
            PAHFL(JL)=ALV*PQHFL(JL)+ALS*PCVS(JL)*PDHFQS(JL)
            PTHFL(JL)=PTHFL(JL)+PAHFL(JL)
            PAHFL(JL)=PAHFLM(JL)+ZDIAGT*PAHFL(JL)
            LO=LALAND(JL).OR.PTSM1M(JL).GT.TMELT

C     LG- replacing the surface temperature by the skin temperature

            IF (LTVEG)
     &           LO=LALAND(JL).OR.TVEGM(JL).GT.TMELT

            ZALVS=ALS*PCVS(JL)+ALV*(1.-PCVS(JL))
            ZALVS=CVMGT(ZALVS,ALS,LO)

C     LG- replacing the surface temperature by the skin temperature

            IF (LTVEG) THEN
               PDHFT(JL)=(1.-(VEGFRAC(JL)+WSFRAC(JL)))*
     &              (-ZTPFAC2*ZCONS15*ZTCOE(JL)*(ZCPTS(JL)/PTSM1M(JL)+
     &              (ZALVS-ZCONS16*PTSM1M(JL))*ZCSAT(JL)*ZDQS(JL))+PDHFT(JL))+
     &              (VEGFRAC(JL)+WSFRAC(JL))*
     &              (-ZTPFAC2*ZCONS15*ZTCOE(JL)*(ZCPTS(JL)/TVEGM(JL)+
     &              (ZALVS-ZCONS16*TVEGM(JL))*ZCSAT(JL)*ZDQS(JL))+PDHFT(JL))

            ELSE
               PDHFT(JL)=-ZTPFAC2*ZCONS15*ZTCOE(JL)*(ZCPTS(JL)/PTSM1M(JL)+
     &              (ZALVS-ZCONS16*PTSM1M(JL))*ZCSAT(JL)*ZDQS(JL))+PDHFT(JL)

            ENDIF

 582     CONTINUE
C     
C     *         5.85     COMPUTATION OF BOUNDARY LAYER HEIGHT
C     
         DO 585 JL=KIDIA,KFDIA
            ZDU2=MAX(ZEPDU2,PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)
            ZRHO=PAPHM1(JL,KLEVP1)/ZRD/
     &           (PTM1(JL,KLEV)*(1.+VTMPC1*PQM1(JL,KLEV)))
            ZSENKF(JL)=-(PAHFS(JL)-PAHFSM(JL))/(ZDIAGT*ZRHO*ZCPD)
            ZLATKF(JL)=-(ALV*PQHFL(JL)+ALS*PCVS(JL)*PDHFQS(JL))/(ZRHO*ALV)
            ZUSTAR1(JL)=((PUSTAR3(JL)-PUSTAR3M(JL))/ZDIAGT)**(1./3.)
            ZCDH(JL)=ZCFH(JL,KLEV)/(ZCONS12*SQRT(ZDU2)*ZRD*ZRHO)
            ZCDM(JL)=ZCFM(JL,KLEV)/(ZCONS12*SQRT(ZDU2)*ZRD*ZRHO)
C     
C     EMPTY SOME WORKSPACE
C     
            ZPBLH(JL)=0.
            ZOBUKL(JL)=0.
            ZBUOYPR(JL)=0.
 585     CONTINUE
C     

         CALL EC_PBLHGHT(KLON,KLEV,KIDIA,KFDIA,
     &        ZTETA1,PQM1,PGEOM1,PUM1,PVM1,ZTESS,
     &        ZSENKF,ZLATKF,ZUSTAR1,
     &        ZCDN,ZCDH,ZCDM,
     &        ZPBLH,ZOBUKL,ZBUOYPR)
C     
         DO 587 JL=KIDIA,KFDIA
            PBLH  (JL) = PBLHM  (JL) + ZDIAGT*ZPBLH  (JL)
 587     CONTINUE

C     LG- see XTTROPO for the assignment of different compartments
C     in the vertical for the 1-D model

C     LG- assignment of PBL height to parameter KPBLHE

         DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
               PBLHGT(JL)=ZPBLH(JL)
               Z(JL,JK)=PGEOM1(JL,JK)/G
               IF (ZPBLH(JL).GE.Z(JL,JK)) THEN
                  IZPBLH(JL)=JK-1
                  GOTO 138 
               ENDIF
            ENDDO
         ENDDO
 138     CONTINUE

C     LG- calculation of level in meters height from the PBL height
C     calculated in the original version of ECHAM (see loop 344)
C     and a recalculation of PBL height, as calculated by subroutine
C     EC_PBLHGHT, to specific level
         
         DO 139 JL=KIDIA,KFDIA
            HPBL(JL)=Z(JL,IHPBL(JL))

C     LG-  KPBLHE has been assigned the level of the PBL height calculated
C     by the subroutine EC_PBLHGT since it appears that the calculated
C     PBL height resembles more realistically the separation between the PBL
C     and the free troposphere compared to the PBL height calculated from the
C     Ekman Layer height and Convective PBL height, KPBLHE is used to
C     distinguish the PBL from the free troposhere for the budget calculations
C     by determining value of the parameter IWHERE as a function of the level. 

            KPBLHE(JL)=IZPBLH(JL)

C     LG- 20031118+ added the calculation of the PBL height which reflects
C     only the increasing height tendency and the residual PBL height

            ZPBLHPLUS=-999.99
	    ZPBLH_RESID=-999.99

	    IF (ZPBLH(JL).GE.ZPBLHOLD) 
     &         ZPBLHPLUS=ZPBLH(1)

	    IF (ZPBLH(1).LT.ZPBLHMAX) THEN
	       ZPBLH_RESID=ZPBLHMAX 
	    ENDIF

C     LG- resetting the maximum boundary layer height for sunrise
            IF (ZRGOLD.LT.1.E-5.AND.RG.GT.0) IRESET=1
            IF (IRESET.AND.ZPBLH(1).GT.ZPBLHOLD) THEN
              ZPBLHMAX=-999.99
              IRESET=0
	    ENDIF
	    IF (RG.GT.0.AND.ZPBLH(1).GT.ZPBLHMAX) ZPBLHMAX=ZPBLH(1)

            ZPBLHOLD=ZPBLH(1)
            ZRGOLD=RG
 139     ENDDO  

C     LG- writing of PBL height to output file

         IF (NSTEP.EQ.0) THEN
            OPEN(UNIT=NUNPBLH,FILE='/data/ganzevl/racmo/output/PBLH.out',
     *           STATUS='UNKNOWN')
            WRITE(NUNPBLH,'(1a)')'The calculated PBL height'
            WRITE(NUNPBLH,'(2a)')'Prognostic PBL height (ECHAM4) and ',
     *           'diagnostic PBL height (EC_PBLHGT) (see ec4_vdiff.f)'     
            WRITE(NUNPBLH,'(a10,a15,4a10)')'istep','time',
     *        'ECHAM4','EC_PBLHGT','EC_dz+','EC_RESID'
         ENDIF

         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) 
     *        WRITE(NUNPBLH,'(1x,i9.9,1x,a14,4f10.2)')
     *        NSTEP,LDATLTIME,HPBL(1),ZPBLH(1),ZPBLHPLUS,ZPBLH_RESID
         
         IF (NSTEP.EQ.NSTOP) CLOSE(NUNPBLH)

C     LG- writing of the Kh at the PBL top to calculate entrainment 
C     fluxes from the free troposphere-PBL gradient. Actually the
C     Kh is written away for the reference height directly above and
C     under the diagnostically derived PBL height. An interpolated 
C     Kh is then estimated using a linear interpolation between the
C     two reference heights 
         
         IF (NSTEP.EQ.0) THEN
            OPEN(UNIT=NUNKHPBL,FILE='/data/ganzevl/racmo/output/Kh_PBLtop.out',
     *           STATUS='UNKNOWN')
            WRITE(NUNKHPBL,'(2a)')'Eddy diffusivity for heat [m-2 s-1]',
     *           ' at the top of the Planetary Boundary layer'
            WRITE(NUNKHPBL,'(4a6,6a15)')
     *           'JDAY','GMT','lon','lat','Kh > PBL top','z[m] > PBL top',
     *           'Kh < PBL top','z[m] < PBL top','PBL height [m]','intpol. Kh'

            NLEV1=1
            NLEV2=NLEVATM+NLEVVEG
            NLEVPR=NLEV2-NLEV1+1

         ENDIF
         
         IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0) THEN
            ZDZ=Z(1,KPBLHE(1))-Z(1,KPBLHE(1)+1)
            WRITE(NUNKHPBL,'(2i6,2f6.2,6f15.2)')
     *           JDAY,GMT,-LONG/API*180.,LAT/API*180.,
     *           ZKH(1,KPBLHE(1)),Z(1,KPBLHE(1)),
     *           ZKH(1,KPBLHE(1)+1),Z(1,KPBLHE(1)+1),ZPBLH(1),
     *           (1.-(Z(1,KPBLHE(1))-ZPBLH(1)))/ZDZ*ZKH(1,KPBLHE(1))+
     *           (1.-(ZPBLH(1)-Z(1,KPBLHE(1))))/ZDZ*ZKH(1,KPBLHE(1)+1)
         ENDIF

         IF (NSTEP.EQ.NSTOP) CLOSE(NUNKHPBL)
         
C     LG- end writing Kh and windspeed profiles

C     - - -
C     Determine localization budget bit-wise:
C     Bit 0 0:Biosphere    1:Planetary Boundary Layer
C     Bit 1 0:Troposphere  1:Stratosphere

         DO JK=1,KLEVEL
            DO JL=KIDIA,KFDIA
               IPBL=0
               ITROP=0
               ISTRAT=0
               IF (JK.LE.KLEV) IPBL=1
               IF (JK.LT.KPBLHE(JL)) ITROP=1
               IF (JK.LT.KTRHE(JL)) ISTRAT=1
               IWHERE(JL,JK)=1+1*IPBL+1*ITROP+1*ISTRAT
            ENDDO
         ENDDO

C     LG- end

C     
C     *       5.95  COMPUTE NEW T2MAX AND MIN
C     
 595     CONTINUE
C     
         
         DO 597 JL=KIDIA,KFDIA
            LO1=ZRICLS(JL).GE.0.
            ZRAT=ZHTQ/PGEOM1(JL,KLEV)
            ZCBN=LOG(1.+(EXP(ZBN(JL))-1.)*ZRAT)
            ZCBS=-(ZBN(JL)-ZBH(JL))*ZRAT

            ! mz_lg_20051216+
            IF ((EXP(ZBN(JL)-ZBH(JL))-1.)*ZRAT.LT.-1) THEN
              print *,'ec4_vdiff: line 4788, conflict calculation T2m'
              print *,'value coeff. in ZCBU',(EXP(ZBN(JL)-ZBH(JL))-1.)*zrat
            ENDIF
            ! mz_lg_20051216-

            ZCBU=-LOG(1.+MAX((EXP(ZBN(JL)-ZBH(JL))-1.)*ZRAT,0.99))
            ZRED=(ZCBN+CVMGT(ZCBS,ZCBU,LO1))/ZBH(JL)
            ZH2M=ZCPTS(JL)+ZRED*(ZCPTGZ(JL,KLEV)-ZCPTS(JL))
            ZT2=(ZH2M-ZHTQ)/(CPD*(1.+VTMPC2*PQM1(JL,KLEV)))
            PTEMP2(JL)=PTEMP2M(JL)+ZDIAGT*ZT2

C     LG- assigning the value of the 2 m meter temperature to T2M(JL)

            T2M(JL)=ZT2

C     LG- end      

            PT2MAX(JL)=MAX(PT2MAXM(JL),ZT2)
            PT2MIN(JL)=MIN(PT2MINM(JL),ZT2)
C     
C     *          5.96   2M DEW POINT
C     
            IT=PTM1(JL,KLEV)*1000.
            ZQS1=TLUCUA(IT)/PAPM1(JL,KLEV)
            ZQS1=ZQS1/(1.-VTMPC1*ZQS1)
            ZRH2M=MAX(ZEPHUM,PQM1(JL,KLEV)/ZQS1)
C     
            LO=ZT2.GT.TMELT
            ZCVM3=CVMGT(C3LES,C3IES,LO)
            ZCVM4=CVMGT(C4LES,C4IES,LO)
            ZAPH2M=PAPHM1(JL,KLEVP1)*
     *           (1.-ZHTQ/(RD*ZT2*(1.+VTMPC1*PQM1(JL,KLEV))))
            IT=ZT2*1000.
            ZQS2=TLUCUA(IT)/ZAPH2M
            ZQS2=ZQS2/(1.-VTMPC1*ZQS2)
            ZQ2M=ZRH2M*ZQS2
            ZFRAC=LOG(ZAPH2M*ZQ2M/(C2ES*(1.+VTMPC1*ZQ2M)))/ZCVM3
            PDEW2(JL)=PDEW2M(JL)+ZDIAGT*MIN(ZT2,(TMELT-ZFRAC*ZCVM4)
     *           /(1.-ZFRAC))

C     
C     *          5.97   10M WIND COMPONENTS, MAX 10M WINDSPEED
C     
            ZRAT=ZHUV/PGEOM1(JL,KLEV)
            ZCBN=LOG(1.+(EXP(ZBN(JL))-1.)*ZRAT)
            ZCBS=-(ZBN(JL)-ZBM(JL))*ZRAT
            ZCBU=-LOG(MAX(1.E10,1.+(EXP(ZBN(JL)-ZBM(JL))-1.)*ZRAT))
            ZRED=(ZCBN+CVMGT(ZCBS,ZCBU,LO1))/ZBM(JL)
            ZU10=ZRED*PUM1(JL,KLEV)
            ZV10=ZRED*PVM1(JL,KLEV)
            ZSPEED=SQRT(ZU10**2+ZV10**2)

C     LG- assigning the value of the 10 meter temperature to W10M(JL)

            WS10M(JL)=ZSPEED

C     LG- end  

            PU10(JL)=PU10M(JL)+ZDIAGT*ZU10
            PV10(JL)=PV10M(JL)+ZDIAGT*ZV10
            PWIMAX(JL)=MAX(PWIMAXM(JL),ZSPEED)
            PWIND10(JL)=PWIND10M(JL)+ZDIAGT*ZSPEED
 597     CONTINUE

C=============================================================
C     EVM  --------------------------------------------------------
C     IF *LSRFFLUX_READ* = .TRUE. READ THE SURFACE FLUXES
C     FOR SENSIBLE AND LATENT HEAT IN *INIFLUX*
C     
         IF (LSRFFLUX_READ) THEN
            IF (LSURF) THEN
C     
C     PUT LSURF TO FALSE
C     
               LSURF = .FALSE.
               WRITE (0,'(A)') ' LSURF PUT TO .FALSE. IN *EC4_VDIFF* '
            END IF
C     
C     ----------------------------------------------
C     PUT DRAG COEFFICIENT REQUIRED TO DETERMINE
C     SURFACE TEMPERATURE AT NEXT TIME STEP INTO *ZCFQSRF*
C     IF  WRIH<=1. (  STABLE)  THEN ZCFHSRF = ZCDN
C     IF  WRIH> 1. (UNSTABLE)  THEN ZCFHSRF = ZCDN*WRIH
C     
            DO 333 JL=KIDIA,KFDIA
               IF (WRIH(JL).GT.1.) THEN
                  ZCFHSRF(JL)=ZCDN(JL)*WRIH(JL)
               ELSE
                  ZCFHSRF(JL)=ZCDN(JL)
               END IF
 333        CONTINUE
C     
            IRFLX1=1
            IRFLX0=1
            IF (NSTEP.EQ.0) IRFLX0=0
            DO IRFL=IRFLX0,IRFLX1
C     
               CALL INIFLUX(
     +              NUNINFLX,   ! nuninflx,
     +              KLON   ,    ! nhor,
     +              KLEV   ,    ! nlev,
     +              KIDIA  ,    ! kstart,
     +              KFDIA  ,    ! kstop,
     +              NSTEP+I,    ! nstep,
     +              NSTOP  ,    ! nstop,
     +              TWODT  ,    ! twodt,
     +              LTS_OBS,    ! ltsobs,
     +              CREFLUX,    ! creflux,
C     
     +              PTM1   ,    ! t,
     +              PQM1   ,    ! q,
     +              PUM1   ,    ! u,
     +              PVM1   ,    ! v,
     +              PGEOM1 ,    ! gpot,
     +              ZDPH   ,    ! dph,
     +              APZERO ,    ! apzero,
     +              PAPHM1 ,    ! pres,
C     
     +              LALAND ,    ! loland,
C     
     +              PTSM1M ,    ! tsold, current surface temp!!
     +              PTSM   ,    ! tsobs, new calculated surf. temp.
     +              PWSM1M ,    ! wsold, current surface moisture
     +              PWSM   ,    ! wsupd,    "
C     
     +              ZSENKF ,    ! senf,
     +              ZLATKF ,    ! latf,
     +              ZCFHSRF     ! cdragh,
C     
     +              )
C     
            ENDDO 
C     
C     DETERMINE AND T_SURF FOR THE NEXT TIMESTEP
C     &
C     DETERMINE BETA*Q_SURF FOR THE NEXT TIMESTEP
C     
            DO 335 JL=KIDIA,KFDIA
               ZCFH(JL,KLEV)=ZCFH(JL,KLEV)/(ZTPFAC1*ZTMST*G)
               LO=LALAND(JL).OR.PTSM1M(JL).GT.ZTMELT

C     LG-   replacing the surface temperature by the skin temperature

               IF (LTVEG)
     &              LO=LALAND(JL).OR.TVEGM(JL).GT.ZTMELT

               ZALVS=ALS*PCVS(JL)+ALV*(1.-PCVS(JL))
               ZALVS=CVMGT(ZALVS,ALS,LO)
               ZCPCOR     = ZLATKF(JL)*PTSM1M(JL)*ZCPD*VTMPC2/ZALVS

C     LG-   replacing the surface temperature by the skin temperature

               IF (LTVEG)
     &              ZCPCOR     = (1.-(VEGFRAC(JL)+WSFRAC(JL)))*
     &              ZLATKF(JL)*PTSM1M(JL)*ZCPD*VTMPC2/ZALVS+
     &              (VEGFRAC(JL)+WSFRAC(JL))*
     &              ZLATKF(JL)*TVEGM(JL)*ZCPD*VTMPC2/ZALVS

               ZSENKF(JL) = ZSENKF(JL) + ZCPCOR
               ZS0 = ZTPFAC1*ZTDIF(JL,KLEV) + (1.-ZTPFAC1)*ZCPTGZ(JL,KLEV)
     +              -  ZSENKF(JL)/ZCFH(JL,KLEV)
               ZQ0 = ZTPFAC1*ZQDIF(JL,KLEV)+(1.-ZTPFAC1)*PQM1(JL,KLEV)
     +              - ZLATKF(JL)/ZCFH(JL,KLEV)/ZALVS
               PTSM(JL)=ZS0/(ZCPD*(1+VTMPC2*ZQ0))

C     
               LO1=(PTSM1M(JL)-ZTMELT).GT.0.

C     LG-   replacing the surface temperature by the skin temperature

               IF (LTVEG)
     &              LO1=(TVEGM(JL)-ZTMELT).GT.0.

               ZCVM3=CVMGT(C3LES,C3IES,LO1)
               ZCVM4=CVMGT(C4LES,C4IES,LO1)
               ZCVM5=CVMGT(C5LES,C5IES,LO1)
               IT=PTSM1M(JL)*1000.

C     LG-   replacing the surface temperature by the skin temperature

               IF (LTVEG)
     &              IT=(1.-(VEGFRAC(JL)+WSFRAC(JL)))*PTSM1M(JL)*1000.+
     &              (VEGFRAC(JL)+WSFRAC(JL))*TVEGM(JL)*1000.

               ZQS(JL)=TLUCUA(IT)/PAPHM1(JL,KLEVP1)
               ZQS(JL)=ZQS(JL)/(1.-VTMPC1*ZQS(JL))
               ZBETA(JL)=
     +              (ZQ0    - ZTPFAC1*ZQDIF(JL,KLEV)-(1.-ZTPFAC1)*PQM1(JL,KLEV))
     +              /(ZQS(JL)- ZTPFAC1*ZQDIF(JL,KLEV)-(1.-ZTPFAC1)*PQM1(JL,KLEV))
 335        CONTINUE
C     
C     EVM-----------------------------------------------------------------
C     
C     *         5.83     REDETERMINE SURFACE FLUXES FROM OBSERVED VALUES
C     
C     DO 583 JL=KIDIA,KFDIA
C     AHFS(JL)=AHFSM(JL)+ZSENKF(JL)*ZDIAGT
C     AHFL(JL)=AHFLM(JL)+ZLATKF(JL)*ZDIAGT
C     DHFT(JL)=0.
C     DHFQ(JL)=0.
C     DHFQW(JL)=0.
C     DHFQS(JL)=0.
C     583    CONTINUE
C     ENDIF (LSRFFLUX_READ)
C     EVM-----------------------------------------------------------------
         ENDIF
C=============================================================
C     EVM----------------------------------------------------------------
C     
C     ------------------------------------------------------------------
C     
C     *         6.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
C     --------- ------------ -- ---------- -- ----------
C     
 600     CONTINUE
C***  
      ELSE
C***  
         DO  601 JL=KIDIA,KFDIA
            PAZ0(JL)=PAZ0M(JL)
            PCVGHL(JL)=0.
            PVDIS(JL)=PVDISM(JL)
            PUSTR(JL)=PUSTRM(JL)
            PVSTR(JL)=PVSTRM(JL)
            PAHFS(JL)=PAHFSM(JL)
            PAHFL(JL)=PAHFLM(JL)
            PEVAP(JL)=PEVAPM(JL)
            PTHFL(JL)=0.
            PDHFT(JL)=0.
            PQHFL(JL)=0.
            PXHFL(JL)=0.
            PDHFQW(JL)=0.
            PDHFQS(JL)=0.
            PRSFL(JL)=0.
            PTEMP2(JL)=PTEMP2M(JL)
            PT2MAX(JL)=PT2MAXM(JL)
            PT2MIN(JL)=PT2MINM(JL)
            PDEW2(JL)=PDEW2M(JL)
            PU10(JL)=PU10M(JL)
            PV10(JL)=PV10M(JL)
            PWIND10(JL)=PWIND10M(JL)
            PUSTAR3(JL)=PUSTAR3M(JL)
            PWIMAX(JL)=PWIMAXM(JL)
 601     CONTINUE
C     DVDISZ(IROW)=0.
C     DHFSZ(IROW)=0.
C     DEVAPZ(IROW)=0.
      END IF

C***  
C     
C     ------------------------------------------------------------------
C     

      RETURN
      END

C
C  ******************************************************
C  *  S U B R O U T I N E    W N D P R O F              *
C  *                                                    *
C  * SUBROUTINE WNDPROF DETERMINES THE MEAN WIND SPEED  *
C  * PROFILE WITHIN THE CANOPY USING A FORM OF THE      *
C  * EXPONENTIAL WIND PROFILE (CIONCO).  THIS IS LATE   *
C  * USED IN THE COMPUTATION OF THE LEAF BOUNDARY LAYER *
C  * RESISTANCE TERM                                    *
C  ******************************************************
C
      SUBROUTINE WNDPROF

      IMPLICIT NONE

      INCLUDE 'parchem.h'
      INCLUDE 'comchem.h'

      REAL BETA,ALPHA,D0,Z0

      INTEGER PTYPE,I,K

      PTYPE=IPROF
      
C LG- for the Tropical Rainforest vertical biomass distribution of
C     the Cuieiras and Jaru site, PTYPE resembles the value of IPROF for
C     the top LAD profile 

      IF (IPROF.GE.5) PTYPE=3     
      
C
C     ************************
C     * SET ALPHA AND BETA   *
C     ************************
C

      BETA = 1-(PTYPE-1)*.25
      IF(HC.GT.10) THEN
        ALPHA = LAI
        IF(ALPHA.GT.4.0) ALPHA = 4.0
      ELSE
        ALPHA = .65*LAI
        IF(ALPHA.GT.3.0) ALPHA = 3.0
      ENDIF

C
C    ******************************************************
C    *  COMPUTATION OF D0 AND Z0 AS A FUNCTION            *
C    *  OF LAI BASED ON MODEL COMPUTATIONS OF             *
C    *  SHAW AND PERIERA (1982) AGRICULTURAL METEOROLOGY  *
C    ******************************************************
C

      D0 = HC*(.05+0.5*LAI**0.20 + (PTYPE-1)*.05)
C
      IF(LAI.LT.1) THEN
        Z0 = HC*0.1
      ELSE
        Z0 = HC*(0.23 - 0.1*LAI**0.25 - (PTYPE-1)*.015)
      ENDIF
C
C LG- The displacement height normally estimated being 2/3 the canopy height.
C     However, a method to derive the displacement height and the surface
C     roughness from the canopy height and LAI has been introduced (12-1999)
C     (see subroutine vegetation.f).
C     USTAR is the friction velocity for the vegetation, 

      IF (DISP.LT.HC) D0=DISP
      Z0=Z0M
      U(NLEVV) = 2.5*USTAR(1)*MAX(1.E-10,ALOG((HC-D0+Z0)/Z0))

C
C   *****************************************
C   *  ESTIMATE WITHIN CANOPY WIND PROFILE  *
C   *  USING A FORM OF CIONCO MODEL         *
C   *****************************************
C
      DO 10 I=1,NLEVV
         K=NLEVV+1-I
         U(K)=U(NLEVV)*EXP(-ALPHA*((1-FLOAT(K)/FLOAT(NLEVV))**BETA))
 10   CONTINUE

      RETURN
      END

